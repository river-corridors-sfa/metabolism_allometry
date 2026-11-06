rm(list = ls())
source("scripts/00_setup.R")

# Matched, complete-year no-uparea Maavara reach predictions. Monthly rates
# are g O2 m^-2 d^-1; ER is a positive respiration magnitude. Each ecoregion
# is partitioned into 100 equal LOCAL river-surface-area bands after ordering
# reaches by upstream drainage area. NEP here is the oxygen-rate proxy GPP - ER;
# no routing or upstream flux is summed.
rates <- readRDS("data/processed/03_filtered_no_uparea_monthly_metabolism.rds") %>%
  transmute(comid = as.character(comid), month = as.integer(month),
            gpp = gpp_g_o2_m2_d, er = er_g_o2_m2_d)
stopifnot(!anyNA(rates[c("comid", "month")]),
          all(rates$month %in% 1:12), !anyDuplicated(rates[c("comid", "month")]))
attrs <- readr::read_csv(
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",
  col_select = c(COMID, lengthkm, merit_width, uparea),
  col_types = cols(COMID = col_character(), .default = col_double()),
  show_col_types = FALSE) %>%
  transmute(comid = COMID, drainage_km2 = uparea,
            river_area_m2 = lengthkm * 1000 * merit_width)
regions <- readRDS("data/processed/04c_no_uparea_comid_huc08_key.rds") %>%
  transmute(comid = as.character(comid), ecoregion = l1_ecoregion)
stopifnot(!anyDuplicated(attrs$comid), !anyDuplicated(regions$comid))

monthly <- rates %>% left_join(attrs, by = "comid") %>%
  left_join(regions, by = "comid") %>%
  mutate(days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[month],
         valid = is.finite(gpp) & gpp >= 0 & is.finite(er) & er >= 0 &
           is.finite(drainage_km2) & drainage_km2 > 0 &
           is.finite(river_area_m2) & river_area_m2 > 0 & !is.na(ecoregion))
complete_ids <- monthly %>% group_by(comid) %>%
  summarise(complete = n_distinct(month) == 12 && all(valid), .groups = "drop") %>%
  filter(complete) %>% pull(comid)
message("Complete GPP and ER reaches: ", length(complete_ids), " / ",
        n_distinct(monthly$comid), "; excluded: ",
        n_distinct(monthly$comid) - length(complete_ids))
if (!length(complete_ids)) stop("No complete matched reaches")

annual <- monthly %>% filter(comid %in% complete_ids) %>%
  group_by(comid, ecoregion, drainage_km2) %>%
  summarise(river_area_m2 = sum(river_area_m2 * days) / sum(days),
            gpp_rate = sum(gpp * days) / sum(days),
            er_rate = sum(er * days) / sum(days), .groups = "drop")
stopifnot(nrow(annual) == length(complete_ids))
ordered <- annual %>% group_by(ecoregion) %>%
  arrange(drainage_km2, comid, .by_group = TRUE) %>%
  mutate(end_pct = 100 * cumsum(river_area_m2) / sum(river_area_m2),
         start_pct = lag(end_pct, default = 0)) %>%
  ungroup()

# Overlap in percentage points is proportional to actual river surface area
# within an ecoregion. A reach crossing band boundaries contributes to each
# band in proportion to its area there. The x position is the area-weighted
# geometric mean drainage area, suitable for a log-scale axis.
make_bins <- function(part) {
  purrr::map_dfr(seq_len(100), function(bin) {
    overlap <- pmax(0, pmin(part$end_pct, bin) -
                       pmax(part$start_pct, bin - 1))
    total <- sum(overlap)
    if (!is.finite(total) || total <= 0) stop("Empty river-area band")
    tibble(ecoregion = part$ecoregion[1], band = bin,
           drainage_km2 = exp(sum(log(part$drainage_km2) * overlap) / total),
           mean_gpp = sum(part$gpp_rate * overlap) / total,
           mean_er = sum(part$er_rate * overlap) / total)
  })
}
bins <- ordered %>% group_split(ecoregion) %>% purrr::map_dfr(make_bins)
stopifnot(nrow(bins) == 100 * n_distinct(annual$ecoregion),
          all(is.finite(bins$drainage_km2)), all(bins$drainage_km2 > 0),
          all(is.finite(bins$mean_gpp)), all(is.finite(bins$mean_er)))

bins <- bins %>% mutate(nep = mean_gpp - mean_er)
stopifnot(all(is.finite(bins$nep)))
summary <- bins %>% group_by(ecoregion) %>%
  summarise(bands = n(), positive_nep_bands = sum(nep > 0),
            fraction_positive_nep = positive_nep_bands / bands,
            .groups = "drop")
message("Fraction of equal-area bands with NEP > 0:")
print(summary, n = Inf)

figure <- ggplot(bins, aes(drainage_km2, nep)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed",
             linewidth = 0.55) +
  geom_line(colour = "#007F5B", linewidth = 0.85) +
  geom_point(colour = "#007F5B", size = 0.8) +
  facet_wrap(vars(ecoregion), ncol = 3, scales = "free_y") +
  scale_x_log10(labels = scales::label_number_auto()) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Upstream drainage area (km²; log scale)",
       y = "NEP proxy: mean GPP − mean ER (g O₂ m⁻² d⁻¹)",
       caption = stringr::str_wrap(paste0(
         "CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea GPP and ER predictions (n = ",
         format(nrow(annual), big.mark = ","), "). Each point represents 1% of an ecoregion's river surface area. ",
         "Positive values indicate modeled net oxygen production; negative values indicate modeled net oxygen consumption. ",
         "This is an oxygen-rate NEP proxy, not a carbon balance."), 120)) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F2F4F5", colour = NA),
        strip.text = element_text(colour = "#253746", face = "bold", size = 9),
        panel.spacing = unit(0.9, "lines"),
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(10, 10, 10, 10))
ggsave("figures/45_nep_by_drainage_ecoregion.png", figure,
       width = 12, height = 8.5, dpi = 300, bg = "white")
