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
# band in proportion to its area there. Each band ends at a whole percentage
# point on the cumulative river-surface-area x-axis.
make_bins <- function(part) {
  purrr::map_dfr(seq_len(100), function(bin) {
    overlap <- pmax(0, pmin(part$end_pct, bin) -
                       pmax(part$start_pct, bin - 1))
    total <- sum(overlap)
    if (!is.finite(total) || total <= 0) stop("Empty river-area band")
    tibble(ecoregion = part$ecoregion[1], band = bin,
           mean_gpp = sum(part$gpp_rate * overlap) / total,
           mean_er = sum(part$er_rate * overlap) / total)
  })
}
bins <- ordered %>% group_split(ecoregion) %>% purrr::map_dfr(make_bins)
stopifnot(nrow(bins) == 100 * n_distinct(annual$ecoregion),
          all(is.finite(bins$mean_gpp)), all(is.finite(bins$mean_er)))

bins <- bins %>% mutate(nep = mean_gpp - mean_er) %>%
  group_by(ecoregion) %>%
  mutate(total_absolute_nep = sum(abs(nep)),
         cumulative_nep_pct = 100 * cumsum(nep) / total_absolute_nep,
         area_pct = band) %>%
  ungroup()
stopifnot(all(is.finite(bins$total_absolute_nep)),
          all(bins$total_absolute_nep > 0),
          all(is.finite(bins$cumulative_nep_pct)))
ends <- bins %>% group_by(ecoregion) %>% slice_tail(n = 1) %>% ungroup() %>%
  select(ecoregion, final_signed_nep_pct = cumulative_nep_pct)
message("Final signed NEP as % of total absolute band NEP:")
print(ends, n = Inf)
points <- bind_rows(
  bins %>% select(ecoregion, area_pct, cumulative_nep_pct),
  bins %>% distinct(ecoregion) %>%
    mutate(area_pct = 0, cumulative_nep_pct = 0)) %>%
  arrange(ecoregion, area_pct)
reference <- points %>% left_join(ends, by = "ecoregion") %>%
  mutate(uniform_pct = area_pct * final_signed_nep_pct / 100)

figure <- ggplot(points, aes(area_pct, cumulative_nep_pct)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.45) +
  geom_line(data = reference, aes(y = uniform_pct),
            colour = "grey60", linetype = "dashed", linewidth = 0.55) +
  geom_line(colour = "#253746", linewidth = 0.9) +
  facet_wrap(vars(ecoregion), ncol = 3) +
  scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = c(-100, -50, 0, 50, 100),
                     expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(-100, 100)) +
  labs(x = "Cumulative river surface area (%)",
       y = "Cumulative signed NEP (% of absolute NEP)",
       caption = stringr::str_wrap(paste0(
         "CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea GPP and ER predictions (n = ",
         format(nrow(annual), big.mark = ","), "). Each 1% river-area band contributes area-weighted mean GPP − mean ER. ",
         "Cumulative signed NEP is divided by the sum of absolute band NEP, so curves need not end at 100%. ",
         "Dashed lines show uniform net balance across river area. This oxygen-based proxy is not a carbon balance."), 120)) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F2F4F5", colour = NA),
        strip.text = element_text(colour = "#253746", face = "bold", size = 9),
        panel.spacing = unit(0.9, "lines"),
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(10, 10, 10, 10))
ggsave("figures/46_cumulative_nep_percent_by_ecoregion.png", figure,
       width = 12, height = 8.5, dpi = 300, bg = "white")
