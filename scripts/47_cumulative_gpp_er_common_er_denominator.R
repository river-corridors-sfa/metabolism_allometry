rm(list = ls())
source("scripts/00_setup.R")

# Matched, complete-year no-uparea Maavara reach predictions. Monthly rates
# are g O2 m^-2 d^-1; ER is a positive respiration magnitude. Each ecoregion
# is partitioned into 100 equal LOCAL river-surface-area bands after ordering
# reaches by upstream drainage area. Both cumulative rate curves use the same
# total ER denominator; no routing or upstream flux is summed.
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

# Each band covers exactly 1% of river surface area. Thus summing its
# area-weighted mean rate is proportional to summing its local oxygen flux,
# while the inputs and calculation remain in annual mean areal-rate units.
# Both curves use the SAME denominator: sum of all 100 ER band means.
bins <- bins %>% group_by(ecoregion) %>%
  mutate(total_er = sum(mean_er),
         area_pct = band,
         gpp_pct_of_er = 100 * cumsum(mean_gpp) / total_er,
         er_pct_of_er = 100 * cumsum(mean_er) / total_er) %>%
  ungroup()
stopifnot(all(is.finite(bins$total_er)), all(bins$total_er > 0),
          all(is.finite(bins$gpp_pct_of_er)),
          all(is.finite(bins$er_pct_of_er)))
ends <- bins %>% group_by(ecoregion) %>% slice_tail(n = 1) %>% ungroup()
stopifnot(all(abs(ends$er_pct_of_er - 100) < 1e-8))
message("Final cumulative GPP and ER as % of total ER:")
print(ends %>% select(ecoregion, gpp_pct_of_er, er_pct_of_er), n = Inf)
points <- bind_rows(
  bins %>% transmute(ecoregion, x = area_pct,
                     gpp_y = gpp_pct_of_er, er_y = er_pct_of_er),
  bins %>% distinct(ecoregion) %>%
    mutate(x = 0, gpp_y = 0, er_y = 0)) %>%
  arrange(ecoregion, x)

# Split the light fill where cumulative GPP crosses cumulative ER.
make_bands <- function(part) {
  gap <- part$gpp_y - part$er_y
  crossing <- which(head(gap, -1) * tail(gap, -1) < 0)
  if (length(crossing)) {
    fraction <- -gap[crossing] / (gap[crossing + 1] - gap[crossing])
    extra <- tibble(
      ecoregion = part$ecoregion[1],
      x = part$x[crossing] + fraction * (part$x[crossing + 1] - part$x[crossing]),
      gpp_y = part$gpp_y[crossing] + fraction *
        (part$gpp_y[crossing + 1] - part$gpp_y[crossing]),
      er_y = part$er_y[crossing] + fraction *
        (part$er_y[crossing + 1] - part$er_y[crossing]))
    part <- bind_rows(part, extra) %>% arrange(x)
  }
  gap <- part$gpp_y - part$er_y
  higher <- if_else(head(gap, -1) + tail(gap, -1) >= 0, "GPP", "ER")
  run <- cumsum(c(1L, as.integer(diff(match(higher, c("GPP", "ER"))) != 0)))
  map_dfr(unique(run), function(id) {
    interval <- which(run == id)
    part[min(interval):(max(interval) + 1), ] %>%
      mutate(higher = higher[min(interval)],
             ribbon = paste(ecoregion, id, sep = "_"))
  })
}
ribbons <- points %>% group_split(ecoregion) %>% map_dfr(make_bands)
lines <- points %>%
  pivot_longer(c(gpp_y, er_y), names_to = "process", values_to = "y") %>%
  mutate(process = recode(process, gpp_y = "GPP", er_y = "ER"))
upper <- max(105, ceiling(max(bins$gpp_pct_of_er) / 10) * 10)
figure <- ggplot() +
  geom_ribbon(data = ribbons,
              aes(x = x, ymin = pmin(gpp_y, er_y), ymax = pmax(gpp_y, er_y),
                  fill = higher, group = ribbon), alpha = 0.16) +
  geom_abline(slope = 1, intercept = 0, colour = "grey65",
              linetype = "dashed", linewidth = 0.5) +
  geom_line(data = lines, aes(x, y, colour = process, group = process),
            linewidth = 0.85) +
  facet_wrap(vars(ecoregion), ncol = 3) +
  scale_colour_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_fill_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = c(0, 50, 100), expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, upper)) +
  labs(x = "Cumulative river surface area (%)",
       y = "Cumulative annual areal rate (% of total ER)",
       colour = "Process", fill = "Higher cumulative total",
       caption = stringr::str_wrap(paste0(
         "CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea GPP and ER predictions (n = ",
         format(nrow(annual), big.mark = ","), "). Within each Level I ecoregion, reaches are ordered by upstream drainage area. ",
         "Each 1% river-area band contributes its area-weighted mean annual areal rate; both curves use the sum of ER band means as denominator. ",
         "ER ends at 100%; GPP may end above or below 100%."), 120)) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F2F4F5", colour = NA),
        strip.text = element_text(colour = "#253746", face = "bold", size = 9),
        panel.spacing = unit(0.9, "lines"),
        legend.position = "top", legend.box = "horizontal",
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(10, 10, 10, 10))
ggsave("figures/47_cumulative_gpp_er_common_er_denominator.png", figure,
       width = 12, height = 8.5, dpi = 300, bg = "white")
