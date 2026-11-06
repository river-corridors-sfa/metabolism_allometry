rm(list = ls())
source("scripts/00_setup.R")

# Median-based percentage curves. Monthly rates are g O2 m^-2 d^-1.
# The y-axis is a robust median-rate contribution proxy, not summed local flux.
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
        n_distinct(monthly$comid))

annual <- monthly %>% filter(comid %in% complete_ids) %>%
  group_by(comid, ecoregion, drainage_km2) %>%
  summarise(river_area_m2 = sum(river_area_m2 * days) / sum(days),
            GPP = sum(gpp * days) / sum(days),
            ER = sum(er * days) / sum(days), .groups = "drop")
stopifnot(nrow(annual) == length(complete_ids))

ordered <- annual %>% group_by(ecoregion) %>%
  arrange(drainage_km2, comid, .by_group = TRUE) %>%
  mutate(end_pct = 100 * cumsum(river_area_m2) / sum(river_area_m2),
         start_pct = lag(end_pct, default = 0)) %>%
  ungroup()

# For each 1%-wide surface-area band, use the area-weighted median annual
# areal rate of all reaches overlapping that band. A reach crossing a band
# boundary contributes only its overlapping river area to the median weight.
# Every band has the same area, so accumulating band median rates gives the
# median-based share of the total. This does not equal actual oxygen flux.
weighted_median <- function(values, weights) {
  keep <- weights > 0
  if (!any(keep)) stop("Empty surface-area bin")
  values <- values[keep]
  weights <- weights[keep]
  order <- order(values)
  values[order][which(cumsum(weights[order]) >= sum(weights) / 2)[1]]
}
make_bins <- function(part) {
  map_dfr(seq_len(100), function(bin) {
    overlap <- pmax(0, pmin(part$end_pct, bin) -
                       pmax(part$start_pct, bin - 1))
    tibble(ecoregion = part$ecoregion[1], area_pct = bin,
           GPP = weighted_median(part$GPP, overlap),
           ER = weighted_median(part$ER, overlap))
  })
}
bins <- ordered %>% group_split(ecoregion) %>% map_dfr(make_bins)
stopifnot(nrow(bins) == 100 * n_distinct(annual$ecoregion),
          all(is.finite(bins$GPP)), all(is.finite(bins$ER)))
curves <- bins %>% group_by(ecoregion) %>%
  mutate(gpp_pct = 100 * cumsum(GPP) / sum(GPP),
         er_pct = 100 * cumsum(ER) / sum(ER)) %>%
  ungroup()
ends <- curves %>% group_by(ecoregion) %>% slice_tail(n = 1) %>% ungroup()
stopifnot(all(abs(ends$gpp_pct - 100) < 1e-8),
          all(abs(ends$er_pct - 100) < 1e-8))
points <- curves %>% transmute(ecoregion, x = area_pct,
                                gpp_y = gpp_pct, er_y = er_pct) %>%
  bind_rows(curves %>% distinct(ecoregion) %>%
              transmute(ecoregion, x = 0, gpp_y = 0, er_y = 0)) %>%
  arrange(ecoregion, x)

# Split shaded bands at interpolated crossings.
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
             band = paste(ecoregion, id, sep = "_"))
  })
}
bands <- points %>% group_split(ecoregion) %>% map_dfr(make_bands)
line_data <- points %>%
  pivot_longer(c(gpp_y, er_y), names_to = "process", values_to = "y") %>%
  mutate(process = recode(process, gpp_y = "GPP", er_y = "ER"))

figure <- ggplot() +
  geom_ribbon(data = bands,
              aes(x = x, ymin = pmin(gpp_y, er_y), ymax = pmax(gpp_y, er_y),
                  fill = higher, group = band), alpha = 0.16) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60",
              linetype = "dashed", linewidth = 0.5) +
  geom_line(data = line_data, aes(x, y, colour = process, group = process),
            linewidth = 0.85) +
  facet_wrap(vars(ecoregion), ncol = 3) +
  scale_colour_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_fill_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = c(0, 50, 100), expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "Cumulative river surface area (%)",
       y = "Cumulative median-based metabolism proxy (%)",
       colour = "Process", fill = "Higher proxy",
       caption = stringr::str_wrap(paste0(
         "CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea GPP and ER predictions (n = ",
         format(nrow(annual), big.mark = ","), "). Reaches ordered by upstream drainage area within each Level I ecoregion. ",
         "Each 1% surface-area band uses its area-weighted median annual areal rate; cumulative shares are normalized within each process and ecoregion. ",
         "This median-based proxy is not observed or predicted total oxygen flux."), 120)) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F2F4F5", colour = NA),
        strip.text = element_text(colour = "#253746", face = "bold", size = 9),
        panel.spacing = unit(0.9, "lines"),
        legend.position = "top", legend.box = "horizontal",
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(10, 10, 10, 10))
ggsave("figures/43_median_rate_percent_curves_ecoregion.png", figure,
       width = 12, height = 8.5, dpi = 300, bg = "white")
