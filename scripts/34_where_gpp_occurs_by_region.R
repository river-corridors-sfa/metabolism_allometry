rm(list = ls())
source("scripts/00_setup.R")

# Same CONUS-filtered Maavara no-uparea predictions as the main figures.
# Local GPP rate is g O2 m^-2 d^-1; lengthkm * 1000 * merit_width is
# local river surface area (m2), distinct from upstream drainage area (km2).
rates <- readRDS("data/processed/03_filtered_no_uparea_monthly_metabolism.rds") %>%
  transmute(comid = as.character(comid), month = as.integer(month),
            gpp = gpp_g_o2_m2_d)
stopifnot(!anyNA(rates[c("comid", "month")]),
          all(rates$month %in% 1:12),
          !anyDuplicated(rates[c("comid", "month")]))

attrs <- readr::read_csv(
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",
  col_select = c(COMID, lengthkm, merit_width, uparea),
  col_types = cols(COMID = col_character(), .default = col_double()),
  show_col_types = FALSE) %>%
  transmute(comid = COMID, river_area_m2 = lengthkm * 1000 * merit_width,
            drainage_km2 = uparea)
key <- readRDS("data/processed/04c_no_uparea_comid_huc08_key.rds") %>%
  transmute(comid = as.character(comid), ecoregion = l1_ecoregion,
            huc2 = substr(huc8, 1, 2))
stopifnot(!anyDuplicated(attrs$comid), !anyDuplicated(key$comid))

monthly <- rates %>%
  left_join(attrs, by = "comid") %>%
  left_join(key, by = "comid") %>%
  mutate(days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[month],
         valid = is.finite(gpp) & gpp >= 0 & is.finite(river_area_m2) &
           river_area_m2 > 0 & is.finite(drainage_km2) & drainage_km2 > 0 &
           !is.na(ecoregion) & !is.na(huc2))
complete_ids <- monthly %>% group_by(comid) %>%
  summarise(complete = n_distinct(month) == 12 && all(valid), .groups = "drop") %>%
  filter(complete) %>% pull(comid)
message("GPP full-year complete reaches: ", length(complete_ids), " / ",
        n_distinct(monthly$comid), "; excluded: ",
        n_distinct(monthly$comid) - length(complete_ids),
        "; invalid GPP reach-months: ", sum(!is.finite(monthly$gpp) | monthly$gpp < 0))
if (!length(complete_ids)) stop("No complete GPP reaches")

annual <- monthly %>% filter(comid %in% complete_ids) %>%
  group_by(comid, ecoregion, huc2, drainage_km2) %>%
  summarise(river_area_m2 = sum(river_area_m2 * days) / sum(days),
            gpp_g_o2_y = sum(gpp * river_area_m2 * days), .groups = "drop")
stopifnot(nrow(annual) == length(complete_ids),
          all(annual$gpp_g_o2_y >= 0), all(annual$river_area_m2 > 0))

make_figure <- function(group_col, label, path, columns) {
  curves <- annual %>%
    mutate(region = .data[[group_col]]) %>%
    group_by(region) %>%
    arrange(drainage_km2, comid, .by_group = TRUE) %>%
    mutate(area_pct = 100 * cumsum(river_area_m2) / sum(river_area_m2),
           gpp_pct = 100 * cumsum(gpp_g_o2_y) / sum(gpp_g_o2_y)) %>%
    ungroup()
  stopifnot(all(is.finite(curves$gpp_pct)))
  ends <- curves %>% group_by(region) %>% slice_tail(n = 1) %>% ungroup()
  stopifnot(all(abs(ends$area_pct - 100) < 1e-8),
            all(abs(ends$gpp_pct - 100) < 1e-8))
  starts <- curves %>% distinct(region) %>%
    mutate(area_pct = 0, gpp_pct = 0)
  plot_data <- bind_rows(starts, curves %>% select(region, area_pct, gpp_pct)) %>%
    arrange(region, area_pct)
  contributions <- plot_data %>% group_by(region) %>%
    summarise(final20 = 100 - approx(area_pct, gpp_pct, xout = 80)$y,
              .groups = "drop") %>%
    mutate(annotation = sprintf("Final 20%%: %.1f%% GPP", final20))
  message(label, " final-20% GPP contributions:")
  print(contributions, n = Inf)

  # Each reach contributes one observation at the midpoint of its surface-area
  # interval. Counts per 5 percentage-point bin reveal whether many small
  # reaches cluster at a particular position along the cumulative x-axis.
  # Scale each region's tallest bin to 28% of panel height: histogram heights
  # are comparable within a region, not between regions.
  histogram <- curves %>% group_by(region) %>%
    mutate(midpoint_pct = area_pct - 50 * river_area_m2 / sum(river_area_m2),
           bin = pmin(19L, floor(midpoint_pct / 5))) %>%
    count(region, bin, name = "n_reaches") %>%
    group_by(region) %>%
    mutate(height = 28 * n_reaches / max(n_reaches), x = 2.5 + 5 * bin) %>%
    ungroup()
  region_colours <- setNames(viridisLite::viridis(n_distinct(curves$region),
                                                   option = "D", end = 0.85),
                             sort(unique(curves$region)))

  figure <- ggplot(plot_data, aes(area_pct, gpp_pct)) +
    {if (group_col == "ecoregion")
      geom_col(data = histogram, aes(x = x, y = height, fill = region),
               inherit.aes = FALSE, width = 5, alpha = 0.28)
      else NULL} +
    geom_abline(slope = 1, intercept = 0, colour = "grey65",
                linetype = "dashed", linewidth = 0.4) +
    geom_line(colour = "#009E73", linewidth = 0.75) +
    geom_vline(xintercept = 80, colour = "#253746", linewidth = 1) +
    {if (group_col == "huc2") geom_text(data = contributions,
              aes(x = 98, y = 10, label = annotation),
              inherit.aes = FALSE, hjust = 1, colour = "#007F5B", size = 2.7)
      else NULL} +
    facet_wrap(vars(region), ncol = columns) +
    scale_fill_manual(values = region_colours, guide = "none") +
    scale_x_continuous(breaks = c(0, 50, 100),
                       expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = c(0, 50, 100),
                       expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(x = "Cumulative river surface area (%)",
         y = "Cumulative annual GPP (%)",
         caption = stringr::str_wrap(paste0(
           "CONUS-filtered MERIT reaches; complete 12-month Maavara no-uparea GPP predictions (n = ",
           format(nrow(annual), big.mark = ","), "). Grouped by ", label,
           "; reaches ordered by upstream drainage area within each group.",
           if (group_col == "ecoregion")
             " Bars show reach counts per 5% x-axis bin, scaled to each panel's tallest bin."
           else ""), 115)) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_rect(fill = "#F2F4F5", colour = NA),
          strip.text = element_text(colour = "#253746", face = "bold", size = 9),
          panel.spacing = unit(0.9, "lines"),
          plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
          plot.margin = margin(10, 10, 10, 10))
  ggsave(path, figure, width = 12, height = if (group_col == "huc2") 11 else 8,
         dpi = 300, bg = "white")
}

make_figure("ecoregion", "Level I ecoregion",
            "figures/34_gpp_cumulative_by_ecoregion.png", 3)
make_figure("huc2", "HUC2",
            "figures/34_gpp_cumulative_by_huc2.png", 4)
