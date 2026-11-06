rm(list = ls())
source("scripts/00_setup.R")

# CONUS-filtered Maavara no-uparea predictions; local reach contributions only.
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
  transmute(comid = as.character(comid), ecoregion = l1_ecoregion,
            watershed_km2 = unitarea_km2)
stopifnot(!anyDuplicated(attrs$comid), !anyDuplicated(regions$comid))

monthly <- rates %>% left_join(attrs, by = "comid") %>%
  left_join(regions, by = "comid") %>%
  mutate(days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[month],
         valid = is.finite(gpp) & gpp >= 0 & is.finite(er) & er >= 0 &
           is.finite(drainage_km2) & drainage_km2 > 0 &
           is.finite(river_area_m2) & river_area_m2 > 0 &
           is.finite(watershed_km2) & watershed_km2 > 0 & !is.na(ecoregion))
complete_ids <- monthly %>% group_by(comid) %>%
  summarise(complete = n_distinct(month) == 12 && all(valid), .groups = "drop") %>%
  filter(complete) %>% pull(comid)
message("Complete GPP and ER reaches: ", length(complete_ids), " / ",
        n_distinct(monthly$comid), "; excluded: ",
        n_distinct(monthly$comid) - length(complete_ids))
if (!length(complete_ids)) stop("No complete reaches")

annual <- monthly %>% filter(comid %in% complete_ids) %>%
  group_by(comid, ecoregion, drainage_km2, watershed_km2) %>%
  summarise(river_area_m2 = sum(river_area_m2 * days) / sum(days),
            GPP = sum(gpp * river_area_m2 * days),
            ER = sum(er * river_area_m2 * days), .groups = "drop")
stopifnot(nrow(annual) == length(complete_ids))

# Incremental MERIT catchment area (unitarea), in km2, is additive across
# reaches. Summing upstream drainage area would repeatedly count shared land.
# Reaches are still ordered by upstream drainage area, as in the prior plots.
ordered <- annual %>% group_by(ecoregion) %>%
  arrange(drainage_km2, comid, .by_group = TRUE) %>%
  mutate(cumulative_watershed_km2 = cumsum(watershed_km2),
         cumulative_gpp_gg_y = cumsum(GPP) / 1e9,
         cumulative_er_gg_y = cumsum(ER) / 1e9,
         cumulative_gpp_er = if_else(cumsum(ER) > 0,
                                      cumsum(GPP) / cumsum(ER), NA_real_)) %>%
  ungroup()
stopifnot(all(is.finite(ordered$cumulative_watershed_km2)),
          all(is.finite(ordered$cumulative_gpp_gg_y)),
          all(is.finite(ordered$cumulative_er_gg_y)),
          all(is.finite(ordered$cumulative_gpp_er)))
totals <- ordered %>% group_by(ecoregion) %>% slice_tail(n = 1) %>% ungroup()
message("Ecoregion totals: incremental watershed km2, annual GPP/ER Gg O2, and cumulative ratio:")
print(totals %>% select(ecoregion, cumulative_watershed_km2,
                        cumulative_gpp_gg_y, cumulative_er_gg_y,
                        cumulative_gpp_er), n = Inf)

flux <- ordered %>%
  select(ecoregion, cumulative_watershed_km2,
         GPP = cumulative_gpp_gg_y, ER = cumulative_er_gg_y) %>%
  pivot_longer(c(GPP, ER), names_to = "process", values_to = "value")
flux <- bind_rows(
  flux %>% distinct(ecoregion, process) %>%
    mutate(cumulative_watershed_km2 = 0, value = 0), flux) %>%
  mutate(process = factor(process, levels = c("GPP", "ER"))) %>%
  arrange(process, ecoregion, cumulative_watershed_km2)
ratio <- ordered %>%
  transmute(ecoregion, cumulative_watershed_km2,
            value = cumulative_gpp_er)
palette <- setNames(viridisLite::viridis(n_distinct(ordered$ecoregion), end = 0.9),
                    sort(unique(ordered$ecoregion)))
base_theme <- theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "#F2F4F5", colour = NA),
        strip.text = element_text(colour = "#253746", face = "bold", size = 12),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "bottom", legend.box = "vertical",
        plot.margin = margin(10, 10, 10, 10))
flux_plot <- ggplot(flux, aes(cumulative_watershed_km2, value,
                              colour = ecoregion, group = ecoregion)) +
  geom_line(linewidth = 0.85) +
  facet_wrap(vars(process), nrow = 1, scales = "fixed") +
  scale_colour_manual(values = palette) +
  scale_x_continuous(labels = scales::label_number(big.mark = ","),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(labels = scales::label_number(big.mark = ","),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cumulative incremental watershed area (km²)",
       y = "Cumulative annual metabolism (Gg O₂ yr⁻¹)",
       colour = "Level I ecoregion") + base_theme +
  theme(legend.position = "none")
ratio_plot <- ggplot(ratio, aes(cumulative_watershed_km2, value,
                                colour = ecoregion, group = ecoregion)) +
  geom_hline(yintercept = 1, colour = "grey75", linetype = "dashed") +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = palette) +
  scale_x_continuous(labels = scales::label_number(big.mark = ","),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cumulative incremental watershed area (km²)",
       y = "Cumulative GPP:ER ratio", colour = "Level I ecoregion") +
  base_theme +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none") +
  ggtitle("GPP:ER")
caption <- stringr::str_wrap(paste0(
  "CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea GPP and ER predictions (n = ",
  format(nrow(annual), big.mark = ","), "). Reaches ordered by upstream drainage area within each Level I ecoregion. ",
  "X accumulates non-overlapping MERIT reach catchments (unitarea); ratio is cumulative GPP divided by cumulative ER."), 170)
legend_data <- tibble(ecoregion = names(palette),
                      x = rep(1:5, times = 2),
                      y = rep(c(2, 1), each = 5))
legend <- ggplot(legend_data, aes(x, y, colour = ecoregion)) +
  geom_point(size = 2.3) +
  geom_text(aes(label = ecoregion), hjust = 0, nudge_x = 0.08,
            size = 3.0, show.legend = FALSE) +
  scale_colour_manual(values = palette, guide = "none") +
  scale_x_continuous(limits = c(0.9, 6.1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.7, 2.3), expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") + theme_void()
panels <- cowplot::plot_grid(flux_plot, ratio_plot, nrow = 1,
                            rel_widths = c(2, 1), align = "h", axis = "tb")
figure <- cowplot::plot_grid(
  panels, legend,
  cowplot::ggdraw() + cowplot::draw_label(caption, x = 0.015, hjust = 0,
                                          size = 9, colour = "grey35"),
  ncol = 1, rel_heights = c(1, 0.14, 0.06))
ggsave("figures/38_watershed_area_metabolism_ratio_by_ecoregion.png", figure,
       width = 19, height = 7.5, dpi = 300, bg = "white")
