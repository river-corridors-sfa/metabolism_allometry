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
if (!length(complete_ids)) stop("No complete reaches")

annual <- monthly %>% filter(comid %in% complete_ids) %>%
  group_by(comid, ecoregion, drainage_km2) %>%
  summarise(river_area_m2 = sum(river_area_m2 * days) / sum(days),
            GPP = sum(gpp * river_area_m2 * days),
            ER = sum(er * river_area_m2 * days), .groups = "drop")
stopifnot(nrow(annual) == length(complete_ids))

# Keep each ecoregion's own running totals, with identical axes in both panels.
# 1 km2 = 1e6 m2; 1 Gg O2 = 1e9 g O2.
ordered <- annual %>% group_by(ecoregion) %>%
  arrange(drainage_km2, comid, .by_group = TRUE) %>%
  mutate(cumulative_river_km2 = cumsum(river_area_m2) / 1e6) %>%
  ungroup()
curves <- ordered %>%
  pivot_longer(c(GPP, ER), names_to = "process", values_to = "annual_g_o2") %>%
  group_by(ecoregion, process) %>%
  mutate(cumulative_g_o2_y = cumsum(annual_g_o2),
         cumulative_gg_o2_y = cumulative_g_o2_y / 1e9) %>%
  ungroup()
stopifnot(all(is.finite(curves$cumulative_river_km2)),
          all(is.finite(curves$cumulative_gg_o2_y)),
          all(curves$cumulative_gg_o2_y >= 0))
starts <- curves %>% distinct(ecoregion, process) %>%
  mutate(cumulative_river_km2 = 0, cumulative_gg_o2_y = 0)
plot_data <- bind_rows(
  starts %>% select(ecoregion, process, cumulative_river_km2, cumulative_gg_o2_y),
  curves %>% select(ecoregion, process, cumulative_river_km2, cumulative_gg_o2_y)) %>%
  mutate(process = factor(process, levels = c("GPP", "ER"))) %>%
  arrange(process, ecoregion, cumulative_river_km2)
totals <- curves %>% group_by(ecoregion, process) %>%
  slice_tail(n = 1) %>% ungroup()
message("Absolute annual ecoregion totals (km2 river surface; Gg O2/year):")
print(totals %>% select(ecoregion, process, cumulative_river_km2,
                        cumulative_gg_o2_y), n = Inf)

palette <- setNames(viridisLite::viridis(n_distinct(curves$ecoregion), end = 0.9),
                    sort(unique(curves$ecoregion)))
figure <- ggplot(plot_data,
                 aes(cumulative_river_km2, cumulative_gg_o2_y,
                     colour = ecoregion, group = ecoregion)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(vars(process), nrow = 1, scales = "fixed") +
  scale_colour_manual(values = palette) +
  scale_x_continuous(labels = scales::label_number(big.mark = ","),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(labels = scales::label_number(big.mark = ","),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cumulative river surface area (km²)",
       y = "Cumulative annual metabolism (Gg O₂ yr⁻¹)", colour = "Level I ecoregion",
       caption = stringr::str_wrap(paste0(
         "CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea GPP and ER predictions (n = ",
         format(nrow(annual), big.mark = ","), "). Reaches ordered by upstream drainage area within each ecoregion. ",
         "Each curve ends at that ecoregion's absolute river surface area and annual oxygen flux."), 135)) +
  theme_classic(base_size = 12) +
  theme(strip.background = element_rect(fill = "#F2F4F5", colour = NA),
        strip.text = element_text(colour = "#253746", face = "bold", size = 12),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "bottom", legend.box = "vertical",
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(10, 10, 10, 10)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
ggsave("figures/37_absolute_metabolism_by_ecoregion.png", figure,
       width = 14, height = 6.5, dpi = 300, bg = "white")
