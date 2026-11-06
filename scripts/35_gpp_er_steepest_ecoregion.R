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

ordered <- annual %>% group_by(ecoregion) %>%
  arrange(drainage_km2, comid, .by_group = TRUE) %>%
  mutate(area_pct = 100 * cumsum(river_area_m2) / sum(river_area_m2),
         area_mid_pct = area_pct - 50 * river_area_m2 / sum(river_area_m2)) %>%
  ungroup()
curves <- ordered %>%
  pivot_longer(c(GPP, ER), names_to = "process", values_to = "annual_g_o2") %>%
  group_by(ecoregion, process) %>%
  mutate(metabolism_pct = 100 * cumsum(annual_g_o2) / sum(annual_g_o2),
         # Delta metabolism % / delta river area % for this single reach.
         slope = (annual_g_o2 / sum(annual_g_o2)) /
           (river_area_m2 / sum(river_area_m2)),
         metabolism_mid_pct = metabolism_pct -
           50 * annual_g_o2 / sum(annual_g_o2)) %>%
  ungroup()
stopifnot(all(is.finite(curves$slope)), all(curves$slope >= 0))
ends <- curves %>% group_by(ecoregion, process) %>% slice_tail(n = 1) %>% ungroup()
stopifnot(all(abs(ends$area_pct - 100) < 1e-8),
          all(abs(ends$metabolism_pct - 100) < 1e-8))

# Ties are broken by smaller drainage area, then COMID. Circles are placed
# halfway along the selected empirical reach increment.
peaks <- curves %>% group_by(ecoregion, process) %>%
  arrange(desc(slope), drainage_km2, comid, .by_group = TRUE) %>%
  slice_head(n = 1) %>% ungroup()
message("Steepest reach increments by ecoregion and process:")
print(peaks %>% select(ecoregion, process, comid, drainage_km2, slope), n = Inf)

starts <- curves %>% distinct(ecoregion, process) %>%
  mutate(area_pct = 0, metabolism_pct = 0)
plot_data <- bind_rows(starts,
                       curves %>% select(ecoregion, process, area_pct, metabolism_pct)) %>%
  arrange(ecoregion, process, area_pct)
figure <- ggplot(plot_data, aes(area_pct, metabolism_pct, colour = process)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey68",
              linetype = "dashed", linewidth = 0.45) +
  geom_line(linewidth = 0.78) +
  geom_point(data = peaks,
             aes(x = area_mid_pct, y = metabolism_mid_pct, colour = process),
             inherit.aes = FALSE, shape = 21, fill = "white",
             stroke = 1.25, size = 3.1, show.legend = FALSE) +
  facet_wrap(vars(ecoregion), ncol = 3) +
  scale_colour_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = c(0, 50, 100), expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(x = "Cumulative river surface area (%)",
       y = "Cumulative annual metabolism (%)", colour = NULL,
       caption = stringr::str_wrap(paste0(
         "CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea GPP and ER predictions (n = ",
         format(nrow(annual), big.mark = ","), "). Reaches ordered by upstream drainage area within each Level I ecoregion. ",
         "Circles mark the largest single-reach rise in metabolism percentage per rise in river surface-area percentage."), 115)) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.justification = "left",
        strip.background = element_rect(fill = "#F2F4F5", colour = NA),
        strip.text = element_text(colour = "#253746", face = "bold", size = 9),
        panel.spacing = unit(0.9, "lines"),
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(10, 10, 10, 10))
ggsave("figures/35_gpp_er_steepest_ecoregion.png", figure,
       width = 12, height = 8.5, dpi = 300, bg = "white")
