rm(list = ls())
source("scripts/00_setup.R")

# Annual local reach GPP - ER, grouped by HUC8. This is net oxygen production
# from predicted river metabolism (an oxygen-based NEP proxy), not a watershed
# carbon budget. Only reaches with valid GPP and ER in all 12 months are used.
rates <- readRDS("data/processed/03_filtered_no_uparea_monthly_metabolism.rds") %>%
  transmute(comid = as.character(comid), month = as.integer(month),
            gpp = gpp_g_o2_m2_d, er = er_g_o2_m2_d)
stopifnot(!anyNA(rates[c("comid", "month")]),
          all(rates$month %in% 1:12), !anyDuplicated(rates[c("comid", "month")]))
attrs <- readr::read_csv(
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",
  col_select = c(COMID, lengthkm, merit_width),
  col_types = cols(COMID = col_character(), .default = col_double()),
  show_col_types = FALSE) %>%
  transmute(comid = COMID, river_area_m2 = lengthkm * 1000 * merit_width)
key <- readRDS("data/processed/04c_no_uparea_comid_huc08_key.rds") %>%
  transmute(comid = as.character(comid), huc8 = as.character(huc8))
stopifnot(!anyDuplicated(attrs$comid), !anyDuplicated(key$comid))

monthly <- rates %>% left_join(attrs, by = "comid") %>%
  left_join(key, by = "comid") %>%
  mutate(days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[month],
         valid = is.finite(gpp) & gpp >= 0 & is.finite(er) & er >= 0 &
           is.finite(river_area_m2) & river_area_m2 > 0 & !is.na(huc8))
complete_ids <- monthly %>% group_by(comid) %>%
  summarise(complete = n_distinct(month) == 12 && all(valid), .groups = "drop") %>%
  filter(complete) %>% pull(comid)
message("Complete GPP and ER reaches: ", length(complete_ids), " / ",
        n_distinct(monthly$comid), "; excluded: ",
        n_distinct(monthly$comid) - length(complete_ids))

huc_totals <- monthly %>% filter(comid %in% complete_ids) %>%
  group_by(huc8) %>%
  summarise(n_reaches = n_distinct(comid),
            gpp_gg_o2_y = sum(gpp * river_area_m2 * days) / 1e9,
            er_gg_o2_y = sum(er * river_area_m2 * days) / 1e9,
            net_gg_o2_y = gpp_gg_o2_y - er_gg_o2_y, .groups = "drop")
stopifnot(all(is.finite(huc_totals$net_gg_o2_y)),
          sum(huc_totals$n_reaches) == length(complete_ids))
message("HUC8s with GPP > ER: ", sum(huc_totals$net_gg_o2_y > 0),
        " / ", nrow(huc_totals))
print(huc_totals %>% summarise(min_net = min(net_gg_o2_y),
                                median_net = median(net_gg_o2_y),
                                max_net = max(net_gg_o2_y)))

boundaries <- sf::read_sf("data/processed/04c_huc08_boundaries.gpkg") %>%
  transmute(huc8 = as.character(huc8)) %>% st_transform(5070)
stopifnot(!anyDuplicated(boundaries$huc8),
          all(huc_totals$huc8 %in% boundaries$huc8))
map_data <- boundaries %>% left_join(huc_totals, by = "huc8")

conus <- rnaturalearth::ne_states(
  country = "United States of America", returnclass = "sf") %>%
  filter(!postal %in% c("AK", "HI", "PR")) %>%
  st_transform(5070) %>%
  summarise(geometry = st_union(geometry))
extent <- st_bbox(conus)

# A symmetric pseudo-log color scale retains sign and displays both modest
# and extreme HUC8 totals without saturating nearly every polygon.
limit <- max(abs(huc_totals$net_gg_o2_y))
ticks <- c(-100, -10, 0, 10, 100)
ticks <- ticks[abs(ticks) <= limit]
figure <- ggplot() +
  geom_sf(data = conus, fill = "white", color = "grey65", linewidth = 0.2) +
  geom_sf(data = map_data, aes(fill = net_gg_o2_y),
          color = "white", linewidth = 0.025) +
  geom_sf(data = conus, fill = NA, color = "#253746", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, trans = scales::pseudo_log_trans(sigma = 1),
                       limits = c(-limit, limit), breaks = ticks,
                       labels = scales::label_number(),
                       na.value = "grey85", name = "GPP − ER\n(Gg O₂ yr⁻¹)") +
  coord_sf(crs = st_crs(5070), xlim = extent[c("xmin", "xmax")],
           ylim = extent[c("ymin", "ymax")], expand = FALSE, datum = NA) +
  labs(caption = stringr::str_wrap(paste0(
    "CONUS HUC8s; Maavara no-uparea monthly predictions summed to annual local reach totals (n = ",
    format(length(complete_ids), big.mark = ","), " complete GPP/ER reaches). ",
    "Blue: ER exceeds GPP; red: GPP exceeds ER. Symmetric pseudo-log color scale; gray: no valid reaches."), 135)) +
  theme_void(base_size = 12) +
  theme(legend.position = "right",
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(8, 8, 8, 8)) +
  guides(fill = guide_colorbar(barheight = grid::unit(5.5, "cm")))
ggsave("figures/40_huc8_annual_net_metabolism_map.png", figure,
       width = 13, height = 8, dpi = 300, bg = "white")
