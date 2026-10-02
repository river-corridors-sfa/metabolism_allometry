# One-off: local metabolism rate vs local drainage area ------------------------

source("scripts/00_setup.R")

merit_cat <- read_sf(
  "data/processed/merit_data/03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
) %>%
  clean_names() %>%
  st_drop_geometry() %>%
  transmute(
    comid = as.character(comid),
    drainage_area_km2 = unitarea
  )

stream_area <- read_csv(
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",
  show_col_types = FALSE
) %>%
  clean_names() %>%
  transmute(
    comid = as.character(comid),
    stream_area_m2 = lengthkm * 1000 * merit_width)

metabolism <- readRDS("data/processed/01_assembled_maavara_monthly_metabolism.rds") %>%
  clean_names() %>%
  mutate(comid = as.character(comid)) %>%
  select(comid, month, gpp_g_o2_m2_d, er_g_o2_m2_d) %>%
  pivot_longer(
    c(gpp_g_o2_m2_d, er_g_o2_m2_d),
    names_to = "metabolism",
    values_to = "rate_g_o2_m2_d") %>%
  inner_join(stream_area, by = "comid") %>%
  inner_join(merit_cat, by = "comid") %>%
  mutate(rate_g_o2_d = rate_g_o2_m2_d * stream_area_m2) %>% 
  filter(
    rate_g_o2_d > 0,
    drainage_area_km2 > 0) 
  

metabolism %>% 
  slice(1:500000) %>% 
  ggplot(aes(log10(drainage_area_km2), 
             log10(rate_g_o2_m2_d),
             color = metabolism)) + 
  geom_point(alpha = 0.1, show.legend = F) + 
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed") +
  facet_wrap(~metabolism, nrow = 1, scales = "free") + 
  theme_bw()



# Add upstream drainage area
upstream_area <- read_sf(
  "data/processed/merit_data/03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
) %>%
  clean_names() %>%
  st_drop_geometry() %>%
  transmute(
    comid = as.character(comid),
    upstream_area_km2 = uparea
  )

metabolism_network <- metabolism %>%
  select(-drainage_area_km2) %>%
  inner_join(upstream_area, by = "comid") %>%
  # inner_join(YOUR_HUC8_LOOKUP, by = "comid") %>%
  filter(
    rate_g_o2_m2_d > 0,
    upstream_area_km2 > 0
  )

# Then, within each HUC8 x month x metabolism type:
network_slopes <- metabolism_network %>%
  group_by(huc8, month, metabolism) %>%
  filter(n() >= 20) %>%
  summarise(
    slope = coef(
      lm(log10(rate_g_o2_m2_d) ~ log10(upstream_area_km2))
    )[2],
    r2 = summary(
      lm(log10(rate_g_o2_m2_d) ~ log10(upstream_area_km2))
    )$r.squared,
    n = n(),
    .groups = "drop"
  )

# Are within-network slopes generally positive?
ggplot(network_slopes, aes(slope)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~metabolism, nrow = 1) +
  labs(
    x = "Within-HUC8 slope: per-area metabolism vs upstream drainage area",
    y = "HUC8-month count"
  ) +
  theme_bw()

