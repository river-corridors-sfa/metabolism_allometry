# 05b_screen_centroid_mismatch_reaches.R

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(
  furrr,
  future,
  ggplot2,
  nhdplusTools,
  scales,
  sf
)


# 0. Settings ------------------------------------------------------------------

centroid_threshold_km <- 1

n_workers <- max(1, parallel::detectCores() - 2)

path_flowlines <- "data/processed/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp"
path_metabolism <- "data/processed/03_maavara_monthly_metabolism.rds"
path_centroid_check <- "data/processed/centroid_difference_MERIT_Hydro_CONUS_catchments_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01.csv"

dir_05b <- "data/processed/05b_centroid_mismatch"
dir.create(dir_05b, recursive = TRUE, showWarnings = FALSE)


# 1. Read inputs ----------------------------------------------------------------

flowlines_sf <- read_sf(path_flowlines) |>
  clean_names() |>
  mutate(
    comid = as.character(comid),
    nxt_dw = as.character(nxt_dw)
  )

metabolism <- readRDS(path_metabolism) |>
  mutate(comid = as.character(comid))

centroid_check <- read_csv(path_centroid_check, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    comid = as.character(comid),
    dist_km = centroid_distance_m / 1000
  )


# 2. Create HUC8 lookup ---------------------------------------------------------

huc_level <- "huc08"
huc_col <- str_replace(huc_level, "huc0", "huc")

huc_boundaries <- nhdplusTools::get_huc(
  AOI = st_as_sfc(st_bbox(flowlines_sf)),
  type = huc_level,
  t_srs = st_crs(flowlines_sf)
) |>
  clean_names() |>
  select(all_of(huc_col))

flowline_hucs <- flowlines_sf |>
  st_point_on_surface() |>
  st_join(huc_boundaries, join = st_within) |>
  st_drop_geometry() |>
  mutate(comid = as.character(comid)) |>
  select(comid, huc = all_of(huc_col), uparea)


# 3. Identify removed COMIDs ----------------------------------------------------

removed_comids <- centroid_check |>
  filter(dist_km > centroid_threshold_km) |>
  distinct(comid, dist_km) |>
  arrange(desc(dist_km))

write_rds(
  removed_comids,
  file.path(dir_05b, "05b_removed_comids.rds")
)

write_csv(
  removed_comids,
  file.path(dir_05b, "05b_removed_comids.csv")
)


# 4. Save filtered metabolism for 06 -------------------------------------------

metabolism_filtered <- metabolism |>
  anti_join(removed_comids, by = "comid")

write_rds(
  metabolism_filtered,
  file.path(dir_05b, "05b_maavara_monthly_metabolism_centroid_filtered.rds")
)


# 5. Save removed spatial features ---------------------------------------------

removed_flowlines_sf <- flowlines_sf |>
  inner_join(removed_comids, by = "comid") |>
  mutate(
    dist_class = case_when(
      dist_km <= 5 ~ "1–5 km",
      dist_km <= 25 ~ "5–25 km",
      TRUE ~ ">25 km"
    )
  )

write_sf(
  removed_flowlines_sf,
  file.path(dir_05b, "05b_removed_centroid_mismatch_reaches.gpkg"),
  delete_dsn = TRUE
)


# 6. Map removed reaches --------------------------------------------------------

p_removed_map <- ggplot() +
  geom_sf(
    data = st_geometry(flowlines_sf),
    color = "grey88",
    linewidth = 0.05
  ) +
  geom_sf(
    data = removed_flowlines_sf,
    aes(color = dist_class),
    linewidth = 0.25
  ) +
  coord_sf(datum = NA) +
  scale_color_viridis_d(option = "C", end = 0.9) +
  labs(
    title = "MERIT reaches removed due to centroid mismatch",
    color = "Centroid offset"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(dir_05b, "05b_removed_centroid_mismatch_map.png"),
  p_removed_map,
  width = 10,
  height = 6,
  dpi = 300
)


# 7. HUC8 reach-loss diagnostics ------------------------------------------------

flowline_lengths <- flowlines_sf |>
  mutate(length_km = as.numeric(st_length(geometry)) / 1000) |>
  st_drop_geometry() |>
  select(comid, length_km)

huc08_reach_impact <- flowline_hucs |>
  filter(!is.na(huc)) |>
  left_join(flowline_lengths, by = "comid") |>
  distinct(huc, comid, length_km) |>
  mutate(removed = comid %in% removed_comids$comid) |>
  group_by(huc) |>
  summarise(
    n_reaches_total = n(),
    n_reaches_removed = sum(removed),
    pct_reaches_removed = n_reaches_removed / n_reaches_total * 100,
    length_total_km = sum(length_km, na.rm = TRUE),
    length_removed_km = sum(length_km[removed], na.rm = TRUE),
    pct_length_removed = length_removed_km / length_total_km * 100,
    .groups = "drop"
  )


# 8. HUC8 monthly metabolic leverage -------------------------------------------

plan(multisession, workers = n_workers)

metabolism_long <- metabolism |>
  inner_join(flowline_hucs |> select(comid, huc), by = "comid") |>
  select(huc, comid, month, gpp, er) |>
  pivot_longer(
    cols = c(gpp, er),
    names_to = "flux",
    values_to = "value"
  ) |>
  mutate(
    removed = comid %in% removed_comids$comid
  )

huc08_flux_impact <- metabolism_long |>
  group_split(huc, month, flux) |>
  future_map_dfr(\(x) {
    x |>
      summarise(
        huc = first(huc),
        month = first(month),
        flux = first(flux),
        n_reaches_total = n_distinct(comid),
        n_reaches_removed = n_distinct(comid[removed]),
        flux_total = sum(value, na.rm = TRUE),
        flux_removed = sum(value[removed], na.rm = TRUE),
        pct_flux_removed = flux_removed / flux_total * 100,
        removed_flux_rank_pctl = mean(percent_rank(value)[removed], na.rm = TRUE),
        .groups = "drop"
      )
  })

plan(sequential)


# 9. Combined HUC8 impact screen ------------------------------------------------

huc08_flux_summary <- huc08_flux_impact |>
  group_by(huc, flux) |>
  summarise(
    max_pct_flux_removed = max(pct_flux_removed, na.rm = TRUE),
    median_pct_flux_removed = median(pct_flux_removed, na.rm = TRUE),
    max_removed_flux_rank_pctl = max(removed_flux_rank_pctl, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = flux,
    values_from = c(
      max_pct_flux_removed,
      median_pct_flux_removed,
      max_removed_flux_rank_pctl
    )
  )

huc08_centroid_filter_impact <- huc08_reach_impact |>
  left_join(huc08_flux_summary, by = "huc") |>
  mutate(
    affected_by_count = pct_reaches_removed > 5,
    affected_by_length = pct_length_removed > 5,
    affected_by_gpp = max_pct_flux_removed_gpp > 5,
    affected_by_er = max_pct_flux_removed_er > 5,
    likely_allometry_sensitive = if_any(
      c(affected_by_count, affected_by_length, affected_by_gpp, affected_by_er),
      ~ replace_na(.x, FALSE)
    )
  ) |>
  arrange(desc(likely_allometry_sensitive), desc(pct_reaches_removed))

write_rds(
  huc08_centroid_filter_impact,
  file.path(dir_05b, "05b_huc08_centroid_filter_impact.rds")
)

write_csv(
  huc08_centroid_filter_impact,
  file.path(dir_05b, "05b_huc08_centroid_filter_impact.csv")
)

write_rds(
  huc08_flux_impact,
  file.path(dir_05b, "05b_huc08_monthly_flux_impact.rds")
)

write_csv(
  huc08_flux_impact,
  file.path(dir_05b, "05b_huc08_monthly_flux_impact.csv")
)


# 10. Quick QA figure -----------------------------------------------------------

p_huc08_impact <- huc08_centroid_filter_impact |>
  filter(n_reaches_removed > 0) |>
  ggplot(aes(x = pct_reaches_removed, y = max_pct_flux_removed_gpp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = n_reaches_removed, color = likely_allometry_sensitive), alpha = 0.75) +
  scale_color_viridis_d(option = "C", end = 0.85) +
  scale_size_continuous(range = c(1.5, 7)) +
  labs(
    title = "Potential HUC8 allometry sensitivity to centroid filtering",
    x = "Removed reaches (%)",
    y = "Maximum monthly GPP removed (%)",
    color = "Sensitive?",
    size = "Removed reaches"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(
  file.path(dir_05b, "05b_huc08_centroid_filter_impact.png"),
  p_huc08_impact,
  width = 7,
  height = 5,
  dpi = 300
)