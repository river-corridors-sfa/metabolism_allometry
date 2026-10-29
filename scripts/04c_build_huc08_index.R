## 04c_build_huc08_index.R
##
## PURPOSE
## Build small, durable HUC8 resources for downstream filtering and mapping.
## The spatial assignment already exists on the filtered MERIT data (`huc`),
## so this script does not repeat a COMID-to-HUC spatial join.


# 0. Setup ---------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(nhdplusTools)

tic("run 04c")

refresh_huc08_boundaries <- FALSE
boundary_simplification_m <- 500

dir_processed <- "data/processed"
dir_merit <- file.path(dir_processed, "merit_data")

path_flowlines <- file.path(
  dir_merit,
  "03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
)
path_catchments <- file.path(
  dir_merit,
  "03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
)
path_metabolism <- file.path(
  dir_processed,
  "03_filtered_no_uparea_monthly_metabolism.rds"
)
path_ecoregions <- file.path(
  dir_processed,
  "04b_huc08_l1_ecoregion_lookup.csv"
)
path_spatial_coverage <- file.path(
  dir_processed,
  "03_no_uparea_huc_spatial_coverage.rds"
)

path_boundaries <- file.path(dir_processed, "04c_huc08_boundaries.gpkg")
path_comid_key <- file.path(dir_processed, "04c_no_uparea_comid_huc08_key.rds")
path_summary <- file.path(dir_processed, "04c_no_uparea_huc08_summary.csv")


# 1. Read the retained network -------------------------------------------------

flowlines <- read_sf(path_flowlines) %>%
  clean_names() %>%
  mutate(
    comid = as.character(comid),
    huc8 = str_pad(as.character(huc), 8, pad = "0")
  )

catchments <- read_sf(path_catchments) %>%
  clean_names() %>%
  mutate(
    comid = as.character(comid),
    huc8 = str_pad(as.character(huc), 8, pad = "0")
  )

retained_huc8 <- flowlines %>%
  st_drop_geometry() %>%
  filter(!is.na(huc8)) %>%
  distinct(huc8)


# 2. Cache HUC8 polygons -------------------------------------------------------
##
## The USGS request is only needed for the first build (or an explicit refresh).
## All downstream scripts should read the saved GeoPackage.

if (!file.exists(path_boundaries) || refresh_huc08_boundaries) {
  conus_aoi <- st_as_sfc(st_bbox(flowlines)) %>%
    st_as_sf()

  huc08_boundaries <- get_huc(
    AOI = conus_aoi,
    type = "huc08",
    t_srs = st_crs(flowlines),
    buffer = 0
  ) %>%
    clean_names() %>%
    transmute(
      huc8 = str_pad(as.character(huc8), 8, pad = "0"),
      huc8_name = name,
      states = states,
      wbd_area_km2 = as.numeric(areasqkm),
      geometry
    ) %>%
    semi_join(retained_huc8, by = "huc8") %>%
    arrange(huc8)

  # The source WBD geometry is unnecessarily detailed for CONUS and HUC8 maps.
  # Simplification keeps this reusable cache small; authoritative area remains
  # available in `wbd_area_km2` and is not recalculated from simplified shapes.
  huc08_boundaries <- huc08_boundaries %>%
    st_transform(5070) %>%
    st_simplify(
      dTolerance = boundary_simplification_m,
      preserveTopology = TRUE
    ) %>%
    st_transform(st_crs(flowlines))

  st_write(
    huc08_boundaries,
    path_boundaries,
    delete_dsn = TRUE,
    quiet = TRUE
  )
} else {
  huc08_boundaries <- read_sf(path_boundaries)
}


# 3. Summarize reach-level metabolism availability ----------------------------

metabolism_coverage <- readRDS(path_metabolism) %>%
  clean_names() %>%
  mutate(comid = as.character(comid)) %>%
  group_by(comid) %>%
  summarise(
    n_months_present = n_distinct(month),
    n_good_gpp_months = sum(
      is.finite(gpp_g_o2_m2_d) & gpp_g_o2_m2_d > 0
    ),
    n_good_er_months = sum(
      is.finite(er_g_o2_m2_d) & er_g_o2_m2_d > 0
    ),
    n_good_both_months = sum(
      is.finite(gpp_g_o2_m2_d) & gpp_g_o2_m2_d > 0 &
        is.finite(er_g_o2_m2_d) & er_g_o2_m2_d > 0
    ),
    .groups = "drop"
  )


# 4. Build the compact COMID key ----------------------------------------------
##
## `is_headwater_retained` answers the analysis-specific question: does this
## reach have another reach upstream within the retained, filtered network?
## `is_headwater_merit` preserves the source-network indicator for comparison.

# Same retained-edge rule as 05, computed from 03 flowlines so 04c runs before 05.
edges_downstream <- flowlines %>%
  st_drop_geometry() %>%
  transmute(from_comid = as.character(comid), to_comid = as.character(nxt_dw)) %>%
  filter(!is.na(to_comid), to_comid != "0", to_comid %in% flowlines$comid)

ecoregions <- read_csv(path_ecoregions, show_col_types = FALSE) %>%
  transmute(
    huc8 = str_pad(as.character(huc), 8, pad = "0"),
    l1_code,
    l1_ecoregion
  )

huc08_attributes <- huc08_boundaries %>%
  st_drop_geometry() %>%
  select(huc8, huc8_name, states)

comid_huc08_key <- flowlines %>%
  st_drop_geometry() %>%
  transmute(
    comid,
    huc8,
    uparea_km2 = uparea,
    stream_order = order,
    is_headwater_merit = maxup == 0,
    is_headwater_retained = !comid %in% edges_downstream$to_comid
  ) %>%
  left_join(
    catchments %>%
      st_drop_geometry() %>%
      transmute(comid, unitarea_km2 = unitarea),
    by = "comid"
  ) %>%
  left_join(metabolism_coverage, by = "comid") %>%
  left_join(huc08_attributes, by = "huc8") %>%
  left_join(ecoregions, by = "huc8") %>%
  mutate(
    across(
      c(n_months_present, n_good_gpp_months, n_good_er_months,
        n_good_both_months),
      ~ replace_na(.x, 0L)
    ),
    headwater_agreement = is_headwater_merit == is_headwater_retained
  ) %>%
  arrange(huc8, comid)


# 5. Build a HUC8 screening summary ------------------------------------------

spatial_coverage <- readRDS(path_spatial_coverage) %>%
  transmute(
    huc8 = str_pad(as.character(huc), 8, pad = "0"),
    huc8_area_km2 = huc_area_km2,
    retained_catchment_area_km2 = catchment_area_km2,
    pct_area_covered = pct_huc_covered_by_catchments,
    n_catchments
  )

huc08_centroids <- huc08_boundaries %>%
  select(huc8) %>%
  st_transform(5070) %>%
  st_point_on_surface() %>%
  st_transform(4326) %>%
  mutate(
    centroid_longitude = st_coordinates(.)[, 1],
    centroid_latitude = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(huc8, centroid_longitude, centroid_latitude)

huc08_summary <- comid_huc08_key %>%
  filter(!is.na(huc8)) %>%
  group_by(huc8) %>%
  summarise(
    n_flowlines = n(),
    n_flowlines_any_metabolism = sum(n_good_both_months > 0),
    n_complete_flowlines = sum(n_good_both_months == 12),
    n_headwaters = sum(is_headwater_retained),
    has_headwater = any(is_headwater_retained),
    n_headwater_disagreements = sum(!headwater_agreement, na.rm = TRUE),
    min_uparea_km2 = min(uparea_km2, na.rm = TRUE),
    max_uparea_km2 = max(uparea_km2, na.rm = TRUE),
    log10_uparea_range = log10(max_uparea_km2 / min_uparea_km2),
    .groups = "drop"
  ) %>%
  left_join(huc08_attributes, by = "huc8") %>%
  left_join(ecoregions, by = "huc8") %>%
  left_join(spatial_coverage, by = "huc8") %>%
  left_join(huc08_centroids, by = "huc8") %>%
  select(
    huc8,
    huc8_name,
    states,
    l1_code,
    l1_ecoregion,
    centroid_longitude,
    centroid_latitude,
    huc8_area_km2,
    retained_catchment_area_km2,
    pct_area_covered,
    n_catchments,
    n_flowlines,
    n_flowlines_any_metabolism,
    n_complete_flowlines,
    n_headwaters,
    has_headwater,
    n_headwater_disagreements,
    min_uparea_km2,
    max_uparea_km2,
    log10_uparea_range
  ) %>%
  arrange(huc8)


# 6. Checks and outputs --------------------------------------------------------

stopifnot(
  nrow(comid_huc08_key) == nrow(flowlines),
  n_distinct(comid_huc08_key$comid) == nrow(comid_huc08_key),
  n_distinct(huc08_summary$huc8) == nrow(huc08_summary)
)

saveRDS(comid_huc08_key, path_comid_key)
write_csv(huc08_summary, path_summary)

tibble(
  n_flowlines = nrow(comid_huc08_key),
  n_huc8 = nrow(huc08_summary),
  n_missing_huc8 = sum(is.na(comid_huc08_key$huc8)),
  n_missing_ecoregion = sum(is.na(huc08_summary$l1_ecoregion)),
  n_huc8_with_20_complete_reaches = sum(
    huc08_summary$n_complete_flowlines >= 20
  ),
  n_huc8_with_headwater = sum(huc08_summary$has_headwater)
) %>%
  print()

message("HUC8 boundaries saved: ", path_boundaries)
message("COMID-HUC8 key saved: ", path_comid_key)
message("HUC8 summary saved: ", path_summary)

toc()
