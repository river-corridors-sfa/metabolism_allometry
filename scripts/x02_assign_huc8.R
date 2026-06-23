# 02_assign_huc8.R
# Assign MERIT-Hydro COMIDs to HUC8s using representative catchment points

rm(list = ls())
source("scripts/00_setup.R")


# 1. Read and prep datasets --------------------------------------------------

pacman::p_load(nhdplusTools,
               tictoc)

test_huc8 <- "10100004"

merit <- st_read(
  file.path(dir_dryad,
            "MERIT_Hydro_CONUS_catchments",
            "MERIT_Hydro_CONUS_catchments.shp"),
  quiet = TRUE
) |>
  clean_names() |>
  select(comid, unitarea, lengthkm)

merit_riv <- st_read(
  "data/merit_data/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp",
  quiet = TRUE
) |>
  clean_names() |>
  select(comid, next_down_id, order, uparea, lengthkm)

merit_aoi <- merit |>
  st_transform(4326) |>
  st_bbox() |>
  st_as_sfc() |>
  st_sf(geometry = _)

tic("pull huc8")

huc8 <- nhdplusTools::get_huc(
  AOI = merit_aoi,
  type = "huc08",
  t_srs = common_crs
) |>
  clean_names() |>
  select(huc8, name) |>
  st_make_valid()

toc()

merit_pts <- merit |>
  st_transform(common_crs) |>
  st_make_valid() |>
  st_point_on_surface()


# 2.1 Test case before full assignment ----------------------------------------

test_boundary <- huc8 |>
  filter(huc8 == test_huc8)

test_pts <- merit_pts |>
  st_join(test_boundary |> select(huc8), join = st_within) |>
  filter(huc8 == test_huc8)

test_riv <- merit_riv |>
  st_transform(common_crs) |>
  st_filter(test_boundary)

test_summary <- tibble(
  test_huc8 = test_huc8,
  n_test_points = nrow(test_pts),
  n_test_flowlines = nrow(test_riv),
  huc8_area_km2 = as.numeric(st_area(test_boundary)) / 1e6
)

test_summary

ggplot() +
  geom_sf(data = test_boundary, fill = "grey95", color = "black", linewidth = 0.8) +
  geom_sf(data = test_riv, aes(linewidth = order), alpha = 0.7) +
  geom_sf(data = test_pts, size = 0.5, alpha = 0.8) +
  scale_linewidth_continuous(range = c(0.15, 1.2)) +
  labs(
    title = paste("HUC8 assignment diagnostic:", test_huc8),
    linewidth = "Order"
  ) +
  theme_void()


# 2.2 additional checks ------------------------------------------------------

# Test-case diagnostics ------------------------------------------------------

test_summary

test_pts |>
  st_drop_geometry() |>
  summarize(
    n_points = n(),
    n_comid = n_distinct(comid),
    n_duplicate_comid = n() - n_distinct(comid)
  )

test_riv |>
  st_drop_geometry() |>
  summarize(
    n_flowlines = n(),
    n_flowline_comid = n_distinct(comid),
    pct_flowline_comids_in_points = mean(comid %in% test_pts$comid) * 100,
    pct_point_comids_in_flowlines = mean(test_pts$comid %in% comid) * 100
  )


# 3. Run full assignment -----------------------------------------------------

tic("assign huc8")

merit_huc8 <- merit_pts |>
  st_join(huc8 |> select(huc8), join = st_within) |>
  st_drop_geometry() |>
  select(comid, huc8)

toc()


# 4. Diagnostics before writing ---------------------------------------------

merit_huc8 |>
  summarize(
    n_rows = n(),
    n_comid = n_distinct(comid),
    n_missing_huc8 = sum(is.na(huc8)),
    pct_missing_huc8 = mean(is.na(huc8)) * 100
  )

merit_huc8 |>
  count(comid, name = "n_huc8") |>
  count(n_huc8, name = "n_comid")

merit_huc8 |>
  count(huc8, sort = TRUE) |>
  print(n = 25)

merit_huc8 |>
  filter(is.na(huc8)) |>
  slice_head(n = 10)

test_assigned <- merit_pts |>
  inner_join(
    merit_huc8 |> filter(huc8 == test_huc8),
    by = "comid"
  )

ggplot() +
  geom_sf(data = test_boundary, fill = "grey95", color = "black", linewidth = 0.8) +
  geom_sf(data = test_riv, aes(linewidth = order), alpha = 0.7) +
  geom_sf(data = test_assigned, size = 0.5, alpha = 0.8) +
  scale_linewidth_continuous(range = c(0.15, 1.2)) +
  labs(
    title = paste("Post-assignment check:", test_huc8),
    linewidth = "Order"
  ) +
  theme_void()


# 5. Write outputs -----------------------------------------------------------
# Uncomment only after diagnostics look correct.

# saveRDS(huc8, "data/processed/huc8_merit_bbox.rds", compress = "xz")
# saveRDS(merit_huc8, "data/processed/merit_huc8_crosswalk.rds", compress = "xz")
