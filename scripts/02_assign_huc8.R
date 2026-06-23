# 02_assign_huc8.R
# Pull HUC8 boundaries and assign MERIT-Hydro catchments

rm(list = ls())
source("scripts/00_setup.R")


# 1. Packages ---------------------------------------------------------------

pacman::p_load(nhdplusTools,
               tictoc)


# 2. Paths ------------------------------------------------------------------

path_merit <- file.path(
  dir_dryad,
  "MERIT_Hydro_CONUS_catchments",
  "MERIT_Hydro_CONUS_catchments.shp"
)

dir_processed <- file.path("data/processed")
dir.create(dir_processed, recursive = TRUE, showWarnings = FALSE)

path_huc8 <- file.path(dir_processed, "huc8_merit_bbox.rds")
path_merit_huc8 <- file.path(dir_processed, "merit_huc8_crosswalk.rds")


# 3. Read MERIT catchments ---------------------------------------------------

tic("read merit") #3.5s

merit <- st_read(path_merit, quiet = TRUE) |>
  clean_names() %>% 
  st_transform(common_crs)

toc()


# 4. Pull HUC8 boundaries ----------------------------------------------------

tic("pull huc8")

huc8 <- nhdplusTools::get_huc(
  AOI = st_as_sfc(st_bbox(merit)),
  type = "huc08") |>
  clean_names() |>
  select(huc8, name)

toc()

# 5. Test case --------------------------------------------------------------

i = 1

huc = huc8 %>% slice(1)
  
ggplot() + geom_sf(data = huc)



# 5. Assign MERIT catchments to HUC8 ----------------------------------------

tic("assign huc8")

merit_huc8 <- merit |>
  st_point_on_surface() |>
  st_join(huc8 |> select(huc8), join = st_within) |>
  st_drop_geometry() |>
  transmute(comid, huc8)

toc()


# 6. Diagnostics ------------------------------------------------------------

merit_huc8 |>
  count(is.na(huc8), name = "n")

merit_huc8 |>
  count(huc8, sort = TRUE)


# 7. Save outputs ------------------------------------------------------------

saveRDS(huc8, path_huc8, compress = "xz")
saveRDS(merit_huc8, path_merit_huc8, compress = "xz")
