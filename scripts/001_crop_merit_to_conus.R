## This script crops the downloaded MERIT datasets down to the spatial extent of 
## CONUS. We do this step first since this dataset is one of the major bottlenecks
## due to its size and the number of different computations we need to do on it
## We do need to compute for both the catchments ("cat") and flowine ("riv") files, 
## the catchment files are used to verify we're matching Maavara, and ID any 
## problematic catchments. The flowline file is why we actually need this dataset


# 0. Setup ---------------------------------------------------------------------

rm(list = ls())
source("scripts/00_setup.R")

pacman::p_load(
  tidyverse,
  ggthemes,
  sf,
  janitor,
  tictoc)


# 0. Settings ------------------------------------------------------------------

dir_merit_raw <- "data/raw/merit_data"
dir_merit_processed <- "data/processed/merit_data"

path_cat_raw <- file.path(dir_merit_raw, "cat_pfaf_7_MERIT_Hydro_v07_Basins_v01.shp")
path_riv_raw <- file.path(dir_merit_raw, "riv_pfaf_7_MERIT_Hydro_v07_Basins_v01.shp")

path_cat_conus <- file.path(dir_merit_processed, "cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp")
path_riv_conus <- file.path(dir_merit_processed, "riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp")
path_log <- file.path(dir_merit_processed, "001_crop_merit_to_conus_log.csv")

dir.create(dir_merit_processed, recursive = TRUE, showWarnings = FALSE)


# 1. Define CONUS boundary -----------------------------------------------------

# Decision: Use built-in Natural Earth country/state geometries to define CONUS,
# excluding AK, HI, and territories.
conus <- rnaturalearth::ne_states(
  country = "United States of America",
  returnclass = "sf"
) |>
  clean_names() |>
  filter(!postal %in% c("AK", "HI", "PR")) |>
  summarise(geometry = st_union(geometry)) |>
  st_transform(common_crs)

ggplot() + 
  geom_sf(data = conus)


# 2. Read MERIT inputs ---------------------------------------------------------

tic("Read MERIT inputs") # 33s

cat_raw <- read_sf(path_cat_raw) |>
  clean_names() |>
  st_transform(common_crs)

riv_raw <- read_sf(path_riv_raw) |>
  clean_names() |>
  st_transform(common_crs)

toc()


# 3. Crop to CONUS -------------------------------------------------------------

# Decision: Crop catchments using an internal representative point rather than
# polygon intersection to avoid expensive validity repair on all MERIT polygons.
# Flowlines are cropped by intersection because line geometries do not trigger the
# same invalid-polygon overlay issue.

tic("Crop MERIT to CONUS") # 39s!! So much faster than st_make_valid()!

cat_conus_assigned <- cat_raw %>%
  mutate(in_conus = lengths(st_intersects(st_point_on_surface(geometry), conus)) > 0) 

riv_conus_assigned <- riv_raw %>%
  mutate(in_conus = lengths(st_intersects(geometry, conus)) > 0)

toc()


# 4. Visual checks -------------------------------------------------------------

## Calculate rows for 
cats_kept <- nrow(cat_conus_assigned %>% filter(in_conus == TRUE))
cats_total <- nrow(cat_raw)
rivs_kept <- nrow(riv_conus_assigned %>% filter(in_conus == TRUE))
rivs_total <- nrow(riv_raw)


tic("save plots")
p_cat_crop <- ggplot() +
  geom_sf(data = conus, color = "black", fill = "white", linewidth = 0.3) +
  geom_sf(data = cat_conus_assigned, linewidth = 0.1, color = "blue") +
  facet_wrap(~in_conus, ncol = 1) +
  ggtitle(paste0(cats_kept, "catchments of ", cats_total, " kept")) + 
  theme_map() +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(dir_merit_processed, "001_cat_conus_crop_check.png"),
  p_cat_crop, width = 7, height = 9, dpi = 300)

p_riv_crop <- ggplot() +
  geom_sf(data = conus, color = "black", fill = "white", linewidth = 0.3) +
  geom_sf(data = riv_conus_assigned, linewidth = 0.1, color = "blue") +
  facet_wrap(~in_conus, ncol = 1) +
  ggtitle(paste0(rivs_kept, "flowlines of ", rivs_total, " kept")) + 
  theme_void() +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(dir_merit_processed, "001_riv_conus_crop_check.png"),
       p_riv_crop, width = 7, height = 9, dpi = 300)
toc()


# 5. Remove non-CONUS parts ----------------------------------------------------

tic("Filter out non-CONUS") 
cat_conus <- cat_conus_assigned %>% 
  filter(in_conus == TRUE)

riv_conus <- riv_conus_assigned %>% 
  filter(in_conus == TRUE)
toc()


# 6. Basic provenance checks ---------------------------------------------------

log_tbl <- tibble(
  dataset = c("catchments", "flowlines"),
  input_path = c(path_cat_raw, path_riv_raw),
  output_path = c(path_cat_conus, path_riv_conus),
  n_input_rows = c(nrow(cat_raw), nrow(riv_raw)),
  n_output_rows = c(nrow(cat_conus), nrow(riv_conus)),
  n_input_comids = c(n_distinct(cat_raw$comid), n_distinct(riv_raw$comid)),
  n_output_comids = c(n_distinct(cat_conus$comid), n_distinct(riv_conus$comid)),
  n_output_missing_comid = c(sum(is.na(cat_conus$comid)), sum(is.na(riv_conus$comid))),
  crs_epsg = common_crs,
  conus_rule = "Natural Earth US states; exclude AK, HI, PR; retain features intersecting CONUS"
)

write_csv(log_tbl, path_log)

log_tbl


# 7. Write outputs -------------------------------------------------------------

tic("Write CONUS MERIT shapefiles")

write_sf(cat_conus %>% filter(), path_cat_conus, delete_dsn = TRUE)
write_sf(riv_conus, path_riv_conus, delete_dsn = TRUE)

toc()
