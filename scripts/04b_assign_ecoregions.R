## 04b_assign_ecoregions.R
##
## PURPOSE
## Assign each retained HUC8 to a North American Level I terrestrial ecoregion.
##
## We retrieve Level I ecoregions directly from the ArcGIS REST FeatureServer,
## assign HUC8s using an internal representative point, and save a lightweight
## HUC8-to-ecoregion lookup for downstream analyses.


# 0. Setup ---------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")

tic("run 04b")

dir_processed <- "data/processed"
dir_merit <- file.path(dir_processed, "merit_data")

path_flowlines <- file.path(
  dir_merit,
  "03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
)

path_ecoregion_lookup <- file.path(
  dir_processed,
  "04b_huc08_l1_ecoregion_lookup.csv"
)


# ArcGIS REST FeatureServer query for North American Level I ecoregions
ecoregion_l1_url <- paste0(
  "https://services7.arcgis.com/oF9CDB4lUYF7Um9q/arcgis/rest/services/",
  "NA_Terrestrial_Ecoregions_Level_1/FeatureServer/5/query?",
  "where=1%3D1",
  "&outFields=LEVEL1%2CNameL1_En",
  "&returnGeometry=true",
  "&outSR=4326",
  "&f=geojson"
)


# 1. Read retained MERIT catchments --------------------------------------------

path_cat <- file.path(
  dir_merit,
  "03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
)

merit_cat <- read_sf(path_cat) %>%
  clean_names() %>%
  mutate(
    huc = str_pad(as.character(huc), 8, pad = "0")
  )


# 2. Create one representative point per HUC8 ----------------------------------
#
# Use catchment interior points rather than downloading HUC polygons.
# The mean location of retained catchment points provides one representative
# location for assigning each HUC8 to a Level I ecoregion.

huc_points <- merit_cat %>%
  st_point_on_surface() %>%
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  group_by(huc) %>%
  summarise(
    x = mean(x),
    y = mean(y),
    .groups = "drop"
  ) %>%
  st_as_sf(
    coords = c("x", "y"),
    crs = st_crs(merit_cat)
  )


# 3. Get Level I ecoregions -----------------------------------------------------

ecoregions_l1 <- read_sf(ecoregion_l1_url) %>%
  clean_names() %>%
  transmute(
    l1_code = level1,
    l1_ecoregion = name_l1_en
  ) %>%
  st_transform(
    st_crs(huc_points)
  )


# 4. Assign HUC8s to ecoregions -------------------------------------------------

huc_ecoregions <- huc_points %>%
  st_join(
    ecoregions_l1,
    join = st_within
  ) %>%
  st_drop_geometry() %>%
  select(
    huc,
    l1_code,
    l1_ecoregion
  ) %>%
  distinct() %>%
  arrange(huc)


# 5. Checks --------------------------------------------------------------------

stopifnot(
  nrow(huc_ecoregions) == n_distinct(merit_cat$huc),
  n_distinct(huc_ecoregions$huc) == nrow(huc_ecoregions)
)

huc_ecoregions %>%
  summarise(
    n_hucs = n(),
    n_ecoregions = n_distinct(l1_ecoregion),
    n_missing_ecoregion = sum(is.na(l1_ecoregion))
  ) %>%
  print()


# 6. Write output ---------------------------------------------------------------

write_csv(
  huc_ecoregions,
  path_ecoregion_lookup
)

message(
  "HUC8 ecoregion lookup saved: ",
  path_ecoregion_lookup
)

toc()