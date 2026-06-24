# 02_compare_maavara_and_merit.R
# Test whether Maavara CONUS catchments match downloaded MERIT v00 cat/riv files

rm(list = ls())

source("scripts/00_setup.R")

p_load(tictoc)


# 1. Paths ---------------------------------------------------------------------

path_maavara_cat <- file.path("data/dryad_data/MERIT_Hydro_CONUS_catchments/MERIT_Hydro_CONUS_catchments.shp")
#path_merit_cat   <- file.path("data/merit_data/pfaf_07_cat_3sMERIT.shp")
#path_merit_riv   <- file.path("data/merit_data/pfaf_07_riv_3sMERIT.shp")
path_merit_cat   <- file.path("data/merit_data/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01.shp")
path_merit_riv   <- file.path("data/merit_data/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01.shp")

# 2. Read spatial data ---------------------------------------------------------

message("Reading Maavara catchments...")
maavara_cat <- read_sf(path_maavara_cat) |>
  clean_names() 

message("Reading MERIT catchments...")
merit_cat <- read_sf(path_merit_cat) |>
  clean_names()

message("Reading MERIT rivers...")
merit_riv <- read_sf(path_merit_riv) |>
  clean_names()

conus_map <- rnaturalearth::ne_states(country = "United States of America") %>% 
  filter(!grepl("Alaska|Hawaii|Puerto Rico", name))


# 2. Trim MERIT to Maavara -----------------------------------------------------

tic("filter merit_cat")
merit_cat_filtered <- merit_cat %>% 
  filter(comid %in% unique(maavara_cat$comid))
toc()

tic("filter merit_riv")
merit_riv_filtered <- merit_riv %>% 
  filter(comid %in% unique(maavara_cat$comid))
toc()


excluded_maavara <- maavara_cat %>% 
  filter(!comid %in% unique(merit_cat$comid)) 

conus_map <- rnaturalearth::ne_states(country = "United States of America") %>% 
  filter(!grepl("Alaska|Hawaii|Puerto Rico", name))

ggplot() + 
  geom_sf(data = conus_map, fill = NA) + 
  geom_sf(data = excluded_maavara) + 
  theme_map()

ggplot() + 
  geom_sf(data = conus_map %>% filter(grepl("Dakota|Minnesota", name)), fill = NA) + 
  geom_sf(data = excluded_maavara) + 
  theme_map()


# excluded_maavara_export <- excluded_maavara |>
#   st_make_valid() |>
#   mutate(
#     centroid = st_centroid(geometry),
#     centroid_lon = st_coordinates(centroid)[, 1],
#     centroid_lat = st_coordinates(centroid)[, 2],
#     geometry_wkt = st_as_text(geometry)
#   ) |>
#   st_drop_geometry() |>
#   select(
#     comid,
#     unitarea,
#     lengthkm,
#     centroid_lon,
#     centroid_lat,
#     geometry_wkt
#   )
# 
# write_csv(
#   excluded_maavara_export,
#   file.path("data/maavara_missing_414_comids.csv")
# )


# 3. 

sf_use_s2(FALSE)

tic()
maavara_cat_centroid <- maavara_cat |>
  st_transform(5070) |>
  mutate(centroid = st_centroid(geometry)) |>
  st_transform(4326)

merit_cat_centroid <- merit_cat_filtered |>
  st_transform(5070) |>
  mutate(centroid = st_centroid(geometry)) |>
  st_transform(4326)
toc()

sf_use_s2(TRUE)

ggplot() + 
  geom_sf(data = conus_map, fill = NA) + 
  geom_sf(data = maavara_cat_centroid, alpha = 0.2) + 
  theme_map()
ggsave("figures/260623_maavara_centroid_map_07.png", 
#ggsave("figures/260623_maavara_centroid_map.png", 
       width = 10, height = 6)

ggplot() + 
  geom_sf(data = conus_map, fill = NA) + 
  geom_sf(data = merit_cat_centroid, alpha = 0.2) + 
  theme_map()
ggsave("figures/260623_merit_centroid_map_07.png", 
#ggsave("figures/260623_merit_centroid_map.png", 
       width = 10, height = 6)


# 4. Compare distance between centroids ----------------------------------------


centroid_distance <- maavara_cat_centroid |>
  #st_drop_geometry() |>
  select(comid, maavara_centroid = centroid) |>
  inner_join(
    merit_cat_centroid |>
      st_drop_geometry() |>
      select(comid, merit_centroid = centroid),
    by = "comid") |>
  mutate(
    centroid_distance_m = as.numeric(st_distance(maavara_centroid, merit_centroid, by_element = TRUE)),
    centroid_distance_km = centroid_distance_m / 1000)

ggplot() + 
  geom_sf(data = conus_map, fill = NA) + 
  geom_sf(data = centroid_distance, 
          aes(fill = centroid_distance_km),
          color = NA, alpha = 0.8) + 
  scale_fill_viridis_c(option = "turbo") + 
  theme_map() + 
  theme(legend.position = "bottom")
ggsave("figures/260623_centroid_difference_map_07.png", 
#ggsave("figures/260623_centroid_difference_map.png", 
       width = 10, height = 7)


ggplot(centroid_distance, aes(centroid_distance_km)) + 
  geom_density()








