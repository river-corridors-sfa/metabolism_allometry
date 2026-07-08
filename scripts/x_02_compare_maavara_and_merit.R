# 02_compare_maavara_and_merit.R
# Test whether Maavara CONUS catchments match downloaded MERIT v00 cat/riv files

rm(list = ls())

source("scripts/00_setup.R")

p_load(tictoc)

#merit_version = "07_bugfix"
merit_version = "07"

# 1. Paths ---------------------------------------------------------------------

path_maavara_cat <- file.path("data/dryad_data/MERIT_Hydro_CONUS_catchments/MERIT_Hydro_CONUS_catchments.shp")
#path_merit_cat   <- file.path("data/merit_data/pfaf_07_cat_3sMERIT.shp")
#path_merit_riv   <- file.path("data/merit_data/pfaf_07_riv_3sMERIT.shp")
path_merit_cat   <- file.path("data/merit_data/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01.shp")
path_merit_riv   <- file.path("data/merit_data/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01.shp")
#path_merit_cat   <- file.path("data/merit_data/cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp")
#path_merit_riv   <- file.path("data/merit_data/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_bugfix1.shp")


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

rm(merit_cat)
rm(merit_riv)

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


# 3. Calculate centroids -------------------------------------------------------

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
ggsave(paste0("figures/260623_maavara_centroid_map_", merit_version, ".png"), 
#ggsave("figures/260623_maavara_centroid_map.png", 
       width = 10, height = 6)

ggplot() + 
  geom_sf(data = conus_map, fill = NA) + 
  geom_sf(data = merit_cat_centroid, alpha = 0.2) + 
  theme_map()
ggsave(paste0("figures/260623_merit_centroid_map_", merit_version, ".png"), 
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
ggsave(paste0("figures/260623_centroid_difference_map_", merit_version, ".png"), 
#ggsave("figures/260623_centroid_difference_map.png", 
       width = 10, height = 7)

ggplot(centroid_distance, aes(centroid_distance_km)) + 
  geom_density()

write_sf(merit_cat_filtered, paste0("data/processed/cat_pfaf_7_MERIT_Hydro_", merit_version, "_Basins_CONUS.shp"))
write_sf(merit_riv_filtered, paste0("data/processed/riv_pfaf_7_MERIT_Hydro_", merit_version, "_Basins_CONUS.shp"))

# Percent >1 km
(nrow(centroid_distance %>% filter(centroid_distance_m > 1000))/
  nrow(centroid_distance))*100


