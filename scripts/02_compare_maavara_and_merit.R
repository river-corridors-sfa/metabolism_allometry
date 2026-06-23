# 02_compare_maavara_and_merit.R
# Test whether Maavara CONUS catchments match downloaded MERIT v00 cat/riv files

rm(list = ls())

source("scripts/00_setup.R")

p_load(tictoc)


# 1. Paths ---------------------------------------------------------------------

path_maavara_cat <- file.path("data/dryad_data/MERIT_Hydro_CONUS_catchments/MERIT_Hydro_CONUS_catchments.shp")
path_merit_cat   <- file.path("data/merit_data/pfaf_07_cat_3sMERIT.shp")
path_merit_riv   <- file.path("data/merit_data/pfaf_07_riv_3sMERIT.shp")


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


# 2. Trim MERIT to CONUS -------------------------------------------------------

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









  





