# 01_interrogate_science_adu9843_data.R

rm(list = ls())

require(pacman)
pacman::p_load(
  tidyverse,
  sf,
  arrow,
  fs,
  janitor)


# 1. Paths ----------------------------------------------------------------

dir_parent = "/Users/regi350/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/RC/rcsfa-ai4do/maavara_nep"
dir_data  <- file.path(dir_parent, "dryad_data")


# 2. File inventory -------------------------------------------------------

file_inventory <- dir_info(dir_data, recurse = TRUE, type = "file") |>
  mutate(
    ext = str_to_lower(path_ext(path)),
    size_mb = size / 1024^2, 
    file_name = path_file(path)
  ) |>
  dplyr::select(path, file_name, ext, size_mb) |>
  arrange(desc(size_mb))

write_csv(file_inventory, file.path(dir_out, "file_inventory.csv"))

file_inventory |>
  count(ext, sort = TRUE)

file_inventory |>
  select(file_name, ext, size_mb) |>
  print(n = 100)

#
# a bunch of non-recorded scans of files to figure out what we want
#


# 5. Inspect prediction and geometry schemas ------------------------------

gpp_jan <- read_csv(
  file.path(dir_data, "scaling_GPP_pred", "scaling_GPP_pred_Jan.csv"),
  n_max = 10,
  show_col_types = FALSE
) |>
  clean_names()

er_jan <- read_csv(
  file.path(dir_data, "scaling_ER_pred", "scaling_ER_pred_Jan.csv"),
  n_max = 10,
  show_col_types = FALSE
) |>
  clean_names()


catchments_head <- read_sf(file.path(
  dir_data,
  "MERIT_Hydro_CONUS_catchments",
  "MERIT_Hydro_CONUS_catchments.shp"
)) %>% 
  slice(1:10) %>% 
  clean_names()


glimpse(gpp_jan)
names(gpp_jan)

glimpse(er_jan)
names(er_jan)

glimpse(catchments_head)
names(catchments_head)
st_crs(catchments_head)

