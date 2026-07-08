# 01_inspect_inputs.R
# Inspect Dryad inputs for HUC8-scale NEP allometry

rm(list = ls())
source("scripts/00_setup.R")


# 1. Define input paths ------------------------------------------------------

path_merit <- file.path(
  dir_dryad,
  "MERIT_Hydro_CONUS_catchments",
  "MERIT_Hydro_CONUS_catchments.shp"
)

path_qmonth <- file.path(dir_dryad, "scaling_dataset_Qmonth_06082024.csv")

path_gpp_jan <- file.path(dir_dryad, "scaling_GPP_pred", "scaling_GPP_pred_Jan.csv")
path_er_jan <- file.path(dir_dryad, "scaling_ER_pred", "scaling_ER_pred_Jan.csv")


# 2. Inspect tabular inputs --------------------------------------------------

qmonth <- read_csv(path_qmonth, n_max = 10) |> clean_names()
gpp_jan <- read_csv(path_gpp_jan, n_max = 10) |> clean_names()
er_jan <- read_csv(path_er_jan, n_max = 10) |> clean_names()

names(qmonth)
names(gpp_jan)
names(er_jan)

glimpse(qmonth)
glimpse(gpp_jan)
glimpse(er_jan)


# 3. Inspect spatial input ---------------------------------------------------

merit <- read_sf(path_merit, n_max = 10) |> clean_names()

names(merit)
glimpse(merit)

st_crs(merit)
st_geometry_type(merit)


# 4. Check likely join fields ------------------------------------------------

intersect(names(qmonth), names(gpp_jan))
intersect(names(qmonth), names(er_jan))
intersect(names(qmonth), names(merit))
intersect(names(gpp_jan), names(er_jan))

