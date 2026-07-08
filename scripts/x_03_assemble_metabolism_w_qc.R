# 03_assemble_metabolism_merit_inputs.R

rm(list = ls())

source("scripts/00_setup.R")


# 2. Longify Maavara ER/GPP predictions ----------------------------------------

months <- month.abb

er <- file.path("data/dryad_data/scaling_ER_pred", paste0("scaling_ER_pred_", months, ".csv")) |>
  set_names(months) |>
  map_dfr(
    \(x) read_csv(x, show_col_types = FALSE),
    .id = "month") |>
  clean_names() |>
  mutate(comid = as.character(comid)) |>
  pivot_longer(
    starts_with("er_pred"),
    names_to = "variable",
    values_to = "er"
  ) |>
  select(comid, month, er) %>% 
  drop_na(er)

gpp <- file.path("data/dryad_data/scaling_GPP_pred", paste0("scaling_GPP_pred_", months, ".csv")) |>
  set_names(months) |>
  map_dfr(
    \(x) read_csv(x, show_col_types = FALSE),
    .id = "month") |>
  clean_names() |>
  mutate(comid = as.character(comid)) |>
  pivot_longer(
    starts_with("gpp_pred"),
    names_to = "variable",
    values_to = "gpp"
  ) |>
  select(comid, month, gpp) %>% 
  drop_na(gpp)

nep <- full_join(gpp, er, by = c("comid", "month")) %>% 
  mutate(nep = gpp - er)

saveRDS(nep, "data/processed/03_maavara_monthly_metabolism.rds")


# Join with centroid QC information

centroid_qc <- read_csv(
  "data/processed/centroid_difference_MERIT_Hydro_CONUS_catchments_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01.csv",
  show_col_types = FALSE) %>% 
  mutate(comid = as.character(comid))


metabolism_allometry <- nep |>
  left_join(centroid_qc, by = "comid")








