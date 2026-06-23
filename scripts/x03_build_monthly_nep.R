# 03_build_monthly_nep.R
# Build monthly GPP, ER, and NEP table

rm(list = ls())
source("scripts/00_setup.R")


# 1. Read attributes ---------------------------------------------------------

qmonth <- read_csv(file.path(dir_dryad, "scaling_dataset_Qmonth_06082024.csv"),
                   show_col_types = FALSE) |>
  clean_names() |>
  select(comid, unitarea, lengthkm, uparea, order)

merit_huc8 <- readRDS("data/processed/merit_huc8_crosswalk.rds")


# 2. Read monthly GPP and ER -------------------------------------------------

gpp <- fs::dir_ls(file.path(dir_dryad, "scaling_GPP_pred"), glob = "*.csv") |>
  map_dfr(\(x) {
    read_csv(x, show_col_types = FALSE) |>
      clean_names() |>
      pivot_longer(-comid, names_to = "month", values_to = "gpp")
  }) |>
  mutate(month = str_remove(month, "gpp_pred_mean_q_"))

er <- fs::dir_ls(file.path(dir_dryad, "scaling_ER_pred"), glob = "*.csv") |>
  map_dfr(\(x) {
    read_csv(x, show_col_types = FALSE) |>
      clean_names() |>
      pivot_longer(-comid, names_to = "month", values_to = "er")
  }) |>
  mutate(month = str_remove(month, "er_pred_mean_q_"))


# 3. Build monthly NEP table -------------------------------------------------

month_levels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

nep_monthly <- gpp |>
  left_join(er, by = c("comid", "month")) |>
  left_join(qmonth, by = "comid") |>
  left_join(merit_huc8, by = "comid") |>
  mutate(
    month = str_to_title(month),
    month = factor(month, levels = month_levels),
    nep = gpp - er
  ) |>
  select(comid, huc8, month, gpp, er, nep, unitarea, lengthkm, uparea, order)


# 4. Diagnostics -------------------------------------------------------------

nep_monthly |>
  summarize(
    n_rows = n(),
    n_comid = n_distinct(comid),
    n_huc8 = n_distinct(huc8, na.rm = TRUE),
    n_missing_huc8 = sum(is.na(huc8)),
    n_missing_nep = sum(is.na(nep))
  )

nep_monthly |>
  count(month)

nep_monthly |>
  count(is.na(huc8), name = "n")


# 5. Save output -------------------------------------------------------------

saveRDS(nep_monthly, "data/processed/monthly_nep.rds", compress = "xz")