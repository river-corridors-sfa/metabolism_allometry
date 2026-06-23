# 04_explore_poc_huc8.R
# Select and inspect a proof-of-concept HUC8

rm(list = ls())
source("scripts/00_setup.R")


# 1. Read data ---------------------------------------------------------------

nep_monthly <- readRDS("data/processed/monthly_nep.rds")

merit <- st_read(
  file.path(
    dir_dryad,
    "MERIT_Hydro_CONUS_catchments",
    "MERIT_Hydro_CONUS_catchments.shp"
  ),
  quiet = TRUE
) |>
  clean_names()


# 2. Summarize HUC8 candidates ----------------------------------------------

huc8_candidates <- nep_monthly |>
  filter(!is.na(huc8)) |>
  group_by(huc8) |>
  summarize(
    n_comid = n_distinct(comid),
    min_order = min(order, na.rm = TRUE),
    max_order = max(order, na.rm = TRUE),
    n_month = n_distinct(month),
    min_uparea = min(uparea, na.rm = TRUE),
    max_uparea = max(uparea, na.rm = TRUE),
    log10_area_span = log10(max_uparea / min_uparea),
    mean_nep = mean(nep, na.rm = TRUE),
    pct_autotrophic = mean(nep > 0, na.rm = TRUE) * 100,
    pct_heterotrophic = mean(nep < 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) |>
  arrange(desc(max_order), desc(n_comid))


# 3. Select POC HUC8 ---------------------------------------------------------

# poc_huc8 <- huc8_candidates |>
#   filter(n_month == 12,
#          max_order >= 5,
#          n_comid >= 200,
#          log10_area_span >= 2,
#          pct_autotrophic >= 20,
#          pct_autotrophic <= 80) |>
#   arrange(abs(pct_autotrophic - 50), desc(max_order), desc(n_comid)) |>
#   slice(1) |>
#   pull(huc8)

poc_huc8 <- "10100004"


# 4. Build POC spatial data --------------------------------------------------

poc_nep <- nep_monthly |>
  filter(huc8 == poc_huc8)

poc_summary <- poc_nep |>
  group_by(comid) |>
  summarize(
    huc8 = first(huc8),
    order = first(order),
    unitarea = first(unitarea),
    lengthkm = first(lengthkm),
    uparea = first(uparea),
    mean_gpp = mean(gpp, na.rm = TRUE),
    mean_er = mean(er, na.rm = TRUE),
    mean_nep = mean(nep, na.rm = TRUE),
    pct_autotrophic = mean(nep > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

poc_map <- merit |>
  left_join(poc_summary, by = "comid") |>
  filter(huc8 == poc_huc8) |>
  st_transform(common_crs)


# 5. Diagnostics -------------------------------------------------------------

poc_huc8

huc8_candidates |>
  filter(max_order >= 5,
         n_comid >= 200,
         log10_area_span >= 2) |>
  print(n = 25)

poc_nep |>
  count(month)

poc_nep |>
  distinct(comid, order) |>
  count(order)

poc_summary |>
  summarize(
    n_comid = n(),
    min_order = min(order, na.rm = TRUE),
    max_order = max(order, na.rm = TRUE),
    min_uparea = min(uparea, na.rm = TRUE),
    max_uparea = max(uparea, na.rm = TRUE),
    log10_area_span = log10(max_uparea / min_uparea),
    mean_nep = mean(mean_nep, na.rm = TRUE),
    pct_comids_autotrophic = mean(mean_nep > 0, na.rm = TRUE) * 100
  )


# 6. Maps --------------------------------------------------------------------

ggplot(poc_map) +
  geom_sf(aes(fill = order), linewidth = 0) +
  scale_fill_viridis_c(option = "viridis") +
  labs(
    title = paste("POC HUC8", poc_huc8),
    fill = "Stream order"
  )

ggplot(poc_map) +
  geom_sf(aes(fill = log10(uparea)), linewidth = 0) +
  scale_fill_viridis_c(option = "viridis") +
  labs(
    title = paste("POC HUC8", poc_huc8),
    fill = "log10 upstream area"
  )

ggplot(poc_map) +
  geom_sf(aes(fill = mean_nep), linewidth = 0) +
  scale_fill_viridis_c(option = "magma") +
  labs(
    title = paste("POC HUC8", poc_huc8),
    fill = "Mean NEP"
  )


# 7. Save outputs ------------------------------------------------------------

saveRDS(huc8_candidates, "data/processed/huc8_candidates.rds", compress = "xz")
saveRDS(poc_summary, "data/processed/poc_nep_summary.rds", compress = "xz")
saveRDS(poc_map, "data/processed/poc_map.rds", compress = "xz")