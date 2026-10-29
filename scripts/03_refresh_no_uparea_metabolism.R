# Reuse script 03's approved catchments; spatial decisions are independent of rates.
# Run from the repository root after script 01. Do not rerender 03 for this migration.
library(dplyr)
library(sf)
library(readr)
library(janitor)

input <- readRDS("data/processed/01_assembled_no_uparea_monthly_metabolism.rds")
catchments <- read_sf("data/processed/merit_data/03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp") %>%
  clean_names() %>% st_drop_geometry() %>% transmute(comid = as.character(comid))
if (anyDuplicated(catchments$comid)) stop("Duplicate approved catchment COMIDs")
filtered <- input %>% mutate(comid = as.character(comid)) %>%
  inner_join(catchments, by = "comid")
baseline <- readRDS("data/processed/03_filtered_maavara_monthly_metabolism.rds") %>%
  mutate(comid = as.character(comid)) %>% select(comid, month)
if (anyDuplicated(filtered[c("comid", "month")]) ||
    nrow(anti_join(filtered, baseline, by = c("comid", "month"))) ||
    nrow(anti_join(baseline, filtered, by = c("comid", "month")))) {
  stop("Filtered COMID-month coverage differs from original script 03 output")
}
saveRDS(filtered, "data/processed/03_filtered_no_uparea_monthly_metabolism.rds")
write_csv(tibble(stage = c("assembled", "approved_catchments"),
                 n_rows = c(nrow(input), nrow(filtered)),
                 n_comids = c(n_distinct(input$comid), n_distinct(filtered$comid))),
          "data/processed/03_no_uparea_attrition.csv")
