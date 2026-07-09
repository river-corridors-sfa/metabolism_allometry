# Test - plot monthly NEP for NM

rm(list = ls())

require(pacman)
pacman::p_load(
  tidyverse,
  sf,
  fs,
  janitor,
  tigris, 
  ggthemes
)

green <- "#2E7D32"
gray  <- "#D0D0C8"
red   <- "#B23A48"


# 1. Paths ----------------------------------------------------------------

dir_parent = "/Users/regi350/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/RC/rcsfa-ai4do/maavara_nep"
dir_data  <- file.path(dir_parent, "dryad_data")
dir_figs  <- file.path(dir_parent, "figures")

dir_create(dir_figs)


# 2. New Mexico outline ---------------------------------------------------

options(tigris_use_cache = TRUE)

nm <- states(cb = TRUE, year = 2023) |>
  clean_names() |>
  filter(name == "New Mexico") |>
  st_transform(4326)

ggplot() + 
  geom_sf(data = nm, fill = "forestgreen", alpha = 0.1) + 
  theme_map()

# 3. New Mexico MERIT catchments ------------------------------------------

path_catchments <- file.path(
  dir_data,
  "MERIT_Hydro_CONUS_catchments",
  "MERIT_Hydro_CONUS_catchments.shp"
)

catchments <- read_sf(
  path_catchments,
  wkt_filter = st_as_text(st_geometry(nm))
) |>
  clean_names() |>
  st_make_valid()

nm_catchments <- catchments |>
  st_filter(nm)

nm_comids <- nm_catchments |>
  st_drop_geometry() |>
  distinct(comid)

ggplot() + 
  geom_sf(data = nm, fill = "forestgreen", alpha = 0.1) + 
  geom_sf(data = nm_catchments, fill = NA) + 
  theme_map()


# 4. Monthly NEP ----------------------------------------------------------

read_month_nep <- function(month_abbr) {
  
  gpp <- read_csv(
    file.path(dir_data, "scaling_GPP_pred", paste0("scaling_GPP_pred_", month_abbr, ".csv")),
    show_col_types = FALSE
  ) |>
    clean_names() |>
    semi_join(nm_comids, by = "comid") |>
    rename(gpp = 2)
  
  er <- read_csv(
    file.path(dir_data, "scaling_ER_pred", paste0("scaling_ER_pred_", month_abbr, ".csv")),
    show_col_types = FALSE
  ) |>
    clean_names() |>
    semi_join(nm_comids, by = "comid") |>
    rename(er = 2)
  
  gpp |>
    left_join(er, by = "comid") |>
    mutate(
      month = month_abbr,
      nep = gpp - er
    )
}

monthly_nep <- month.abb |>
  map(read_month_nep) |>
  list_rbind() |>
  mutate(month = factor(month, levels = month.abb))


# 5. Join predictions to catchments ---------------------------------------

nm_nep <- nm_catchments |>
  select(comid, geometry) |>
  left_join(monthly_nep, by = "comid") |>
  filter(!is.na(nep))


# 6. Quick checks ---------------------------------------------------------

monthly_nep |>
  summarise(
    n_comid = n_distinct(comid),
    n_month = n_distinct(month),
    nep_min = min(nep, na.rm = TRUE),
    nep_max = max(nep, na.rm = TRUE)
  ) |>
  print()

monthly_nep |>
  group_by(month) |>
  summarise(
    n = n(),
    nep_median = median(nep, na.rm = TRUE),
    .groups = "drop"
  ) |>
  print(n = 12)


# 7. Plot -----------------------------------------------------------------

nm_map <- nm |>
  st_transform(5070)

nm_nep_map <- nm_nep |>
  st_transform(5070)

p <- ggplot() +
  geom_sf(
    data = nm_nep_map,
    aes(fill = nep),
    color = NA
  ) +
  geom_sf(
    data = nm_map,
    fill = NA,
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_viridis_c(option = "magma") +
  facet_wrap(~ month, ncol = 4) +
  coord_sf(datum = NA) +
  theme_void() +
  labs(fill = "NEP")


# 8. Classify NEP regimes -------------------------------------------------

nep_cols <- c(
  "Heterotrophic"  = "#B23A48",
  "Circumneutral" = "#D0D0C8",
  "Autotrophic"   = "#2E7D32"
)

nm_nep_map <- nm_nep_map |>
  mutate(
    nep_class = case_when(
      nep > 1  ~ "Autotrophic",
      nep < -1 ~ "Heterotrophic",
      TRUE     ~ "Circumneutral"
    ),
    nep_class = factor(
      nep_class,
      levels = c("Heterotrophic", "Circumneutral", "Autotrophic")
    )
  )

# 10. Plot monthly NEP classes --------------------------------------------

p_nep_class <- ggplot() +
  geom_sf(
    data = nm_nep_map,
    aes(fill = nep_class),
    color = NA
  ) +
  geom_sf(
    data = nm_map,
    fill = NA,
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = nep_cols, drop = FALSE) +
  facet_wrap(~ month, ncol = 4) +
  coord_sf(datum = NA) +
  theme_void() +
  labs(fill = "NEP class") +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10)
  )

p_nep_class


# 9. Save plots -------------------------------------------------

ggsave(
  file.path(dir_figs, "new_mexico_monthly_nep_catchments.png"),
  plot = p,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  file.path(dir_figs, "new_mexico_monthly_nep_classes.png"),
  plot = p_nep_class,
  width = 12,
  height = 9,
  dpi = 300
)





