# One-off CONUS maps of the reach-level cumulative values used in scaling.

rm(list = ls())
source("scripts/00_setup.R")

pacman::p_load(
  tidyverse,
  sf,
  janitor,
  rnaturalearth,
  ggplot2,
  ggthemes,
  viridis
)

dir_project <- "/Users/regi350/offline_work/GitHub/RC/metabolism_allometry"
dir_processed <- file.path(dir_project, "data/processed")
dir_merit <- file.path(dir_processed, "merit_data")
dir_figs <- file.path(dir_project, "figures")

path_allometry_data <- file.path(dir_processed, "005_huc08_allometry_input.rds")
path_riv <- file.path(
  dir_merit,
  "002_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
)

map_month <- 7L

dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)

# This is the filtered reach-month table passed directly to the scaling models.
scaling_values <- readRDS(path_allometry_data) %>%
  mutate(downstream_comid = as.character(downstream_comid))

july_gpp_values <- scaling_values %>%
  filter(month == map_month, metabolism_metric == "GPP") %>%
  select(
    downstream_comid,
    cumulative_gpp_g_o2_d = cumulative_flux_g_o2_d,
    cumulative_area_km2
  ) %>%
  distinct(downstream_comid, .keep_all = TRUE)

flowlines <- read_sf(path_riv, quiet = TRUE) %>%
  clean_names() %>%
  transmute(downstream_comid = as.character(comid), geometry) %>%
  inner_join(july_gpp_values, by = "downstream_comid") %>%
  st_transform(common_crs)

conus_states <- rnaturalearth::ne_states(
  country = "United States of America",
  returnclass = "sf"
) %>%
  clean_names() %>%
  filter(!postal %in% c("AK", "HI", "PR")) %>%
  st_transform(common_crs)

flowline_map_theme <- theme_map() +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    plot.margin = margin(8, 8, 8, 8)
  )

map_base <- ggplot() +
  geom_sf(data = conus_states, fill = "gray97", color = "white", linewidth = 0.25) +
  coord_sf(datum = NA) +
  flowline_map_theme

p_conus_cumulative_gpp <- map_base +
  geom_sf(
    data = flowlines,
    aes(color = cumulative_gpp_g_o2_d),
    linewidth = 0.22,
    alpha = 0.9
  ) +
  scale_color_viridis_c(
    option = "magma",
    trans = "log10",
    name = "Cumulative GPP\n(g O₂ d⁻¹)",
    labels = scales::label_scientific(digits = 2),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(7, "cm"),
      barheight = grid::unit(0.35, "cm")
    )
  )

p_conus_cumulative_area <- map_base +
  geom_sf(
    data = flowlines,
    aes(color = cumulative_area_km2),
    linewidth = 0.22,
    alpha = 0.9
  ) +
  scale_color_viridis_c(
    option = "magma",
    trans = "log10",
    name = "Cumulative area\n(km²)",
    labels = scales::label_number(
      scale_cut = scales::cut_short_scale(),
      accuracy = 1
    ),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(7, "cm"),
      barheight = grid::unit(0.35, "cm")
    )
  )

ggsave(
  file.path(dir_figs, "010a_conus_flowlines_cumulative_gpp_july.png"),
  p_conus_cumulative_gpp,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(dir_figs, "010a_conus_flowlines_cumulative_area.png"),
  p_conus_cumulative_area,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

message(
  "Mapped ", scales::comma(nrow(flowlines)),
  " flowlines using the exact filtered scaling-input values."
)
