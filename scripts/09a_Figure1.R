
# 1. Setup ---------------------------------------------------------------------

rm(list = ls())

start = Sys.time()

source("scripts/00_setup.R")

pacman::p_load(
  forcats,
  glue,
  janitor,
  ggthemes,
  nhdplusTools,
  patchwork,
  scales,
  sf,
  rnaturalearth,
  tigris,
  tidyverse)

theme_set(theme_bw())

dir_data <- "data/processed"
dir_figs <- "figures"

huc_level <- "huc08"
huc_col <- str_replace(huc_level, "huc0", "huc")

month_levels <- month.abb
flux_levels <- c("GPP", "ER")
class_levels <- c("Sublinear", "Linear", "Superlinear", "Unclear")

# Annual class is assigned from the median slope across coherent months.
# HUC8-flux combinations with fewer than this many coherent months are labelled Unclear.
min_valid_months <- 6
min_reaches <- 10
coherence_r2 <- 0.8

slope_breaks <- c(0.8, 1.2)

class_colors <- c(
  "Sublinear" = "#2166ac",
  "Linear" = "#fdae61",
  "Superlinear" = "#b2182b",
  "Unclear" = "#cfcfcf"
)


# 2. Read in data --------------------------------------------------------------

path_flowlines <- file.path(dir_data, "riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp")
path_metrics <- file.path(dir_data, paste0("06_", huc_level, "_allometry_metrics.csv"))
path_l1_lookup <- file.path(dir_data, "08_huc08_l1_ecoregion_lookup.csv")

metrics <- read_csv(path_metrics, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    month = factor(month, levels = month_levels),
    metabolism = recode(
      metabolism_var,
      c_gpp_flux_g_o2_d = "GPP",
      c_er_flux_g_o2_d = "ER",
      c_gpp_flux_gO2_d = "GPP",
      c_er_flux_gO2_d = "ER"
    ),
    valid_fit = n >= min_reaches & r2 >= coherence_r2 & is.finite(slope)
  )


# 3. Establish allometry classes -----------------------------------------------

annual_classes <- metrics |>
  filter(metabolism %in% flux_levels) |>
  group_by(huc, metabolism) |>
  summarise(
    n_months = n_distinct(month),
    n_valid_months = sum(valid_fit, na.rm = TRUE),
    annual_median_slope = if_else(
      n_valid_months >= min_valid_months,
      median(slope[valid_fit], na.rm = TRUE),
      NA_real_
    ),
    annual_median_r2 = if_else(
      n_valid_months >= min_valid_months,
      median(r2[valid_fit], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) |>
  mutate(
    slope_class = case_when(
      is.na(annual_median_slope) ~ "Unclear",
      annual_median_slope < slope_breaks[1] ~ "Sublinear",
      annual_median_slope <= slope_breaks[2] ~ "Linear",
      annual_median_slope > slope_breaks[2] ~ "Superlinear"
    ),
    slope_class = factor(slope_class, levels = class_levels),
    metabolism = factor(metabolism, levels = flux_levels)
  )

annual_classes |>
  count(metabolism, slope_class) |>
  group_by(metabolism) |>
  mutate(pct = n / sum(n) * 100)


# 4. Prep other map elements ---------------------------------------------------

## HUCs
flowlines_sf <- read_sf(path_flowlines) |>
  clean_names()

huc_boundaries <- nhdplusTools::get_huc(
  AOI = st_as_sfc(st_bbox(flowlines_sf)),
  type = huc_level,
  t_srs = 4326
) |>
  clean_names() |>
  mutate(huc = str_pad(as.character(.data[[huc_col]]), width = 8, pad = "0")) |>
  select(huc, geometry) |>
  st_make_valid()

## CONUS boundary
conus_border <- tigris::states(cb = TRUE, resolution = "20m", class = "sf") |>
  clean_names() |>
  filter(!stusps %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS")) |>
  st_transform(st_crs(huc_boundaries)) |>
  summarise(geometry = st_union(geometry))

## Great Lakes
great_lakes <- rnaturalearth::ne_download(
  scale = 50,
  type = "lakes",
  category = "physical",
  returnclass = "sf"
) |>
  clean_names() |>
  filter(name %in% c("Lake Superior", "Lake Michigan", "Lake Huron", "Lake Erie", "Lake Ontario")) |>
  st_transform(st_crs(huc_boundaries))

## Quick test
# ggplot() + 
#   geom_sf(data = great_lakes, fill = "lightblue", color = "white") + 
#   geom_sf(data = conus_border, fill = NA, color = "black", lwd = 0.5) + 
#   theme_map()


## Another less quick but more important test
# ggplot() + 
#   geom_sf(data = huc_boundaries) + 
#   geom_sf(data = great_lakes, fill = "lightblue", color = "white") + 
#   geom_sf(data = conus_border, fill = NA, color = "black", lwd = 0.5) + 
#   theme_map()


huc_flux_grid <- tidyr::expand_grid(
  huc = huc_boundaries$huc,
  metabolism = factor(flux_levels, levels = flux_levels)
) |>
  left_join(annual_classes, by = c("huc", "metabolism")) |>
  mutate(
    slope_class = replace_na(as.character(slope_class), "Unclear"),
    slope_class = factor(slope_class, levels = class_levels),
    n_valid_months = replace_na(n_valid_months, 0L)
  )

map_data <- huc_boundaries |>
  left_join(huc_flux_grid, by = "huc")


# 5. Pull in EPA Level 1 ecoregions --------------------------------------------

build_huc_l1_lookup <- function() {
  
  l1_url <- paste0(
    "https://services7.arcgis.com/oF9CDB4lUYF7Um9q/ArcGIS/rest/services/",
    "NA_Terrestrial_Ecoregions_Level_1/FeatureServer/5/query?",
    "where=1%3D1",
    "&outFields=LEVEL1,NameL1_En",
    "&returnGeometry=true",
    "&f=geojson"
  )
  
  ecoregions_l1 <- read_sf(l1_url) |>
    clean_names() |>
    transmute(
      l1_code = as.character(level1),
      l1_name = name_l1_en,
      geometry
    ) |>
    st_make_valid() |>
    st_transform(5070)
  
  huc_boundaries |>
    st_transform(5070) |>
    st_point_on_surface() |>
    st_join(ecoregions_l1, join = st_within) |>
    st_drop_geometry() |>
    select(huc, l1_code, l1_name) |>
    distinct() |>
    filter(!is.na(l1_name)) |>
    mutate(ecoregion_l1 = paste(l1_code, l1_name))
}

huc_l1_lookup <- if (file.exists(path_l1_lookup)) {
  read_csv(path_l1_lookup, show_col_types = FALSE)
} else {
  lookup <- build_huc_l1_lookup()
  write_csv(lookup, path_l1_lookup)
  lookup
}

huc_l1_lookup <- huc_l1_lookup |>
  clean_names() |>
  mutate(huc = str_pad(as.character(huc), width = 8, pad = "0"))

if (!"l1_name" %in% names(huc_l1_lookup) & "l1_ecoregion" %in% names(huc_l1_lookup)) {
  huc_l1_lookup <- huc_l1_lookup |>
    rename(l1_name = l1_ecoregion)
}

if (!"ecoregion_l1" %in% names(huc_l1_lookup)) {
  huc_l1_lookup <- huc_l1_lookup |>
    mutate(ecoregion_l1 = paste(l1_code, l1_name))
}

huc_l1_lookup <- huc_l1_lookup |>
  mutate(
    l1_code = as.character(l1_code),
    l1_name = str_replace_all(l1_name, "–", "-"),
    ecoregion_l1 = str_replace_all(ecoregion_l1, "–", "-")
  )


ecoregion_classes <- huc_l1_lookup |>
  select(huc, l1_name) |>
  inner_join(huc_flux_grid, by = "huc") |>
  filter(!is.na(l1_name), metabolism %in% flux_levels) |>
  mutate(
    slope_class = factor(slope_class, levels = class_levels),
    metabolism = factor(metabolism, levels = flux_levels)
  )

ecoregion_order <- ecoregion_classes |>
  filter(metabolism == "GPP") |>
  group_by(l1_name) |>
  summarise(
    pct_superlinear = mean(slope_class == "Superlinear", na.rm = TRUE) * 100,
    n_huc8 = n_distinct(huc),
    .groups = "drop"
  ) |>
  arrange(pct_superlinear) |>
  pull(l1_name)

bar_data <- ecoregion_classes |>
  mutate(l1_name = factor(l1_name, levels = ecoregion_order)) |>
  count(metabolism, l1_name, slope_class, name = "n") |>
  group_by(metabolism, l1_name) |>
  mutate(
    n_total = sum(n),
    pct = n / n_total * 100
  ) |>
  ungroup()

n_labels <- bar_data |>
  distinct(metabolism, l1_name, n_total)

write_csv(
  bar_data,
  file.path(dir_data, "09_l1_ecoregion_annual_scaling_class_summary.csv")
)



plot_map <- function(selected_flux){
  map_data %>% 
    filter(metabolism == selected_flux) %>% 
    ggplot() +
    geom_sf(
      aes(fill = slope_class),
      color = alpha("white", 0.18),
      linewidth = 0.015) +
    scale_fill_manual(
      values = class_colors,
      breaks = class_levels,
      drop = FALSE,
      name = NULL) +
    coord_sf(datum = NA) +
    labs(title = paste0(selected_flux, " annual allometry slopes")) +
    geom_sf(data = great_lakes, fill = "lightblue", color = "white") + 
    geom_sf(data = conus_border, fill = NA, color = "black", lwd = 0.5) + 
    theme_map() + 
    theme(legend.position = "bottom")
}

p1_gpp_map <- plot_map("GPP")
p2_er_map <- plot_map("ER")

p_maps <- plot_grid(p1_gpp_map, p2_er_map, nrow = 1)
ggsave("figures/09_01_maps.png", 
       width = 10, height = 4)


bar_order <- bar_data |>
  filter(slope_class == "Unclear") |>
  group_by(l1_name) |>
  summarise(pct_unclear = mean(pct, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(pct_unclear)) |>
  pull(l1_name)

plot_data <- bar_data |>
  mutate(
    l1_name = factor(l1_name, levels = bar_order),
    slope_class = factor(slope_class, c("Superlinear", "Linear", "Sublinear", "Unclear"))
  )

label_data <- n_labels |>
  mutate(l1_name = factor(l1_name, levels = bar_order))

p_bar <- ggplot(plot_data, aes(pct, l1_name, fill = slope_class)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(
    data = label_data,
    aes(x = 103, y = l1_name, label = paste0("n=", n_total)),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.7
  ) +
  facet_wrap(~ metabolism, nrow = 1) +
  scale_fill_manual(
    values = class_colors,
    breaks = c("Unclear", "Sublinear", "Linear", "Superlinear"),
    name = NULL
  ) +
  scale_x_continuous(
    labels = label_percent(scale = 1),
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 115),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14)
  )

p_bar
ggsave("figures/09_02_barcharts.png", 
       width = 8, height = 5)

final_plot <- plot_grid(p_maps, p_bar, ncol = 1)
ggsave("figures/09_Figure1_v1.png", 
       plot = final_plot, 
       width = 12, height = 12)

Sys.time() - start

