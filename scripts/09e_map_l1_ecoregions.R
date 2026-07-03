# 09e_map_l1_ecoregions.R

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(
  cowplot,
  janitor,
  sf,
  tigris,
  tidyverse
)


# 0. Settings ------------------------------------------------------------------

dir_figs <- "figures"

path_l1_lookup <- "data/processed/08_huc08_l1_ecoregion_lookup.csv"

l1_url <- paste0(
  "https://services7.arcgis.com/oF9CDB4lUYF7Um9q/ArcGIS/rest/services/",
  "NA_Terrestrial_Ecoregions_Level_1/FeatureServer/5/query?",
  "where=1%3D1",
  "&outFields=LEVEL1,NameL1_En",
  "&returnGeometry=true",
  "&f=geojson"
)

options(tigris_use_cache = TRUE)


# 1. Ecoregion labels ----------------------------------------------------------

l1_labels <- tribble(
  ~l1_ecoregion,                    ~l1_abbrev,
  "Northern Forests",               "NF",
  "Eastern Temperate Forests",      "ETF",
  "Great Plains",                   "GP",
  "North American Deserts",         "NAD",
  "Northwestern Forested Mountains","NFM",
  "Temperate Sierras",              "TS",
  "Southern Semi-Arid Highlands",   "SSH",
  "Marine West Coast Forests",      "MWCF",
  "Mediterranean California",       "MC",
  "Tropical Humid Forests",         "THF"
)


# 2. Read HUC8 lookup ----------------------------------------------------------

huc_l1_lookup <- read_csv(path_l1_lookup, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    l1_code = as.character(l1_code)
  )

l1_counts <- huc_l1_lookup |>
  count(l1_ecoregion, name = "n_huc8") |>
  left_join(l1_labels, by = "l1_ecoregion")


# 3. Read and clip Level I ecoregions -----------------------------------------

conus_states <- states(cb = TRUE, resolution = "20m", class = "sf") |>
  clean_names() |>
  filter(!stusps %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS")) |>
  st_transform(5070)

conus_outline <- conus_states |>
  summarise(geometry = st_union(geometry))

ecoregions_l1 <- read_sf(l1_url) |>
  clean_names() |>
  transmute(
    l1_code = as.character(level1),
    l1_ecoregion = name_l1_en,
    geometry
  ) |>
  st_make_valid() |>
  st_transform(5070) |>
  filter(l1_ecoregion %in% l1_labels$l1_ecoregion) |>
  st_intersection(conus_outline) |>
  st_collection_extract("POLYGON") |>
  left_join(l1_labels, by = "l1_ecoregion") |>
  left_join(l1_counts, by = c("l1_ecoregion", "l1_abbrev")) |>
  mutate(
    legend_label = paste0(l1_ecoregion, " (", l1_abbrev, ")"),
    facet_label = paste0(l1_abbrev, "\n(n = ", n_huc8, ")")
  )


# 4. Colors --------------------------------------------------------------------

legend_levels <- ecoregions_l1 |>
  st_drop_geometry() |>
  arrange(l1_ecoregion) |>
  distinct(legend_label, .keep_all = TRUE) |>
  pull(legend_label)

l1_colors <- hcl.colors(length(legend_levels), palette = "Dark 3")
names(l1_colors) <- legend_levels

ecoregions_l1 <- ecoregions_l1 |>
  mutate(legend_label = factor(legend_label, levels = legend_levels))


# 5. Big CONUS map -------------------------------------------------------------

p_big <- ggplot() +
  geom_sf(
    data = ecoregions_l1,
    aes(fill = legend_label),
    color = "white",
    linewidth = 0.15
  ) +
  geom_sf(
    data = conus_outline,
    fill = NA,
    color = "grey20",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = l1_colors,
    name = NULL,
    guide = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(linewidth = 0)
    )
  ) +
  coord_sf(datum = NA) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )


# 6. Small multiples -----------------------------------------------------------

facet_order <- ecoregions_l1 |>
  st_drop_geometry() |>
  arrange(l1_ecoregion) |>
  distinct(l1_abbrev, facet_label) |>
  pull(facet_label)

facet_data <- ecoregions_l1 |>
  st_drop_geometry() |>
  distinct(l1_abbrev, facet_label) |>
  pull(l1_abbrev) |>
  map_dfr(\(selected_l1) {
    ecoregions_l1 |>
      mutate(
        selected_l1 = selected_l1,
        is_selected = l1_abbrev == selected_l1
      )
  }) |>
  mutate(facet_label = factor(facet_label, levels = facet_order))

p_facets <- ggplot(facet_data) +
  geom_sf(
    aes(fill = is_selected),
    color = "white",
    linewidth = 0.05
  ) +
  geom_sf(
    data = conus_outline,
    fill = NA,
    color = "grey50",
    linewidth = 0.12
  ) +
  facet_wrap(~facet_label, ncol = 5) +
  scale_fill_manual(
    values = c(
      "FALSE" = "grey92",
      "TRUE" = "#08306B"
    ),
    guide = "none"
  ) +
  coord_sf(datum = NA) +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = unit(0.45, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )


# 7. Combine and save ----------------------------------------------------------

p_09e_l1_map <- plot_grid(
  p_big,
  p_facets,
  ncol = 1,
  rel_heights = c(1.25, 1)
)

p_09e_l1_map

ggsave(
  file.path(dir_figs, "09e_l1_ecoregion_map.png"),
  p_09e_l1_map,
  width = 11,
  height = 11,
  dpi = 400
)

ggsave(
  file.path(dir_figs, "09e_l1_ecoregion_map_big_only.png"),
  p_big,
  width = 11,
  height = 6,
  dpi = 400
)

ggsave(
  file.path(dir_figs, "09e_l1_ecoregion_map_facets_only.png"),
  p_facets,
  width = 11,
  height = 6,
  dpi = 400
)



