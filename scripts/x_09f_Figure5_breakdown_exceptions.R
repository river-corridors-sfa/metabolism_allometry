

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(
  cowplot,
  janitor,
  nhdplusTools,
  scales,
  sf,
  tigris,
  tidyverse
)


# 0. Settings ------------------------------------------------------------------

dir_data <- "data/processed"
dir_figs <- "figures"

path_metrics <- file.path(dir_data, "06_huc08_allometry_metrics.csv")
path_l1_lookup <- file.path(dir_data, "08_huc08_l1_ecoregion_lookup.csv")

min_reaches <- 10
coherence_r2 <- 0.8
min_valid_months <- 6
seasonal_consistency_threshold <- 0.6

flux_levels <- c("ER", "GPP")

exception_levels <- c(
  "Coherent superlinear",
  "Coherent linear/sublinear",
  "Seasonally inconsistent",
  "Low coherence",
  "Insufficient data"
)

exception_colors <- c(
  "Coherent superlinear" = "#b2182b",
  "Coherent linear/sublinear" = "#fdae61",
  "Seasonally inconsistent" = "#762a83",
  "Low coherence" = "#2166ac",
  "Insufficient data" = "#cfcfcf"
)

options(tigris_use_cache = TRUE)


# 1. Read data -----------------------------------------------------------------

metrics <- read_csv(path_metrics, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    metabolism = recode(
      metabolism_var,
      c_gpp_flux_g_o2_d = "GPP",
      c_er_flux_g_o2_d = "ER",
      c_gpp_flux_gO2_d = "GPP",
      c_er_flux_gO2_d = "ER"
    )
  ) |>
  filter(metabolism %in% flux_levels)

huc_l1_lookup <- read_csv(path_l1_lookup, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    l1_ecoregion = str_replace_all(l1_ecoregion, "–", "-")
  )


# 2. Classify monthly fits -----------------------------------------------------

monthly_classes <- metrics |>
  mutate(
    has_enough_data = n >= min_reaches & is.finite(slope) & is.finite(r2),
    is_coherent = has_enough_data & r2 >= coherence_r2,
    slope_class = case_when(
      !is_coherent ~ NA_character_,
      slope < 0.8 ~ "Sublinear",
      slope <= 1.2 ~ "Linear",
      slope > 1.2 ~ "Superlinear"
    )
  )


# 3. Summarize annual HUC8 exception classes ----------------------------------

dominant_class <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

huc_exception_summary <- monthly_classes |>
  group_by(huc, metabolism) |>
  summarise(
    n_months = n(),
    n_data_months = sum(has_enough_data, na.rm = TRUE),
    n_coherent_months = sum(is_coherent, na.rm = TRUE),
    prop_coherent = n_coherent_months / n_months,
    median_r2 = median(r2[has_enough_data], na.rm = TRUE),
    median_slope = median(slope[is_coherent], na.rm = TRUE),
    dominant_slope_class = dominant_class(slope_class),
    dominant_slope_fraction = max(prop.table(table(slope_class)), na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    dominant_slope_fraction = if_else(
      is.infinite(dominant_slope_fraction),
      NA_real_,
      dominant_slope_fraction
    ),
    exception_class = case_when(
      n_data_months < min_valid_months ~ "Insufficient data",
      n_coherent_months < min_valid_months ~ "Low coherence",
      dominant_slope_fraction < seasonal_consistency_threshold ~ "Seasonally inconsistent",
      dominant_slope_class == "Superlinear" ~ "Coherent superlinear",
      dominant_slope_class %in% c("Linear", "Sublinear") ~ "Coherent linear/sublinear",
      TRUE ~ "Low coherence"
    ),
    exception_class = factor(exception_class, levels = exception_levels),
    metabolism = factor(metabolism, levels = flux_levels)
  ) |>
  left_join(huc_l1_lookup, by = "huc")

write_csv(
  huc_exception_summary,
  file.path(dir_data, "09f_figure5_breakdown_exception_huc8_summary.csv")
)


# 4. Read HUC8 boundaries ------------------------------------------------------

conus_border <- tigris::states(cb = TRUE, resolution = "20m", class = "sf") |>
  clean_names() |>
  filter(!stusps %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS")) |>
  st_transform(5070) |>
  summarise(geometry = st_union(geometry))

huc08_boundaries <- nhdplusTools::get_huc(
  AOI = conus_border,
  type = "huc08"
) |>
  st_transform(5070) |>
  clean_names() |>
  mutate(huc = str_pad(as.character(huc8), width = 8, pad = "0")) |>
  select(huc, geometry)


# 5. Prepare map and ecoregion summary data -----------------------------------

map_data <- huc08_boundaries |>
  left_join(huc_exception_summary, by = "huc") |>
  mutate(
    exception_class = replace_na(
      as.character(exception_class),
      "Insufficient data"
    ),
    exception_class = factor(exception_class, levels = exception_levels)
  )

bar_data <- huc_exception_summary |>
  filter(!is.na(l1_ecoregion)) |>
  count(l1_ecoregion, metabolism, exception_class, name = "n_huc8") |>
  group_by(l1_ecoregion, metabolism) |>
  mutate(
    prop_huc8 = n_huc8 / sum(n_huc8),
    label_n = sum(n_huc8)
  ) |>
  ungroup()

ecoregion_order <- bar_data |>
  filter(exception_class %in% c("Low coherence", "Insufficient data", "Seasonally inconsistent")) |>
  group_by(l1_ecoregion) |>
  summarise(exception_fraction = mean(prop_huc8, na.rm = TRUE), .groups = "drop") |>
  arrange(exception_fraction) |>
  pull(l1_ecoregion)

bar_data <- bar_data |>
  mutate(
    l1_ecoregion = factor(l1_ecoregion, levels = ecoregion_order),
    exception_class = factor(exception_class, levels = exception_levels)
  )

write_csv(
  bar_data,
  file.path(dir_data, "09f_figure5_breakdown_exception_ecoregion_summary.csv")
)


# 6. Plot ----------------------------------------------------------------------

p_map <- map_data |>
  ggplot() +
  geom_sf(
    aes(fill = exception_class),
    color = alpha("white", 0.2),
    linewidth = 0.015
  ) +
  geom_sf(
    data = conus_border,
    fill = NA,
    color = "black",
    linewidth = 0.2
  ) +
  facet_wrap(~ metabolism, nrow = 1) +
  scale_fill_manual(
    values = exception_colors,
    name = NULL,
    drop = FALSE
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(0.6, "lines")
  )

p_bars <- bar_data |>
  ggplot(aes(
    x = l1_ecoregion,
    y = prop_huc8,
    fill = exception_class
  )) +
  geom_col(width = 0.75) +
  facet_wrap(~ metabolism, nrow = 1) +
  coord_flip() +
  scale_y_continuous(
    labels = label_percent(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_fill_manual(
    values = exception_colors,
    name = NULL,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Fraction of HUC8s"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10)
  )

p_figure5 <- cowplot::plot_grid(
  p_map,
  p_bars,
  ncol = 1,
  rel_heights = c(1.1, 1),
  align = "v"
)

p_figure5


# 7. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figs, "09f_figure5_breakdown_exceptions.png"),
  p_figure5,
  width = 10,
  height = 9,
  dpi = 400
)





