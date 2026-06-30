# 06b_plot_scaling_outputs.R

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(
  nhdplusTools,
  sf,
  tidyverse,
  viridis
)

theme_set(theme_bw())


# 0. Settings ------------------------------------------------------------------

dir_data <- "data/processed"
dir_figs <- "figures"

huc_level <- "huc08"
huc_col <- str_replace(huc_level, "huc0", "huc")

compromised_threshold <- 5

path_flowlines <- "data/processed/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp"

month_levels <- month.abb

huc2_lookup <- tribble(
  ~huc2, ~huc2_region,
  "01", "New England",
  "02", "Mid-Atlantic",
  "03", "South Atlantic-Gulf",
  "04", "Great Lakes",
  "05", "Ohio",
  "06", "Tennessee",
  "07", "Upper Mississippi",
  "08", "Lower Mississippi",
  "09", "Souris-Red-Rainy",
  "10", "Missouri",
  "11", "Arkansas-White-Red",
  "12", "Texas-Gulf",
  "13", "Rio Grande",
  "14", "Upper Colorado",
  "15", "Lower Colorado",
  "16", "Great Basin",
  "17", "Pacific Northwest",
  "18", "California"
)


# 1. Read data -----------------------------------------------------------------

metrics <- read_csv(
  file.path(dir_data, paste0("06_", huc_level, "_allometry_metrics.csv")),
  show_col_types = FALSE
) |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    month = factor(month, levels = month_levels),
    metabolism = recode(
      metabolism_var,
      c_gpp_flux_gO2_d = "GPP",
      c_er_flux_gO2_d = "ER"
    ),
    coherent = r2 >= 0.8
  )

diagnostics <- read_csv(
  file.path(dir_data, paste0("06_", huc_level, "_diagnostics.csv")),
  show_col_types = FALSE
) |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    huc2 = str_sub(huc, 1, 2),
    compromised = pct_reaches_removed >= compromised_threshold |
      pct_length_removed >= compromised_threshold |
      likely_allometry_sensitive
  ) |>
  left_join(huc2_lookup, by = "huc2")

allometry <- metrics |>
  inner_join(diagnostics, by = c("huc_level", "huc"))


# 2. Summary plot: HUCs and reaches dropped ------------------------------------

drop_summary <- diagnostics |>
  summarise(
    hucs_total = n(),
    hucs_compromised = sum(compromised, na.rm = TRUE),
    reaches_retained = sum(n_reaches, na.rm = TRUE),
    reaches_removed = sum(n_reaches_removed, na.rm = TRUE),
    reaches_prefilter = sum(n_reaches + n_reaches_removed, na.rm = TRUE)
  ) |>
  pivot_longer(everything(), names_to = "metric", values_to = "value")

p_drop_summary <- drop_summary |>
  ggplot(aes(x = metric, y = value)) +
  geom_col(fill = "grey35") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Centroid mismatch filtering summary",
    x = NULL,
    y = "Count"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(dir_figs, "06b_centroid_filter_count_summary.png"),
  p_drop_summary,
  width = 7,
  height = 5,
  dpi = 300
)


# 3. Summary plot: regional compromised HUCs -----------------------------------

p_compromised_region <- diagnostics |>
  group_by(huc2, huc2_region) |>
  summarise(
    n_huc8 = n(),
    n_compromised = sum(compromised, na.rm = TRUE),
    pct_compromised = n_compromised / n_huc8 * 100,
    reaches_removed = sum(n_reaches_removed, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(huc2_label = paste(huc2, huc2_region)) |>
  arrange(pct_compromised) |>
  mutate(huc2_label = factor(huc2_label, levels = huc2_label)) |>
  ggplot(aes(x = pct_compromised, y = huc2_label, fill = reaches_removed)) +
  geom_col() +
  scale_fill_viridis_c(option = "C", labels = comma) +
  labs(
    title = "Regional distribution of compromised HUC8s",
    x = "Compromised HUC8s (%)",
    y = "HUC2 region",
    fill = "Removed\nreaches"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(dir_figs, "06b_compromised_hucs_by_region.png"),
  p_compromised_region,
  width = 8,
  height = 6,
  dpi = 300
)


# 4. CONUS map: compromised HUCs ------------------------------------------------

flowlines_sf <- read_sf(path_flowlines) |>
  clean_names()

huc_boundaries <- nhdplusTools::get_huc(
  AOI = st_as_sfc(st_bbox(flowlines_sf)),
  type = huc_level,
  t_srs = st_crs(flowlines_sf)
) |>
  clean_names() |>
  mutate(huc = str_pad(as.character(.data[[huc_col]]), width = 8, pad = "0")) |>
  left_join(
    diagnostics |>
      select(huc, compromised, pct_reaches_removed, pct_length_removed),
    by = "huc"
  )

p_compromised_map <- ggplot() +
  geom_sf(data = huc_boundaries, fill = "grey93", color = NA) +
  geom_sf(
    data = huc_boundaries |> filter(compromised),
    aes(fill = pmax(pct_reaches_removed, pct_length_removed, na.rm = TRUE)),
    color = NA
  ) +
  scale_fill_viridis_c(
    option = "C",
    end = 0.9,
    name = "Max removed\nreaches/length (%)"
  ) +
  coord_sf(datum = NA) +
  labs(title = "HUC8s compromised by centroid mismatch filtering") +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(dir_figs, "06b_compromised_huc08_conus_map.png"),
  p_compromised_map,
  width = 10,
  height = 6,
  dpi = 300
)


# 5. Regional differences in GPP and ER scaling slopes -------------------------

huc_summary <- allometry |>
  group_by(huc, huc2, huc2_region, metabolism) |>
  summarise(
    pct_coherent = mean(coherent, na.rm = TRUE) * 100,
    mean_r2 = mean(r2, na.rm = TRUE),
    median_slope = median(slope, na.rm = TRUE),
    slope_iqr = IQR(slope, na.rm = TRUE),
    .groups = "drop"
  )

huc2_slopes <- huc_summary |>
  group_by(huc2, huc2_region, metabolism) |>
  summarise(
    median_slope = median(median_slope, na.rm = TRUE),
    n_huc8 = n_distinct(huc),
    .groups = "drop"
  ) |>
  mutate(huc2_label = paste(huc2, huc2_region))

huc2_order <- huc2_slopes |>
  filter(metabolism == "GPP") |>
  arrange(median_slope) |>
  pull(huc2_label)

p_huc2_slopes <- huc2_slopes |>
  mutate(huc2_label = factor(huc2_label, levels = huc2_order)) |>
  ggplot(aes(x = median_slope, y = huc2_label, color = metabolism)) +
  geom_line(aes(group = huc2_label), color = "grey75", linewidth = 0.6) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_viridis_d(option = "D", end = 0.85) +
  labs(
    title = "Regional differences in median GPP and ER scaling slopes",
    x = "Median HUC8 allometric slope",
    y = "HUC2 region",
    color = "Flux"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(dir_figs, "06b_regional_gpp_er_scaling_slopes.png"),
  p_huc2_slopes,
  width = 8,
  height = 6,
  dpi = 300
)


# 6. Do coherent HUC8s share GPP and ER scaling rules? -------------------------

gpp_er_pairs <- allometry |>
  select(huc, month, metabolism, slope, r2, coherent) |>
  pivot_wider(
    names_from = metabolism,
    values_from = c(slope, r2, coherent)
  ) |>
  filter(coherent_GPP, coherent_ER)

p_gpp_er_coupling <- gpp_er_pairs |>
  ggplot(aes(x = slope_GPP, y = slope_ER, color = month)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey35") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey45") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey45") +
  geom_point(alpha = 0.35, size = 1.2) +
  coord_equal() +
  scale_color_viridis_d(option = "C", end = 0.95) +
  labs(
    title = "Do coherent HUC8s share GPP and ER scaling rules?",
    x = "GPP allometric slope",
    y = "ER allometric slope",
    color = "Month"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(dir_figs, "06b_coherent_huc8_gpp_er_slope_coupling.png"),
  p_gpp_er_coupling,
  width = 6.5,
  height = 6,
  dpi = 300
)