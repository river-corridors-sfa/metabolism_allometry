# 06c_plot_null_coherence_outputs.R

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

path_flowlines <- "data/processed/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp"
path_coherence <- file.path(
  dir_data,
  "06a_null_coherence",
  paste0("06a_", huc_level, "_observed_vs_null_coherence.csv")
)

month_levels <- month.abb


# 1. Read data -----------------------------------------------------------------

coherence <- read_csv(path_coherence, show_col_types = FALSE) |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    month = factor(month, levels = month_levels),
    metabolism = recode(
      metabolism_var,
      c_gpp_flux_gO2_d = "GPP",
      c_er_flux_gO2_d = "ER"
    )
  ) |>
  filter(!is.na(r2), !is.na(r2_null_95))


# 2. Observed R2 versus null R2 -------------------------------------------------

p_r2_observed_vs_null <- coherence |>
  ggplot(aes(x = r2_null_95, y = r2)) +
  geom_bin2d(bins = 80) +
  geom_abline(slope = 1, intercept = 0, color = "grey30", linetype = "dashed") +
  geom_hline(yintercept = 0.8, color = "grey45", linetype = "dotted") +
  facet_wrap(~ metabolism) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", labels = scales::comma) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Observed allometric coherence exceeds network-preserving null expectations",
    x = "Null 95th percentile R²",
    y = "Observed R²",
    fill = "Count"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(dir_figs, "06c_observed_vs_null_r2.png"),
  p_r2_observed_vs_null,
  width = 8,
  height = 5.5,
  dpi = 300
)


# 3. R2 margin above null -------------------------------------------------------

p_r2_margin <- coherence |>
  ggplot(aes(x = r2_minus_null_95, fill = metabolism, color = metabolism)) +
  geom_density(alpha = 0.35, linewidth = 0.9) +
  geom_vline(xintercept = 0, color = "grey35", linetype = "dashed") +
  scale_fill_viridis_d(option = "D", end = 0.85) +
  scale_color_viridis_d(option = "D", end = 0.85) +
  labs(
    title = "Observed R² margin above shuffled-network null",
    x = "Observed R² - null 95th percentile R²",
    y = "Density",
    fill = "Flux",
    color = "Flux"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(dir_figs, "06c_r2_margin_density.png"),
  p_r2_margin,
  width = 7,
  height = 5,
  dpi = 300
)


# 4. Monthly coherence summary -------------------------------------------------

monthly_coherence <- coherence |>
  group_by(month, metabolism) |>
  summarise(
    n_tests = n(),
    pct_observed_coherent = mean(observed_coherent, na.rm = TRUE) * 100,
    pct_stronger_than_null = mean(stronger_than_null, na.rm = TRUE) * 100,
    pct_rcc_like_coherent = mean(rcc_like_coherent, na.rm = TRUE) * 100,
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = starts_with("pct_"),
    names_to = "metric",
    values_to = "pct_tests"
  ) |>
  mutate(
    metric = recode(
      metric,
      pct_observed_coherent = "Observed R² ≥ 0.8",
      pct_stronger_than_null = "Observed R² > null 95th",
      pct_rcc_like_coherent = "RCC-like coherent"
    ),
    metric = factor(
      metric,
      levels = c(
        "Observed R² ≥ 0.8",
        "Observed R² > null 95th",
        "RCC-like coherent"
      )
    )
  )

p_monthly_coherence <- monthly_coherence |>
  ggplot(aes(x = month, y = pct_tests, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ metabolism, nrow = 1) +
  scale_fill_viridis_d(option = "C", end = 0.9) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Monthly frequency of coherent HUC8 allometries",
    x = NULL,
    y = "HUC8-month tests (%)",
    fill = NULL
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(dir_figs, "06c_monthly_coherence_summary.png"),
  p_monthly_coherence,
  width = 9,
  height = 5,
  dpi = 300
)


# 5. HUC8 summary map -----------------------------------------------------------

huc_coherence <- coherence |>
  group_by(huc, metabolism) |>
  summarise(
    n_months = n_distinct(month),
    pct_observed_coherent = mean(observed_coherent, na.rm = TRUE) * 100,
    pct_rcc_like_coherent = mean(rcc_like_coherent, na.rm = TRUE) * 100,
    median_r2_margin = median(r2_minus_null_95, na.rm = TRUE),
    .groups = "drop"
  )

flowlines_sf <- read_sf(path_flowlines) |>
  clean_names()

huc_boundaries <- nhdplusTools::get_huc(
  AOI = st_as_sfc(st_bbox(flowlines_sf)),
  type = huc_level,
  t_srs = 4326
) |>
  clean_names() |>
  mutate(huc = str_pad(as.character(.data[[huc_col]]), width = 8, pad = "0"))

huc_coherence_map <- huc_boundaries |>
  left_join(huc_coherence, by = "huc")

p_rcc_like_map <- huc_coherence_map |>
  filter(!is.na(metabolism)) |>
  ggplot() +
  geom_sf(aes(fill = pct_rcc_like_coherent), color = NA) +
  facet_wrap(~ metabolism, nrow = 1) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 100),
    na.value = "grey92",
    name = "RCC-like\nmonths (%)"
  ) +
  coord_sf(datum = NA) +
  labs(title = "Frequency of RCC-like coherent allometry by HUC8") +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(dir_figs, "06c_rcc_like_coherence_huc08_map.png"),
  p_rcc_like_map,
  width = 10,
  height = 6,
  dpi = 300
)


# 6. Spatial margin above null --------------------------------------------------

p_r2_margin_map <- huc_coherence_map |>
  filter(!is.na(metabolism)) |>
  ggplot() +
  geom_sf(aes(fill = median_r2_margin), color = NA) +
  facet_wrap(~ metabolism, nrow = 1) +
  scale_fill_viridis_c(
    option = "C",
    na.value = "grey92",
    name = "Median R²\nmargin"
  ) +
  coord_sf(datum = NA) +
  labs(title = "Median observed R² margin above null by HUC8") +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(dir_figs, "06c_r2_margin_huc08_map.png"),
  p_r2_margin_map,
  width = 10,
  height = 6,
  dpi = 300
)


# 7. Compact table --------------------------------------------------------------

coherence_summary <- coherence |>
  group_by(metabolism) |>
  summarise(
    n_tests = n(),
    n_observed_coherent = sum(observed_coherent, na.rm = TRUE),
    n_stronger_than_null = sum(stronger_than_null, na.rm = TRUE),
    n_rcc_like_coherent = sum(rcc_like_coherent, na.rm = TRUE),
    pct_observed_coherent = n_observed_coherent / n_tests * 100,
    pct_rcc_like_coherent = n_rcc_like_coherent / n_tests * 100,
    median_r2_margin = median(r2_minus_null_95, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  coherence_summary,
  file.path(dir_data, "06a_null_coherence", paste0("06c_", huc_level, "_coherence_plot_summary.csv"))
)

coherence_summary