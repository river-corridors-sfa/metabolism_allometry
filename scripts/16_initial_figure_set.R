## 16_initial_figure_set.R
##
## PURPOSE
## Create an initial manuscript figure set testing:
## Q1. How much metabolic allometry is reproduced by network structure?
## Q2. Where and when does metabolism depart from that expectation?


# 0. Setup ---------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(patchwork)

dir_processed <- "data/processed"
dir_merit <- file.path(dir_processed, "merit_data")
dir_figs <- "figures"

cols <- c(
  "ER" = "#D95F02",
  "GPP" = "#1B9E77"
)


# 1. Read existing results ------------------------------------------------------

results <- read_csv(
  file.path(dir_processed, "10_huc08_null_comparison.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    huc = as.character(huc),
    delta_beta = observed_slope - null_slope_median
  )

variance_summary <- read_csv(
  file.path(dir_processed, "15_network_variance_summary.csv"),
  show_col_types = FALSE
)

area_sensitivity <- readRDS(
  file.path(dir_processed, "14_null_area_bin_sensitivity.rds")
)

uparea_sensitivity <- read_csv(
  file.path(dir_processed, "12_uparea_allometry_july.csv"),
  show_col_types = FALSE
)

sf_use_s2(FALSE)

huc_shapes <- read_sf(
  file.path(
    dir_merit,
    "03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
  )
) %>%
  clean_names() %>%
  st_make_valid() %>%
  mutate(huc = as.character(huc)) %>%
  group_by(huc) %>%
  summarise() %>%
  st_transform(5070)


# 2. Figure 1: observed metabolic allometry ------------------------------------

observed_huc <- results %>%
  group_by(huc, metabolism_metric) %>%
  summarise(
    observed_slope = median(observed_slope, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    metabolism_metric %in% c("ER", "GPP"),
    is.finite(observed_slope)
  )

fig1_data <- huc_shapes %>%
  inner_join(observed_huc, by = "huc")

fig1 <- ggplot(fig1_data) +
  geom_sf(
    aes(fill = observed_slope),
    linewidth = 0
  ) +
  scale_fill_viridis_c(
    option = "C",
    name = expression(beta[observed])
  ) +
  facet_wrap(~metabolism_metric) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12)
  )

ggsave(
  file.path(dir_figs, "16_fig1_observed_allometry.png"),
  fig1,
  width = 9,
  height = 5,
  dpi = 300
)


# 3. Figure 2: observed vs network-constrained expectation ----------------------

fig2_labels <- variance_summary %>%
  transmute(
    metabolism_metric,
    label = paste0(
      "Variance reproduced = ",
      scales::percent(
        variance_reproduced,
        accuracy = 1
      )
    )
  )

fig2 <- results %>%
  filter(
    metabolism_metric %in% c("ER", "GPP")
  ) %>%
  ggplot(
    aes(
      x = null_slope_median,
      y = observed_slope,
      color = metabolism_metric
    )
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    alpha = 0.15,
    size = 0.7
  ) +
  geom_text(
    data = fig2_labels,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.4,
    size = 3.5
  ) +
  scale_color_manual(
    values = cols,
    guide = "none"
  ) +
  facet_wrap(~metabolism_metric) +
  coord_equal() +
  labs(
    x = expression(beta[null]),
    y = expression(beta[observed])
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(dir_figs, "16_fig2_observed_vs_null.png"),
  fig2,
  width = 7,
  height = 5,
  dpi = 300
)


# 4. Figure 3: robustness --------------------------------------------------------

uparea_plot_data <- uparea_sensitivity %>%
  select(
    huc,
    metabolism,
    model,
    slope
  ) %>%
  pivot_wider(
    names_from = model,
    values_from = slope
  )

p_uparea <- uparea_plot_data %>%
  ggplot(
    aes(
      x = Full,
      y = `No uparea`,
      color = metabolism
    )
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    alpha = 0.25,
    size = 0.8
  ) +
  scale_color_manual(
    values = cols,
    guide = "none"
  ) +
  facet_wrap(~metabolism) +
  coord_equal() +
  labs(
    x = expression(beta["full RF"]),
    y = expression(beta["RF without upstream area"])
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )


null_sensitivity <- area_sensitivity %>%
  group_by(
    area_bins,
    huc,
    month,
    metabolism_metric
  ) %>%
  summarise(
    null_slope = median(slope, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  inner_join(
    results %>%
      select(
        huc,
        month,
        metabolism_metric,
        observed_slope
      ),
    by = c(
      "huc",
      "month",
      "metabolism_metric"
    )
  ) %>%
  mutate(
    delta_beta = observed_slope - null_slope,
    area_bins = factor(area_bins)
  )

p_null <- null_sensitivity %>%
  ggplot(
    aes(
      x = area_bins,
      y = delta_beta,
      color = metabolism_metric
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  scale_color_manual(
    values = cols,
    guide = "none"
  ) +
  facet_wrap(~metabolism_metric) +
  labs(
    x = "Upstream-area bins preserved by null",
    y = expression(Delta * beta)
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid.minor = element_blank()
  )


fig3 <- p_uparea / p_null +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  file.path(dir_figs, "16_fig3_robustness.png"),
  fig3,
  width = 8,
  height = 9,
  dpi = 300
)


# 5. Figure 4: where observed scaling departs -----------------------------------

huc_delta <- results %>%
  group_by(
    huc,
    metabolism_metric
  ) %>%
  summarise(
    delta_beta = median(
      delta_beta,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    metabolism_metric %in% c("ER", "GPP"),
    is.finite(delta_beta)
  )

fig4_data <- huc_shapes %>%
  inner_join(
    huc_delta,
    by = "huc"
  )

delta_limit <- quantile(
  abs(fig4_data$delta_beta),
  0.98,
  na.rm = TRUE
)

fig4 <- ggplot(fig4_data) +
  geom_sf(
    aes(fill = delta_beta),
    linewidth = 0
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -delta_limit,
      delta_limit
    ),
    oob = scales::squish,
    name = expression(Delta * beta)
  ) +
  facet_wrap(~metabolism_metric) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12)
  )

ggsave(
  file.path(dir_figs, "16_fig4_delta_beta_map.png"),
  fig4,
  width = 9,
  height = 5,
  dpi = 300
)


# 6. Figure 5: ecoregion x month departures ------------------------------------

huc_ecoregions <- read_csv(
  file.path(
    dir_processed,
    "04b_huc08_l1_ecoregion_lookup.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    huc = str_pad(
      as.character(huc),
      8,
      pad = "0"
    )
  )

eco_month <- results %>%
  mutate(
    huc = str_pad(
      as.character(huc),
      8,
      pad = "0"
    )
  ) %>%
  left_join(
    huc_ecoregions,
    by = "huc"
  ) %>%
  filter(
    metabolism_metric %in% c("ER", "GPP"),
    !is.na(l1_ecoregion),
    is.finite(delta_beta)
  ) %>%
  group_by(
    metabolism_metric,
    l1_ecoregion,
    month
  ) %>%
  summarise(
    n_hucs = n_distinct(huc),
    median_delta_beta = median(
      delta_beta,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(
    n_hucs >= 5
  ) %>%
  mutate(
    month = factor(
      month,
      levels = 1:12,
      labels = month.abb
    )
  )

delta_limit <- quantile(
  abs(eco_month$median_delta_beta),
  0.98,
  na.rm = TRUE
)

fig5 <- eco_month %>%
  ggplot(
    aes(
      x = month,
      y = fct_reorder(
        l1_ecoregion,
        median_delta_beta,
        .fun = median
      ),
      fill = median_delta_beta
    )
  ) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(
    ~metabolism_metric,
    nrow = 1
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -delta_limit,
      delta_limit
    ),
    oob = scales::squish,
    name = expression(Delta * beta)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(
    dir_figs,
    "16_fig5_delta_beta_ecoregion_month.png"
  ),
  fig5,
  width = 10,
  height = 7,
  dpi = 300
)