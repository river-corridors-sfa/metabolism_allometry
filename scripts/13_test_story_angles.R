## 13_test_story_angles.R
##
## PURPOSE
## Test two candidate paper directions using existing null-model outputs:
##
## 1. Network structure sets the baseline; metabolism controls departures.
## 4. Where and when does metabolism depart from the network template?
##
## This script only analyzes existing script-10 outputs.


# 0. Setup ---------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")

tic("run 13")

dir_processed <- "data/processed"
dir_figs <- "figures"


# 1. Read null-model comparison ------------------------------------------------

results <- read_csv(
  file.path(
    dir_processed,
    "10_huc08_null_comparison.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    huc = as.character(huc),
    delta_beta = observed_slope - null_slope_median
  )


# ==============================================================================
# ANGLE 1
# Network structure sets the baseline; metabolism controls departures
# ==============================================================================


# 2. Summary diagnostics --------------------------------------------------------

angle1_summary <- results %>%
  group_by(metabolism_metric) %>%
  summarise(
    n = n(),
    observed_null_cor = cor(
      observed_slope,
      null_slope_median,
      use = "complete.obs"
    ),
    median_observed = median(
      observed_slope,
      na.rm = TRUE
    ),
    median_null = median(
      null_slope_median,
      na.rm = TRUE
    ),
    median_delta_beta = median(
      delta_beta,
      na.rm = TRUE
    ),
    median_abs_delta_beta = median(
      abs(delta_beta),
      na.rm = TRUE
    ),
    pct_above_null = mean(
      delta_beta > 0,
      na.rm = TRUE
    ) * 100,
    pct_below_null = mean(
      delta_beta < 0,
      na.rm = TRUE
    ) * 100,
    pct_slope_exceeds_null = mean(
      slope_exceeds_null,
      na.rm = TRUE
    ) * 100,
    pct_null_exceeding = mean(
      null_exceeding,
      na.rm = TRUE
    ) * 100,
    .groups = "drop"
  )

print(angle1_summary)


# 3. Figure: observed vs null slope --------------------------------------------

p_angle1_baseline <- results %>%
  ggplot(
    aes(
      x = null_slope_median,
      y = observed_slope
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
  facet_wrap(
    ~metabolism_metric
  ) +
  coord_equal() +
  labs(
    x = "Median null-model slope",
    y = "Observed slope",
    title =
      "Observed metabolic allometry relative to the network-constrained baseline"
  ) +
  theme_bw()

ggsave(
  file.path(
    dir_figs,
    "13_angle1_observed_vs_null.png"
  ),
  p_angle1_baseline,
  width = 7,
  height = 5,
  dpi = 300
)


# 4. Figure: distribution of departures ----------------------------------------

p_angle1_delta <- results %>%
  ggplot(
    aes(
      x = delta_beta
    )
  ) +
  geom_histogram(
    bins = 60
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~metabolism_metric,
    scales = "free_y"
  ) +
  labs(
    x = expression(
      Delta * beta ==
        beta[observed] - beta[null]
    ),
    y = "HUC8-month fits",
    title =
      "Departures from network-constrained metabolic allometry"
  ) +
  theme_bw()

ggsave(
  file.path(
    dir_figs,
    "13_angle1_delta_beta.png"
  ),
  p_angle1_delta,
  width = 7,
  height = 5,
  dpi = 300
)


# ==============================================================================
# ANGLE 4
# Where and when does metabolism depart from the network template?
# ==============================================================================


# 5. Seasonal structure ---------------------------------------------------------

seasonal_summary <- results %>%
  group_by(
    metabolism_metric,
    month
  ) %>%
  summarise(
    n = n(),
    median_delta_beta = median(
      delta_beta,
      na.rm = TRUE
    ),
    q25 = quantile(
      delta_beta,
      0.25,
      na.rm = TRUE
    ),
    q75 = quantile(
      delta_beta,
      0.75,
      na.rm = TRUE
    ),
    pct_positive = mean(
      delta_beta > 0,
      na.rm = TRUE
    ) * 100,
    pct_slope_exceeds_null = mean(
      slope_exceeds_null,
      na.rm = TRUE
    ) * 100,
    .groups = "drop"
  )

print(seasonal_summary)


p_angle4_season <- results %>%
  mutate(
    month = factor(
      month,
      levels = 1:12,
      labels = month.abb
    )
  ) %>%
  ggplot(
    aes(
      x = month,
      y = delta_beta
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  facet_wrap(
    ~metabolism_metric,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = expression(Delta * beta),
    title =
      "Seasonality of departures from network-constrained allometry"
  ) +
  theme_bw()

ggsave(
  file.path(
    dir_figs,
    "13_angle4_seasonality.png"
  ),
  p_angle4_season,
  width = 9,
  height = 5,
  dpi = 300
)


# 6. Persistence within HUC8s ---------------------------------------------------

huc_summary <- results %>%
  group_by(
    metabolism_metric,
    huc
  ) %>%
  summarise(
    n_months = n(),
    median_delta_beta = median(
      delta_beta,
      na.rm = TRUE
    ),
    mean_delta_beta = mean(
      delta_beta,
      na.rm = TRUE
    ),
    sd_delta_beta = sd(
      delta_beta,
      na.rm = TRUE
    ),
    pct_positive = mean(
      delta_beta > 0,
      na.rm = TRUE
    ) * 100,
    pct_negative = mean(
      delta_beta < 0,
      na.rm = TRUE
    ) * 100,
    sign_consistency = pmax(
      pct_positive,
      pct_negative
    ),
    .groups = "drop"
  ) %>%
  filter(
    n_months >= 6
  )


p_angle4_persistence <- huc_summary %>%
  ggplot(
    aes(
      x = median_delta_beta,
      y = sign_consistency
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 75,
    linetype = "dashed"
  ) +
  geom_point(
    alpha = 0.35
  ) +
  facet_wrap(
    ~metabolism_metric
  ) +
  labs(
    x = expression(
      "Median HUC8 " * Delta * beta
    ),
    y = "% months in dominant direction",
    title =
      "Persistence of metabolic departures within HUC8s"
  ) +
  theme_bw()

ggsave(
  file.path(
    dir_figs,
    "13_angle4_huc_persistence.png"
  ),
  p_angle4_persistence,
  width = 7,
  height = 5,
  dpi = 300
)


# 7. Month-to-month spatial consistency ----------------------------------------

month_cor <- results %>%
  select(metabolism_metric, huc, month, delta_beta) %>%
  pivot_wider(
    names_from = month,
    values_from = delta_beta
  ) %>%
  group_by(metabolism_metric) %>%
  group_modify(~ {
    
    .x %>%
      select(-huc) %>%
      cor(use = "pairwise.complete.obs") %>%
      as.table() %>%
      as.data.frame() %>%
      as_tibble() %>%
      rename(
        month_1 = Var1,
        month_2 = Var2,
        correlation = Freq
      )
  }) %>%
  ungroup()


p_angle4_month_cor <- month_cor %>%
  ggplot(
    aes(
      x = month_1,
      y = month_2,
      fill = correlation
    )
  ) +
  geom_tile() +
  facet_wrap(
    ~metabolism_metric
  ) +
  scale_fill_gradient2(
    midpoint = 0
  ) +
  labs(
    x = "Month",
    y = "Month",
    fill = "Correlation",
    title =
      "Month-to-month persistence of spatial departures"
  ) +
  theme_bw()

ggsave(
  file.path(
    dir_figs,
    "13_angle4_month_correlation.png"
  ),
  p_angle4_month_cor,
  width = 8,
  height = 5,
  dpi = 300
)


# 8. Story diagnostics ----------------------------------------------------------

story_diagnostics <- huc_summary %>%
  group_by(metabolism_metric) %>%
  summarise(
    n_hucs = n(),
    median_abs_huc_delta = median(
      abs(median_delta_beta),
      na.rm = TRUE
    ),
    pct_hucs_consistent_75 = mean(
      sign_consistency >= 75,
      na.rm = TRUE
    ) * 100,
    pct_hucs_consistent_90 = mean(
      sign_consistency >= 90,
      na.rm = TRUE
    ) * 100,
    .groups = "drop"
  ) %>%
  left_join(
    month_cor %>%
      filter(
        month_1 != month_2
      ) %>%
      group_by(metabolism_metric) %>%
      summarise(
        median_month_correlation = median(
          correlation,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "metabolism_metric"
  )

print(story_diagnostics)


# 9. Save summaries -------------------------------------------------------------

write_csv(
  angle1_summary,
  file.path(
    dir_processed,
    "13_angle1_summary.csv"
  )
)

write_csv(
  seasonal_summary,
  file.path(
    dir_processed,
    "13_angle4_seasonal_summary.csv"
  )
)

write_csv(
  huc_summary,
  file.path(
    dir_processed,
    "13_angle4_huc_summary.csv"
  )
)

write_csv(
  story_diagnostics,
  file.path(
    dir_processed,
    "13_story_diagnostics.csv"
  )
)


toc()