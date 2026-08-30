# 1. Setup ----------------------------------------------------------------

require(pacman)
pacman::p_load(tidyverse)

DATA_PATH <- "/Users/regi350/offline_work/GitHub/RC/metabolism_allometry/data/processed/005_huc08_allometry_results.csv"
OUTPUT_DIR <- "/Users/regi350/offline_work/GitHub/RC/metabolism_allometry/figures/260818_test"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

theme_set(
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
)

month_labels <- month.abb


# 2. Read and verify -------------------------------------------------------

raw <- readr::read_csv(
  DATA_PATH,
  col_types = cols(huc = col_character()),
  show_col_types = FALSE
)

expected_columns <- c(
  "huc", "month", "metabolism_metric", "n_reaches",
  "log10_area_range", "median_pct_area_missing", "intercept", "slope",
  "r2", "slope_class_simple"
)

duplicate_count <- raw |>
  count(huc, month, metabolism_metric) |>
  filter(n > 1) |>
  nrow()

verification <- tibble(
  rows = nrow(raw),
  hucs = n_distinct(raw$huc),
  months = n_distinct(raw$month),
  gpp_records = sum(raw$metabolism_metric == "GPP"),
  er_records = sum(raw$metabolism_metric == "ER"),
  duplicates = duplicate_count
)

stopifnot(
  identical(names(raw), expected_columns),
  nrow(raw) == 44860,
  n_distinct(raw$huc) == 1871,
  all(sort(unique(raw$month)) == 1:12),
  sum(raw$metabolism_metric == "GPP") == 22452,
  sum(raw$metabolism_metric == "ER") == 22408,
  duplicate_count == 0,
  !anyNA(raw)
)


# 3. Build paired HUC8-month data -----------------------------------------

paired <- raw |>
  select(
    huc, month, metabolism_metric, slope, r2, n_reaches,
    log10_area_range, median_pct_area_missing, intercept
  ) |>
  pivot_wider(
    names_from = metabolism_metric,
    values_from = c(
      slope, r2, n_reaches, log10_area_range,
      median_pct_area_missing, intercept
    ),
    names_glue = "{.value}_{tolower(metabolism_metric)}"
  ) |>
  filter(!is.na(slope_gpp), !is.na(slope_er)) |>
  rename(beta_gpp = slope_gpp, beta_er = slope_er) |>
  mutate(
    delta_beta = beta_gpp - beta_er,
    ratio_scaling_10x = 10^delta_beta,
    month_factor = factor(month, levels = 1:12, labels = month_labels)
  )

missing_pairs <- tibble(
  paired_huc_months = nrow(paired),
  gpp_without_er = verification$gpp_records - nrow(paired),
  er_without_gpp = verification$er_records - nrow(paired)
)

stopifnot(
  nrow(paired) == 22408,
  missing_pairs$gpp_without_er == 44,
  missing_pairs$er_without_gpp == 0
)

print(verification)
print(missing_pairs)


# 4. Summer-winter primary comparison ------------------------------------

seasonal <- paired |>
  filter(month %in% c(12, 1, 2, 6, 7, 8)) |>
  mutate(season = if_else(month %in% 6:8, "summer", "winter")) |>
  group_by(huc, season) |>
  summarise(
    n_months = n(),
    delta_beta = mean(delta_beta),
    r2_support = min(pmin(r2_gpp, r2_er)),
    reaches_support = min(pmin(n_reaches_gpp, n_reaches_er)),
    area_support = min(pmin(log10_area_range_gpp, log10_area_range_er)),
    .groups = "drop"
  ) |>
  filter(n_months == 3) |>
  select(-n_months) |>
  pivot_wider(
    names_from = season,
    values_from = c(delta_beta, r2_support, reaches_support, area_support)
  ) |>
  drop_na() |>
  transmute(
    huc,
    summer_delta_beta = delta_beta_summer,
    winter_delta_beta = delta_beta_winter,
    seasonal_shift = summer_delta_beta - winter_delta_beta,
    min_r2 = pmin(r2_support_summer, r2_support_winter),
    min_reaches = pmin(reaches_support_summer, reaches_support_winter),
    min_area_range = pmin(area_support_summer, area_support_winter)
  )

seasonal_stats <- seasonal |>
  summarise(
    n_hucs = n(),
    median_summer_delta_beta = median(summer_delta_beta),
    median_winter_delta_beta = median(winter_delta_beta),
    median_seasonal_shift = median(seasonal_shift),
    q25_seasonal_shift = quantile(seasonal_shift, 0.25),
    q75_seasonal_shift = quantile(seasonal_shift, 0.75),
    fraction_positive_shift = mean(seasonal_shift > 0)
  )

p1 <- ggplot(seasonal, aes(seasonal_shift)) +
  geom_histogram(bins = 70, fill = "#3B8BC2", color = "white", linewidth = 0.15) +
  geom_rug(alpha = 0.08, linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey25", linewidth = 0.7) +
  geom_vline(
    xintercept = seasonal_stats$median_seasonal_shift,
    color = "#D1495B", linewidth = 1
  ) +
  annotate(
    "label", x = Inf, y = Inf, hjust = 1.03, vjust = 1.15,
    label = sprintf(
      "Median shift = %.3f\nPositive shift = %.1f%%\nn = %s HUC8s",
      seasonal_stats$median_seasonal_shift,
      100 * seasonal_stats$fraction_positive_shift,
      scales::comma(seasonal_stats$n_hucs)
    ),
    size = 3.6
  ) +
  labs(
    x = expression(paste("Summer ", Delta, beta, " − winter ", Delta, beta)),
    y = "Number of HUC8s"
  )

ggsave(
  file.path(OUTPUT_DIR, "01_summer_winter_delta_beta.png"),
  p1, width = 8.5, height = 5.5, dpi = 320
)


# 5. Full annual cycle -----------------------------------------------------

p2 <- ggplot(paired, aes(month_factor, delta_beta)) +
  geom_violin(fill = "#6BAED6", color = "#24618A", linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  geom_hline(yintercept = 0, color = "grey25", linewidth = 0.6) +
  labs(x = NULL, y = expression(Delta * beta == beta[GPP] - beta[ER]))

ggsave(
  file.path(OUTPUT_DIR, "02_monthly_delta_beta.png"),
  p2, width = 10, height = 5.8, dpi = 320
)


# 6. GPP and ER slopes -----------------------------------------------------

slope_long <- paired |>
  select(huc, month_factor, beta_gpp, beta_er) |>
  pivot_longer(
    c(beta_gpp, beta_er), names_to = "process", values_to = "slope"
  ) |>
  mutate(process = recode(process, beta_gpp = "GPP", beta_er = "ER"))

p3 <- ggplot(slope_long, aes(month_factor, slope, fill = process)) +
  geom_violin(color = NA, alpha = 0.78, scale = "width") +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  facet_wrap(vars(process), ncol = 1) +
  scale_fill_manual(values = c(GPP = "#1B9E77", ER = "#7570B3"), guide = "none") +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.5) +
  labs(x = NULL, y = "Allometric slope")

ggsave(
  file.path(OUTPUT_DIR, "03_monthly_gpp_er_slopes.png"),
  p3, width = 10, height = 8, dpi = 320
)


# 7. Seasonal amplitudes ---------------------------------------------------

amplitudes <- paired |>
  group_by(huc) |>
  filter(n() == 12) |>
  summarise(
    gpp_amplitude = max(beta_gpp) - min(beta_gpp),
    er_amplitude = max(beta_er) - min(beta_er),
    amplitude_difference = gpp_amplitude - er_amplitude,
    .groups = "drop"
  )

amplitude_stats <- amplitudes |>
  summarise(
    n_hucs = n(),
    median_gpp_amplitude = median(gpp_amplitude),
    median_er_amplitude = median(er_amplitude),
    median_paired_difference = median(amplitude_difference),
    fraction_gpp_greater = mean(gpp_amplitude > er_amplitude)
  )

amp_limit <- quantile(c(amplitudes$gpp_amplitude, amplitudes$er_amplitude), 0.995)

p4 <- ggplot(amplitudes, aes(er_amplitude, gpp_amplitude)) +
  geom_abline(slope = 1, intercept = 0, color = "grey25", linewidth = 0.7) +
  geom_point(alpha = 0.22, size = 1.3, color = "#2C7FB8") +
  coord_equal(xlim = c(0, amp_limit), ylim = c(0, amp_limit)) +
  annotate(
    "label", x = Inf, y = -Inf, hjust = 1.03, vjust = -0.15,
    label = sprintf(
      "Median GPP = %.3f\nMedian ER = %.3f\nMedian paired difference = %.3f\nGPP > ER = %.1f%%",
      amplitude_stats$median_gpp_amplitude,
      amplitude_stats$median_er_amplitude,
      amplitude_stats$median_paired_difference,
      100 * amplitude_stats$fraction_gpp_greater
    ),
    size = 3.3
  ) +
  labs(x = "ER seasonal amplitude", y = "GPP seasonal amplitude")

ggsave(
  file.path(OUTPUT_DIR, "04_gpp_er_seasonal_amplitude.png"),
  p4, width = 7.5, height = 7, dpi = 320
)


# 8. Relative scaling effect ----------------------------------------------

p5 <- ggplot(paired, aes(month_factor, ratio_scaling_10x)) +
  geom_violin(fill = "#66C2A4", color = "#238B45", linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  geom_hline(yintercept = 1, color = "grey25", linewidth = 0.6) +
  coord_cartesian(ylim = quantile(paired$ratio_scaling_10x, c(0.01, 0.99))) +
  labs(
    x = NULL,
    y = expression("Relative scaling across a 10× size increase  " * (10^Delta*beta))
  )

ggsave(
  file.path(OUTPUT_DIR, "05_monthly_relative_scaling_effect.png"),
  p5, width = 10, height = 5.8, dpi = 320
)


# 9. Fit-quality diagnostic ------------------------------------------------

monthly_qc <- paired |>
  transmute(
    scope = "Monthly |delta beta|",
    response = abs(delta_beta),
    `Minimum R²` = pmin(r2_gpp, r2_er),
    `Minimum reaches` = pmin(n_reaches_gpp, n_reaches_er),
    `Minimum log10 area range` = pmin(log10_area_range_gpp, log10_area_range_er)
  )

seasonal_qc <- seasonal |>
  transmute(
    scope = "Seasonal |shift|",
    response = abs(seasonal_shift),
    `Minimum R²` = min_r2,
    `Minimum reaches` = min_reaches,
    `Minimum log10 area range` = min_area_range
  )

qc_long <- bind_rows(monthly_qc, seasonal_qc) |>
  pivot_longer(-c(scope, response), names_to = "support_metric", values_to = "support") |>
  group_by(scope, support_metric) |>
  mutate(support_percentile = percent_rank(support)) |>
  ungroup()

qc_correlations <- qc_long |>
  group_by(scope, support_metric) |>
  summarise(
    spearman_rho = cor(response, support, method = "spearman"),
    median_response_lowest_support_quartile = median(response[support_percentile <= 0.25]),
    median_response_other_three_quartiles = median(response[support_percentile > 0.25]),
    .groups = "drop"
  )

seasonal_quality_sensitivity <- seasonal |>
  select(huc, seasonal_shift, min_r2, min_reaches, min_area_range) |>
  pivot_longer(
    c(min_r2, min_reaches, min_area_range),
    names_to = "support_metric", values_to = "support"
  ) |>
  group_by(support_metric) |>
  mutate(support_quartile = ntile(support, 4)) |>
  group_by(support_metric, support_quartile) |>
  summarise(
    n_hucs = n(),
    median_seasonal_shift = median(seasonal_shift),
    fraction_positive_shift = mean(seasonal_shift > 0),
    .groups = "drop"
  )

p6 <- ggplot(qc_long, aes(support_percentile, response)) +
  geom_point(alpha = 0.045, size = 0.7, color = "#2C7FB8") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#D1495B", linewidth = 0.8) +
  facet_grid(vars(scope), vars(support_metric), scales = "free_y") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Fit-support percentile", y = "Absolute effect size")

ggsave(
  file.path(OUTPUT_DIR, "06_quality_diagnostic.png"),
  p6, width = 12, height = 7, dpi = 320
)


# 10. Numerical summaries -------------------------------------------------

monthly_summary <- paired |>
  group_by(month) |>
  summarise(
    n_paired_hucs = n(),
    median_beta_gpp = median(beta_gpp),
    q25_beta_gpp = quantile(beta_gpp, 0.25),
    q75_beta_gpp = quantile(beta_gpp, 0.75),
    median_beta_er = median(beta_er),
    q25_beta_er = quantile(beta_er, 0.25),
    q75_beta_er = quantile(beta_er, 0.75),
    median_delta_beta = median(delta_beta),
    q25_delta_beta = quantile(delta_beta, 0.25),
    q75_delta_beta = quantile(delta_beta, 0.75),
    fraction_delta_beta_positive = mean(delta_beta > 0),
    median_ratio_scaling_10x = median(ratio_scaling_10x),
    median_r2_gpp = median(r2_gpp),
    median_r2_er = median(r2_er),
    median_n_reaches = median(c(n_reaches_gpp, n_reaches_er)),
    .groups = "drop"
  )

readr::write_csv(
  monthly_summary,
  file.path(OUTPUT_DIR, "initial_headline_diagnostics.csv")
)

cat("\nSeasonal comparison:\n")
print(seasonal_stats)
cat("\nAmplitude comparison:\n")
print(amplitude_stats)
cat("\nQuality associations:\n")
print(qc_correlations, n = Inf)
cat("\nSigned seasonal result by fit-support quartile (4 = strongest support):\n")
print(seasonal_quality_sensitivity, n = Inf)
cat("\nMonthly summary:\n")
print(monthly_summary, n = 12, width = Inf)
