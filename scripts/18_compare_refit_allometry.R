## 18_compare_refit_allometry.R
##
## PURPOSE
## Compare HUC8 allometric scaling produced by the paired Maavara random
## forests fit in script 17: one with upstream watershed area (`uparea`) and
## one without it.
##
## The routing, cumulative-flux calculation, and eligibility rules match the
## canonical workflow in 07_calculate_huc8_allometric_scaling.Rmd. Each fit is:
##
##   log10(cumulative metabolism) ~ log10(cumulative drainage area)
##
## for one HUC8 x month x metabolism process x RF variant.


# 0. Setup --------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")
pacman::p_load(Matrix)

tic("run 18")

dir_processed <- "data/processed"
dir_merit <- file.path(dir_processed, "merit_data")
dir_predictions <- file.path(dir_processed, "17_maavara_refits")
dir_figures <- "figures"

dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

path_catchments <- file.path(
  dir_merit,
  "03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
)
path_flowlines <- file.path(
  dir_merit,
  "03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
)
path_upstream <- file.path(dir_merit, "05_upstream_lookup.rds")
path_characteristics <-
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv"

path_results <- file.path(dir_processed, "18_refit_allometry_results.csv")
path_pairs <- file.path(dir_processed, "18_refit_allometry_paired.csv")
path_summary <- file.path(dir_processed, "18_refit_allometry_summary.csv")
path_paired_summary <-
  file.path(dir_processed, "18_refit_allometry_paired_summary.csv")
path_scatter_input <-
  file.path(dir_processed, "18_refit_no_uparea_scatter_input.rds")

# These are the same screening values used in script 07.
missing_area_threshold <- 10
min_reaches <- 20
min_area_range_log10 <- 1

# Terminal modeled reaches used as named reference points in Figure 29.
featured_mouth_reaches <- tibble(
  river = c("Mississippi", "Ohio", "Columbia", "Rio Grande"),
  comid = c("74072169", "74050625", "78018922", "75025589")
)

month_table <- tibble(
  month = 1:12,
  month_abbr = month.abb,
  month_name = month.name
)

scenario_table <- tidyr::crossing(
  model = c("with_uparea", "without_uparea"),
  metabolism = c("GPP", "ER"),
  month_table
) %>%
  mutate(
    model_label = recode(
      model,
      with_uparea = "With watershed area",
      without_uparea = "Without watershed area"
    ),
    prediction_path = file.path(
      dir_predictions,
      model,
      metabolism,
      paste0("scaling_", metabolism, "_pred_", month_abbr, ".csv")
    )
  )

stopifnot(all(file.exists(scenario_table$prediction_path)))


# 1. Read the durable routing and reach resources ------------------------------

message("1/6 Reading reach attributes and routing lookup...")

catchments <- read_sf(path_catchments, quiet = TRUE) %>%
  clean_names() %>%
  st_drop_geometry() %>%
  transmute(
    comid = as.character(comid),
    huc = as.character(huc),
    unitarea = as.numeric(unitarea)
  )

flowlines <- read_sf(path_flowlines, quiet = TRUE) %>%
  clean_names() %>%
  st_drop_geometry() %>%
  transmute(
    comid = as.character(comid),
    uparea = as.numeric(uparea)
  )

characteristics <- read_csv(path_characteristics, show_col_types = FALSE) %>%
  clean_names() %>%
  transmute(
    comid = as.character(comid),
    lengthkm = as.numeric(lengthkm),
    merit_width = as.numeric(merit_width)
  )

reach_index <- catchments %>%
  left_join(flowlines, by = "comid") %>%
  left_join(characteristics, by = "comid") %>%
  mutate(
    idx = row_number(),
    stream_area_m2 = lengthkm * 1000 * merit_width
  )

stopifnot(
  !anyDuplicated(reach_index$comid),
  !anyNA(reach_index$unitarea),
  !anyNA(reach_index$stream_area_m2)
)

upstream_lookup <- readRDS(path_upstream) %>%
  transmute(
    downstream_comid = as.character(downstream_comid),
    upstream_comid = as.character(upstream_comid)
  )


# 2. Build one sparse routing matrix ------------------------------------------

message("2/6 Building sparse routing matrix...")

routing_lookup <- upstream_lookup %>%
  inner_join(
    reach_index %>%
      select(upstream_comid = comid, upstream_idx = idx),
    by = "upstream_comid"
  ) %>%
  inner_join(
    reach_index %>%
      select(downstream_comid = comid, downstream_idx = idx),
    by = "downstream_comid"
  )

routing_matrix <- sparseMatrix(
  i = routing_lookup$downstream_idx,
  j = routing_lookup$upstream_idx,
  x = 1,
  dims = rep(nrow(reach_index), 2)
)

if (length(routing_matrix@x) != nrow(routing_lookup)) {
  stop("The upstream lookup contains duplicated downstream-upstream pairs.")
}

rm(upstream_lookup, routing_lookup)
invisible(gc())


# 3. Calculate cumulative drainage area ---------------------------------------

message("3/6 Calculating cumulative drainage area...")

cumulative_area_km2 <- as.numeric(
  routing_matrix %*% reach_index$unitarea
)

routing_base <- reach_index %>%
  transmute(
    comid,
    huc,
    uparea,
    cumulative_area_km2,
    area_ratio_retained = cumulative_area_km2 / uparea,
    pct_area_missing = pmax(0, 1 - area_ratio_retained) * 100
  )


# 4. Read paired predictions and route all scenarios together -----------------

message("4/6 Reading and routing all 48 prediction scenarios...")

n_reaches <- nrow(reach_index)
n_scenarios <- nrow(scenario_table)

local_flux <- matrix(
  0,
  nrow = n_reaches,
  ncol = n_scenarios,
  dimnames = list(NULL, paste0("scenario_", seq_len(n_scenarios)))
)

scenario_table$n_missing_rate <- NA_integer_
scenario_table$n_matched_comid <- NA_integer_

for (scenario_id in seq_len(n_scenarios)) {
  prediction <- read_csv(
    scenario_table$prediction_path[scenario_id],
    show_col_types = FALSE
  )

  names(prediction)[1] <- "comid"
  names(prediction)[2] <- "rate"

  prediction <- prediction %>%
    transmute(
      comid = as.character(comid),
      rate = as.numeric(rate)
    )

  prediction_row <- match(reach_index$comid, prediction$comid)
  rate <- prediction$rate[prediction_row]
  flux <- rate * reach_index$stream_area_m2

  scenario_table$n_matched_comid[scenario_id] <- sum(!is.na(prediction_row))
  scenario_table$n_missing_rate[scenario_id] <- sum(!is.finite(rate))

  # This matches sum(..., na.rm = TRUE) in script 07. Missing local predictions
  # contribute no flux. The same five archived covariate gaps occur in both
  # paired RF variants, so missingness is held constant across the comparison.
  flux[!is.finite(flux)] <- 0
  local_flux[, scenario_id] <- flux
}

cumulative_flux <- as.matrix(routing_matrix %*% local_flux)

# Save a compact, reproducible subset for Figure 1's pooled scatterplots.
# The complete routed matrix is intentionally not written because it is large.
set.seed(18)
n_scatter_reaches_per_month <- 20000

scatter_input <- map_dfr(
  which(scenario_table$model == "without_uparea"),
  function(scenario_id) {
    scenario <- scenario_table[scenario_id, ]
    paired_scenario_ids <- which(
      scenario_table$metabolism == scenario$metabolism &
        scenario_table$month == scenario$month
    )

    eligible_reaches <- routing_base %>%
      mutate(
        cumulative_flux_g_o2_d = cumulative_flux[, scenario_id],
        paired_flux_is_positive = rowSums(
          is.finite(cumulative_flux[, paired_scenario_ids, drop = FALSE]) &
            cumulative_flux[, paired_scenario_ids, drop = FALSE] > 0
        ) == length(paired_scenario_ids)
      ) %>%
      filter(
        paired_flux_is_positive,
        !is.na(huc),
        cumulative_area_km2 > 0,
        cumulative_flux_g_o2_d > 0,
        pct_area_missing < missing_area_threshold
      ) %>%
      transmute(
        comid,
        month = scenario$month,
        metabolism = scenario$metabolism,
        cumulative_area_km2,
        cumulative_flux_g_o2_d
      )

    sampled_reaches <- slice_sample(
      eligible_reaches,
      n = min(n_scatter_reaches_per_month, nrow(eligible_reaches))
    )

    featured_reaches <- eligible_reaches %>%
      semi_join(featured_mouth_reaches, by = "comid")

    bind_rows(sampled_reaches, featured_reaches) %>%
      left_join(featured_mouth_reaches, by = "comid")
  }
)

saveRDS(scatter_input, path_scatter_input)

rm(local_flux, routing_matrix)
invisible(gc())


# 5. Fit HUC8 x month allometric relationships --------------------------------

message("5/6 Fitting HUC8 x month allometric relationships...")

fit_one_scenario <- function(scenario_id) {
  scenario <- scenario_table[scenario_id, ]

  paired_scenario_ids <- which(
    scenario_table$metabolism == scenario$metabolism &
      scenario_table$month == scenario$month
  )

  # Require the same usable reaches in both RF variants. A few cumulative
  # predictions cross zero between variants; allowing variant-specific filters
  # would confound model differences with a change in the regression sample.
  paired_flux <- cumulative_flux[, paired_scenario_ids, drop = FALSE]
  paired_flux_is_positive <- rowSums(
    is.finite(paired_flux) & paired_flux > 0
  ) == length(paired_scenario_ids)

  routing_base %>%
    mutate(
      cumulative_flux_g_o2_d = cumulative_flux[, scenario_id],
      paired_flux_is_positive = paired_flux_is_positive
    ) %>%
    filter(
      paired_flux_is_positive,
      !is.na(huc),
      is.finite(cumulative_area_km2),
      is.finite(cumulative_flux_g_o2_d),
      cumulative_area_km2 > 0,
      pct_area_missing < missing_area_threshold
    ) %>%
    mutate(
      log_area = log10(cumulative_area_km2),
      log_flux = log10(cumulative_flux_g_o2_d)
    ) %>%
    group_by(huc) %>%
    summarise(
      n_reaches = n(),
      n_distinct_area = n_distinct(cumulative_area_km2),
      log10_area_range = max(log_area) - min(log_area),
      median_pct_area_missing = median(pct_area_missing),
      intercept = mean(log_flux) - cov(log_area, log_flux) /
        var(log_area) * mean(log_area),
      slope = cov(log_area, log_flux) / var(log_area),
      r2 = cor(log_area, log_flux)^2,
      .groups = "drop"
    ) %>%
    filter(
      n_reaches >= min_reaches,
      n_distinct_area >= min_reaches,
      log10_area_range >= min_area_range_log10
    ) %>%
    mutate(
      month = scenario$month,
      month_abbr = scenario$month_abbr,
      metabolism = scenario$metabolism,
      model = scenario$model,
      model_label = scenario$model_label,
      .before = 1
    )
}

allometry_results <- map_dfr(
  seq_len(n_scenarios),
  fit_one_scenario,
  .progress = TRUE
) %>%
  mutate(
    model_label = factor(
      model_label,
      levels = c("With watershed area", "Without watershed area")
    ),
    metabolism = factor(metabolism, levels = c("GPP", "ER")),
    slope_class = case_when(
      r2 < 0.8 ~ "Low R²",
      slope < 0.95 ~ "Sublinear",
      slope <= 1.05 ~ "Near-linear",
      slope > 1.05 ~ "Superlinear",
      TRUE ~ NA_character_
    )
  )

paired_results <- allometry_results %>%
  select(huc, month, month_abbr, metabolism, model, slope, r2, n_reaches) %>%
  pivot_wider(
    names_from = model,
    values_from = c(slope, r2, n_reaches)
  ) %>%
  filter(!is.na(slope_with_uparea), !is.na(slope_without_uparea)) %>%
  mutate(
    delta_slope = slope_without_uparea - slope_with_uparea,
    delta_r2 = r2_without_uparea - r2_with_uparea,
    abs_delta_slope = abs(delta_slope),
    abs_delta_r2 = abs(delta_r2)
  )

if (any(paired_results$n_reaches_with_uparea !=
        paired_results$n_reaches_without_uparea)) {
  stop("Paired fits used different numbers of reaches.")
}

distribution_summary <- allometry_results %>%
  group_by(metabolism, model, model_label) %>%
  summarise(
    n_fits = n(),
    n_hucs = n_distinct(huc),
    slope_median = median(slope),
    slope_q025 = quantile(slope, 0.025),
    slope_q25 = quantile(slope, 0.25),
    slope_q75 = quantile(slope, 0.75),
    slope_q975 = quantile(slope, 0.975),
    r2_median = median(r2),
    r2_q025 = quantile(r2, 0.025),
    r2_q25 = quantile(r2, 0.25),
    r2_q75 = quantile(r2, 0.75),
    r2_q975 = quantile(r2, 0.975),
    pct_r2_ge_08 = mean(r2 >= 0.8) * 100,
    .groups = "drop"
  )

paired_summary <- paired_results %>%
  group_by(metabolism) %>%
  summarise(
    n_pairs = n(),
    n_hucs = n_distinct(huc),
    slope_correlation = cor(slope_with_uparea, slope_without_uparea),
    median_delta_slope = median(delta_slope),
    q025_delta_slope = quantile(delta_slope, 0.025),
    q975_delta_slope = quantile(delta_slope, 0.975),
    median_abs_delta_slope = median(abs_delta_slope),
    pct_slope_increases_without_area = mean(delta_slope > 0) * 100,
    r2_correlation = cor(r2_with_uparea, r2_without_uparea),
    median_delta_r2 = median(delta_r2),
    q025_delta_r2 = quantile(delta_r2, 0.025),
    q975_delta_r2 = quantile(delta_r2, 0.975),
    median_abs_delta_r2 = median(abs_delta_r2),
    pct_r2_increases_without_area = mean(delta_r2 > 0) * 100,
    .groups = "drop"
  )

write_csv(allometry_results, path_results)
write_csv(paired_results, path_pairs)
write_csv(distribution_summary, path_summary)
write_csv(paired_summary, path_paired_summary)


# 6. Visualize distributions and paired changes -------------------------------

message("6/6 Drawing comparison figures...")

model_colors <- c(
  "With watershed area" = "#0072B2",
  "Without watershed area" = "#D55E00"
)

metric_long <- allometry_results %>%
  select(huc, month, metabolism, model_label, slope, r2) %>%
  pivot_longer(
    c(slope, r2),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      slope = "Allometric slope",
      r2 = "R²"
    ),
    metric = factor(metric, levels = c("Allometric slope", "R²"))
  )

reference_lines <- tibble(
  metric = factor(
    c("Allometric slope", "R²"),
    levels = levels(metric_long$metric)
  ),
  reference = c(1, 0.8)
)

p_distributions <- ggplot(
  metric_long,
  aes(x = value, color = model_label, fill = model_label)
) +
  geom_density(alpha = 0.16, linewidth = 0.8, adjust = 1.1) +
  geom_vline(
    data = reference_lines,
    aes(xintercept = reference),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "grey35"
  ) +
  facet_wrap(
    vars(metric, metabolism),
    ncol = 2,
    scales = "free"
  ) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  labs(
    x = NULL,
    y = "Density across HUC8-month fits",
    color = NULL,
    fill = NULL,
    title = "HUC8 allometry with and without watershed area in the RF"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(dir_figures, "18_refit_allometry_distributions.png"),
  p_distributions,
  width = 9,
  height = 7,
  dpi = 300
)

paired_long <- bind_rows(
  paired_results %>%
    transmute(
      huc, month, metabolism,
      metric = "Allometric slope",
      with_area = slope_with_uparea,
      without_area = slope_without_uparea
    ),
  paired_results %>%
    transmute(
      huc, month, metabolism,
      metric = "R²",
      with_area = r2_with_uparea,
      without_area = r2_without_uparea
    )
) %>%
  mutate(
    metric = factor(metric, levels = c("Allometric slope", "R²"))
  )

p_paired <- ggplot(
  paired_long,
  aes(x = with_area, y = without_area)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_hex(bins = 45) +
  scale_fill_viridis_c(option = "magma", trans = "sqrt") +
  facet_wrap(
    vars(metric, metabolism),
    ncol = 2,
    scales = "free"
  ) +
  labs(
    x = "With watershed area",
    y = "Without watershed area",
    fill = "HUC8-month\nfits",
    title = "Paired changes in HUC8 allometric fits"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  file.path(dir_figures, "18_refit_allometry_paired.png"),
  p_paired,
  width = 9,
  height = 7,
  dpi = 300
)

p_slope_month <- ggplot(
  allometry_results,
  aes(x = factor(month), y = slope, fill = model_label)
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35") +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.72,
    position = position_dodge(width = 0.78)
  ) +
  facet_wrap(~ metabolism, ncol = 1) +
  scale_fill_manual(values = model_colors) +
  coord_cartesian(ylim = quantile(allometry_results$slope, c(0.005, 0.995))) +
  labs(
    x = "Month",
    y = "Allometric slope",
    fill = NULL,
    title = "Seasonal distribution of HUC8 slopes"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(dir_figures, "18_refit_allometry_slopes_by_month.png"),
  p_slope_month,
  width = 10,
  height = 7,
  dpi = 300
)

print(distribution_summary, n = Inf)
print(paired_summary, n = Inf)

toc()
