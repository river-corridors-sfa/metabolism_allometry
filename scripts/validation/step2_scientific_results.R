# Read-only comparison of saved refit metrics, 07 fits, and 08-derived 10/11 classes.
# No model fitting, routing, or figure rendering.
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(sf)

out <- "reports/no_uparea_migration/step2_scientific_results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
emit <- function(x, name) write_csv(x, file.path(out, paste0(name, ".csv")), na = "")

# Script 01a uses the same training rows, response, 70/30 split, and seed 24
# for both refits. Test record IDs/predictions were not saved.
metrics <- read_csv("data/processed/17_maavara_refits/model_metrics.csv", show_col_types = FALSE) %>%
  filter(dataset == "test") %>%
  select(process, variant, n, r2_score, rmse, mae) %>%
  pivot_wider(names_from = variant, values_from = c(n, r2_score, rmse, mae)) %>%
  mutate(same_test_n = n_with_uparea == n_without_uparea,
         delta_r2 = r2_score_without_uparea - r2_score_with_uparea,
         delta_rmse = rmse_without_uparea - rmse_with_uparea,
         relative_rmse_pct = 100 * delta_rmse / rmse_with_uparea,
         delta_mae = mae_without_uparea - mae_with_uparea,
         relative_mae_pct = 100 * delta_mae / mae_with_uparea)
stopifnot(nrow(metrics) == 2, all(metrics$same_test_n))
emit(metrics, "paired_refit_test_skill")

fit_paths <- c(original = "data/processed/07_huc08_allometry_results.csv",
               no_uparea = "data/processed/07_no_uparea_huc08_allometry_results.csv")
fits <- bind_rows(lapply(names(fit_paths), function(variant) {
  read_csv(fit_paths[[variant]], show_col_types = FALSE,
           col_types = cols(huc = col_character(), .default = col_guess())) %>%
    mutate(variant = variant)
}))
stopifnot(!anyDuplicated(fits[c("huc", "month", "metabolism_metric", "variant")]))
fit_counts <- fits %>% count(variant, metabolism_metric, name = "eligible_fits")
emit(fit_counts, "eligible_fit_counts")
paired_fits <- fits %>%
  select(huc, month, metabolism_metric, variant, slope, r2, n_reaches) %>%
  pivot_wider(names_from = variant, values_from = c(slope, r2, n_reaches)) %>%
  mutate(composition = case_when(is.na(slope_original) ~ "no_uparea_only",
                                 is.na(slope_no_uparea) ~ "original_only",
                                 TRUE ~ "matched"),
         delta_slope = slope_no_uparea - slope_original,
         delta_r2 = r2_no_uparea - r2_original,
         delta_reaches = n_reaches_no_uparea - n_reaches_original)
support_columns <- c("huc", "month", "metabolism_metric", "downstream_comid")
support_old <- readRDS("data/processed/07_huc08_allometry_input.rds") %>%
  select(all_of(support_columns))
support_new <- readRDS("data/processed/07_no_uparea_huc08_allometry_input.rds") %>%
  select(all_of(support_columns))
if (anyDuplicated(support_old) || anyDuplicated(support_new)) stop("Duplicate 07 support rows")
support_changes <- full_join(
  anti_join(support_old, support_new, by = support_columns) %>%
    count(huc, month, metabolism_metric, name = "rows_lost"),
  anti_join(support_new, support_old, by = support_columns) %>%
    count(huc, month, metabolism_metric, name = "rows_gained"),
  by = c("huc", "month", "metabolism_metric")) %>%
  mutate(rows_lost = replace_na(rows_lost, 0L), rows_gained = replace_na(rows_gained, 0L))
paired_fits <- paired_fits %>%
  left_join(support_changes, by = c("huc", "month", "metabolism_metric")) %>%
  mutate(rows_lost = replace_na(rows_lost, 0L), rows_gained = replace_na(rows_gained, 0L),
         sample_changed = rows_lost + rows_gained > 0)
emit(paired_fits, "paired_fits")
emit(paired_fits %>% count(metabolism_metric, composition), "fit_composition")
emit(paired_fits %>% filter(composition == "matched") %>% group_by(metabolism_metric) %>%
  summarise(n = n(), median_slope_original = median(slope_original),
            median_slope_no_uparea = median(slope_no_uparea),
            median_delta_slope = median(delta_slope),
            p05_delta_slope = quantile(delta_slope, .05),
            p95_delta_slope = quantile(delta_slope, .95),
            median_r2_original = median(r2_original),
            median_r2_no_uparea = median(r2_no_uparea),
            median_delta_r2 = median(delta_r2),
            n_sample_changed = sum(sample_changed),
            .groups = "drop"), "matched_fit_summary")
emit(paired_fits %>% filter(composition == "matched") %>%
  group_by(metabolism_metric, sample_changed) %>%
  summarise(n = n(), median_delta_slope = median(delta_slope),
            median_delta_r2 = median(delta_r2), .groups = "drop"),
  "fit_changes_by_sample_composition")

# Exact 10/11 classification logic, evaluated on both saved 08 files.
classify <- function(path, variant, boundaries) {
  stability <- read_csv(path, show_col_types = FALSE,
                        col_types = cols(huc = col_character(), .default = col_guess())) %>%
    mutate(huc = str_pad(huc, 8, pad = "0"))
  required <- c("huc", "month", "metabolism_metric", "n_reaches", "n_unique_areas",
                "log10_area_range", "n_valid", "baseline_slope", "baseline_r2",
                "slope_min", "slope_max", "r2_min",
                "median_pct_metabolism_missing", "max_pct_metabolism_missing")
  stopifnot(all(required %in% names(stability)),
            !anyDuplicated(stability[c("huc", "month", "metabolism_metric")]),
            all(stability$max_pct_metabolism_missing < 10))
  monthly <- stability %>%
    mutate(supported = if_all(all_of(setdiff(required, c("huc", "month", "metabolism_metric"))),
                              ~ is.finite(.x)) &
             n_reaches >= 20 & n_unique_areas >= 20 & log10_area_range >= 1 & n_valid > 0,
           slope_low = pmin(baseline_slope, slope_min),
           slope_high = pmax(baseline_slope, slope_max),
           weakest_r2 = pmin(baseline_r2, r2_min),
           monthly_class = case_when(
             !supported ~ "Insufficient data",
             weakest_r2 < .8 ~ "Uncertain",
             slope_high < .95 ~ "Sub-linear",
             slope_low >= .95 & slope_high <= 1.05 ~ "Linear",
             slope_low > 1.05 ~ "Super-linear",
             TRUE ~ "Uncertain"),
           available = TRUE) %>%
    select(huc, month, metabolism_metric, monthly_class, available)
  annual <- expand_grid(huc = boundaries, metabolism_metric = c("GPP", "ER"), month = 1:12) %>%
    left_join(monthly, by = c("huc", "metabolism_metric", "month")) %>%
    group_by(huc, metabolism_metric) %>%
    summarise(scaling_class = case_when(
      sum(!is.na(available)) < 12 ~ "Insufficient data",
      n_distinct(monthly_class) > 1 ~ "Seasonally variable",
      TRUE ~ first(monthly_class)), .groups = "drop")
  list(monthly = mutate(monthly, variant = variant),
       annual = mutate(annual, variant = variant))
}
boundaries <- read_sf("data/processed/04c_huc08_boundaries.gpkg") %>%
  st_drop_geometry() %>% transmute(huc = str_pad(as.character(huc8), 8, pad = "0")) %>%
  pull(huc)
stopifnot(!anyDuplicated(boundaries))
classified <- list(
  classify("data/processed/08_huc08_fit_stability.csv", "original", boundaries),
  classify("data/processed/08_no_uparea_huc08_fit_stability.csv", "no_uparea", boundaries))
monthly <- bind_rows(lapply(classified, `[[`, "monthly"))
annual <- bind_rows(lapply(classified, `[[`, "annual"))
emit(monthly, "monthly_classes_from_08_using_10_11_rules")
emit(annual, "annual_classes_from_08_using_10_11_rules")
emit(monthly %>% count(variant, metabolism_metric, monthly_class), "monthly_class_counts")
emit(annual %>% count(variant, metabolism_metric, scaling_class), "annual_class_counts")

monthly_pairs <- monthly %>%
  select(huc, month, metabolism_metric, variant, monthly_class) %>%
  pivot_wider(names_from = variant, values_from = monthly_class) %>%
  mutate(changed = is.na(original) | is.na(no_uparea) | original != no_uparea)
annual_pairs <- annual %>%
  select(huc, metabolism_metric, variant, scaling_class) %>%
  pivot_wider(names_from = variant, values_from = scaling_class) %>%
  mutate(changed = original != no_uparea)
emit(monthly_pairs, "monthly_class_pairs")
emit(annual_pairs, "annual_class_pairs")
emit(monthly_pairs %>% count(metabolism_metric, original, no_uparea, changed),
     "monthly_class_transitions")
emit(annual_pairs %>% count(metabolism_metric, original, no_uparea, changed),
     "annual_class_transitions")
emit(annual_pairs %>% filter(changed), "changed_annual_huc8s")

verify_saved <- function(stage, expected) {
  keys <- if (stage == "monthly") c("huc", "month", "metabolism_metric") else c("huc", "metabolism_metric")
  value <- if (stage == "monthly") "monthly_class" else "scaling_class"
  bind_rows(lapply(c("10", "11"), function(script) {
    path <- file.path("data/processed", paste0(script, "_no_uparea_", stage, "_classes.csv"))
    if (!file.exists(path)) return(tibble(script, stage, file_exists = FALSE,
                                          n_rows = NA_integer_, key_mismatches = NA_integer_,
                                          class_mismatches = NA_integer_))
    saved <- read_csv(path, show_col_types = FALSE,
                      col_types = cols(huc = col_character(), .default = col_guess()))
    compared <- full_join(saved, expected, by = keys, suffix = c("_saved", "_audit"))
    tibble(script, stage, file_exists = TRUE, n_rows = nrow(saved),
           key_mismatches = nrow(compared) - nrow(inner_join(saved, expected, by = keys)),
           class_mismatches = sum(compared[[paste0(value, "_saved")]] !=
                                    compared[[paste0(value, "_audit")]], na.rm = TRUE))
  }))
}
emit(bind_rows(verify_saved("monthly", monthly %>% filter(variant == "no_uparea") %>% select(-variant)),
               verify_saved("annual", annual %>% filter(variant == "no_uparea") %>% select(-variant))),
     "saved_10_11_class_verification")

cat("Step 2 saved-output comparison complete.\n")
