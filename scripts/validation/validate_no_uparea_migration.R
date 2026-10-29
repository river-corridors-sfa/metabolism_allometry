# Run from repository root after the runbook's production steps.
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
out <- "reports/no_uparea_migration"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
emit <- function(x, name) write_csv(x, file.path(out, paste0(name, ".csv")), na = "")
summarize_rates <- function(x, stage) {
  x %>% pivot_longer(c(gpp_gO2_m2_d, er_gO2_m2_d), names_to = "process", values_to = "rate") %>%
    group_by(process) %>% summarise(stage = stage, rows = n(), comids = n_distinct(comid),
      comid_months = n_distinct(paste(comid, month)), months = n_distinct(month),
      missing = sum(is.na(rate)), negative = sum(rate < 0, na.rm = TRUE),
      zero = sum(rate == 0, na.rm = TRUE), nonfinite = sum(!is.finite(rate)),
      mean = mean(rate[is.finite(rate)]), median = median(rate[is.finite(rate)]), .groups = "drop")
}
old <- readRDS("data/processed/01_assembled_maavara_monthly_metabolism.rds") %>% mutate(comid = as.character(comid))
new <- readRDS("data/processed/01_assembled_no_uparea_monthly_metabolism.rds") %>% mutate(comid = as.character(comid))
for (x in list(old, new)) {
  stopifnot(!anyDuplicated(x[c("comid", "month")]), setequal(x$month, 1:12))
}
original_only_keys <- nrow(anti_join(old, new, by = c("comid", "month")))
no_uparea_only_keys <- nrow(anti_join(new, old, by = c("comid", "month")))
if (original_only_keys || no_uparea_only_keys) stop("Source COMID-month key mismatch")
emit(bind_rows(summarize_rates(old, "original"), summarize_rates(new, "no_uparea")), "source_rate_counts")
long <- function(x, label) x %>% pivot_longer(c(gpp_gO2_m2_d, er_gO2_m2_d), names_to = "process", values_to = label)
paired <- full_join(long(old, "original"), long(new, "no_uparea"), by = c("comid", "month", "process"))
coverage <- paired %>% group_by(process) %>% summarise(original_keys = n(),
  no_uparea_keys = n(), matched_finite = sum(is.finite(original) & is.finite(no_uparea)),
  original_only_keys = original_only_keys, no_uparea_only_keys = no_uparea_only_keys,
  original_missing_value = sum(is.na(original)), no_uparea_missing_value = sum(is.na(no_uparea)),
  pct_matched_finite = 100 * matched_finite / n(), .groups = "drop")
emit(coverage, "source_coverage")
agreement <- paired %>% filter(is.finite(original), is.finite(no_uparea)) %>% group_by(process) %>%
  summarise(n = n(), correlation = cor(original, no_uparea),
    bias = mean(no_uparea - original), rmse = sqrt(mean((no_uparea - original)^2)),
    mean_difference = mean(no_uparea - original), .groups = "drop")
emit(agreement, "prediction_agreement")
set.seed(24)
plot_rows <- paired %>% filter(is.finite(original), is.finite(no_uparea)) %>%
  group_by(process) %>% group_modify(~ slice_sample(.x, n = min(nrow(.x), 10000))) %>% ungroup()
ggsave(file.path(out, "paired_predictions.png"),
  ggplot(plot_rows, aes(original, no_uparea)) + geom_point(alpha = .1, size = .3) +
    geom_abline(slope = 1, intercept = 0, color = "red") + facet_wrap(~process, scales = "free") + theme_bw(),
  width = 8, height = 4, dpi = 180)

# These are paired refits on the same split, not a test against the published forest.
metrics <- read_csv("data/processed/17_maavara_refits/model_metrics.csv", show_col_types = FALSE) %>%
  filter(dataset == "test") %>% select(process, variant, n, r2_score, rmse, mae)
skill <- metrics %>% pivot_wider(names_from = variant, values_from = c(n, r2_score, rmse, mae)) %>%
  mutate(identical_evaluation_n = n_with_uparea == n_without_uparea,
    delta_r2 = r2_score_without_uparea - r2_score_with_uparea,
    delta_rmse = rmse_without_uparea - rmse_with_uparea,
    relative_rmse_pct = 100 * delta_rmse / rmse_with_uparea,
    delta_mae = mae_without_uparea - mae_with_uparea,
    relative_mae_pct = 100 * delta_mae / mae_with_uparea)
emit(skill, "paired_refit_test_skill")

stage_counts <- function(path, stage) {
  x <- readRDS(path)
  tibble(stage = stage, rows = nrow(x), comids = n_distinct(x$comid),
    comid_months = n_distinct(paste(x$comid, x$month)), months = n_distinct(x$month))
}
emit(bind_rows(stage_counts("data/processed/01_assembled_maavara_monthly_metabolism.rds", "01_original"),
  stage_counts("data/processed/01_assembled_no_uparea_monthly_metabolism.rds", "01_no_uparea"),
  stage_counts("data/processed/03_filtered_maavara_monthly_metabolism.rds", "03_original"),
  stage_counts("data/processed/03_filtered_no_uparea_monthly_metabolism.rds", "03_no_uparea")), "stage_counts")
screened <- bind_rows(
  readRDS("data/processed/07_huc08_allometry_input.rds") %>% mutate(variant = "original"),
  readRDS("data/processed/07_no_uparea_huc08_allometry_input.rds") %>% mutate(variant = "no_uparea"))
emit(screened %>% group_by(variant, metabolism_metric, month) %>%
  summarise(rows = n(), reaches = n_distinct(downstream_comid), huc8s = n_distinct(huc), .groups = "drop"),
  "screened_reach_month_counts")
emit(screened %>% select(variant, downstream_comid, month, metabolism_metric, huc) %>%
  pivot_wider(names_from = variant, values_from = huc) %>%
  mutate(composition = case_when(is.na(original) ~ "new_only", is.na(no_uparea) ~ "original_only", TRUE ~ "matched")) %>%
  count(metabolism_metric, month, composition), "screened_composition")

fit_paths <- c(original = "data/processed/07_huc08_allometry_results.csv",
  no_uparea = "data/processed/07_no_uparea_huc08_allometry_results.csv")
fits <- lapply(names(fit_paths), function(variant) read_csv(fit_paths[[variant]], show_col_types = FALSE,
  col_types = cols(huc = col_character(), .default = col_guess())) %>%
  mutate(variant = variant)) %>% bind_rows()
emit(fits %>% count(variant, metabolism_metric, month), "fit_counts")
paired_fits <- fits %>% select(huc, month, metabolism_metric, variant, slope, r2, n_reaches) %>%
  pivot_wider(names_from = variant, values_from = c(slope, r2, n_reaches)) %>%
  mutate(delta_slope = slope_no_uparea - slope_original,
    delta_r2 = r2_no_uparea - r2_original,
    delta_reaches = n_reaches_no_uparea - n_reaches_original)
emit(paired_fits, "paired_allometry_fits")
ggsave(file.path(out, "paired_exponents.png"),
  ggplot(paired_fits, aes(slope_original, slope_no_uparea)) + geom_point(alpha = .12, size = .3) +
    geom_abline(slope = 1, intercept = 0, color = "red") + facet_wrap(~metabolism_metric) + theme_bw(),
  width = 8, height = 4, dpi = 180)
stability <- bind_rows(
  read_csv("data/processed/08_huc08_fit_stability.csv", show_col_types = FALSE,
    col_types = cols(huc = col_character(), .default = col_guess())) %>% mutate(variant = "original"),
  read_csv("data/processed/08_no_uparea_huc08_fit_stability.csv", show_col_types = FALSE,
    col_types = cols(huc = col_character(), .default = col_guess())) %>% mutate(variant = "no_uparea"))
emit(stability %>% count(variant, metabolism_metric, month), "stability_counts")
classes <- stability %>% select(huc, month, metabolism_metric, variant, baseline_class) %>%
  pivot_wider(names_from = variant, values_from = baseline_class) %>%
  mutate(changed = !is.na(original) & !is.na(no_uparea) & original != no_uparea)
emit(classes, "paired_stability_classes")
emit(classes %>% count(metabolism_metric, original, no_uparea, changed), "classification_transitions")
segmented <- bind_rows(
  read_csv("figures/12_er_gpp_area_segmented_fit.csv", show_col_types = FALSE) %>% mutate(variant = "original"),
  read_csv("figures/12_no_uparea_er_gpp_area_segmented_fit.csv", show_col_types = FALSE) %>% mutate(variant = "no_uparea"))
emit(segmented, "segmented_comparison")
emit(read_csv("data/processed/07_no_uparea_screen_attrition.csv", show_col_types = FALSE), "no_uparea_screen_attrition")
writeLines(c("Migration validation summary", paste("Source matched finite records:", paste(coverage$process, coverage$matched_finite, collapse = "; ")),
  paste("Paired refit test R2 changes:", paste(skill$process, signif(skill$delta_r2, 4), collapse = "; ")),
  paste("Allometry fit counts:", paste(fits %>% count(variant) %>% transmute(x = paste(variant, n)) %>% pull(x), collapse = "; ")),
  "Inspect all CSV ledgers before scientific acceptance."), file.path(out, "summary.txt"))
