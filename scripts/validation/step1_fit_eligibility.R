# Attribute changes in 07's eligible HUC8-month fits using its saved screened rows.
# This only repeats the three fit eligibility predicates; it does not fit models.
library(data.table)

out <- "reports/no_uparea_migration/step1_filter_attribution"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
screen_paths <- c(original = "data/processed/07_huc08_allometry_input.rds",
                  no_uparea = "data/processed/07_no_uparea_huc08_allometry_input.rds")
fit_paths <- c(original = "data/processed/07_huc08_allometry_results.csv",
               no_uparea = "data/processed/07_no_uparea_huc08_allometry_results.csv")

summarise_groups <- function(path) {
  x <- as.data.table(readRDS(path))
  x[, huc := as.character(huc)]
  x[, .(n_reaches = .N,
        n_distinct_areas = uniqueN(cumulative_area_km2),
        log10_area_range = log10(max(cumulative_area_km2) / min(cumulative_area_km2))),
    by = .(huc, month, metabolism_metric)]
}

groups <- lapply(screen_paths, summarise_groups)
fits <- lapply(fit_paths, function(path) fread(path, colClasses = list(character = "huc")))
for (variant in names(groups)) {
  x <- groups[[variant]]
  x[, eligibility := fcase(
    n_reaches < 20L, "fewer_than_20_reaches",
    n_distinct_areas < 20L, "fewer_than_20_distinct_areas",
    log10_area_range < 1, "area_range_under_one_decade",
    default = "eligible")]
  actual <- fits[[variant]][, .(huc = as.character(huc), month, metabolism_metric)]
  predicted <- x[eligibility == "eligible", .(huc, month, metabolism_metric)]
  setkeyv(actual, c("huc", "month", "metabolism_metric"))
  setkeyv(predicted, c("huc", "month", "metabolism_metric"))
  if (nrow(actual[!predicted]) || nrow(predicted[!actual])) {
    stop("Eligibility replay does not match saved fit keys for ", variant)
  }
}

paired <- merge(groups$original, groups$no_uparea,
                by = c("huc", "month", "metabolism_metric"), all = TRUE,
                suffixes = c("_original", "_no_uparea"))
for (variant in names(groups)) {
  column <- paste0("eligibility_", variant)
  paired[is.na(get(column)), (column) := "no_screened_reaches"]
}
paired[, change := fcase(
  eligibility_original == "eligible" & eligibility_no_uparea != "eligible", "original_only",
  eligibility_original != "eligible" & eligibility_no_uparea == "eligible", "no_uparea_only",
  default = "unchanged")]
changed <- paired[change != "unchanged"]
transitions <- changed[, .N, by = .(metabolism_metric, change,
                                    eligibility_original, eligibility_no_uparea)]
fwrite(changed, file.path(out, "changed_fit_eligibility.csv"))
fwrite(transitions, file.path(out, "fit_eligibility_transitions.csv"))
cat("Fit eligibility replay matched saved 07 fit keys in both variants.\n")
