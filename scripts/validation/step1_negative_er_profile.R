# Saved-output profile of original and no-uparea ER predictions.
# Read-only analysis; does not alter production inputs or filtering.
library(data.table)
library(sf)
library(janitor)

out <- "reports/no_uparea_migration/step1_negative_er_profile"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
old <- as.data.table(readRDS("data/processed/01_assembled_maavara_monthly_metabolism.rds"))[
  , .(comid = as.character(comid), month, old_er = er_gO2_m2_d)]
new <- as.data.table(readRDS("data/processed/01_assembled_no_uparea_monthly_metabolism.rds"))[
  , .(comid = as.character(comid), month, new_er = er_gO2_m2_d)]
x <- merge(old, new, by = c("comid", "month"), all = TRUE)
stopifnot(nrow(x) == 2101392L, !anyDuplicated(x, by = c("comid", "month")))
x[, transition := fcase(
  is.na(new_er), "new_missing",
  old_er < 0 & new_er < 0, "negative_both",
  old_er >= 0 & new_er < 0, "new_negative_only",
  old_er < 0 & new_er >= 0, "old_negative_only",
  default = "nonnegative_both")]

counts <- x[, .(n = .N), by = .(month, transition)]
fwrite(counts, file.path(out, "sign_transitions_by_month.csv"))
fwrite(x[, .(n = .N), by = transition], file.path(out, "sign_transitions_total.csv"))

# Like-for-like comparison of the two saved refits isolates the predictor-set
# change better than comparing either refit with downloaded original predictions.
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
                 "Sep", "Oct", "Nov", "Dec")
refit_signs <- rbindlist(lapply(month_names, function(label) {
  root <- "data/processed/17_maavara_refits"
  with <- fread(file.path(root, "with_uparea", "ER", paste0("scaling_ER_pred_", label, ".csv")),
                colClasses = list(character = "COMID"))
  without <- fread(file.path(root, "without_uparea", "ER", paste0("scaling_ER_pred_", label, ".csv")),
                   colClasses = list(character = "COMID"))
  setnames(with, 2L, "with_er")
  setnames(without, 2L, "without_er")
  paired <- merge(with, without, by = "COMID", all = TRUE)
  stopifnot(nrow(paired) == 175116L, !anyDuplicated(paired$COMID))
  paired[, .(month = label, n = .N,
             with_negative = sum(!is.na(with_er) & with_er < 0),
             without_negative = sum(!is.na(without_er) & without_er < 0),
             new_negative_only = sum(!is.na(with_er) & with_er >= 0 &
                                       !is.na(without_er) & without_er < 0),
             with_negative_only = sum(!is.na(with_er) & with_er < 0 &
                                        !is.na(without_er) & without_er >= 0),
             negative_both = sum(!is.na(with_er) & with_er < 0 &
                                   !is.na(without_er) & without_er < 0),
             without_missing = sum(is.na(without_er)))]
}))
fwrite(refit_signs, file.path(out, "paired_refit_signs_by_month.csv"))

# The refit code multiplies raw ER by -1, with no response-sign exclusion.
training_er <- fread("data/raw/dryad_data/appling_data_06202024.csv",
                     select = "ER")$ER
fwrite(data.table(training_rows = length(training_er),
                  raw_er_negative = sum(training_er < 0),
                  raw_er_positive = sum(training_er > 0),
                  raw_er_zero = sum(training_er == 0),
                  magnitude_target_negative = sum(-training_er < 0)),
       file.path(out, "training_er_signs.csv"))

neg <- x[!is.na(new_er) & new_er < 0]
quant <- function(z) as.list(setNames(as.numeric(quantile(z, c(0, .01, .05, .1, .25, .5,
                                                 .75, .9, .95, .99, 1), na.rm = TRUE)),
                                  c("min", "p01", "p05", "p10", "p25", "median",
                                    "p75", "p90", "p95", "p99", "max")))
overall <- neg[, c(list(n = .N, mean = mean(new_er)), quant(new_er))]
fwrite(overall, file.path(out, "negative_magnitude_quantiles.csv"))
fwrite(neg[, .(n = .N, mean = mean(new_er), median = median(new_er),
               p10 = unname(quantile(new_er, .1)), p90 = unname(quantile(new_er, .9))),
           by = month], file.path(out, "negative_magnitude_by_month.csv"))
fwrite(neg[, .(n = .N), by = .(transition,
    magnitude = cut(-new_er, c(0, .01, .1, .5, 1, 2, 5, 10, Inf), right = FALSE))],
  file.path(out, "negative_magnitude_bins.csv"))

# Positive original rates where the new model predicts negative.
newly_negative <- x[transition == "new_negative_only"]
fwrite(newly_negative[, c(list(n = .N, mean = mean(old_er)), quant(old_er))],
       file.path(out, "old_rates_for_new_negatives.csv"))
fwrite(newly_negative[, .(n = .N), by = .(month,
  old_rate_bin = cut(old_er, c(0, .01, .1, .5, 1, 2, 5, 10, Inf), right = FALSE))],
  file.path(out, "old_rate_bins_for_new_negatives.csv"))

# Persistence by COMID and HUC for reaches in the filtered 03 catchment network.
network <- as.data.table(st_drop_geometry(clean_names(read_sf(
  "data/processed/merit_data/03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"))))[
  , .(comid = as.character(comid), huc = as.character(huc))]
stopifnot(!anyDuplicated(network$comid))
network_x <- merge(x, network, by = "comid", all = FALSE)
per_reach <- network_x[, .(n_negative = sum(!is.na(new_er) & new_er < 0),
                           n_newly_negative = sum(transition == "new_negative_only")),
                       by = .(comid, huc)]
fwrite(per_reach[, .(reaches = .N), by = n_negative],
       file.path(out, "negative_months_per_reach.csv"))
fwrite(per_reach[, .(reaches = .N, reaches_ever_negative = sum(n_negative > 0),
                     negative_reach_months = sum(n_negative)),
                 by = .(huc2 = substr(huc, 1, 2))],
       file.path(out, "negative_by_huc2.csv"))
fwrite(network_x[, .(reach_months = .N,
                     negative_reach_months = sum(!is.na(new_er) & new_er < 0)),
                 by = huc][order(-negative_reach_months)][1:30],
       file.path(out, "top_30_huc8_negative_counts.csv"))
cat("Profile complete.\n")
