# Saved-output audit only. Does not refit, route, or rerun script 07.
library(data.table)
out <- "reports/no_uparea_migration/step1_sample_size"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
emit <- function(x, name) fwrite(x, file.path(out, paste0(name, ".csv")))

source_old <- as.data.table(readRDS("data/processed/01_assembled_maavara_monthly_metabolism.rds"))
source_new <- as.data.table(readRDS("data/processed/01_assembled_no_uparea_monthly_metabolism.rds"))
filtered_old <- as.data.table(readRDS("data/processed/03_filtered_maavara_monthly_metabolism.rds"))
filtered_new <- as.data.table(readRDS("data/processed/03_filtered_no_uparea_monthly_metabolism.rds"))
screen_old <- as.data.table(readRDS("data/processed/07_huc08_allometry_input.rds"))
screen_new <- as.data.table(readRDS("data/processed/07_no_uparea_huc08_allometry_input.rds"))
fit_old <- fread("data/processed/07_huc08_allometry_results.csv", colClasses = list(character = "huc"))
fit_new <- fread("data/processed/07_no_uparea_huc08_allometry_results.csv", colClasses = list(character = "huc"))

key_check <- function(x, columns, label) {
  if (anyDuplicated(x[, ..columns])) stop("Duplicate keys in ", label)
}
key_check(source_old, c("comid", "month"), "original 01")
key_check(source_new, c("comid", "month"), "no-uparea 01")
key_check(screen_old, c("downstream_comid", "month", "metabolism_metric"), "original 07")
key_check(screen_new, c("downstream_comid", "month", "metabolism_metric"), "no-uparea 07")
key_check(fit_old, c("huc", "month", "metabolism_metric"), "original fits")
key_check(fit_new, c("huc", "month", "metabolism_metric"), "no-uparea fits")

stage <- rbindlist(lapply(list(original_01 = source_old, no_uparea_01 = source_new,
                           original_03 = filtered_old, no_uparea_03 = filtered_new),
  function(x) data.table(rows = nrow(x), comids = uniqueN(x$comid),
                         months = uniqueN(x$month), comid_months = uniqueN(x, by = c("comid", "month")))),
  idcol = "stage")
emit(stage, "stage_counts")

rates <- rbindlist(lapply(c("GPP", "ER"), function(process) {
  old_col <- if (process == "GPP") "gpp_gO2_m2_d" else "er_gO2_m2_d"
  filtered_col <- if (process == "GPP") "gpp_g_o2_m2_d" else "er_g_o2_m2_d"
  make <- function(x, variant, stage, col) {
    v <- x[[col]]
    data.table(process, variant, stage, rows = length(v), missing = sum(is.na(v)),
               negative = sum(v < 0, na.rm = TRUE), zero = sum(v == 0, na.rm = TRUE),
               nonfinite = sum(!is.finite(v)), finite = sum(is.finite(v)))
  }
  rbind(make(source_old, "original", "01", old_col),
        make(source_new, "no_uparea", "01", old_col),
        make(filtered_old, "original", "03", filtered_col),
        make(filtered_new, "no_uparea", "03", filtered_col))
}))
emit(rates, "rate_counts")

source_transitions <- rbindlist(lapply(c("GPP", "ER"), function(process) {
  col <- if (process == "GPP") "gpp_gO2_m2_d" else "er_gO2_m2_d"
  old <- source_old[, .(comid, month, original = get(col))]
  new <- source_new[, .(comid, month, no_uparea = get(col))]
  setkey(old, comid, month); setkey(new, comid, month)
  if (nrow(old[!new]) || nrow(new[!old])) stop("01 key mismatch: ", process)
  x <- merge(old, new, by = c("comid", "month"))
  x[, .(process = process, matched_keys = .N,
        new_negative_only = sum(is.finite(no_uparea) & no_uparea < 0 &
                                (!is.finite(original) | original >= 0)),
        old_negative_only = sum(is.finite(original) & original < 0 &
                                (!is.finite(no_uparea) | no_uparea >= 0)),
        both_negative = sum(is.finite(original) & original < 0 &
                            is.finite(no_uparea) & no_uparea < 0),
        newly_missing = sum(is.finite(original) & !is.finite(no_uparea)),
        recovered_missing = sum(!is.finite(original) & is.finite(no_uparea)))]
}))
emit(source_transitions, "source_rate_transitions")

screen_keys <- c("downstream_comid", "month", "metabolism_metric")
old_keys <- unique(screen_old[, ..screen_keys])
new_keys <- unique(screen_new[, ..screen_keys])
setkeyv(old_keys, screen_keys); setkeyv(new_keys, screen_keys)
lost <- old_keys[!new_keys]
gained <- new_keys[!old_keys]
lost[, change := "original_only"]
gained[, change := "no_uparea_only"]
changes <- rbindlist(list(lost, gained))
setcolorder(changes, c("change", screen_keys))
emit(changes, "screened_key_changes")
screen_counts <- rbindlist(lapply(c("GPP", "ER"), function(process) {
  a <- screen_old[metabolism_metric == process]
  b <- screen_new[metabolism_metric == process]
  l <- lost[metabolism_metric == process]
  g <- gained[metabolism_metric == process]
  a_reaches <- unique(a$downstream_comid)
  b_reaches <- unique(b$downstream_comid)
  data.table(process, original_rows = nrow(a), no_uparea_rows = nrow(b),
             matched_rows = nrow(a) - nrow(l), original_only = nrow(l),
             no_uparea_only = nrow(g), net_change = nrow(b) - nrow(a),
             pct_change = 100 * (nrow(b) - nrow(a)) / nrow(a),
             original_reaches = uniqueN(a$downstream_comid),
             no_uparea_reaches = uniqueN(b$downstream_comid),
             reaches_lost_completely = length(setdiff(a_reaches, b_reaches)),
             reaches_gained_completely = length(setdiff(b_reaches, a_reaches)),
             original_huc8s = uniqueN(a$huc), no_uparea_huc8s = uniqueN(b$huc))
}))
emit(screen_counts, "screened_counts")
reach_changes <- rbindlist(lapply(c("GPP", "ER"), function(process) {
  a <- unique(screen_old[metabolism_metric == process, downstream_comid])
  b <- unique(screen_new[metabolism_metric == process, downstream_comid])
  rbind(data.table(process, change = "original_only", downstream_comid = setdiff(a, b)),
        data.table(process, change = "no_uparea_only", downstream_comid = setdiff(b, a)))
}))
emit(reach_changes, "reach_changes")
emit(changes[, .N, by = .(change, metabolism_metric, month)][order(metabolism_metric, month, change)],
     "screened_changes_by_month")

setkeyv(screen_old, screen_keys); setkeyv(screen_new, screen_keys)
lost_baseline <- screen_old[lost, on = screen_keys, nomatch = 0L]
gained_new <- screen_new[gained, on = screen_keys, nomatch = 0L]
coverage_margin <- rbindlist(list(
  lost_baseline[, .(change = "original_only", .N, median_pct_metabolism_missing = median(pct_metabolism_missing),
                    n_8_to_10_pct = sum(pct_metabolism_missing >= 8),
                    n_9_to_10_pct = sum(pct_metabolism_missing >= 9)), by = metabolism_metric],
  gained_new[, .(change = "no_uparea_only", .N, median_pct_metabolism_missing = median(pct_metabolism_missing),
                 n_8_to_10_pct = sum(pct_metabolism_missing >= 8),
                 n_9_to_10_pct = sum(pct_metabolism_missing >= 9)), by = metabolism_metric]))
emit(coverage_margin, "changed_row_coverage_margins")

# A focal reach's own source-rate transition is one possible explanation;
# upstream transitions and coverage can also affect it.
local_status <- rbindlist(lapply(c("GPP", "ER"), function(process) {
  col <- if (process == "GPP") "gpp_gO2_m2_d" else "er_gO2_m2_d"
  x <- merge(source_old[, .(downstream_comid = comid, month, original = get(col))],
             source_new[, .(downstream_comid = comid, month, no_uparea = get(col))],
             by = c("downstream_comid", "month"))
  x[, metabolism_metric := process]
  x
}))
setkeyv(local_status, screen_keys)
changed_status <- merge(changes, local_status, by = screen_keys, all.x = TRUE)
changed_status[, focal_status := fcase(
  !is.finite(no_uparea), "new_missing",
  no_uparea < 0 & is.finite(original) & original >= 0, "new_negative_only",
  is.finite(original) & original < 0 & no_uparea >= 0, "old_negative_only",
  no_uparea < 0 & original < 0, "both_negative",
  default = "neither_negative")]
emit(changed_status[, .N, by = .(change, metabolism_metric, focal_status)],
     "screened_change_focal_rate_status")

fit_keys <- c("huc", "month", "metabolism_metric")
setkeyv(fit_old, fit_keys); setkeyv(fit_new, fit_keys)
fit_changes <- data.table(original_only = nrow(fit_old[!fit_new]),
                          no_uparea_only = nrow(fit_new[!fit_old]),
                          matched = nrow(fit_old) - nrow(fit_old[!fit_new]))
emit(fit_changes, "fit_key_changes")
emit(rbindlist(list(fit_old[, .(variant = "original", .N), by = metabolism_metric],
                    fit_new[, .(variant = "no_uparea", .N), by = metabolism_metric])),
     "fit_counts_by_process")
emit(fread("data/processed/07_no_uparea_screen_attrition.csv"), "no_uparea_screen_rules")

cat("Step 1 saved-output audit complete. Exact original per-rule exclusions are unavailable: original 07 saved only final screened rows.\n")
