# Replay only script 07's cumulative screening for the two saved prediction
# variants. No regression fits, figures, routing rebuild, or baseline overwrite.
# Run from the repository root. This is deliberately expensive; it expands
# the saved upstream lookup against one month at a time.
library(data.table)
library(sf)
library(janitor)
library(readr)

out <- "reports/no_uparea_migration/step1_filter_attribution"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
path_changes <- file.path(out, "changed_screened_rows.csv")
if (file.exists(path_changes)) stop("Refusing to overwrite an existing audit: ", path_changes)

cat_path <- "data/processed/merit_data/03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
riv_path <- "data/processed/merit_data/03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
characteristics_path <- "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv"

cat("Reading saved network and reach attributes...\n")
catchments <- as.data.table(st_drop_geometry(clean_names(read_sf(cat_path))))[, .(
  comid = as.character(comid), huc = as.character(huc), unitarea)]
flowlines <- as.data.table(st_drop_geometry(clean_names(read_sf(riv_path))))[, .(
  comid = as.character(comid), uparea)]
characteristics <- as.data.table(clean_names(read_csv(characteristics_path, show_col_types = FALSE)))[,
  .(comid = as.character(comid), lengthkm, merit_width)]
for (x in list(catchments, flowlines, characteristics)) {
  if (anyDuplicated(x$comid)) stop("Duplicate COMID in saved reach attributes")
}
reach <- merge(catchments, flowlines, by = "comid", all.x = TRUE, sort = FALSE)
reach <- merge(reach, characteristics, by = "comid", all.x = TRUE, sort = FALSE)
reach[, stream_area_m2 := lengthkm * 1000 * merit_width]
rm(catchments, flowlines, characteristics)

route <- as.data.table(readRDS("data/processed/merit_data/05_upstream_lookup.rds"))
route[, `:=`(downstream_comid = as.character(downstream_comid),
             upstream_comid = as.character(upstream_comid))]
if (anyDuplicated(route, by = c("downstream_comid", "upstream_comid"))) {
  stop("Duplicate saved routing pairs")
}
setkey(route, downstream_comid)
downstream <- unique(route$downstream_comid)
chunks <- split(downstream, ceiling(seq_along(downstream) / 5000L))

# This is the same area sum as 07, independent of process and prediction variant.
cat("Reconstructing cumulative upstream area once...\n")
unitarea <- reach[, .(upstream_comid = comid, upstream_unitarea = unitarea)]
area_parts <- vector("list", length(chunks))
for (i in seq_along(chunks)) {
  pairs <- route[.(chunks[[i]]), on = "downstream_comid", nomatch = 0L]
  joined <- merge(pairs, unitarea, by = "upstream_comid", all = FALSE, sort = FALSE)
  area_parts[[i]] <- joined[, .(cumulative_area_km2 = sum(upstream_unitarea, na.rm = TRUE)),
                            by = downstream_comid]
  rm(pairs, joined)
}
area <- rbindlist(area_parts)
rm(area_parts, unitarea)
area <- merge(area, reach[, .(downstream_comid = comid, huc, uparea)],
              by = "downstream_comid", all.x = TRUE, sort = FALSE)
area[, pct_area_missing := pmax(0, 1 - cumulative_area_km2 / uparea) * 100]
setkey(area, downstream_comid)

saved_paths <- c(original = "data/processed/07_huc08_allometry_input.rds",
                 no_uparea = "data/processed/07_no_uparea_huc08_allometry_input.rds")
source_paths <- c(original = "data/processed/01_assembled_maavara_monthly_metabolism.rds",
                  no_uparea = "data/processed/01_assembled_no_uparea_monthly_metabolism.rds")
sources <- lapply(source_paths, function(path) as.data.table(readRDS(path)))
saved <- lapply(saved_paths, function(path) as.data.table(readRDS(path))[, .(
  downstream_comid = as.character(downstream_comid), month,
  metabolism_metric = as.character(metabolism_metric))])
for (variant in names(saved)) {
  if (anyDuplicated(saved[[variant]], by = c("downstream_comid", "month", "metabolism_metric"))) {
    stop("Duplicate saved 07 keys for ", variant)
  }
}

rule_order <- c("huc", "month", "positive_area", "positive_flux", "area_coverage",
                "finite_metabolism_coverage", "metabolism_coverage")
ledger_parts <- list()
verification_parts <- list()
transition_parts <- list()
append_changes <- FALSE

replay_month <- function(source, month_i) {
  # Exact 07 local rule: negative areal rates become NA, then are multiplied by
  # stream area. A nonfinite stream area propagates to the local flux.
  local <- merge(source[month == month_i], reach, by = "comid", all = FALSE, sort = FALSE)
  local[, gpp_flux := fifelse(!is.na(gpp_gO2_m2_d) & gpp_gO2_m2_d < 0,
                             NA_real_, gpp_gO2_m2_d) * stream_area_m2]
  local[, er_flux := fifelse(!is.na(er_gO2_m2_d) & er_gO2_m2_d < 0,
                            NA_real_, er_gO2_m2_d) * stream_area_m2]
  local <- local[, .(upstream_comid = comid, unitarea, gpp_flux, er_flux)]
  parts <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    pairs <- route[.(chunks[[i]]), on = "downstream_comid", nomatch = 0L]
    joined <- merge(pairs, local, by = "upstream_comid", all = FALSE, sort = FALSE)
    parts[[i]] <- joined[, .(
      cumulative_gpp_g_o2_d = if (any(!is.na(gpp_flux))) sum(gpp_flux, na.rm = TRUE) else NA_real_,
      cumulative_er_g_o2_d = if (any(!is.na(er_flux))) sum(er_flux, na.rm = TRUE) else NA_real_,
      gpp_covered_area_km2 = sum(unitarea[!is.na(gpp_flux)], na.rm = TRUE),
      er_covered_area_km2 = sum(unitarea[!is.na(er_flux)], na.rm = TRUE)),
      by = downstream_comid]
    rm(pairs, joined)
  }
  cumulative <- rbindlist(parts)
  rm(parts, local)
  cumulative <- merge(cumulative, area, by = "downstream_comid", all.x = TRUE, sort = FALSE)
  result <- rbindlist(lapply(c("GPP", "ER"), function(process) {
    flux <- if (process == "GPP") cumulative$cumulative_gpp_g_o2_d else cumulative$cumulative_er_g_o2_d
    covered <- if (process == "GPP") cumulative$gpp_covered_area_km2 else cumulative$er_covered_area_km2
    data.table(downstream_comid = cumulative$downstream_comid, month = month_i,
               metabolism_metric = process, huc = cumulative$huc,
               cumulative_area_km2 = cumulative$cumulative_area_km2,
               cumulative_flux_g_o2_d = flux,
               pct_area_missing = cumulative$pct_area_missing,
               pct_metabolism_missing = (1 - covered / cumulative$cumulative_area_km2) * 100)
  }))
  rm(cumulative)
  result[, reason := fcase(
    is.na(huc), "huc",
    is.na(month), "month",
    is.na(cumulative_area_km2) | !(cumulative_area_km2 > 0), "positive_area",
    is.na(cumulative_flux_g_o2_d) | !(cumulative_flux_g_o2_d > 0), "positive_flux",
    is.na(pct_area_missing) | !(pct_area_missing < 10), "area_coverage",
    !is.finite(pct_metabolism_missing), "finite_metabolism_coverage",
    is.na(pct_metabolism_missing) | !(pct_metabolism_missing < 10), "metabolism_coverage",
    default = "retained")]
  result
}

cat("Replaying original and no-uparea screens by month...\n")
for (month_i in 1:12) {
  replay <- lapply(names(sources), function(variant) {
    x <- replay_month(sources[[variant]], month_i)
    actual <- saved[[variant]][month == month_i]
    predicted <- x[reason == "retained", .(downstream_comid, month, metabolism_metric)]
    setkey(predicted, downstream_comid, month, metabolism_metric)
    setkey(actual, downstream_comid, month, metabolism_metric)
    missing_from_replay <- nrow(actual[!predicted])
    extra_in_replay <- nrow(predicted[!actual])
    verification_parts[[length(verification_parts) + 1L]] <<- data.table(
      variant, month = month_i, saved_rows = nrow(actual), replay_rows = nrow(predicted),
      missing_from_replay, extra_in_replay)
    if (missing_from_replay || extra_in_replay) {
      stop("Replay does not match saved 07 for ", variant, " month ", month_i,
           ": missing=", missing_from_replay, ", extra=", extra_in_replay)
    }
    ledger_parts[[length(ledger_parts) + 1L]] <<- rbindlist(lapply(c("GPP", "ER"), function(process) {
      y <- x[metabolism_metric == process]
      removed <- tabulate(match(y$reason, rule_order), nbins = length(rule_order))
      before <- nrow(y) - c(0L, cumsum(removed)[-length(removed)])
      data.table(variant, month = month_i, metabolism_metric = process,
                 rule = rule_order, before, removed, after = before - removed)
    }))
    x[, c("huc", "cumulative_area_km2", "cumulative_flux_g_o2_d", "pct_area_missing") := NULL]
    setkey(x, downstream_comid, month, metabolism_metric)
    x
  })
  names(replay) <- names(sources)
  paired <- merge(replay$original, replay$no_uparea,
                  by = c("downstream_comid", "month", "metabolism_metric"),
                  all = TRUE, suffixes = c("_original", "_no_uparea"))
  paired[is.na(reason_original), reason_original := "absent_from_cumulative"]
  paired[is.na(reason_no_uparea), reason_no_uparea := "absent_from_cumulative"]
  paired[, change := fcase(
    reason_original == "retained" & reason_no_uparea != "retained", "original_only",
    reason_original != "retained" & reason_no_uparea == "retained", "no_uparea_only",
    default = "unchanged")]
  changed <- paired[change != "unchanged", .(downstream_comid, month, metabolism_metric, change,
    reason_original, reason_no_uparea,
    pct_metabolism_missing_original, pct_metabolism_missing_no_uparea)]
  fwrite(changed, path_changes, append = append_changes, col.names = !append_changes)
  append_changes <- TRUE
  transition_parts[[length(transition_parts) + 1L]] <- changed[, .N,
    by = .(metabolism_metric, change, reason_original, reason_no_uparea)]
  cat("  month ", month_i, ": ", nrow(changed), " changed screened keys; both replays match saved 07\n", sep = "")
  rm(replay, paired, changed)
  gc(FALSE)
}

verification <- rbindlist(verification_parts)
ledger <- rbindlist(ledger_parts)
transitions <- rbindlist(transition_parts)[,
  .(n = sum(N)), by = .(metabolism_metric, change, reason_original, reason_no_uparea)]
fwrite(verification, file.path(out, "replay_verification.csv"))
fwrite(ledger, file.path(out, "rule_ledger_by_month_process.csv"))
fwrite(ledger[, .(before = sum(before), removed = sum(removed), after = sum(after)),
              by = .(variant, metabolism_metric, rule)],
       file.path(out, "rule_ledger_by_process.csv"))
fwrite(transitions, file.path(out, "changed_key_rule_transitions.csv"))
cat("Complete: both variants exactly reproduced saved 07 retained keys for all 12 months.\n")
