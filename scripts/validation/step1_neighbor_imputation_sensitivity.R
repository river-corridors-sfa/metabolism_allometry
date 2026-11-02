# Sensitivity only: replace a negative no-uparea areal rate with the mean of
# (a) all direct upstream areal rates and (b) the direct downstream areal rate,
# when every direct neighbor is positive. One simultaneous, noniterative pass.
# Replays 07's existing screens using saved network inputs; no production edits.
library(data.table)
library(sf)
library(janitor)
library(readr)

out <- "reports/no_uparea_migration/step1_neighbor_imputation_sensitivity"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out, "screened_recovery.csv"))) stop("Refusing to overwrite completed sensitivity output")

cat_path <- "data/processed/merit_data/03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
riv_path <- "data/processed/merit_data/03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
cat("Reading saved inputs...\n")
cat_dt <- as.data.table(st_drop_geometry(clean_names(read_sf(cat_path))))[, .(
  comid = as.character(comid), huc = as.character(huc), unitarea)]
riv_dt <- as.data.table(st_drop_geometry(clean_names(read_sf(riv_path))))[, .(
  comid = as.character(comid), uparea)]
chars <- as.data.table(clean_names(read_csv(
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv", show_col_types = FALSE)))[,
  .(comid = as.character(comid), lengthkm, merit_width)]
reach <- merge(merge(cat_dt, riv_dt, by = "comid", all.x = TRUE, sort = FALSE),
               chars, by = "comid", all.x = TRUE, sort = FALSE)
reach[, stream_area_m2 := lengthkm * 1000 * merit_width]
rm(cat_dt, riv_dt, chars)
edges <- as.data.table(readRDS("data/processed/merit_data/05_edges_downstream.rds"))
edges[, `:=`(from_comid = as.character(from_comid), to_comid = as.character(to_comid))]
if (anyDuplicated(edges, by = c("from_comid", "to_comid"))) stop("Duplicate direct edges")
route <- as.data.table(readRDS("data/processed/merit_data/05_upstream_lookup.rds"))
route[, `:=`(downstream_comid = as.character(downstream_comid),
             upstream_comid = as.character(upstream_comid))]
setkey(route, downstream_comid)
chunks <- split(unique(route$downstream_comid),
                ceiling(seq_along(unique(route$downstream_comid)) / 5000L))
unitarea <- reach[, .(upstream_comid = comid, upstream_unitarea = unitarea)]
area_parts <- lapply(chunks, function(ids) {
  pairs <- route[.(ids), on = "downstream_comid", nomatch = 0L]
  merge(pairs, unitarea, by = "upstream_comid", all = FALSE, sort = FALSE)[,
    .(cumulative_area_km2 = sum(upstream_unitarea, na.rm = TRUE)), by = downstream_comid]
})
area <- rbindlist(area_parts)
area <- merge(area, reach[, .(downstream_comid = comid, huc, uparea)],
              by = "downstream_comid", all.x = TRUE, sort = FALSE)
area[, pct_area_missing := pmax(0, 1 - cumulative_area_km2 / uparea) * 100]
rm(area_parts, unitarea)

source <- as.data.table(readRDS("data/processed/01_assembled_no_uparea_monthly_metabolism.rds"))
source[, comid := as.character(comid)]
saved <- as.data.table(readRDS("data/processed/07_no_uparea_huc08_allometry_input.rds"))[, .(
  downstream_comid = as.character(downstream_comid), month, metabolism_metric)]
setkey(saved, downstream_comid, month, metabolism_metric)

candidate_parts <- list()
recovery_parts <- list()
count_parts <- list()
for (month_i in 1:12) {
  local <- copy(source[month == month_i])
  for (process in c("GPP", "ER")) {
    rate_col <- if (process == "GPP") "gpp_gO2_m2_d" else "er_gO2_m2_d"
    rates <- local[, .(comid, rate = get(rate_col))]
    # Every immediate upstream neighbor must have a positive observed prediction.
    upstream <- merge(edges, rates, by.x = "from_comid", by.y = "comid",
                      all.x = TRUE, sort = FALSE)[,
      .(n_upstream = .N, n_positive_upstream = sum(!is.na(rate) & rate > 0),
        upstream_mean = if (all(!is.na(rate) & rate > 0)) mean(rate) else NA_real_),
      by = to_comid]
    downstream <- merge(edges, rates, by.x = "to_comid", by.y = "comid",
                        all.x = TRUE, sort = FALSE)[,
      .(comid = from_comid, downstream_rate = rate)]
    candidates <- merge(rates[!is.na(rate) & rate < 0], upstream,
                        by.x = "comid", by.y = "to_comid", all.x = TRUE, sort = FALSE)
    candidates <- merge(candidates, downstream, by = "comid", all.x = TRUE, sort = FALSE)
    candidates[, eligible := !is.na(upstream_mean) & !is.na(downstream_rate) &
                downstream_rate > 0]
    candidates[, imputed_rate := fifelse(eligible,
                                       (upstream_mean + downstream_rate) / 2, NA_real_)]
    candidate_parts[[length(candidate_parts) + 1L]] <- candidates[, .(
      comid, month = month_i, metabolism_metric = process, original_rate = rate,
      n_upstream, n_positive_upstream, upstream_mean, downstream_rate,
      eligible, imputed_rate)]
    replacements <- candidates[eligible == TRUE, .(comid, imputed_rate)]
    if (nrow(replacements)) {
      setkey(rates, comid)
      rates[replacements, rate := i.imputed_rate, on = "comid"]
    }
    local[, (rate_col) := rates[local, rate, on = "comid"]]
  }

  local <- merge(local, reach, by = "comid", all = FALSE, sort = FALSE)
  local[, gpp_flux := fifelse(!is.na(gpp_gO2_m2_d) & gpp_gO2_m2_d < 0,
                             NA_real_, gpp_gO2_m2_d) * stream_area_m2]
  local[, er_flux := fifelse(!is.na(er_gO2_m2_d) & er_gO2_m2_d < 0,
                            NA_real_, er_gO2_m2_d) * stream_area_m2]
  local <- local[, .(upstream_comid = comid, unitarea, gpp_flux, er_flux)]
  parts <- lapply(chunks, function(ids) {
    pairs <- route[.(ids), on = "downstream_comid", nomatch = 0L]
    merge(pairs, local, by = "upstream_comid", all = FALSE, sort = FALSE)[,
      .(gpp_flux = if (any(!is.na(gpp_flux))) sum(gpp_flux, na.rm = TRUE) else NA_real_,
        er_flux = if (any(!is.na(er_flux))) sum(er_flux, na.rm = TRUE) else NA_real_,
        gpp_covered = sum(unitarea[!is.na(gpp_flux)], na.rm = TRUE),
        er_covered = sum(unitarea[!is.na(er_flux)], na.rm = TRUE)),
      by = downstream_comid]
  })
  cumulative <- merge(rbindlist(parts), area, by = "downstream_comid", all.x = TRUE,
                      sort = FALSE)
  screened <- rbindlist(lapply(c("GPP", "ER"), function(process) {
    flux <- if (process == "GPP") cumulative$gpp_flux else cumulative$er_flux
    covered <- if (process == "GPP") cumulative$gpp_covered else cumulative$er_covered
    data.table(downstream_comid = cumulative$downstream_comid, month = month_i,
               metabolism_metric = process, huc = cumulative$huc,
               cumulative_area_km2 = cumulative$cumulative_area_km2,
               cumulative_flux_g_o2_d = flux,
               pct_area_missing = cumulative$pct_area_missing,
               pct_metabolism_missing = (1 - covered / cumulative$cumulative_area_km2) * 100)
  }))
  screened <- screened[!is.na(huc) & !is.na(month) & cumulative_area_km2 > 0 &
                       !is.na(cumulative_flux_g_o2_d) & cumulative_flux_g_o2_d > 0 &
                       !is.na(pct_area_missing) & pct_area_missing < 10 &
                       is.finite(pct_metabolism_missing) & pct_metabolism_missing < 10,
                       .(downstream_comid, month, metabolism_metric)]
  setkey(screened, downstream_comid, month, metabolism_metric)
  baseline <- saved[month == month_i]
  gained <- screened[!baseline]
  lost <- baseline[!screened]
  if (nrow(lost)) stop("Positive-only imputation unexpectedly lost screened rows in month ", month_i)
  recovery_parts[[length(recovery_parts) + 1L]] <- gained
  count_parts[[length(count_parts) + 1L]] <- data.table(month = month_i,
    metabolism_metric = c("GPP", "ER"),
    baseline_rows = as.integer(table(factor(baseline$metabolism_metric, levels = c("GPP", "ER")))),
    scenario_rows = as.integer(table(factor(screened$metabolism_metric, levels = c("GPP", "ER")))),
    gained_rows = as.integer(table(factor(gained$metabolism_metric, levels = c("GPP", "ER")))))
  cat("month ", month_i, ": ", nrow(gained), " recovered screened rows\n", sep = "")
  gc(FALSE)
}

candidates <- rbindlist(candidate_parts, fill = TRUE)
recovery <- rbindlist(recovery_parts)
counts <- rbindlist(count_parts)
fwrite(candidates, file.path(out, "negative_reach_candidates.csv"))
fwrite(candidates[, .(negative_rows = .N, eligible = sum(eligible)),
                  by = .(metabolism_metric)], file.path(out, "candidate_summary.csv"))
fwrite(recovery, file.path(out, "recovered_screened_rows.csv"))
fwrite(counts, file.path(out, "screened_recovery.csv"))
cat("Sensitivity complete; original outputs unchanged.\n")
