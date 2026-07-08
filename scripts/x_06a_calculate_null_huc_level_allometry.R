# 06a_test_huc08_null_coherence.R

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(
  furrr,
  future,
  nhdplusTools,
  progressr
)


# 0. Settings ------------------------------------------------------------------

set.seed(123)

huc_level <- "huc08"
huc_col <- str_replace(huc_level, "huc0", "huc")

n_iter <- 100          # increase to 1000 for manuscript-grade null tests
min_reaches <- 30
coherence_r2 <- 0.8
n_workers <- max(1, parallel::detectCores() - 2)

path_flowlines <- "data/processed/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp"
path_metabolism <- "data/processed/05b_centroid_mismatch/05b_maavara_monthly_metabolism_centroid_filtered.rds"
path_removed_comids <- "data/processed/05b_centroid_mismatch/05b_removed_comids.rds"
path_upstream_lookup <- "data/processed/04_upstream_lookup.rds"
path_observed_metrics <- paste0("data/processed/06_", huc_level, "_allometry_metrics.csv")

dir_06a <- "data/processed/06a_null_coherence"
dir.create(dir_06a, recursive = TRUE, showWarnings = FALSE)


# 1. Read base data -------------------------------------------------------------

flowlines_sf <- read_sf(path_flowlines) |>
  clean_names() |>
  mutate(
    comid = as.character(comid),
    nxt_dw = as.character(nxt_dw)
  )

removed_comids <- readRDS(path_removed_comids) |>
  mutate(comid = as.character(comid)) |>
  distinct(comid)

metabolism <- readRDS(path_metabolism) |>
  mutate(comid = as.character(comid))

comid_characteristics <- read_csv(
  "data/dryad_data/scaling_dataset_Qmonth_06082024.csv",
  show_col_types = FALSE
) |>
  clean_names() |>
  select(comid, unitarea, lengthkm, merit_width, order, merit_qmean) |>
  mutate(
    comid = as.character(comid),
    stream_area_m2 = lengthkm * 1000 * merit_width
  ) |>
  anti_join(removed_comids, by = "comid")

metabolism_scaled <- metabolism |>
  left_join(comid_characteristics, by = "comid") |>
  transmute(
    comid,
    month,
    gpp_flux_gO2_d = gpp * stream_area_m2,
    er_flux_gO2_d = er * stream_area_m2
  )


# 2. Create filtered HUC8 lookup ------------------------------------------------

huc_boundaries <- nhdplusTools::get_huc(
  AOI = st_as_sfc(st_bbox(flowlines_sf)),
  type = huc_level,
  t_srs = st_crs(flowlines_sf)
) |>
  clean_names() |>
  select(all_of(huc_col))

flowline_hucs <- flowlines_sf |>
  anti_join(removed_comids, by = "comid") |>
  st_point_on_surface() |>
  st_join(huc_boundaries, join = st_within) |>
  st_drop_geometry() |>
  mutate(comid = as.character(comid)) |>
  select(comid, huc = all_of(huc_col), uparea) |>
  filter(!is.na(huc))


# 3. Prepare within-HUC upstream lookup ----------------------------------------

huc_comid_lookup <- flowline_hucs |>
  select(comid, huc)

upstream_lookup_huc <- readRDS(path_upstream_lookup) |>
  mutate(
    downstream_comid = as.character(downstream_comid),
    upstream_comid = as.character(upstream_comid)
  ) |>
  anti_join(removed_comids, by = c("downstream_comid" = "comid")) |>
  anti_join(removed_comids, by = c("upstream_comid" = "comid")) |>
  inner_join(
    huc_comid_lookup,
    by = c("downstream_comid" = "comid")
  ) |>
  rename(huc_downstream = huc) |>
  inner_join(
    huc_comid_lookup,
    by = c("upstream_comid" = "comid")
  ) |>
  rename(huc_upstream = huc) |>
  filter(huc_downstream == huc_upstream) |>
  transmute(
    huc = huc_downstream,
    downstream_comid,
    upstream_comid
  )

hucs <- upstream_lookup_huc |>
  distinct(huc) |>
  arrange(huc) |>
  pull(huc)


# 4. Define null allometry function --------------------------------------------

fit_allometry <- function(cumulative_metabolism) {
  
  cumulative_metabolism |>
    select(downstream_comid, month, uparea, c_gpp_flux_gO2_d, c_er_flux_gO2_d) |>
    pivot_longer(
      c(c_gpp_flux_gO2_d, c_er_flux_gO2_d),
      names_to = "metabolism_var",
      values_to = "cumulative_metabolism"
    ) |>
    filter(
      uparea > 0,
      cumulative_metabolism > 0,
      !is.na(uparea),
      !is.na(cumulative_metabolism)
    ) |>
    group_by(month, metabolism_var) |>
    summarise(
      n = n(),
      fit = list(if (n() >= 3) lm(log10(cumulative_metabolism) ~ log10(uparea)) else NULL),
      .groups = "drop"
    ) |>
    mutate(
      intercept = map_dbl(fit, \(x) if (is.null(x)) NA_real_ else coef(x)[1]),
      slope = map_dbl(fit, \(x) if (is.null(x)) NA_real_ else coef(x)[2]),
      r2 = map_dbl(fit, \(x) if (is.null(x)) NA_real_ else summary(x)$r.squared)
    ) |>
    select(-fit)
}


calculate_null_allometry_by_huc <- function(selected_huc) {
  
  huc_comids <- flowline_hucs |>
    filter(huc == selected_huc) |>
    select(comid, uparea)
  
  huc_upstream_lookup <- upstream_lookup_huc |>
    filter(huc == selected_huc) |>
    select(downstream_comid, upstream_comid)
  
  local_metabolism <- metabolism_scaled |>
    semi_join(huc_comids, by = "comid")
  
  if (n_distinct(local_metabolism$comid) < min_reaches) {
    return(tibble())
  }
  
  map_dfr(seq_len(n_iter), \(i) {
    
    null_metabolism <- local_metabolism |>
      group_by(month) |>
      mutate(
        gpp_flux_gO2_d = sample(gpp_flux_gO2_d),
        er_flux_gO2_d = sample(er_flux_gO2_d)
      ) |>
      ungroup()
    
    cumulative_null <- huc_upstream_lookup |>
      left_join(
        null_metabolism |>
          select(
            upstream_comid = comid,
            month,
            gpp_flux_gO2_d,
            er_flux_gO2_d
          ),
        by = "upstream_comid",
        relationship = "many-to-many"
      ) |>
      group_by(downstream_comid, month) |>
      summarise(
        c_gpp_flux_gO2_d = sum(gpp_flux_gO2_d, na.rm = TRUE),
        c_er_flux_gO2_d = sum(er_flux_gO2_d, na.rm = TRUE),
        .groups = "drop"
      ) |>
      left_join(
        huc_comids,
        by = c("downstream_comid" = "comid")
      )
    
    fit_allometry(cumulative_null) |>
      mutate(
        huc_level = huc_level,
        huc = selected_huc,
        iter = i
      )
  })
}


# 5. Run null allometries -------------------------------------------------------

plan(multisession, workers = n_workers)
handlers(global = TRUE)

with_progress({
  
  p <- progressor(along = hucs)
  
  null_allometry <- hucs |>
    future_map_dfr(\(selected_huc) {
      
      p(message = selected_huc)
      calculate_null_allometry_by_huc(selected_huc)
      
    }, .options = furrr_options(seed = TRUE))
})

plan(sequential)

write_rds(
  null_allometry,
  file.path(dir_06a, paste0("06a_", huc_level, "_null_allometry_metrics.rds"))
)


# 6. Summarise null distribution ------------------------------------------------

quantile_or_na <- function(x, probs) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    quantile(x, probs = probs, na.rm = TRUE)
  }
}

null_summary <- null_allometry |>
  group_by(huc_level, huc, month, metabolism_var) |>
  summarise(
    n_null = sum(!is.na(r2)),
    r2_null_median = median(r2, na.rm = TRUE),
    r2_null_95 = quantile_or_na(r2, 0.95),
    slope_null_median = median(slope, na.rm = TRUE),
    slope_null_025 = quantile_or_na(slope, 0.025),
    slope_null_975 = quantile_or_na(slope, 0.975),
    .groups = "drop"
  )

write_rds(
  null_summary,
  file.path(dir_06a, paste0("06a_", huc_level, "_null_summary.rds"))
)

write_csv(
  null_summary,
  file.path(dir_06a, paste0("06a_", huc_level, "_null_summary.csv"))
)


# 7. Compare observed coherence against null -----------------------------------

observed_metrics <- read_csv(path_observed_metrics, show_col_types = FALSE) |>
  mutate(huc = as.character(huc))

coherence_test <- observed_metrics |>
  left_join(
    null_summary,
    by = c("huc_level", "huc", "month", "metabolism_var")
  ) |>
  mutate(
    observed_coherent = n >= min_reaches & r2 >= coherence_r2,
    stronger_than_null = n >= min_reaches & r2 > r2_null_95,
    rcc_like_coherent = observed_coherent & stronger_than_null,
    r2_minus_null_95 = r2 - r2_null_95
  )

write_rds(
  coherence_test,
  file.path(dir_06a, paste0("06a_", huc_level, "_observed_vs_null_coherence.rds"))
)

write_csv(
  coherence_test,
  file.path(dir_06a, paste0("06a_", huc_level, "_observed_vs_null_coherence.csv"))
)


# 8. Quick summary --------------------------------------------------------------

coherence_summary <- coherence_test |>
  mutate(
    metabolism = recode(
      metabolism_var,
      c_gpp_flux_gO2_d = "GPP",
      c_er_flux_gO2_d = "ER"
    )
  ) |>
  group_by(metabolism) |>
  summarise(
    n_tests = n(),
    n_observed_coherent = sum(observed_coherent, na.rm = TRUE),
    n_stronger_than_null = sum(stronger_than_null, na.rm = TRUE),
    n_rcc_like_coherent = sum(rcc_like_coherent, na.rm = TRUE),
    pct_observed_coherent = n_observed_coherent / n_tests * 100,
    pct_rcc_like_coherent = n_rcc_like_coherent / n_tests * 100,
    .groups = "drop"
  )

write_csv(
  coherence_summary,
  file.path(dir_06a, paste0("06a_", huc_level, "_coherence_summary.csv"))
)

coherence_summary