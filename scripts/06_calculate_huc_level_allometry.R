# 06_calculate_huc_level_allometry.R

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(
  furrr,
  future,
  nhdplusTools,
  progressr
)


# 0. Settings ------------------------------------------------------------------

huc_level <- "huc08"   # options: "huc08", "huc10", "huc12"
huc_col <- str_replace(huc_level, "huc0", "huc")

limiter <- NULL        # use e.g., 10 for testing; NULL runs all HUCs
n_workers <- max(1, parallel::detectCores() - 2)

path_flowlines <- "data/processed/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp"
path_metabolism <- "data/processed/05b_centroid_mismatch/05b_maavara_monthly_metabolism_centroid_filtered.rds"
path_removed_comids <- "data/processed/05b_centroid_mismatch/05b_removed_comids.rds"
path_filter_impact <- "data/processed/05b_centroid_mismatch/05b_huc08_centroid_filter_impact.rds"


# 1. Read base datasets --------------------------------------------------------

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

comid_characteristics <- read_csv("data/dryad_data/scaling_dataset_Qmonth_06082024.csv") |>
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

upstream_lookup <- readRDS("data/processed/04_upstream_lookup.rds") |>
  mutate(
    downstream_comid = as.character(downstream_comid),
    upstream_comid = as.character(upstream_comid)
  ) |>
  anti_join(removed_comids, by = c("downstream_comid" = "comid")) |>
  anti_join(removed_comids, by = c("upstream_comid" = "comid"))


# 2. Create filtered HUC lookup ------------------------------------------------

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
  select(comid, all_of(huc_col), uparea)

hucs <- flowline_hucs |>
  filter(!is.na(.data[[huc_col]])) |>
  distinct(.data[[huc_col]]) |>
  arrange(.data[[huc_col]]) |>
  pull(.data[[huc_col]])

if (!is.null(limiter)) {
  hucs <- hucs[1:limiter]
}


# 3. Split upstream lookup by HUC ---------------------------------------------

dir_upstream_huc <- paste0("data/processed/tmp_", huc_level, "_upstream_lookup")

dir.create(
  dir_upstream_huc,
  showWarnings = FALSE,
  recursive = TRUE
)

huc_comid_lookup <- flowline_hucs |>
  filter(!is.na(.data[[huc_col]])) |>
  transmute(
    comid,
    huc = .data[[huc_col]]
  )

upstream_lookup |>
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
  ) |>
  group_by(huc) |>
  group_walk(\(x, y) {
    
    saveRDS(
      x |>
        ungroup() |>
        select(downstream_comid, upstream_comid),
      file.path(dir_upstream_huc, paste0(y$huc, ".rds"))
    )
    
  })

rm(upstream_lookup)
gc()


# 4. Calculate allometry by HUC ------------------------------------------------

calculate_allometric_scaling_by_month_and_huc <- function(selected_huc) {
  
  huc_comids <- flowline_hucs |>
    filter(.data[[huc_col]] == selected_huc) |>
    select(comid, uparea)
  
  huc_lookup_path <- file.path(dir_upstream_huc, paste0(selected_huc, ".rds"))
  
  if (!file.exists(huc_lookup_path)) {
    return(tibble())
  }
  
  huc_upstream_lookup <- readRDS(huc_lookup_path)
  
  local_metabolism <- metabolism_scaled |>
    semi_join(huc_comids |> select(comid), by = "comid")
  
  cumulative_metabolism <- huc_upstream_lookup |>
    left_join(
      local_metabolism |>
        select(upstream_comid = comid, month, gpp_flux_gO2_d, er_flux_gO2_d),
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
      comid_characteristics |> select(comid, unitarea, order, lengthkm),
      by = c("downstream_comid" = "comid")
    ) |>
    left_join(
      huc_comids |> select(comid, uparea),
      by = c("downstream_comid" = "comid")
    )
  
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
      r2 = map_dbl(fit, \(x) if (is.null(x)) NA_real_ else summary(x)$r.squared),
      huc_level = huc_level,
      huc = selected_huc
    ) |>
    select(-fit)
}

plan(multicore, workers = n_workers)
handlers(global = TRUE)

with_progress({
  
  p <- progressor(along = hucs)
  
  huc_allometry <- hucs |>
    future_map_dfr(\(selected_huc) {
      
      p(message = selected_huc)
      calculate_allometric_scaling_by_month_and_huc(selected_huc)
      
    }, .options = furrr_options(seed = TRUE))
})

plan(sequential)

write_csv(
  huc_allometry,
  paste0("data/processed/06_", huc_level, "_allometry_metrics.csv")
)


# 5. HUC diagnostics -----------------------------------------------------------

huc_areas <- huc_boundaries |>
  mutate(huc_area_km2 = as.numeric(st_area(geometry)) / 1e6) |>
  st_drop_geometry() |>
  rename(huc = all_of(huc_col)) |>
  select(huc, huc_area_km2)

filter_impact <- if (file.exists(path_filter_impact)) {
  readRDS(path_filter_impact) |>
    mutate(huc = as.character(huc))
} else {
  tibble(huc = character())
}

calculate_huc_diagnostics <- function(selected_huc) {
  
  huc_comids <- flowline_hucs |>
    filter(.data[[huc_col]] == selected_huc) |>
    select(comid)
  
  huc_lengths <- comid_characteristics |>
    semi_join(huc_comids, by = "comid")
  
  tibble(
    huc_level = huc_level,
    huc = selected_huc,
    n_reaches = n_distinct(huc_comids$comid),
    total_lengthkm = sum(huc_lengths$lengthkm, na.rm = TRUE)
  ) |>
    left_join(huc_areas, by = "huc")
}

plan(multicore, workers = n_workers)

with_progress({
  
  p <- progressor(along = hucs)
  
  huc_diagnostics <- hucs |>
    future_map_dfr(\(selected_huc) {
      
      p(message = selected_huc)
      calculate_huc_diagnostics(selected_huc)
      
    }, .options = furrr_options(seed = TRUE))
})

plan(sequential)

huc_diagnostics <- huc_diagnostics |>
  left_join(
    filter_impact |>
      select(
        huc,
        n_reaches_total_prefilter = n_reaches_total,
        n_reaches_removed,
        pct_reaches_removed,
        length_removed_km,
        pct_length_removed,
        likely_allometry_sensitive
      ),
    by = "huc"
  ) |>
  mutate(
    n_reaches_removed = replace_na(n_reaches_removed, 0),
    pct_reaches_removed = replace_na(pct_reaches_removed, 0),
    length_removed_km = replace_na(length_removed_km, 0),
    pct_length_removed = replace_na(pct_length_removed, 0),
    likely_allometry_sensitive = replace_na(likely_allometry_sensitive, FALSE)
  )

write_csv(
  huc_diagnostics,
  paste0("data/processed/06_", huc_level, "_diagnostics.csv")
)

