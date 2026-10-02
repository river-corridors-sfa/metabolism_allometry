## 14_null_sensitivity.R
##
## PURPOSE
## Test whether conclusions from the HUC8 null model depend strongly on the
## choice of upstream-area binning used during shuffling.
##
## The null is rerun using 4, 10, and 20 upstream-area bins. If delta-beta is
## stable across these alternatives, the interpretation is robust to reasonable
## changes in how strongly the null preserves metabolism-size structure.


# 0. Setup ---------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")

tic("run 14") #10989s OR 3 HOURS... that's a NEW state record

pacman::p_load(progress)

dir_processed <- "data/processed"
dir_merit <- file.path(dir_processed, "merit_data")
dir_figs <- "figures"

path_cat <- file.path(
  dir_merit,
  "03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
)

path_riv <- file.path(
  dir_merit,
  "03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
)

path_upstream_lookup <- file.path(
  dir_merit,
  "05_upstream_lookup.rds"
)

path_metabolism <- file.path(
  dir_processed,
  "01_assembled_maavara_monthly_metabolism.rds"
)

path_characteristics <-
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv"

path_observed_results <- file.path(
  dir_processed,
  "07_huc08_allometry_results.csv"
)

path_output <- file.path(
  dir_processed,
  "14_null_area_bin_sensitivity.rds"
)

area_bins <- c(4, 10, 20)
n_permutations <- 20
seed <- 4107

missing_area_threshold <- 10
min_reaches <- 20
min_area_range_log10 <- 1

r2_threshold <- 0.8
linear_lower <- 0.95
linear_upper <- 1.05


# 1. Read inputs ---------------------------------------------------------------

merit_cat <- read_sf(path_cat) %>%
  clean_names() %>%
  mutate(comid = as.character(comid))

merit_riv <- read_sf(path_riv) %>%
  clean_names() %>%
  mutate(
    comid = as.character(comid),
    huc = as.character(huc)
  )

upstream_lookup <- readRDS(path_upstream_lookup) %>%
  mutate(
    downstream_comid = as.character(downstream_comid),
    upstream_comid = as.character(upstream_comid)
  )

metabolism <- readRDS(path_metabolism) %>%
  clean_names() %>%
  mutate(comid = as.character(comid))

comid_characteristics <- read_csv(
  path_characteristics,
  show_col_types = FALSE
) %>%
  clean_names() %>%
  mutate(comid = as.character(comid))

observed_results <- read_csv(
  path_observed_results,
  show_col_types = FALSE
) %>%
  select(
    huc,
    month,
    metabolism_metric,
    observed_slope = slope
  )


# 2. Prepare metabolism and reach attributes -----------------------------------

if ("gpp_g_o2_m2_d" %in% names(metabolism)) {
  
  metabolism <- metabolism %>%
    rename(
      gpp_rate_g_o2_m2_d = gpp_g_o2_m2_d,
      er_rate_g_o2_m2_d = er_g_o2_m2_d
    )
  
} else {
  
  metabolism <- metabolism %>%
    rename(
      gpp_rate_g_o2_m2_d = gpp,
      er_rate_g_o2_m2_d = er
    )
}

metabolism <- metabolism %>%
  select(
    comid,
    month,
    gpp_rate_g_o2_m2_d,
    er_rate_g_o2_m2_d
  )

reach_attributes <- merit_cat %>%
  st_drop_geometry() %>%
  select(comid, huc, unitarea) %>%
  left_join(
    merit_riv %>%
      st_drop_geometry() %>%
      select(comid, uparea),
    by = "comid"
  ) %>%
  left_join(
    comid_characteristics %>%
      select(comid, lengthkm, merit_width),
    by = "comid"
  ) %>%
  mutate(
    stream_area_m2 = lengthkm * 1000 * merit_width
  )

retained_upstream_area <- upstream_lookup %>%
  inner_join(
    reach_attributes %>%
      select(
        upstream_comid = comid,
        upstream_unitarea = unitarea
      ),
    by = "upstream_comid"
  ) %>%
  group_by(downstream_comid) %>%
  summarise(
    cumulative_area_km2 = sum(
      upstream_unitarea,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    reach_attributes %>%
      select(
        downstream_comid = comid,
        huc,
        uparea
      ),
    by = "downstream_comid"
  ) %>%
  mutate(
    pct_area_missing =
      pmax(
        0,
        1 - cumulative_area_km2 / uparea
      ) * 100
  )

local_base <- metabolism %>%
  inner_join(
    reach_attributes %>%
      select(
        comid,
        huc,
        uparea,
        stream_area_m2
      ),
    by = "comid"
  )

upstream_lookup <- upstream_lookup %>%
  semi_join(
    retained_upstream_area %>%
      select(downstream_comid),
    by = "downstream_comid"
  )


# 3. Helper functions -----------------------------------------------------------

shuffle_observed <- function(x) {
  
  keep <- which(!is.na(x))
  
  x[keep] <- sample(
    x[keep],
    length(keep),
    replace = FALSE
  )
  
  x
}


classify_allometry <- function(slope, r2) {
  
  case_when(
    r2 < r2_threshold ~ "low_r2",
    slope < linear_lower ~ "sublinear",
    slope <= linear_upper ~ "near_linear",
    slope > linear_upper ~ "superlinear",
    TRUE ~ NA_character_
  )
}


fit_allometry <- function(data) {
  
  data %>%
    group_by(
      huc,
      month,
      metabolism_metric
    ) %>%
    filter(
      n() >= min_reaches,
      n_distinct(cumulative_area_km2) >= min_reaches,
      log10(
        max(cumulative_area_km2) /
          min(cumulative_area_km2)
      ) >= min_area_range_log10
    ) %>%
    group_modify(~ {
      
      fit <- lm(
        log10(cumulative_flux_g_o2_d) ~
          log10(cumulative_area_km2),
        data = .x
      )
      
      tibble(
        slope = coef(fit)[2],
        r2 = summary(fit)$r.squared
      )
    }) %>%
    ungroup() %>%
    mutate(
      slope_class_simple =
        classify_allometry(
          slope,
          r2
        )
    )
}


run_null_permutation <- function(permutation_id, n_bins) {
  
  set.seed(
    seed +
      permutation_id +
      n_bins * 1000
  )
  
  local_null <- local_base %>%
    group_by(huc) %>%
    mutate(
      area_bin = ntile(
        uparea,
        n_bins
      )
    ) %>%
    group_by(
      huc,
      month,
      area_bin
    ) %>%
    mutate(
      gpp_rate_null =
        shuffle_observed(
          gpp_rate_g_o2_m2_d
        ),
      er_rate_null =
        shuffle_observed(
          er_rate_g_o2_m2_d
        )
    ) %>%
    ungroup() %>%
    mutate(
      gpp_flux_g_o2_d =
        gpp_rate_null *
        stream_area_m2,
      er_flux_g_o2_d =
        er_rate_null *
        stream_area_m2
    )
  
  cumulative_null <- upstream_lookup %>%
    inner_join(
      local_null %>%
        select(
          upstream_comid = comid,
          month,
          gpp_flux_g_o2_d,
          er_flux_g_o2_d
        ),
      by = "upstream_comid",
      relationship = "many-to-many"
    ) %>%
    group_by(
      downstream_comid,
      month
    ) %>%
    summarise(
      cumulative_gpp_g_o2_d =
        if_else(
          any(!is.na(gpp_flux_g_o2_d)),
          sum(
            gpp_flux_g_o2_d,
            na.rm = TRUE
          ),
          NA_real_
        ),
      cumulative_er_g_o2_d =
        if_else(
          any(!is.na(er_flux_g_o2_d)),
          sum(
            er_flux_g_o2_d,
            na.rm = TRUE
          ),
          NA_real_
        ),
      .groups = "drop"
    ) %>%
    left_join(
      retained_upstream_area,
      by = "downstream_comid"
    ) %>%
    pivot_longer(
      c(
        cumulative_gpp_g_o2_d,
        cumulative_er_g_o2_d
      ),
      names_to = "metabolism_metric",
      values_to = "cumulative_flux_g_o2_d"
    ) %>%
    mutate(
      metabolism_metric = recode(
        metabolism_metric,
        cumulative_gpp_g_o2_d = "GPP",
        cumulative_er_g_o2_d = "ER"
      )
    ) %>%
    filter(
      !is.na(huc),
      !is.na(month),
      cumulative_area_km2 > 0,
      cumulative_flux_g_o2_d > 0,
      pct_area_missing <
        missing_area_threshold
    )
  
  fit_allometry(cumulative_null) %>%
    mutate(
      area_bins = n_bins,
      permutation = permutation_id,
      .before = 1
    )
}


# 4. Run sensitivity nulls ------------------------------------------------------

runs <- expand_grid(
  area_bins = area_bins,
  permutation = seq_len(
    n_permutations
  )
)

pb <- progress_bar$new(
  total = nrow(runs),
  format = paste0(
    "  Null sensitivity [:bar] :percent | ",
    ":current/:total | elapsed :elapsed | ETA :eta"
  ),
  clear = FALSE,
  width = 80
)

null_sensitivity <- map2_dfr(
  runs$permutation,
  runs$area_bins,
  ~ {
    
    out <- run_null_permutation(
      permutation_id = .x,
      n_bins = .y
    )
    
    pb$tick()
    
    out
  }
)

saveRDS(
  null_sensitivity,
  path_output
)


# 5. Compare observed slopes with each null -------------------------------------

sensitivity <- null_sensitivity %>%
  group_by(
    area_bins,
    huc,
    month,
    metabolism_metric
  ) %>%
  summarise(
    null_slope = median(
      slope,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  inner_join(
    observed_results,
    by = c(
      "huc",
      "month",
      "metabolism_metric"
    )
  ) %>%
  mutate(
    delta_beta =
      observed_slope -
      null_slope
  )


# 6. Summarize sensitivity ------------------------------------------------------

sensitivity_summary <- sensitivity %>%
  group_by(
    metabolism_metric,
    area_bins
  ) %>%
  summarise(
    n = n(),
    median_delta_beta =
      median(
        delta_beta,
        na.rm = TRUE
      ),
    median_abs_delta_beta =
      median(
        abs(delta_beta),
        na.rm = TRUE
      ),
    pct_positive =
      mean(
        delta_beta > 0,
        na.rm = TRUE
      ) * 100,
    pct_negative =
      mean(
        delta_beta < 0,
        na.rm = TRUE
      ) * 100,
    .groups = "drop"
  )

print(
  sensitivity_summary
)

write_csv(
  sensitivity_summary,
  file.path(
    dir_processed,
    "14_null_area_bin_sensitivity_summary.csv"
  )
)


# 7. Plot sensitivity -----------------------------------------------------------

p <- sensitivity %>%
  mutate(
    area_bins = factor(
      area_bins
    )
  ) %>%
  ggplot(
    aes(
      x = area_bins,
      y = delta_beta
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_boxplot() +
  facet_wrap(
    ~metabolism_metric, scales = "free"
  ) +
  labs(
    x = "Upstream-area bins preserved by null",
    y = expression(
      Delta * beta
    ),
    title =
      "Sensitivity of metabolic departures to null-model structure"
  ) +
  theme_bw()

ggsave(
  file.path(
    dir_figs,
    "14_null_area_bin_sensitivity.png"
  ),
  p,
  width = 7,
  height = 5,
  dpi = 300
)


message("Done.")

toc()