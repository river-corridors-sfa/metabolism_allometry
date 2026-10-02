## 10_null_model_huc8.R
##
## PURPOSE
## Test whether observed HUC8 allometric scaling depends on the spatial organization
## of reach-scale metabolism, rather than arising from network accumulation alone.
##
## NULL MODEL
## Within each HUC8 x month x upstream-area quartile, independently shuffle the
## Maavara areal GPP and ER rates among reaches. We preserve:
##   - the broad relationship between metabolism and upstream watershed area,
##   - the observed distribution of metabolism rates within each HUC8-month,
##   - missing-value locations,
##   - reach geometry / stream surface area,
##   - drainage area,
##   - routing topology,
##   - month and HUC membership.
##
## We destroy finer-scale correspondence between a reach's metabolism rate and its
## position in the river network. After shuffling, each rate is multiplied by the
## destination reach's real stream surface area, cumulative upstream metabolism is
## recalculated, and the exact HUC8 allometry rules from script 07 are reapplied.
##
## PRE-SPECIFIED NULL TEST
## For each HUC8 x month x metabolism fit:
##   - slope: observed value must fall outside the central 95% null interval
##            (two-sided test)
##   - R2:    observed value must exceed the 97.5th percentile of the null
## A fit is "null-exceeding" only when BOTH conditions are met.
##
## Initial interpretation rule, evaluated separately for GPP and ER:
##   GO        >= 75% of observed fits are null-exceeding
##   AMBIGUOUS  > 25% and < 75%
##   NO        <= 25%
##
## NOTE
## Script 07 contains a pct_metabolism_missing reporting field that is not used in
## filtering, regression fitting, slope calculation, R2, or class assignment.
## It is intentionally omitted here because it does not affect the null statistics.


# 0. Setup ---------------------------------------------------------------------

source("scripts/00_setup.R")

tic("run 10")

pacman::p_load(progress)


# Paths
dir_processed <- "data/processed"
dir_merit <- file.path(dir_processed, "merit_data")
dir_figs <- "figures"

dir.create(dir_figs, showWarnings = FALSE, recursive = TRUE)

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

path_null_results <- file.path(
  dir_processed,
  "10_huc08_null_permutations.rds"
)

path_null_comparison <- file.path(
  dir_processed,
  "10_huc08_null_comparison.csv"
)

path_decision_summary <- file.path(
  dir_processed,
  "10_huc08_null_decision_summary.csv"
)

path_category_matrix <- file.path(
  dir_processed,
  "10_huc08_null_category_matrix.csv"
)


# Null-model settings
n_permutations <- 50
seed <- 4107


# Exact allometry rules carried forward from script 07
missing_area_threshold <- 10
min_reaches <- 20
min_area_range_log10 <- 1

r2_threshold <- 0.8
linear_lower <- 0.95
linear_upper <- 1.05


# Pre-specified null thresholds
null_lower <- 0.025
null_upper <- 0.975
go_threshold <- 0.75
no_threshold <- 0.25


# 1. Read static inputs ---------------------------------------------------------

message("1/7 Reading spatial, routing, metabolism, and observed allometry data...")

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
    observed_slope = slope,
    observed_r2 = r2,
    observed_class = slope_class_simple
  )


# Script 01 currently writes these names after clean_names(). This fallback keeps
# script 10 compatible with the alternate raw names handled in script 07.
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


# 2. Rebuild the static reach attributes used in script 07 ---------------------

message("2/7 Building the static reach and drainage-area tables...")


# WHAT: Recreate the same reach surface area used to convert Maavara's areal
# metabolism rates to total metabolism for each reach.
#
# WHY: The null shuffles only the areal rate. Each shuffled rate is then applied
# to the REAL surface area of the reach where it lands, preserving river geometry.
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
      select(comid, lengthkm, merit_width, order),
    by = "comid"
  ) %>%
  mutate(
    stream_area_m2 = lengthkm * 1000 * merit_width
  )


# WHAT: Calculate cumulative retained drainage area once.
#
# WHY: Drainage area and routing never change under the null, so recalculating
# this inside every permutation would waste substantial time.
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
    cumulative_area_km2 = sum(upstream_unitarea, na.rm = TRUE),
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
    area_ratio_retained = cumulative_area_km2 / uparea,
    pct_area_missing = pmax(0, 1 - area_ratio_retained) * 100
  )


# WHAT: Create one row per reach-month and divide reaches into upstream-area
# quartiles within each HUC8.
#
# WHY: Rates will only be shuffled among reaches in the same area quartile.
# This preserves broad metabolism-vs-network-size structure while destroying
# finer spatial organization within those size classes.
local_base <- metabolism %>%
  inner_join(
    reach_attributes %>%
      select(
        comid,
        huc,
        unitarea,
        uparea,
        stream_area_m2
      ),
    by = "comid"
  ) %>%
  group_by(huc) %>%
  mutate(
    area_bin = ntile(uparea, 4)
  ) %>%
  ungroup()


# Optional testing subset.
# Uncomment slice_sample() to run a random subset before running all HUC8s.
test_hucs <- local_base %>%
  distinct(huc) %>%
  # slice_sample(n = 100) %>%   # COMMENT OUT FOR FULL RUN
  pull(huc)


local_base <- local_base %>%
  filter(huc %in% test_hucs)

retained_upstream_area <- retained_upstream_area %>%
  filter(huc %in% test_hucs)

observed_results <- observed_results %>%
  filter(huc %in% test_hucs)


# Restrict the large routing lookup to downstream reaches that are actually
# included in this run. This saves memory and time, especially for test subsets.
upstream_lookup <- upstream_lookup %>%
  semi_join(
    retained_upstream_area %>%
      select(downstream_comid),
    by = "downstream_comid"
  )


# 3. Define small, readable helper functions -----------------------------------


# Shuffle only observed values and leave NA locations unchanged.
#
# WHY: We want to randomize spatial metabolism organization, not spatial patterns
# in data availability. This also keeps every permutation directly comparable.
shuffle_observed <- function(x) {
  
  keep <- which(!is.na(x))
  
  x[keep] <- sample(
    x[keep],
    length(keep),
    replace = FALSE
  )
  
  x
}


# Apply the exact slope / R2 classification used in script 07.
classify_allometry <- function(slope, r2) {
  
  case_when(
    r2 < r2_threshold ~ "low_r2",
    slope < linear_lower ~ "sublinear",
    slope <= linear_upper ~ "near_linear",
    slope > linear_upper ~ "superlinear",
    TRUE ~ NA_character_
  )
}


# Fit the same HUC8 x month allometric regressions used in script 07.
fit_allometry <- function(allometry_data) {
  
  allometry_data %>%
    group_by(huc, month, metabolism_metric) %>%
    filter(
      n() >= min_reaches,
      n_distinct(cumulative_area_km2) >= min_reaches,
      log10(
        max(cumulative_area_km2, na.rm = TRUE) /
          min(cumulative_area_km2, na.rm = TRUE)
      ) >= min_area_range_log10
    ) %>%
    group_modify(~ {
      
      fit <- lm(
        log10(cumulative_flux_g_o2_d) ~ log10(cumulative_area_km2),
        data = .x
      )
      
      tibble(
        slope = coef(fit)[2],
        r2 = summary(fit)$r.squared
      )
      
    }) %>%
    ungroup() %>%
    mutate(
      slope_class_simple = classify_allometry(slope, r2)
    )
}


# Run one complete null realization.
run_null_permutation <- function(permutation_id) {
  
  # Give every permutation its own reproducible random-number sequence.
  set.seed(seed + permutation_id)
  
  
  # WHAT: Shuffle GPP and ER rates independently within
  # HUC8 x month x upstream-area quartile.
  #
  # WHY: This preserves seasonality and broad upstream-area dependence while
  # removing finer spatial placement of high and low metabolic rates.
  local_null <- local_base %>%
    group_by(huc, month, area_bin) %>%
    mutate(
      gpp_rate_null = shuffle_observed(gpp_rate_g_o2_m2_d),
      er_rate_null = shuffle_observed(er_rate_g_o2_m2_d)
    ) %>%
    ungroup() %>%
    mutate(
      gpp_flux_g_o2_d = gpp_rate_null * stream_area_m2,
      er_flux_g_o2_d = er_rate_null * stream_area_m2
    )
  
  
  # WHAT: Route and sum the shuffled local fluxes through the REAL network.
  #
  # WHY: Network topology is part of the null background we want to preserve.
  # If scaling remains after spatial randomization, routing/accumulation may be
  # sufficient to reproduce that portion of the observed scaling pattern.
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
    group_by(downstream_comid, month) %>%
    summarise(
      
      cumulative_gpp_g_o2_d = if_else(
        any(!is.na(gpp_flux_g_o2_d)),
        sum(gpp_flux_g_o2_d, na.rm = TRUE),
        NA_real_
      ),
      
      cumulative_er_g_o2_d = if_else(
        any(!is.na(er_flux_g_o2_d)),
        sum(er_flux_g_o2_d, na.rm = TRUE),
        NA_real_
      ),
      
      .groups = "drop"
    ) %>%
    left_join(
      retained_upstream_area,
      by = "downstream_comid"
    )
  
  
  # Apply the same positivity and retained-area filters as script 07.
  allometry_null <- cumulative_null %>%
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
      pct_area_missing < missing_area_threshold
    )
  
  
  fit_allometry(allometry_null) %>%
    mutate(
      permutation = permutation_id,
      .before = 1
    )
}


# 4. Run the null permutations --------------------------------------------------

message(
  "3/7 Running ",
  n_permutations,
  " null permutations sequentially..."
)

message(
  "This is the computationally expensive step. ",
  "Progress will update after each complete permutation."
)


pb <- progress_bar$new(
  total = n_permutations,
  format = paste0(
    "  Null model [:bar] :percent | ",
    ":current/:total | elapsed :elapsed | ETA :eta"
  ),
  clear = FALSE,
  width = 80
)


null_results <- vector(
  "list",
  n_permutations
)


for (i in seq_len(n_permutations)) {
  
  null_results[[i]] <- run_null_permutation(i)
  
  pb$tick()
}


null_results <- bind_rows(null_results)


saveRDS(
  null_results,
  path_null_results
)


message(
  "Null permutations complete and saved: ",
  path_null_results
)


# 5. Compare each observed fit with its null distribution ----------------------

message(
  "4/7 Comparing observed slope and R2 with their null distributions..."
)


null_summary <- null_results %>%
  group_by(
    huc,
    month,
    metabolism_metric
  ) %>%
  summarise(
    
    n_null = n(),
    
    null_slope_median = median(
      slope,
      na.rm = TRUE
    ),
    
    null_slope_q025 = quantile(
      slope,
      null_lower,
      na.rm = TRUE
    ),
    
    null_slope_q975 = quantile(
      slope,
      null_upper,
      na.rm = TRUE
    ),
    
    null_r2_median = median(
      r2,
      na.rm = TRUE
    ),
    
    null_r2_q025 = quantile(
      r2,
      null_lower,
      na.rm = TRUE
    ),
    
    null_r2_q975 = quantile(
      r2,
      null_upper,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


comparison <- observed_results %>%
  inner_join(
    null_summary,
    by = c(
      "huc",
      "month",
      "metabolism_metric"
    )
  ) %>%
  mutate(
    
    # Slope is deliberately two-sided: either unusually low OR unusually high
    # relative to the shuffled network counts as departure from the null.
    slope_exceeds_null =
      observed_slope < null_slope_q025 |
      observed_slope > null_slope_q975,
    
    # R2 is one-sided because the hypothesis is that observed spatial
    # organization produces a stronger relationship than random placement.
    r2_exceeds_null =
      observed_r2 > null_r2_q975,
    
    # This is the pre-specified primary per-fit criterion.
    null_exceeding =
      slope_exceeds_null &
      r2_exceeds_null
  )


# Add empirical permutation p-values / percentiles for transparent reporting.
# The +1 correction avoids zero p-values with a finite number of permutations.
null_empirical <- null_results %>%
  inner_join(
    observed_results,
    by = c(
      "huc",
      "month",
      "metabolism_metric"
    )
  ) %>%
  group_by(
    huc,
    month,
    metabolism_metric
  ) %>%
  summarise(
    
    slope_null_percentile =
      mean(
        slope <= observed_slope,
        na.rm = TRUE
      ),
    
    r2_null_percentile =
      mean(
        r2 <= observed_r2,
        na.rm = TRUE
      ),
    
    slope_p_two_sided = pmin(
      1,
      2 * pmin(
        (1 + sum(
          slope <= observed_slope,
          na.rm = TRUE
        )) / (n() + 1),
        
        (1 + sum(
          slope >= observed_slope,
          na.rm = TRUE
        )) / (n() + 1)
      )
    ),
    
    r2_p_upper =
      (1 + sum(
        r2 >= observed_r2,
        na.rm = TRUE
      )) / (n() + 1),
    
    same_category_probability =
      mean(
        slope_class_simple == observed_class,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


comparison <- comparison %>%
  left_join(
    null_empirical,
    by = c(
      "huc",
      "month",
      "metabolism_metric"
    )
  )


write_csv(
  comparison,
  path_null_comparison
)


# 6. Apply the pre-specified GO / AMBIGUOUS / NO rule --------------------------

message(
  "5/7 Summarizing the pre-specified decision rule..."
)


decision_summary <- comparison %>%
  group_by(
    metabolism_metric
  ) %>%
  summarise(
    
    n_observed_fits = n(),
    
    n_null_exceeding =
      sum(
        null_exceeding,
        na.rm = TRUE
      ),
    
    proportion_null_exceeding =
      mean(
        null_exceeding,
        na.rm = TRUE
      ),
    
    median_same_category_probability =
      median(
        same_category_probability,
        na.rm = TRUE
      ),
    
    decision = case_when(
      proportion_null_exceeding >= go_threshold ~ "GO",
      proportion_null_exceeding <= no_threshold ~ "NO",
      TRUE ~ "AMBIGUOUS"
    ),
    
    .groups = "drop"
  )


write_csv(
  decision_summary,
  path_decision_summary
)


print(
  decision_summary
)


# 7. Category similarity / dissimilarity matrix --------------------------------

message(
  "6/7 Building category-reproduction matrix..."
)


category_matrix <- null_results %>%
  inner_join(
    observed_results,
    by = c(
      "huc",
      "month",
      "metabolism_metric"
    )
  ) %>%
  count(
    metabolism_metric,
    observed_class,
    null_class = slope_class_simple,
    name = "n"
  ) %>%
  group_by(
    metabolism_metric,
    observed_class
  ) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  ungroup()


write_csv(
  category_matrix,
  path_category_matrix
)


# 8. Diagnostic figures --------------------------------------------------------

message(
  "7/7 Saving diagnostic figures..."
)


# Figure 1: direct summary of the pre-specified decision statistic.
p_decision <- decision_summary %>%
  ggplot(
    aes(
      x = metabolism_metric,
      y = proportion_null_exceeding
    )
  ) +
  geom_col(
    width = 0.65
  ) +
  geom_hline(
    yintercept = c(
      no_threshold,
      go_threshold
    ),
    linetype = "dashed"
  ) +
  geom_text(
    aes(
      label = paste0(
        scales::percent(
          proportion_null_exceeding,
          accuracy = 1
        ),
        "\n",
        decision
      )
    ),
    vjust = -0.3
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = "Observed fits exceeding null",
    title = "Null-model decision summary",
    subtitle = paste0(
      "Fit exceeds null only when slope AND R2 ",
      "meet pre-specified criteria"
    )
  ) +
  theme_bw()


ggsave(
  file.path(
    dir_figs,
    "10_null_decision_summary.png"
  ),
  p_decision,
  width = 7,
  height = 5,
  dpi = 300
)


# Figure 2: how often each observed class is reproduced after spatial shuffling.
# A strong diagonal means the null often recreates the observed classifications;
# weak diagonal agreement means classifications depend on spatial organization.
p_category <- category_matrix %>%
  ggplot(
    aes(
      x = null_class,
      y = observed_class,
      fill = proportion
    )
  ) +
  geom_tile() +
  geom_text(
    aes(
      label = scales::percent(
        proportion,
        accuracy = 1
      )
    )
  ) +
  scale_fill_viridis_c(
    labels = scales::percent,
    limits = c(0, 1),
    name = "Null\nproportion"
  ) +
  facet_wrap(
    ~metabolism_metric
  ) +
  labs(
    x = "Null-model class",
    y = "Observed class",
    title = "Observed vs null allometry classifications"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )


ggsave(
  file.path(
    dir_figs,
    "10_null_category_matrix.png"
  ),
  p_category,
  width = 9,
  height = 5,
  dpi = 300
)


# Figure 3: where observed slope and R2 fall within their null distributions.
#
# For slope, values near 0 or 1 are unusual under the two-sided null.
# For R2, values near 1 indicate stronger fits than expected under the null.
percentile_plot_data <- comparison %>%
  select(
    huc,
    month,
    metabolism_metric,
    slope_null_percentile,
    r2_null_percentile
  ) %>%
  pivot_longer(
    c(
      slope_null_percentile,
      r2_null_percentile
    ),
    names_to = "metric",
    values_to = "null_percentile"
  ) %>%
  mutate(
    metric = recode(
      metric,
      slope_null_percentile = "Slope",
      r2_null_percentile = "R2"
    )
  )


p_percentiles <- percentile_plot_data %>%
  ggplot(
    aes(
      x = null_percentile
    )
  ) +
  geom_histogram(
    binwidth = 0.05,
    boundary = 0
  ) +
  geom_vline(
    data = tibble(
      metric = c(
        "Slope",
        "Slope",
        "R2"
      ),
      cutoff = c(
        null_lower,
        null_upper,
        null_upper
      )
    ),
    aes(
      xintercept = cutoff
    ),
    linetype = "dashed"
  ) +
  facet_grid(
    metric ~ metabolism_metric,
    scales = "free_y"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  labs(
    x = "Observed value percentile within null distribution",
    y = "Number of HUC8-month fits",
    title = "Observed slope and R2 relative to shuffled null distributions"
  ) +
  theme_bw()


ggsave(
  file.path(
    dir_figs,
    "10_null_metric_percentiles.png"
  ),
  p_percentiles,
  width = 9,
  height = 6,
  dpi = 300
)


message("Done.")

message(
  "Full permutations: ",
  path_null_results
)

message(
  "Observed-null comparison: ",
  path_null_comparison
)

message(
  "Decision summary: ",
  path_decision_summary
)

message(
  "Category matrix: ",
  path_category_matrix
)

toc()