# 11_ecoregion_scaling_bars.R
# Run from the project root after 04b and 08. Saves one stacked-bar PNG.
# Reuses 04b's point-based L1 assignments, not dominant polygon coverage.
# Monthly and annual classification rules match script 10.

local({
  library(tidyverse)
  library(sf)

  # 1. Settings ---------------------------------------------------------------
  use_test_output <- FALSE       # TRUE previews 08's three-HUC test only
  min_reaches <- 20L
  min_area_range <- 1
  missing_coverage_cutoff <- 10  # 07 requires <10% for area and metabolism
  slope_bounds <- c(0.95, 1.05)
  r2_cutoff <- 0.8
  dir_processed <- "data/processed"
  dir_figures <- "figures"
  dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)
  path_input <- file.path(dir_processed, if (use_test_output) {
    "08_huc08_fit_stability_test.csv"
  } else "08_huc08_fit_stability.csv")
  if (!file.exists(path_input)) stop("Missing 08 summary: ", path_input,
    ". Run 08 in full mode, or set use_test_output <- TRUE for a preview.")

  classes <- c("Sub-linear", "Linear", "Super-linear", "Uncertain",
               "Insufficient data", "Seasonally variable")
  palette <- c("Sub-linear" = "#0072B2", "Linear" = "#009E73",
               "Super-linear" = "#D55E00", "Uncertain" = "#E5C849",
               "Insufficient data" = "#D6D6D6", "Seasonally variable" = "#9B69AF")

  # 2. Classify each monthly fit ---------------------------------------------
  stability <- read_csv(path_input, show_col_types = FALSE,
                        col_types = cols(huc = col_character())) %>%
    mutate(huc = str_pad(huc, 8, pad = "0"))
  required <- c("huc", "month", "metabolism_metric", "n_reaches", "n_unique_areas",
                "log10_area_range", "n_valid", "baseline_slope", "baseline_r2",
                "slope_min", "slope_max", "r2_min",
                "median_pct_metabolism_missing", "max_pct_metabolism_missing")
  stopifnot(all(required %in% names(stability)),
            all(stability$month %in% 1:12),
            all(stability$metabolism_metric %in% c("GPP", "ER")),
            !anyNA(stability$huc),
            !anyDuplicated(stability[c("huc", "month", "metabolism_metric")]))
  if (any(!is.finite(stability$max_pct_metabolism_missing) |
          stability$max_pct_metabolism_missing >= missing_coverage_cutoff)) {
    stop("08 output violates 07's <10% metabolism coverage rule. Rerun 07 and 08.")
  }

  monthly <- stability %>%
    mutate(
      # Full fit and all VALID deletion refits must meet the same thresholds.
      supported = if_all(all_of(setdiff(required, c("huc", "month", "metabolism_metric"))),
                         ~ is.finite(.x)) &
        n_reaches >= min_reaches & n_unique_areas >= min_reaches &
        log10_area_range >= min_area_range & n_valid > 0,
      slope_low = pmin(baseline_slope, slope_min),
      slope_high = pmax(baseline_slope, slope_max),
      weakest_r2 = pmin(baseline_r2, r2_min),
      monthly_class = case_when(
        !supported ~ "Insufficient data",
        weakest_r2 < r2_cutoff ~ "Uncertain",
        slope_high < slope_bounds[1] ~ "Sub-linear",
        slope_low >= slope_bounds[1] & slope_high <= slope_bounds[2] ~ "Linear",
        slope_low > slope_bounds[2] ~ "Super-linear",
        TRUE ~ "Uncertain"
      ),
      available = TRUE
    ) %>% select(huc, month, metabolism_metric, monthly_class, available)
  rm(stability)

  # 3. Require agreement across all 12 months --------------------------------
  boundaries <- read_sf(file.path(dir_processed, "04c_huc08_boundaries.gpkg")) %>%
    st_drop_geometry() %>%
    transmute(huc = str_pad(as.character(huc8), 8, pad = "0"))
  stopifnot(!anyNA(boundaries$huc), !anyDuplicated(boundaries$huc),
            all(monthly$huc %in% boundaries$huc))
  if (use_test_output) boundaries <- boundaries %>% filter(huc %in% monthly$huc)

  annual <- expand_grid(huc = boundaries$huc, metabolism_metric = c("GPP", "ER"),
                        month = 1:12) %>%
    left_join(monthly, by = c("huc", "metabolism_metric", "month")) %>%
    group_by(huc, metabolism_metric) %>% summarise(
      scaling_class = case_when(
        # A missing month is incomplete coverage, not evidence of seasonality.
        sum(!is.na(available)) < 12 ~ "Insufficient data",
        n_distinct(monthly_class) > 1 ~ "Seasonally variable",
        TRUE ~ first(monthly_class)
      ), .groups = "drop"
    ) %>% mutate(scaling_class = factor(scaling_class, levels = classes))
  rm(monthly, boundaries)


  # 4. Join the existing point-based L1 assignments ----------------------------
  regions <- read_csv(file.path(dir_processed, "04b_huc08_l1_ecoregion_lookup.csv"),
                       col_types = cols(huc = col_character()),
                       show_col_types = FALSE) %>%
    transmute(huc = str_pad(huc, 8, pad = "0"), l1_code, l1_ecoregion)
  stopifnot(!anyNA(regions$huc), !anyDuplicated(regions$huc))
  regional <- annual %>% left_join(regions, by = "huc", relationship = "many-to-one") %>%
    mutate(l1_ecoregion = replace_na(l1_ecoregion, "Unassigned"))
  # Preserve unassigned HUC8s so no study HUC silently leaves the denominator.
  region_order <- regional %>% distinct(l1_code, l1_ecoregion) %>%
    arrange(l1_code, l1_ecoregion) %>% pull(l1_ecoregion) %>% unique()
  bars <- regional %>% count(l1_ecoregion, metabolism_metric, scaling_class,
                             name = "n_hucs") %>%
    group_by(l1_ecoregion, metabolism_metric) %>%
    mutate(total_hucs = sum(n_hucs), fraction = n_hucs / total_hucs) %>%
    ungroup() %>%
    mutate(l1_ecoregion = factor(l1_ecoregion, levels = rev(region_order)))
  totals <- bars %>% distinct(l1_ecoregion, metabolism_metric, total_hucs)

  # 5. Companion stacked bar chart -------------------------------------------
  figure <- ggplot(bars, aes(x = fraction, y = l1_ecoregion, fill = scaling_class)) +
    geom_col(width = 0.75, position = position_stack(reverse = TRUE)) +
    geom_text(data = totals, aes(x = 1.02, y = l1_ecoregion,
                                 label = paste0("n=", total_hucs)),
                inherit.aes = FALSE, hjust = 0, size = 3) +
    facet_wrap(~ metabolism_metric, nrow = 1) +
    scale_fill_manual(values = palette, limits = classes, drop = FALSE,
                       name = "Scaling behavior") +
    scale_x_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.25),
                        expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(labels = function(x) str_wrap(x, 28)) +
    coord_cartesian(xlim = c(0, 1.18)) +
    labs(x = "Percentage of HUC8s", y = NULL,
         caption = paste0(
           "07 sets negative local rates to NA; reach-months require <",
           missing_coverage_cutoff,
           "% missing upstream area and process-specific metabolism coverage.")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
           panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
           plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  path_output <- file.path(dir_figures, if (use_test_output) {
    "figure_2_ecoregion_scaling_bars_test.png"
  } else "figure_2_ecoregion_scaling_bars.png")
  ggsave(path_output, figure, width = 14,
          height = max(5, 2.5 + 0.38 * length(region_order)), dpi = 300, bg = "white")
  message("Saved: ", path_output)
  invisible(NULL)
})
