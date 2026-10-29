# 10_figure1_scaling_maps.R
# Run from the project root after 08. Writes Figure 1 and its criteria table.
# Classes describe seasonal agreement and tested deletion sensitivity, not CIs.

library(tidyverse)
library(sf)
library(gridExtra)
library(rnaturalearth)

# 1. Settings ---------------------------------------------------------------
use_test_output <- FALSE          # TRUE previews 08's three-HUC test only
min_reaches <- 20L
min_area_range <- 1
missing_coverage_cutoff <- 10  # 07 requires <10% for area and metabolism
slope_bounds <- c(0.95, 1.05)
r2_cutoff <- 0.8
dir_processed <- "data/processed"
dir_figures <- "figures"
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)
path_input <- file.path(dir_processed, if (use_test_output) {
  "08_no_uparea_huc08_fit_stability_test.csv"
} else "08_no_uparea_huc08_fit_stability.csv")
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
  mutate(huc = str_pad(as.character(huc8), 8, pad = "0")) %>%
  select(huc) %>%
  st_transform(5070)
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
map_data <- boundaries %>% left_join(annual, by = "huc",
                                     relationship = "one-to-many")
rm(monthly, annual)

# 4. Figure 1 --------------------------------------------------------------
# Low-resolution country outline is bundled data; no HUC boundary download.
conus <- ne_states(
  country = "United States of America",
  returnclass = "sf") %>%
  filter(!postal %in% c("AK", "HI", "PR")) %>%
  st_transform(5070) %>%
  summarise(geometry = st_union(geometry))

extent <- st_bbox(conus)
figure1 <- ggplot() +
  geom_sf(data = conus, fill = "white", color = "grey55", linewidth = 0.25) +
  geom_sf(data = map_data, aes(fill = scaling_class),
          color = "white", linewidth = 0.04) +
  geom_sf(data = conus, fill = NA, color = "black", linewidth = 0.5) +
  facet_wrap(~ metabolism_metric, nrow = 1) +
  scale_fill_manual(values = palette, limits = classes, drop = FALSE,
                     name = "Scaling behavior") +
  coord_sf(crs = st_crs(5070), xlim = extent[c("xmin", "xmax")],
           ylim = extent[c("ymin", "ymax")], expand = FALSE, datum = NA) +
  theme_void(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 15),
         plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
path_map <- file.path(dir_figures, if (use_test_output) {
  "figure_1_no_uparea_scaling_maps_test.png"
} else "figure_1_no_uparea_scaling_maps.png")
ggsave(path_map, figure1, width = 14, height = 6.4, dpi = 300, bg = "white")

# 5. Supplemental criteria table -------------------------------------------
# Generate labels from settings so the saved table matches the map rules.
low <- format(slope_bounds[1], trim = TRUE)
high <- format(slope_bounds[2], trim = TRUE)
r2 <- format(r2_cutoff, trim = TRUE)
support_rule <- paste0("Full fit has <", min_reaches, " reaches, <", min_reaches,
  " distinct drainage areas, <", min_area_range,
  " log10 unit of area range, or undefined slope/R2; or no valid deletion refits. ",
  "Missing/nonfinite required summary values also receive this category.")
strong_rule <- paste0("Full-fit R2 and every valid deletion R2 >= ", r2, ".")
criteria <- tribble(
  ~Category, ~`Monthly criteria`, ~`Annual map criteria`,
  "Sub-linear", paste0("Full-fit slope and every valid deletion slope < ", low,
                        ". ", strong_rule), "All 12 available months are Sub-linear.",
  "Linear", paste0("Full-fit slope and every valid deletion slope lie within [",
                    low, ", ", high, "] (inclusive). ", strong_rule,
                    " This denotes near-linear scaling."), "All 12 available months are Linear.",
  "Super-linear", paste0("Full-fit slope and every valid deletion slope > ", high,
                          ". ", strong_rule), "All 12 available months are Super-linear.",
  "Uncertain", paste0("Support requirements pass, but the full fit or any valid deletion has R2 < ",
                       r2, ", or full/deletion slopes span multiple scaling categories."),
                       "All 12 available months are Uncertain.",
  "Insufficient data", support_rule,
                       "Any monthly result is absent, or all 12 available months are Insufficient data.",
  "Seasonally variable", "Not assigned at the monthly level.",
                       "All 12 monthly results are available, but their categories differ. Includes mixtures with Uncertain or Insufficient data."
) %>% mutate(across(everything(), ~ str_wrap(.x,
  width = if (cur_column() == "Category") 20 else if (cur_column() == "Monthly criteria") 69 else 47)))
table <- tableGrob(criteria, rows = NULL, theme = ttheme_minimal(
  base_size = 11, padding = grid::unit(c(5, 5), "mm"),
  core = list(fg_params = list(hjust = 0, x = 0.02),
              bg_params = list(fill = rep(c("white", "#F1F3F5"), length.out = 6))),
  colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.02),
                 bg_params = list(fill = "#DCE4EA"))
))
table$widths <- grid::unit(c(0.17, 0.48, 0.35), "npc")
supplement <- arrangeGrob(
  grid::textGrob("Supplemental Figure S1  Scaling classification criteria",
                 x = 0, hjust = 0, gp = grid::gpar(fontsize = 16, fontface = "bold")),
  table,
  grid::textGrob(paste(
    "Apply monthly rules in this order: Insufficient data, Uncertain, then the three scaling categories.",
    paste0("07 sets negative local rates to NA and requires <",
           missing_coverage_cutoff,
           "% missing upstream area and process-specific metabolism coverage per reach-month."),
    "Valid deletion refits meet the same reach-count, distinct-area and area-range requirements as the full fit.",
    "Invalid deletions are counted separately in 08. Agreement refers only to tested sections, not all possible deletions.",
    "Thresholds: 07's working slope band; R2 cutoff follows Regier et al. (2025), doi:10.1029/2024JG008344.",
    "Deletion sensitivity is an adaptation; neither its range nor the working near-linear band is a confidence interval.",
    sep = "\n"), x = 0, hjust = 0, gp = grid::gpar(fontsize = 10)),
  ncol = 1, heights = c(0.07, 0.78, 0.15)
)
path_table <- file.path(dir_figures, "figure_S1_no_uparea_scaling_criteria.png")
ggsave(path_table, supplement, width = 15, height = 9, dpi = 300, bg = "white")
message("Saved: ", path_map, "\nSaved: ", path_table)
invisible(NULL)
