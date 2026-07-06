
# 1. Setup ---------------------------------------------------------------------

rm(list = ls())

start = Sys.time()

source("scripts/00_setup.R")

pacman::p_load(
  forcats,
  glue,
  janitor,
  ggthemes,
  nhdplusTools,
  patchwork,
  scales,
  sf,
  rnaturalearth,
  tigris,
  tidyverse)

theme_set(theme_bw())

dir_data <- "data/processed"
dir_figs <- "figures"

huc_level <- "huc08"
huc_col <- str_replace(huc_level, "huc0", "huc")

month_levels <- month.abb
flux_levels <- c("GPP", "ER")
class_levels <- c("Sublinear", "Linear", "Superlinear", "Unclear")

# Annual class is assigned from the median slope across coherent months.
# HUC8-flux combinations with fewer than this many coherent months are labelled Unclear.
min_valid_months <- 6
min_reaches <- 10
coherence_r2 <- 0.8

slope_breaks <- c(0.8, 1.2)

class_colors <- c(
  "Sublinear" = "#2166ac",
  "Linear" = "#fdae61",
  "Superlinear" = "#b2182b",
  "Unclear" = "#cfcfcf"
)


# 2. Read in data --------------------------------------------------------------

path_flowlines <- file.path(dir_data, "riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp")
path_metrics <- file.path(dir_data, paste0("06_", huc_level, "_allometry_metrics.csv"))
path_l1_lookup <- file.path(dir_data, "08_huc08_l1_ecoregion_lookup.csv")

metrics <- read_csv(path_metrics, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    month = factor(month, levels = month_levels),
    metabolism = recode(
      metabolism_var,
      c_gpp_flux_g_o2_d = "GPP",
      c_er_flux_g_o2_d = "ER",
      c_gpp_flux_gO2_d = "GPP",
      c_er_flux_gO2_d = "ER"
    ),
    valid_fit = n >= min_reaches & r2 >= coherence_r2 & is.finite(slope)
  )


# Percent above linear expectation --------------------------------------------

huc_l1_lookup <- read_csv(path_l1_lookup, show_col_types = FALSE) %>% 
  clean_names() |>
  mutate(huc = str_pad(as.character(huc), width = 8, pad = "0"))


excess_summary <- metrics |>
  filter(
    metabolism %in% flux_levels,
    n >= min_reaches,
    r2 >= coherence_r2,
    is.finite(slope)) |>
  left_join(huc_l1_lookup, by = "huc") |>
  mutate(
    excess_vs_linear_pct_10x = (10^(slope - 1) - 1) * 100, 
    excess_vs_linear_pct_2x = (2^(slope - 1) - 1) * 100)

conus_excess_summary <- excess_summary |>
  group_by(metabolism) |>
  summarise(
    n_tests = n(),
    median_pct_10x = median(excess_vs_linear_pct_10x, na.rm = TRUE),
    p05_pct_10x = quantile(excess_vs_linear_pct_10x, 0.05, na.rm = TRUE),
    p95_pct_10x = quantile(excess_vs_linear_pct_10x, 0.95, na.rm = TRUE),
    median_pct_2x = median(excess_vs_linear_pct_2x, na.rm = TRUE),
    p05_pct_2x = quantile(excess_vs_linear_pct_2x, 0.05, na.rm = TRUE),
    p95_pct_2x = quantile(excess_vs_linear_pct_2x, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(l1_ecoregion = "CONUS", .before = metabolism)

l1_excess_summary <- excess_summary |>
  group_by(l1_ecoregion, metabolism) |>
  summarise(
    n_tests = n(),
    median_pct_10x = median(excess_vs_linear_pct_10x, na.rm = TRUE),
    p05_pct_10x = quantile(excess_vs_linear_pct_10x, 0.05, na.rm = TRUE),
    p95_pct_10x = quantile(excess_vs_linear_pct_10x, 0.95, na.rm = TRUE),
    median_pct_2x = median(excess_vs_linear_pct_2x, na.rm = TRUE),
    p05_pct_2x = quantile(excess_vs_linear_pct_2x, 0.05, na.rm = TRUE),
    p95_pct_2x = quantile(excess_vs_linear_pct_2x, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

excess_table <- bind_rows(
  conus_excess_summary,
  l1_excess_summary
) |>
  mutate(
    across(
      c(
        median_pct_10x, p05_pct_10x, p95_pct_10x,
        median_pct_2x, p05_pct_2x, p95_pct_2x
      ),
      ~ round(.x, 1)
    )
  ) |>
  arrange(metabolism, desc(l1_ecoregion == "CONUS"), l1_ecoregion)

write_csv(
  excess_table,
  file.path(dir_data, "09_excess_metabolism_vs_linear_scaling_summary.csv")
)

dodge <- position_dodge(width = 0.9)

excess_table |> 
  drop_na(l1_ecoregion) |> 
  ggplot(aes(x = l1_ecoregion,y = median_pct_2x, 
    color = metabolism,fill = metabolism)) + 
  geom_col(position = dodge, alpha = 0.4) + 
  geom_errorbar(aes(ymin = p05_pct_2x, ymax = p95_pct_2x), 
    position = dodge, width = 0.2) + 
  geom_hline(yintercept = 0) + 
  coord_flip()


