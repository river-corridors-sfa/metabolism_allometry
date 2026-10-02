## 15_quantify_network_baseline.R
##
## PURPOSE
## Quantify how much observed variation in metabolic allometry is reproduced by
## the network-constrained null, and identify where and when departures are largest.


# 0. Setup ---------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")

tic("run 15") # 11.5 minutes

dir_processed <- "data/processed"
dir_merit <- file.path(dir_processed, "merit_data")
dir_figs <- "figures"


# 1. Read data -----------------------------------------------------------------

results <- read_csv(
  file.path(dir_processed, "10_huc08_null_comparison.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    huc = as.character(huc),
    delta_beta = observed_slope - null_slope_median
  )

huc_shapes <- read_sf(
  file.path(
    dir_merit,
    "03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
  )
) %>%
  clean_names() %>%
  st_make_valid() %>%
  mutate(huc = as.character(huc)) %>%
  group_by(huc) %>%
  summarise()


# 2. Quantify variation reproduced by the null ---------------------------------

network_summary <- results %>%
  group_by(metabolism_metric) %>%
  summarise(
    n = n(),
    
    # Fraction of observed slope variation reproduced on the 1:1 scale
    variance_reproduced =
      1 -
      sum(
        (observed_slope - null_slope_median)^2,
        na.rm = TRUE
      ) /
      sum(
        (observed_slope - mean(observed_slope, na.rm = TRUE))^2,
        na.rm = TRUE
      ),
    
    # Fraction of between-fit variation associated with the null slope
    r2_observed_null =
      cor(
        observed_slope,
        null_slope_median,
        use = "complete.obs"
      )^2,
    
    median_delta_beta =
      median(delta_beta, na.rm = TRUE),
    
    median_abs_delta_beta =
      median(abs(delta_beta), na.rm = TRUE),
    
    .groups = "drop"
  )

print(network_summary)

write_csv(
  network_summary,
  file.path(dir_processed, "15_network_variance_summary.csv")
)


# 3. Spatial distribution of departures ----------------------------------------

huc_delta <- results %>%
  group_by(huc, metabolism_metric) %>%
  summarise(
    median_delta_beta = median(delta_beta, na.rm = TRUE),
    median_abs_delta_beta = median(abs(delta_beta), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(median_delta_beta),
    metabolism_metric %in% c("ER", "GPP")
  )


map_data <- huc_shapes %>%
  inner_join(huc_delta, by = "huc")


p_map <- ggplot(map_data) +
  geom_sf(
    aes(fill = median_delta_beta),
    linewidth = 0
  ) +
  scale_fill_gradient2(
    midpoint = 0,
    name = expression(
      "Median " * Delta * beta
    )
  ) +
  facet_wrap(~metabolism_metric) +
  coord_sf(
    xlim = c(-125, -66),
    ylim = c(24, 50)
  ) +
  labs(
    title = "Spatial distribution of departures from network-constrained allometry"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(dir_figs, "15_delta_beta_map.png"),
  p_map,
  width = 10,
  height = 5.5,
  dpi = 300
)


# 4. Seasonal distribution of departures ---------------------------------------

p_month <- results %>%
  filter(
    metabolism_metric %in% c("ER", "GPP"),
    is.finite(delta_beta)
  ) %>%
  mutate(
    month = factor(
      month,
      levels = 1:12,
      labels = month.abb
    )
  ) %>%
  ggplot(
    aes(
      x = month,
      y = delta_beta
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_violin(
    trim = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 1.5
  ) +
  facet_wrap(
    ~metabolism_metric,
    ncol = 1
  ) +
  labs(
    x = NULL,
    y = expression(
      Delta * beta ==
        beta[observed] - beta[null]
    ),
    title = "Seasonal distribution of departures from network-constrained allometry"
  ) +
  theme_bw()

ggsave(
  file.path(dir_figs, "15_delta_beta_by_month.png"),
  p_month,
  width = 9,
  height = 7,
  dpi = 300
)


# 5. Identify strongest candidate control points / moments ----------------------

departure_thresholds <- results %>%
  group_by(metabolism_metric) %>%
  summarise(
    threshold = quantile(
      abs(delta_beta),
      0.90,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

strong_departures <- results %>%
  left_join(
    departure_thresholds,
    by = "metabolism_metric"
  ) %>%
  filter(
    abs(delta_beta) >= threshold
  ) %>%
  arrange(
    metabolism_metric,
    desc(abs(delta_beta))
  )

write_csv(
  strong_departures,
  file.path(dir_processed, "15_strongest_delta_beta.csv")
)


toc()