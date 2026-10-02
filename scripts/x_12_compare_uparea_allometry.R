## 12_compare_uparea_allometry.R
##
## PURPOSE
## Test whether removing explicit upstream watershed area from the Maavara
## RF changes HUC8-scale allometric scaling.
##
## This is intentionally a July-only diagnostic using predictions generated
## in script 11.


# 0. Setup ---------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(Matrix)

tic("run 12")

dir_processed <- "data/processed"
dir_figs <- "figures"
dir_merit <- file.path(dir_processed, "merit_data")


# 1. Read inputs ---------------------------------------------------------------

predictions <- readRDS(
  file.path(
    dir_processed,
    "11_uparea_sensitivity_predictions_test.rds"
  )
) %>%
  mutate(
    COMID = as.character(COMID)
  )


# Reach attributes used in the existing allometry workflow
flowlines <- read_sf(
  file.path(
    dir_merit,
    "03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
  )
) %>%
  clean_names() %>%
  mutate(
    comid = as.character(comid)
  )


catchments <- read_sf(
  file.path(
    dir_merit,
    "03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
  )
) %>%
  clean_names() %>%
  mutate(
    comid = as.character(comid),
    huc = as.character(huc)
  )


upstream_lookup <- readRDS(
  file.path(
    dir_merit,
    "05_upstream_lookup.rds"
  )
) %>%
  mutate(
    downstream_comid = as.character(downstream_comid),
    upstream_comid = as.character(upstream_comid)
  )


# 2. Prepare reach-level metabolism --------------------------------------------

conus <- read_csv(
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",
  show_col_types = FALSE
) %>%
  clean_names() %>%
  mutate(
    comid = as.character(comid)
  )


reach_data <- catchments %>%
  st_drop_geometry() %>%
  select(
    comid,
    huc,
    unitarea
  ) %>%
  left_join(
    conus %>%
      select(
        comid,
        lengthkm,
        merit_width
      ),
    by = "comid"
  ) %>%
  left_join(
    predictions %>%
      select(
        COMID,
        metabolism,
        prediction_full,
        prediction_no_uparea
      ) %>%
      rename(
        comid = COMID
      ),
    by = "comid"
  ) %>%
  pivot_longer(
    cols = c(
      prediction_full,
      prediction_no_uparea
    ),
    names_to = "model",
    values_to = "rate"
  ) %>%
  mutate(
    model = recode(
      model,
      prediction_full = "Full",
      prediction_no_uparea = "No uparea"
    ),
    stream_area_m2 =
      lengthkm * 1000 * merit_width,
    local_flux =
      rate * stream_area_m2
  )
  



# 3. Build sparse routing matrix ------------------------------------------------

reach_index <- catchments %>%
  st_drop_geometry() %>%
  transmute(
    comid,
    idx = row_number()
  )


routing_lookup <- upstream_lookup %>%
  inner_join(
    reach_index %>%
      rename(
        upstream_comid = comid,
        upstream_idx = idx
      ),
    by = "upstream_comid"
  ) %>%
  inner_join(
    reach_index %>%
      rename(
        downstream_comid = comid,
        downstream_idx = idx
      ),
    by = "downstream_comid"
  )


routing_matrix <- sparseMatrix(
  i = routing_lookup$downstream_idx,
  j = routing_lookup$upstream_idx,
  x = 1,
  dims = c(
    nrow(reach_index),
    nrow(reach_index)
  )
)


# Cumulative drainage area
unitarea <- catchments %>%
  st_drop_geometry() %>%
  right_join(
    reach_index,
    by = "comid"
  ) %>%
  arrange(idx) %>%
  pull(unitarea) %>%
  replace_na(0)


cumulative_area <- as.numeric(
  routing_matrix %*% unitarea
)


routing_base <- reach_index %>%
  mutate(
    cumulative_area_km2 =
      cumulative_area
  ) %>%
  left_join(
    catchments %>%
      st_drop_geometry() %>%
      select(
        comid,
        huc
      ),
    by = "comid"
  )


# 4. Route each metabolism/model combination -----------------------------------

route_one <- function(metric, model_name) {
  
  local <- reach_data %>%
    filter(
      metabolism == metric,
      model == model_name
    ) %>%
    right_join(
      reach_index,
      by = "comid"
    ) %>%
    arrange(idx)
  
  
  local_flux <- replace_na(
    local$local_flux,
    0
  )
  
  
  tibble(
    comid = reach_index$comid,
    
    metabolism = metric,
    
    model = model_name,
    
    cumulative_flux =
      as.numeric(
        routing_matrix %*% local_flux
      )
  )
}


routed <- bind_rows(
  route_one("GPP", "Full"),
  route_one("GPP", "No uparea"),
  route_one("ER", "Full"),
  route_one("ER", "No uparea")
) %>%
  left_join(
    routing_base,
    by = "comid"
  ) %>%
  filter(
    cumulative_area_km2 > 0,
    cumulative_flux > 0
  )


# 5. Calculate HUC8 allometry ---------------------------------------------------

allometry <- routed %>%
  group_by(
    huc,
    metabolism,
    model
  ) %>%
  filter(
    n() >= 20,
    n_distinct(cumulative_area_km2) >= 20,
    log10(
      max(cumulative_area_km2) /
        min(cumulative_area_km2)
    ) >= 1
  ) %>%
  group_modify(~ {
    
    fit <- lm(
      log10(cumulative_flux) ~
        log10(cumulative_area_km2),
      data = .x
    )
    
    tibble(
      slope = coef(fit)[2],
      r2 = summary(fit)$r.squared,
      n_reaches = nrow(.x)
    )
    
  }) %>%
  ungroup()


write_csv(
  allometry,
  file.path(
    dir_processed,
    "12_uparea_allometry_july.csv"
  )
)


# 6. Direct full vs no-uparea comparison ---------------------------------------

comparison <- allometry %>%
  select(
    huc,
    metabolism,
    model,
    slope,
    r2
  ) %>%
  pivot_wider(
    names_from = model,
    values_from = c(
      slope,
      r2
    )
  ) %>%
  mutate(
    delta_slope =
      `slope_No uparea` -
      slope_Full,
    
    delta_r2 =
      `r2_No uparea` -
      r2_Full
  )


print(
  comparison %>%
    group_by(metabolism) %>%
    summarise(
      n_hucs = n(),
      
      slope_correlation =
        cor(
          slope_Full,
          `slope_No uparea`,
          use = "complete.obs"
        ),
      
      median_delta_slope =
        median(
          delta_slope,
          na.rm = TRUE
        ),
      
      median_abs_delta_slope =
        median(
          abs(delta_slope),
          na.rm = TRUE
        ),
      
      median_delta_r2 =
        median(
          delta_r2,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
)


# 7. Plot -----------------------------------------------------------------------

p <- comparison %>%
  ggplot(
    aes(
      x = slope_Full,
      y = `slope_No uparea`
    )
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    alpha = 0.4
  ) +
  facet_wrap(
    ~metabolism
  ) +
  coord_equal() +
  labs(
    x = "Full-model slope",
    y = "No-uparea slope",
    title =
      "Effect of upstream-area removal on July HUC8 allometry"
  ) +
  theme_bw()


ggsave(
  file.path(
    dir_figs,
    "12_uparea_allometry_july.png"
  ),
  p,
  width = 7,
  height = 5,
  dpi = 300
)


toc()