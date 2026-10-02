## 11_uparea_sensitivity.R
##
## PURPOSE
## Test how much explicit inclusion of upstream watershed area ("uparea")
## affects Maavara-style GPP and ER predictions.
##
## DESIGN
## 1. Reproduce the published/full RF model as closely as possible.
## 2. Fit a second model with EVERYTHING identical except removal of uparea.
## 3. Compare:
##      - CONUS prediction distributions
##      - full-model minus no-uparea prediction residuals
##      - train/test R2
##
## IMPORTANT
## The RF fitting is done with Python scikit-learn through reticulate because the
## published models were built with sklearn RandomForestRegressor.


# 0. Setup ---------------------------------------------------------------------

Sys.setenv(
  RETICULATE_PYTHON =
    "/Users/regi350/Library/r-miniconda-arm64/envs/r-reticulate/bin/python"
)

source("scripts/00_setup.R")

pacman::p_load(reticulate, yardstick)

tic("run 11")

dir_figs <- "figures"
dir_processed <- "data/processed"

dir.create(dir_figs, showWarnings = FALSE, recursive = TRUE)


# 1. Read training and CONUS prediction data -----------------------------------

training <- read_csv(
  "data/raw/dryad_data/appling_data_06202024.csv",
  show_col_types = FALSE
)

conus <- read_csv(
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",
  show_col_types = FALSE
)


# 2. Define final published predictor sets -------------------------------------

# Exact predictor order used by the final seven-feature GPP notebook model.
gpp_predictors <- c(
  "doy",
  "temp.water",
  "discharge",
  "shortwave",
  "uparea",
  "TCC_mean",
  "uparea_prop_shrubs"
)

# Exact predictor order used by the final seven-feature ER notebook model.
er_predictors <- c(
  "doy",
  "temp.water",
  "discharge",
  "shortwave",
  "uparea",
  "TCC_mean",
  "uparea_prop_wet"
)

gpp_predictors_no_area <- setdiff(gpp_predictors, "uparea")
er_predictors_no_area <- setdiff(er_predictors, "uparea")


# 3. Published-model settings ---------------------------------------------------

# Exact settings used by the final reduced GPP and ER models.
test_size <- 0.30
random_state <- 24

n_estimators <- 500
max_depth <- 30
min_samples_split <- 2
min_samples_leaf <- 1
max_features <- 5
max_samples <- 0.4


# 4. Connect to sklearn ---------------------------------------------------------

sklearn_ensemble <- import("sklearn.ensemble")
sklearn_model_selection <- import("sklearn.model_selection")


# 5. Helper: fit paired full and no-uparea models -------------------------------

fit_rf_pair <- function(data,
                        response,
                        predictors,
                        predictors_no_area) {
  
  model_data <- data %>%
    select(all_of(c(response, predictors))) %>%
    drop_na()
  
  # Convert explicitly to objects sklearn understands correctly.
  x_full <- model_data %>%
    select(all_of(predictors)) %>%
    as.matrix()
  
  y <- model_data[[response]]
  
  
  # Create ONE train/test split shared by both models.
  split <- sklearn_model_selection$train_test_split(
    x_full,
    y,
    test_size = test_size,
    random_state = as.integer(random_state)
  )
  
  x_train_full <- split[[1]]
  x_test_full <- split[[2]]
  y_train <- split[[3]]
  y_test <- split[[4]]
  
  
  # Identify columns retained in the no-uparea model.
  no_area_cols <- match(
    predictors_no_area,
    predictors
  )
  
  x_train_no_area <- x_train_full[, no_area_cols, drop = FALSE]
  x_test_no_area <- x_test_full[, no_area_cols, drop = FALSE]
  
  
  # Full model: exact published RF settings.
  rf_full <- sklearn_ensemble$RandomForestRegressor(
    n_estimators = as.integer(n_estimators),
    max_depth = as.integer(max_depth),
    min_samples_split = as.integer(min_samples_split),
    min_samples_leaf = as.integer(min_samples_leaf),
    max_features = as.integer(max_features),
    max_samples = max_samples,
    random_state = as.integer(random_state),
    n_jobs = -1L,
    oob_score = TRUE
  )
  
  rf_full$fit(
    x_train_full,
    y_train
  )
  
  
  # Identical model with only uparea removed.
  rf_no_area <- sklearn_ensemble$RandomForestRegressor(
    n_estimators = as.integer(n_estimators),
    max_depth = as.integer(max_depth),
    min_samples_split = as.integer(min_samples_split),
    min_samples_leaf = as.integer(min_samples_leaf),
    max_features = as.integer(max_features),
    max_samples = max_samples,
    random_state = as.integer(random_state),
    n_jobs = -1L,
    oob_score = TRUE
  )
  
  rf_no_area$fit(
    x_train_no_area,
    y_train
  )
  
  
  # Generate predictions once for performance calculations.
  pred_train_full <- rf_full$predict(x_train_full)
  pred_test_full <- rf_full$predict(x_test_full)
  
  pred_train_no_area <- rf_no_area$predict(x_train_no_area)
  pred_test_no_area <- rf_no_area$predict(x_test_no_area)
  
  
  # The published notebooks report Pearson correlation squared as R2.
  performance <- tibble(
    model = rep(
      c("Full", "No uparea"),
      each = 2
    ),
    
    dataset = rep(
      c("Train", "Test"),
      times = 2
    ),
    
    r2 = c(
      cor(y_train, pred_train_full)^2,
      cor(y_test, pred_test_full)^2,
      cor(y_train, pred_train_no_area)^2,
      cor(y_test, pred_test_no_area)^2
    )
  )
  
  
  list(
    full = rf_full,
    no_area = rf_no_area,
    performance = performance
  )
}


# 6. Fit GPP and ER paired models ----------------------------------------------

message("Fitting paired GPP models...")

gpp_models <- fit_rf_pair(
  data = training,
  response = "GPP",
  predictors = gpp_predictors,
  predictors_no_area = gpp_predictors_no_area)


message("Fitting paired ER models...")

er_models <- fit_rf_pair(
  data = training %>% mutate(ER = -ER),
  response = "ER",
  predictors = er_predictors,
  predictors_no_area = er_predictors_no_area)


# 7. First validation: reproduced base-model performance -----------------------

rf_performance <- bind_rows(
  gpp_models$performance %>%
    mutate(metabolism = "GPP"),
  
  er_models$performance %>%
    mutate(metabolism = "ER")
)

rf_performance

# Published notebook values for the full seven-feature models should be
# approximately:
# GPP: Train R2 = 0.9271, Test R2 = 0.8278
# ER:  Train R2 = 0.9041, Test R2 = 0.7785


# 8. Prepare one CONUS month ----------------------------------------------------

month_test <- 7
month_name <- month.abb[month_test]

conus_month <- conus %>%
  transmute(
    COMID,
    
    # Static predictors
    uparea,
    TCC_mean,
    uparea_prop_shrubs = up_shrubs,
    uparea_prop_wet = up_wet,
    
    # Monthly predictors
    discharge = .data[[paste0(
      "Q_avg_",
      month.name[month_test]
    )]],
    
    temp.water = .data[[paste0(
      "ML_river_temp_",
      month_test
    )]],
    
    shortwave = .data[[paste0(
      "shortwave_",
      tolower(month.abb[month_test]),
      "_mean"
    )]],
    
    # Notebook predictions use the middle day of each month.
    # For July this is DOY 197.
    doy = 197
  ) %>%
  drop_na()


# 9. Generate paired CONUS predictions -----------------------------------------

# GPP predictor matrices in the exact training-model predictor order.
x_gpp_full <- conus_month %>%
  select(all_of(gpp_predictors)) %>%
  as.matrix()

x_gpp_no_area <- conus_month %>%
  select(all_of(gpp_predictors_no_area)) %>%
  as.matrix()


gpp_predictions <- conus_month %>%
  mutate(
    prediction_full =
      gpp_models$full$predict(x_gpp_full),
    
    prediction_no_uparea =
      gpp_models$no_area$predict(x_gpp_no_area),
    
    residual =
      prediction_full - prediction_no_uparea,
    
    metabolism = "GPP"
  )


# ER predictor matrices in the exact training-model predictor order.
x_er_full <- conus_month %>%
  select(all_of(er_predictors)) %>%
  as.matrix()

x_er_no_area <- conus_month %>%
  select(all_of(er_predictors_no_area)) %>%
  as.matrix()


er_predictions <- conus_month %>%
  mutate(
    prediction_full =
      er_models$full$predict(x_er_full),
    
    prediction_no_uparea =
      er_models$no_area$predict(x_er_no_area),
    
    residual =
      prediction_full - prediction_no_uparea,
    
    metabolism = "ER"
  )


predictions <- bind_rows(
  gpp_predictions,
  er_predictions
)

saveRDS(
  predictions,
  file.path(
    dir_processed,
    "11_uparea_sensitivity_predictions_test.rds"
  )
)


# 10. Figure 1 - prediction distributions --------------------------------------

prediction_long <- predictions %>%
  select(
    metabolism,
    prediction_full,
    prediction_no_uparea
  ) %>%
  pivot_longer(
    starts_with("prediction_"),
    names_to = "model",
    values_to = "prediction"
  ) %>%
  mutate(
    model = recode(
      model,
      prediction_full = "Full",
      prediction_no_uparea = "No uparea"
    )
  )


p_distribution <- prediction_long %>%
  ggplot(
    aes(
      x = model,
      y = prediction
    )
  ) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  facet_wrap(
    ~metabolism,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = expression(
      "Predicted metabolism (g O"[2]~
        m^-2~
        d^-1*")"
    ),
    title =
      "Predicted metabolism with and without upstream area"
  ) +
  theme_bw()


ggsave(
  file.path(
    dir_figs,
    "11_prediction_distributions.png"
  ),
  p_distribution,
  width = 7,
  height = 5,
  dpi = 300
)


# 11. Figure 2 - prediction residuals ------------------------------------------

p_residuals <- predictions %>%
  ggplot(
    aes(
      x = residual
    )
  ) +
  geom_histogram(
    bins = 60
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~metabolism,
    scales = "free"
  ) +
  labs(
    x =
      "Full prediction - no-uparea prediction",
    y = "CONUS reaches",
    title =
      "Effect of upstream area on metabolism predictions"
  ) +
  theme_bw()


ggsave(
  file.path(
    dir_figs,
    "11_prediction_residuals.png"
  ),
  p_residuals,
  width = 8,
  height = 5,
  dpi = 300
)


# 12. Figure 3 - train/test R2 --------------------------------------------------

p_r2 <- rf_performance %>%
  ggplot(
    aes(
      x = dataset,
      y = r2,
      group = model,
      shape = model
    )
  ) +
  geom_point(
    size = 3,
    position = position_dodge(
      width = 0.25
    )
  ) +
  facet_wrap(
    ~metabolism
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    x = NULL,
    y = expression(R^2),
    shape = "Model",
    title =
      "Random-forest performance with and without upstream area"
  ) +
  theme_bw()


ggsave(
  file.path(
    dir_figs,
    "11_train_test_r2.png"
  ),
  p_r2,
  width = 7,
  height = 5,
  dpi = 300
)


toc()