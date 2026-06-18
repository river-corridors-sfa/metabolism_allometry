rm(list = ls())

require(pacman)
pacman::p_load(
  tidyverse,
  ranger,
  fs,
  janitor
)


# 1. Paths ----------------------------------------------------------------

dir_parent = "/Users/regi350/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/RC/rcsfa-ai4do/maavara_nep"
dir_data  <- file.path(dir_parent, "dryad_data")
dir_out   <- file.path(dir_parent, "model_outputs")

dir_create(dir_out)


# 2. Read training data ---------------------------------------------------

appling <- read_csv(
  file.path(dir_data, "appling_data_06202024.csv"),
  show_col_types = FALSE
) |>
  clean_names() |>
  mutate(
    random_numbers = runif(n())
  )


# 3. Prepare GPP model data -----------------------------------------------

gpp_predictors <- c(
  "doy",
  "temp_water",
  "discharge",
  "shortwave",
  "uparea",
  "tcc_mean",
  "uparea_prop_shrubs"
)

gpp_data <- appling |>
  select(
    gpp,
    all_of(gpp_predictors)
  ) |>
  drop_na()

set.seed(24)

test_ids <- sample(
  seq_len(nrow(gpp_data)),
  size = floor(0.3 * nrow(gpp_data))
)

gpp_train <- gpp_data[-test_ids, ]
gpp_test  <- gpp_data[test_ids, ]


# 4. Fit final GPP random forest ------------------------------------------

gpp_rf <- ranger(
  gpp ~ .,
  data = gpp_train,
  num.trees = 500,
  mtry = 5,
  max.depth = 30,
  sample.fraction = 0.4,
  replace = TRUE,
  importance = "impurity",
  seed = 24
)

gpp_rf


# 5. Model performance ----------------------------------------------------

calc_metrics <- function(obs, pred) {
  tibble(
    pearson_r = cor(obs, pred),
    r2_cor = cor(obs, pred)^2,
    rmse = sqrt(mean((obs - pred)^2)),
    bias = mean(pred - obs)
  )
}

gpp_train_pred <- predict(gpp_rf, gpp_train)$predictions
gpp_test_pred  <- predict(gpp_rf, gpp_test)$predictions

gpp_metrics <- bind_rows(
  train = calc_metrics(gpp_train$gpp, gpp_train_pred),
  test  = calc_metrics(gpp_test$gpp, gpp_test_pred),
  .id = "dataset"
)

gpp_metrics


# 6. Feature importance ---------------------------------------------------

gpp_importance <- enframe(
  gpp_rf$variable.importance,
  name = "variable",
  value = "importance"
) |>
  arrange(desc(importance))

gpp_importance

p_importance <- gpp_importance |>
  ggplot(aes(x = importance, y = reorder(variable, importance))) +
  geom_col() +
  theme_minimal() +
  labs(
    x = "Importance",
    y = NULL
  )

p_importance

ggsave(
  file.path(dir_out, "gpp_rf_feature_importance.png"),
  plot = p_importance,
  width = 6,
  height = 4,
  dpi = 300
)


# 7. Observed vs predicted ------------------------------------------------

gpp_test_results <- gpp_test |>
  mutate(
    pred = gpp_test_pred
  )

p_test <- gpp_test_results |>
  ggplot(aes(x = gpp, y = pred)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = expression("Observed GPP (g O"[2]~m^{-2}~day^{-1}*")"),
    y = expression("Predicted GPP (g O"[2]~m^{-2}~day^{-1}*")")
  )

p_test

ggsave(
  file.path(dir_out, "gpp_rf_observed_vs_predicted.png"),
  plot = p_test,
  width = 6,
  height = 5,
  dpi = 300
)





