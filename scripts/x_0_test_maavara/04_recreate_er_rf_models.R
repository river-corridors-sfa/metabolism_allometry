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
  clean_names()


# 3. Prepare ER model data ------------------------------------------------

er_predictors <- c(
  "doy",
  "temp_water",
  "discharge",
  "shortwave",
  "uparea",
  "tcc_mean",
  "uparea_prop_wet"
)

er_data <- appling |>
  mutate(
    er = er * -1
  ) |>
  select(
    er,
    all_of(er_predictors)
  ) |>
  drop_na()

set.seed(24)

test_ids <- sample(
  seq_len(nrow(er_data)),
  size = floor(0.3 * nrow(er_data))
)

er_train <- er_data[-test_ids, ]
er_test  <- er_data[test_ids, ]


# 4. Fit final ER random forest -------------------------------------------

er_rf <- ranger(
  er ~ .,
  data = er_train,
  num.trees = 500,
  mtry = 5,
  max.depth = 30,
  sample.fraction = 0.4,
  replace = TRUE,
  importance = "impurity",
  seed = 24
)

er_rf


# 5. Model performance ----------------------------------------------------

calc_metrics <- function(obs, pred) {
  tibble(
    pearson_r = cor(obs, pred),
    r2_cor = cor(obs, pred)^2,
    rmse = sqrt(mean((obs - pred)^2)),
    bias = mean(pred - obs)
  )
}

er_train_pred <- predict(er_rf, er_train)$predictions
er_test_pred  <- predict(er_rf, er_test)$predictions

er_metrics <- bind_rows(
  train = calc_metrics(er_train$er, er_train_pred),
  test  = calc_metrics(er_test$er, er_test_pred),
  .id = "dataset"
)

er_metrics


# 6. Feature importance ---------------------------------------------------

er_importance <- enframe(
  er_rf$variable.importance,
  name = "variable",
  value = "importance"
) |>
  arrange(desc(importance))

er_importance

p_importance <- er_importance |>
  ggplot(aes(x = importance, y = reorder(variable, importance))) +
  geom_col() +
  theme_minimal() +
  labs(
    x = "Importance",
    y = NULL
  )

p_importance

ggsave(
  file.path(dir_out, "er_rf_feature_importance.png"),
  plot = p_importance,
  width = 6,
  height = 4,
  dpi = 300
)


# 7. Observed vs predicted ------------------------------------------------

er_test_results <- er_test |>
  mutate(
    pred = er_test_pred
  )

p_test <- er_test_results |>
  ggplot(aes(x = er, y = pred)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    x = expression("Observed ER (g O"[2]~m^{-2}~day^{-1}*")"),
    y = expression("Predicted ER (g O"[2]~m^{-2}~day^{-1}*")")
  )

p_test

ggsave(
  file.path(dir_out, "er_rf_observed_vs_predicted.png"),
  plot = p_test,
  width = 6,
  height = 5,
  dpi = 300
)