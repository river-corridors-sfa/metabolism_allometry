# 05_fit_allometry_poc.R
# Fit screening allometries for the proof-of-concept HUC8

rm(list = ls())
source("scripts/00_setup.R")


# 1. Settings ---------------------------------------------------------------

poc_huc8 <- "10100004"


# 2. Read and prepare data --------------------------------------------------

nep_monthly <- readRDS("data/processed/monthly_nep.rds")

poc_nep <- nep_monthly |>
  filter(huc8 == poc_huc8) |>
  mutate(
    nep_pos = if_else(nep > 0, nep, NA_real_),
    nep_neg = if_else(nep < 0, abs(nep), NA_real_)
  )


# 3. Define allometry function ----------------------------------------------

fit_allometry <- function(df) {
  
  if (nrow(df) < 10) {
    return(tibble(
      n = nrow(df),
      slope = NA_real_,
      intercept = NA_real_,
      r2 = NA_real_,
      slope_lwr = NA_real_,
      slope_upr = NA_real_
    ))
  }
  
  fit <- lm(log10(value) ~ log10(uparea), data = df)
  
  tibble(
    n = nrow(df),
    slope = coef(fit)[2],
    intercept = coef(fit)[1],
    r2 = summary(fit)$r.squared,
    slope_lwr = confint(fit)[2, 1],
    slope_upr = confint(fit)[2, 2]
  )
}


# 4. Fit monthly screening allometries --------------------------------------

allometry_poc <- poc_nep |>
  select(comid, month, uparea, gpp, er, nep_pos, nep_neg) |>
  pivot_longer(
    cols = c(gpp, er, nep_pos, nep_neg),
    names_to = "metric",
    values_to = "value"
  ) |>
  filter(uparea > 0,
         value > 0,
         !is.na(value)) |>
  group_by(month, metric) |>
  group_modify(\(x, y) fit_allometry(x)) |>
  ungroup() |>
  mutate(
    allometry_class = case_when(
      is.na(r2) ~ "Too few points",
      r2 < 0.8 ~ "Uncertain",
      slope_lwr < 1 & slope_upr < 1 ~ "Sublinear",
      slope_lwr < 1 & slope_upr > 1 ~ "Linear",
      slope_lwr > 1 & slope_upr > 1 ~ "Super-linear",
      TRUE ~ "Uncertain"
    )
  )


# 5. Diagnostics ------------------------------------------------------------

allometry_poc |>
  arrange(metric, month)

allometry_poc |>
  count(metric, allometry_class)

allometry_poc |>
  group_by(metric) |>
  summarize(
    min_n = min(n, na.rm = TRUE),
    max_n = max(n, na.rm = TRUE),
    min_r2 = min(r2, na.rm = TRUE),
    max_r2 = max(r2, na.rm = TRUE),
    min_slope = min(slope, na.rm = TRUE),
    max_slope = max(slope, na.rm = TRUE),
    .groups = "drop"
  )


# 6. Plot allometry trajectories --------------------------------------------

ggplot(allometry_poc, aes(month, slope, group = metric, color = metric)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(option = "viridis") +
  labs(
    title = paste("Monthly allometry slopes, HUC8", poc_huc8),
    x = NULL,
    y = "Slope",
    color = NULL
  )

ggplot(allometry_poc, aes(month, r2, group = metric, color = metric)) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(option = "viridis") +
  labs(
    title = paste("Monthly allometry strength, HUC8", poc_huc8),
    x = NULL,
    y = expression(R^2),
    color = NULL
  )


# 7. Plot example scaling relationships -------------------------------------

poc_long <- poc_nep |>
  select(comid, month, uparea, gpp, er, nep_pos, nep_neg) |>
  pivot_longer(
    cols = c(gpp, er, nep_pos, nep_neg),
    names_to = "metric",
    values_to = "value"
  ) |>
  filter(uparea > 0,
         value > 0,
         !is.na(value))

ggplot(
  poc_long |> filter(month %in% c("Jan", "Apr", "Jul", "Oct")),
  aes(uparea, value)
) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(metric ~ month, scales = "free_y") +
  labs(
    title = paste("Screening allometries, HUC8", poc_huc8),
    x = "Upstream area",
    y = "Metric value"
  )


# 8. Save output -------------------------------------------------------------

saveRDS(allometry_poc, "data/processed/allometry_poc.rds", compress = "xz")