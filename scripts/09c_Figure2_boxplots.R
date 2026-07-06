
# 1. Setup ---------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(tidyverse, janitor, patchwork, scales)

dir_data <- "data/processed"
dir_figs <- "figures"

min_n <- 30
min_r2 <- 0.8

flux_colors <- c(
  "ER" = "#440154",
  "GPP" = "#7ad151"
)


# 2. Read in data --------------------------------------------------------------

allometry <- read_csv(file.path(dir_data, "06_huc08_allometry_metrics.csv"), show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), 8, pad = "0"),
    metabolism = case_when(
      metabolism_var == "c_gpp_flux_gO2_d" ~ "GPP",
      metabolism_var == "c_er_flux_gO2_d" ~ "ER",
      TRUE ~ NA_character_
    )
  )

huc_l1 <- read_csv(file.path(dir_data, "08_huc08_l1_ecoregion_lookup.csv"), show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), 8, pad = "0"),
    l1_code = as.character(l1_code),
    l1_name = l1_ecoregion,
    l1_name = str_remove(l1_name, "^[0-9]+\\s+")
  )

annual_metrics <- allometry |>
  filter(n >= min_n, r2 >= min_r2) |>
  group_by(huc, metabolism) |>
  summarise(
    slope = median(slope, na.rm = TRUE),
    r2 = median(r2, na.rm = TRUE),
    intercept = median(intercept, na.rm = TRUE),
    n_months = n(),
    .groups = "drop"
  ) |>
  left_join(huc_l1 |> select(huc, l1_name), by = "huc") |>
  filter(!is.na(l1_name), metabolism %in% c("GPP", "ER"))

ecoregion_order <- annual_metrics |>
  group_by(l1_name, metabolism) |>
  summarise(mean_slope = mean(slope, na.rm = TRUE), .groups = "drop") |>
  group_by(l1_name) |>
  summarise(mean_slope = mean(mean_slope, na.rm = TRUE), .groups = "drop") |>
  arrange(mean_slope) |>
  pull(l1_name)

plot_data <- annual_metrics |>
  pivot_longer(
    cols = c(slope, r2, intercept),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    l1_name = factor(l1_name, levels = ecoregion_order),
    metabolism = factor(metabolism, levels = c("ER", "GPP")),
    metric = recode(metric, slope = "Slope", r2 = "R²", intercept = "Intercept"),
    metric = factor(metric, levels = c("Slope", "R²", "Intercept"))
  )


# 3. Plot function -------------------------------------------------------------

plot_metric <- function(selected_metric, x_lab, add_ref = NULL) {
  p <- plot_data |>
    filter(metric == selected_metric) |>
    ggplot(aes(x = value, y = l1_name, fill = metabolism)) +
    geom_boxplot(
      width = 0.55,
      outlier.alpha = 0.22,
      outlier.size = 1.1,
      position = position_dodge(width = 0.7)
    ) +
    scale_fill_manual(values = flux_colors, name = NULL) +
    labs(x = x_lab, y = NULL) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "grey90"),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_text(size = 10)
    )
  
  if (!is.null(add_ref)) {
    p <- p + geom_vline(xintercept = add_ref, linetype = "dashed", linewidth = 0.35)
  }
}


# 4. Slopes figure -------------------------------------------------------------

pacman::p_load(multcompView)

slope_data <- plot_data |>
  filter(metric == "Slope")

get_tukey_letters <- function(df) {
  
  l1_lookup <- df |>
    distinct(l1_name) |>
    arrange(l1_name) |>
    mutate(l1_id = paste0("eco_", row_number()))
  
  df_safe <- df |>
    left_join(l1_lookup, by = "l1_name") |>
    mutate(l1_id = factor(l1_id, levels = l1_lookup$l1_id))
  
  fit <- aov(value ~ l1_id, data = df_safe)
  
  tukey <- TukeyHSD(fit, "l1_id")$l1_id[, "p adj"]
  
  letters <- multcompView::multcompLetters(tukey)$Letters |>
    enframe(name = "l1_id", value = "tukey_letter")
  
  l1_lookup |>
    left_join(letters, by = "l1_id") |>
    transmute(
      l1_name,
      tukey_letter = replace_na(tukey_letter, "a")
    )
}

slope_range <- range(slope_data$value, na.rm = TRUE)
slope_pad <- diff(slope_range) * 0.10

slope_letters <- slope_data |>
  group_by(metabolism) |>
  group_modify(~get_tukey_letters(.x)) |>
  ungroup() |>
  mutate(
    l1_name = factor(l1_name, levels = levels(slope_data$l1_name)),
    x = slope_range[2] + slope_pad
  )

p_09b_slope <- slope_data |>
  ggplot(aes(x = value, y = l1_name, fill = metabolism)) +
  geom_boxplot(
    width = 0.55,
    outlier.alpha = 0.22,
    outlier.size = 1.1,
    position = position_dodge(width = 0.7)
  ) +
  geom_text(
    data = slope_letters,
    aes(
      x = x,
      y = l1_name,
      label = tukey_letter,
      color = metabolism,
      group = metabolism
    ),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.7),
    hjust = 0,
    size = 3
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.35
  ) +
  scale_fill_manual(values = flux_colors, name = NULL) +
  scale_color_manual(values = flux_colors, guide = "none") +
  scale_x_continuous(
    limits = c(slope_range[1], slope_range[2] + slope_pad * 3),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Annual median allometric slope", y = NULL) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey90"),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 10)
  )

p_09b_slope
ggsave(file.path(dir_figs, "09b_slope_ecoregion_boxplot.png"),
  p_09b_slope,
  width = 7,
  height = 6,
  dpi = 400)  







