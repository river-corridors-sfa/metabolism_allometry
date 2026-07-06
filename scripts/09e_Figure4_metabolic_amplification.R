

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(
  janitor,
  scales,
  tidyverse
)


# 0. Settings ------------------------------------------------------------------

dir_data <- "data/processed"
dir_figs <- "figures"

path_metrics <- file.path(dir_data, "06_huc08_allometry_metrics.csv")
path_l1_lookup <- file.path(dir_data, "08_huc08_l1_ecoregion_lookup.csv")

min_reaches <- 10
coherence_r2 <- 0.8

flux_levels <- c("ER", "GPP")

flux_colors <- c(
  ER = "#7B3294",
  GPP = "#1B9E77"
)


# 1. Read data -----------------------------------------------------------------

metrics <- read_csv(path_metrics, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    metabolism = recode(
      metabolism_var,
      c_gpp_flux_g_o2_d = "GPP",
      c_er_flux_g_o2_d = "ER",
      c_gpp_flux_gO2_d = "GPP",
      c_er_flux_gO2_d = "ER"
    )
  )

huc_l1_lookup <- read_csv(path_l1_lookup, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    l1_ecoregion = str_replace_all(l1_ecoregion, "–", "-")
  )


# 2. Calculate doubling-scale metabolic amplification --------------------------

amplification_data <- metrics |>
  filter(
    metabolism %in% flux_levels,
    n >= min_reaches,
    r2 >= coherence_r2,
    is.finite(slope)
  ) |>
  left_join(huc_l1_lookup, by = "huc") |>
  mutate(
    total_pct_2x = (2^slope - 1) * 100,
    excess_pct_2x = (2^(slope - 1) - 1) * 100
  )

summarise_amplification <- function(data) {
  data |>
    summarise(
      n_tests = n(),
      n_huc8 = n_distinct(huc),
      median_total_pct_2x = median(total_pct_2x, na.rm = TRUE),
      p25_total_pct_2x = quantile(total_pct_2x, 0.25, na.rm = TRUE, names = FALSE),
      p75_total_pct_2x = quantile(total_pct_2x, 0.75, na.rm = TRUE, names = FALSE),
      median_excess_pct_2x = median(excess_pct_2x, na.rm = TRUE),
      p25_excess_pct_2x = quantile(excess_pct_2x, 0.25, na.rm = TRUE, names = FALSE),
      p75_excess_pct_2x = quantile(excess_pct_2x, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
}

conus_amplification_summary <- amplification_data |>
  group_by(metabolism) |>
  summarise_amplification() |>
  mutate(l1_ecoregion = "CONUS", .before = metabolism)

l1_amplification_summary <- amplification_data |>
  filter(!is.na(l1_ecoregion)) |>
  group_by(l1_ecoregion, metabolism) |>
  summarise_amplification()

amplification_table <- bind_rows(
  conus_amplification_summary,
  l1_amplification_summary
) |>
  mutate(
    across(
      c(
        median_total_pct_2x, p25_total_pct_2x, p75_total_pct_2x,
        median_excess_pct_2x, p25_excess_pct_2x, p75_excess_pct_2x
      ),
      ~ round(.x, 1)
    )
  )


# 3. Prepare plot data ---------------------------------------------------------

ecoregion_order <- amplification_table |>
  filter(l1_ecoregion != "CONUS") |>
  group_by(l1_ecoregion) |>
  summarise(
    order_value = mean(median_total_pct_2x, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(order_value)) |>
  pull(l1_ecoregion)

plot_data <- amplification_table |>
  mutate(
    l1_ecoregion = factor(
      l1_ecoregion,
      levels = rev(c("CONUS", ecoregion_order))
    ),
    metabolism = factor(metabolism, levels = flux_levels)
  )

write_csv(
  amplification_table,
  file.path(dir_data, "09e_figure4_metabolic_amplification_summary.csv")
)


# 4. Plot ----------------------------------------------------------------------

dodge <- position_dodge(width = 0.75)

p_figure4 <- plot_data |>
  ggplot(aes(
    x = l1_ecoregion,
    y = median_total_pct_2x,
    fill = metabolism,
    color = metabolism
  )) +
  geom_col(
    position = dodge,
    width = 0.65,
    alpha = 0.45
  ) +
  geom_errorbar(
    aes(
      ymin = p25_total_pct_2x,
      ymax = p75_total_pct_2x
    ),
    position = dodge,
    width = 0.18,
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = 100,
    linetype = "dashed",
    linewidth = 0.35
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  scale_fill_manual(values = flux_colors, name = NULL) +
  scale_color_manual(values = flux_colors, guide = "none") +
  labs(
    x = NULL,
    y = "Increase in cumulative metabolism when watershed area doubles"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 10)
  )

p_figure4


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figs, "09e_figure4_metabolic_amplification.png"),
  p_figure4,
  width = 8,
  height = 6,
  dpi = 400
)