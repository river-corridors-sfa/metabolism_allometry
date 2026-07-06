
# 1. Setup ---------------------------------------------------------------------

rm(list = ls())

source("scripts/00_setup.R")

pacman::p_load(
  janitor,
  scales,
  tidyverse,
  viridis
)

theme_set(theme_bw())

dir_data <- "data/processed"
dir_figs <- "figures"

path_metrics <- file.path(dir_data, "06_huc08_allometry_metrics.csv")
path_l1_lookup <- file.path(dir_data, "08_huc08_l1_ecoregion_lookup.csv")

month_levels <- month.abb
flux_levels <- c("ER", "GPP")


# 2. Read data -----------------------------------------------------------------

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
    )
  ) |>
  filter(!is.na(metabolism), n >= 10, r2 >= 0.8, is.finite(slope))

huc_l1_lookup <- read_csv(path_l1_lookup, show_col_types = FALSE) |>
  clean_names() |>
  mutate(
    huc = str_pad(as.character(huc), width = 8, pad = "0"),
    l1_ecoregion = str_replace_all(l1_ecoregion, "–", "-")
  )


# 3. Summarize monthly median slopes ------------------------------------------

heatmap_data <- metrics |>
  inner_join(huc_l1_lookup, by = "huc") |>
  group_by(l1_ecoregion, metabolism, month) |>
  summarise(
    median_slope = median(slope, na.rm = TRUE),
    n_huc8 = n_distinct(huc),
    .groups = "drop"
  )

ecoregion_order <- heatmap_data |>
  filter(metabolism == "GPP") |>
  group_by(l1_ecoregion) |>
  summarise(mean_slope = mean(median_slope, na.rm = TRUE), .groups = "drop") |>
  arrange(mean_slope) |>
  pull(l1_ecoregion)

heatmap_data <- heatmap_data |>
  mutate(
    l1_ecoregion = factor(l1_ecoregion, levels = ecoregion_order),
    metabolism = factor(metabolism, levels = flux_levels),
    label = round(median_slope, 2)
  )


# 4. Plot ----------------------------------------------------------------------

p_monthly_l1_heatmap <- heatmap_data |>
  drop_na(l1_ecoregion) %>% 
  ggplot(aes(x = month, y = l1_ecoregion, fill = median_slope)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 3) +
  facet_wrap(~ metabolism, nrow = 1) +
  scale_fill_viridis_c(
    option = "C",
    name = "Median\nslope"
  ) +
  labs(
    title = "Monthly median HUC8 allometric slope by Level I ecoregion",
    x = NULL,
    y = NULL
  ) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey85", color = "grey30"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 15)
  )

p_monthly_l1_heatmap


# 5. Save ----------------------------------------------------------------------

ggsave(
  file.path(dir_figs, "09c_Figure2_heatmap.png"),
  p_monthly_l1_heatmap,
  width = 10,
  height = 7,
  dpi = 400
)