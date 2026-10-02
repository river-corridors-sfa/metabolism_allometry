## 25_figure2_ecoregion_slopes.R
##
## Four-panel Figure 2:
##   Top: HUC8 slope distributions by L1 ecoregion
##   Bottom: mean monthly slopes by L1 ecoregion
##
## Run from the project root:
## Rscript scripts/25_figure2_ecoregion_slopes.R


# 0. Setup ---------------------------------------------------------------------

source("scripts/00_setup.R")

pacman::p_load(patchwork)

dir_processed <- "data/processed"
dir_figs <- "figures"

path_allometry_results <- file.path(dir_processed, "18_refit_allometry_results.csv")
path_ecoregions <- file.path(dir_processed, "04b_huc08_l1_ecoregion_lookup.csv")

output_png <- file.path(dir_figs, "25_figure2_ecoregion_slopes.png")
output_pdf <- file.path(dir_figs, "25_figure2_ecoregion_slopes.pdf")

# Easy-to-edit figure settings
figure_width <- 12
figure_height <- 10


# 1. Read completed HUC8 slopes and ecoregions ---------------------------------

ecoregions <- read_csv(path_ecoregions, show_col_types = FALSE) %>%
  mutate(huc = str_pad(as.character(huc), width = 8, pad = "0"))

slopes <- read_csv(path_allometry_results, show_col_types = FALSE) %>%
  mutate(huc = str_pad(as.character(huc), width = 8, pad = "0")) %>%
  left_join(ecoregions, by = "huc") %>%
  filter(
    model == "without_uparea",
    metabolism %in% c("GPP", "ER"),
    !is.na(l1_ecoregion),
    is.finite(slope)
  )

# Use one ecoregion order in all four panels, from lower to higher median slope.
ecoregion_order <- slopes %>%
  group_by(l1_ecoregion) %>%
  summarise(median_slope = median(slope), .groups = "drop") %>%
  arrange(median_slope) %>%
  pull(l1_ecoregion)

slopes <- slopes %>%
  mutate(
    l1_ecoregion = factor(l1_ecoregion, levels = ecoregion_order),
    metabolism = factor(metabolism, levels = c("GPP", "ER"))
  )

monthly_slopes <- slopes %>%
  group_by(metabolism, l1_ecoregion, month) %>%
  summarise(mean_slope = mean(slope), .groups = "drop") %>%
  mutate(month = factor(month, levels = 1:12, labels = month.abb))

# Use the range of tile values so seasonal and regional differences are visible.
tile_limits <- c(1, max(monthly_slopes$mean_slope))


# 2. Top row: HUC8 slope distributions -----------------------------------------

p_box <- ggplot(slopes, aes(x = slope, y = l1_ecoregion)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_violin(aes(fill = l1_ecoregion), colour = NA, alpha = 0.85) +
  facet_wrap(
    ~metabolism,
    nrow = 1,
    labeller = as_labeller(c(GPP = "GPP", ER = "ER"))
  ) +
  scale_fill_viridis_d(guide = "none") +
  labs(x = expression("HUC8 allometric slope (" * beta * ")"), y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )


# 3. Bottom row: mean monthly slopes -------------------------------------------

p_tile <- ggplot(monthly_slopes, aes(x = month, y = l1_ecoregion, fill = mean_slope)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(
    option = "C",
    limits = tile_limits,
    oob = scales::squish,
    name = expression("Mean slope (" * beta * ")")
  ) +
  facet_wrap(
    ~metabolism,
    nrow = 1,
    labeller = as_labeller(c(GPP = "GPP", ER = "ER"))
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


# 4. Combine and save ----------------------------------------------------------

figure_2 <- p_box / p_tile +
  plot_layout(heights = c(1, 0.75), guides = "collect")

ggsave(output_png, figure_2,
       width = figure_width, height = figure_height, dpi = 320, bg = "white")
ggsave(output_pdf, figure_2,
       width = figure_width, height = figure_height, bg = "white")
