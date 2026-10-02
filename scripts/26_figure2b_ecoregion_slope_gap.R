## 26_figure2b_ecoregion_slope_gap.R
##
## Figure 2b: monthly mean HUC8 slopes by L1 ecoregion, with a third panel for
## the GPP minus ER slope difference.
##
## Run from the project root:
## Rscript scripts/26_figure2b_ecoregion_slope_gap.R


# 0. Setup ---------------------------------------------------------------------

source("scripts/00_setup.R")

pacman::p_load(patchwork)

dir_processed <- "data/processed"
dir_figs <- "figures"

path_allometry_results <- file.path(dir_processed, "18_refit_allometry_results.csv")
path_ecoregions <- file.path(dir_processed, "04b_huc08_l1_ecoregion_lookup.csv")

output_png <- file.path(dir_figs, "26_figure2b_ecoregion_slope_gap.png")
output_pdf <- file.path(dir_figs, "26_figure2b_ecoregion_slope_gap.pdf")

# Easy-to-edit figure settings
figure_width <- 15
figure_height <- 6
label_size <- 2


# 1. Read completed no-watershed-area results ----------------------------------

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

# Keep the same ecoregion order in all three panels.
ecoregion_order <- slopes %>%
  group_by(l1_ecoregion) %>%
  summarise(median_slope = median(slope), .groups = "drop") %>%
  arrange(median_slope) %>%
  pull(l1_ecoregion)

monthly_slopes <- slopes %>%
  group_by(metabolism, l1_ecoregion, month) %>%
  summarise(mean_slope = mean(slope), .groups = "drop") %>%
  mutate(
    l1_ecoregion = factor(l1_ecoregion, levels = ecoregion_order),
    month = factor(month, levels = 1:12, labels = month.abb)
  )

slope_limits <- range(monthly_slopes$mean_slope)

slope_gap <- monthly_slopes %>%
  select(metabolism, l1_ecoregion, month, mean_slope) %>%
  pivot_wider(names_from = metabolism, values_from = mean_slope) %>%
  mutate(gpp_minus_er = GPP - ER)

# Symmetric limits make red and blue directly comparable. The 98th percentile
# retains contrast while preventing one extreme cell from washing out the map.
gap_limit <- quantile(abs(slope_gap$gpp_minus_er), 0.98, na.rm = TRUE)


# 2. Tile-plot helper -----------------------------------------------------------

tile_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    legend.position = "bottom"
  )

plot_slope_tiles <- function(metric) {
  ggplot(filter(monthly_slopes, metabolism == metric),
         aes(x = month, y = l1_ecoregion, fill = mean_slope)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", mean_slope)), size = label_size) +
    scale_fill_viridis_c(
      option = "C",
      limits = slope_limits,
      name = expression("Mean slope (" * beta * ")")
    ) +
    labs(title = metric) +
    tile_theme
}


# 3. Three panels --------------------------------------------------------------

p_er <- plot_slope_tiles("ER")
p_gpp <- plot_slope_tiles("GPP")

p_gap <- ggplot(slope_gap, aes(x = month, y = l1_ecoregion, fill = gpp_minus_er)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.2f", gpp_minus_er)), size = label_size) +
  scale_fill_gradient2(
    low = "#D7191C",
    mid = "#BDBDBD",
    high = "#0571B0",
    midpoint = 0,
    limits = c(-gap_limit, gap_limit),
    oob = scales::squish,
    name = expression("GPP - ER (" * beta * ")")
  ) +
  labs(title = "GPP - ER") +
  tile_theme


# 4. Combine and save ----------------------------------------------------------

figure_2b <- p_er + p_gpp + p_gap +
  plot_layout(nrow = 1, guides = "collect")

ggsave(output_png, figure_2b,
       width = figure_width, height = figure_height, dpi = 320, bg = "white")
ggsave(output_pdf, figure_2b,
       width = figure_width, height = figure_height, bg = "white")
