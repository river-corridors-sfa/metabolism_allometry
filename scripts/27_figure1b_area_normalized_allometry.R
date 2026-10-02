## 27_figure1b_area_normalized_allometry.R
##
## Figure 1b: cumulative no-area-model fluxes normalized by cumulative upstream
## area. The y-axis is therefore g O2 km-2 d-1, and mapped slopes are beta - 1.
##
## Run from the project root:
## Rscript scripts/27_figure1b_area_normalized_allometry.R


# 0. Setup ---------------------------------------------------------------------

source("scripts/00_setup.R")

pacman::p_load(patchwork, hexbin)

dir_processed <- "data/processed"
dir_figs <- "figures"

path_scatter_input <- file.path(dir_processed, "18_refit_no_uparea_scatter_input.rds")
path_allometry_results <- file.path(dir_processed, "18_refit_allometry_results.csv")
path_huc8 <- file.path(dir_processed, "04c_huc08_boundaries.gpkg")

output_png <- file.path(dir_figs, "27_figure1b_area_normalized_allometry.png")
output_pdf <- file.path(dir_figs, "27_figure1b_area_normalized_allometry.pdf")

# Easy-to-edit figure settings
hex_bins <- 75
figure_width <- 10
figure_height <- 10
map_crs <- 5070


# 1. Read completed no-watershed-area outputs ----------------------------------

flux_data <- readRDS(path_scatter_input) %>%
  mutate(
    metabolism = factor(metabolism, levels = c("GPP", "ER")),
    cumulative_rate_g_o2_km2_d = cumulative_flux_g_o2_d /
      cumulative_area_km2
  )

# Dividing cumulative flux by cumulative area subtracts one from the log-log
# slope: log10(flux / area) = intercept + (beta - 1) * log10(area).
huc8_slopes <- read_csv(path_allometry_results, show_col_types = FALSE) %>%
  mutate(huc = str_pad(as.character(huc), width = 8, pad = "0")) %>%
  filter(model == "without_uparea", metabolism %in% c("GPP", "ER"), is.finite(slope)) %>%
  group_by(huc, metabolism) %>%
  summarise(area_normalized_slope = median(slope - 1), .groups = "drop") %>%
  mutate(metabolism = factor(metabolism, levels = c("GPP", "ER")))

huc8 <- read_sf(path_huc8, quiet = TRUE) %>%
  rename(huc = huc8) %>%
  select(huc) %>%
  st_make_valid() %>%
  st_transform(map_crs)

map_data <- huc8 %>%
  inner_join(huc8_slopes, by = "huc")

map_limit <- quantile(abs(map_data$area_normalized_slope), 0.98, na.rm = TRUE)


# 2. Top row: area-normalized cumulative flux ----------------------------------

p_scatter <- ggplot(
  flux_data,
  aes(x = cumulative_area_km2, y = cumulative_rate_g_o2_km2_d)
) +
  #geom_point(alpha = 0.01) + 
  geom_hex(bins = hex_bins) +
  #scale_x_log10(labels = scales::label_number_auto()) +
  #scale_y_log10(labels = scales::label_number_auto()) +
  scale_fill_viridis_c(
    option = "D",
    trans = "log10",
    name = "Reach-month count"
  ) +
  facet_wrap(~metabolism, nrow = 1) +
  #geom_abline(slope = 1, intercept = 0) + 
  geom_smooth(method = "lm") + 
  labs(
    x = expression("Cumulative upstream area (km"^2*")"),
    y = expression("Cumulative flux / area (g O"[2]*" km"^-2*" d"^-1*")")
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

p_scatter

# 3. Bottom row: area-normalized HUC8 slopes -----------------------------------

p_map <- ggplot(map_data) +
  geom_sf(aes(fill = area_normalized_slope), colour = NA) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-map_limit, map_limit),
    oob = scales::squish,
    name = expression("Area-normalized HUC8 slope (" * beta - 1 * ")")
  ) +
  facet_wrap(~metabolism, nrow = 1) +
  coord_sf(datum = NA) +
  theme_void(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right"
  )


# 4. Combine and save ----------------------------------------------------------

figure_1b <- p_scatter / p_map +
  plot_layout(heights = c(1, 0.75), guides = "collect")

ggsave(output_png, figure_1b,
       width = figure_width, height = figure_height, dpi = 320, bg = "white")
ggsave(output_pdf, figure_1b,
       width = figure_width, height = figure_height, bg = "white")
