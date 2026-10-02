## 24_figure1_cumulative_flux_allometry.R
##
## Four-panel Figure 1:
##   Top: pooled reach-month cumulative-area relationships
##   Bottom: annual (median monthly) HUC8 allometry slopes
##
## Run from the project root:
## Rscript scripts/24_figure1_cumulative_flux_allometry.R


# 0. Setup ---------------------------------------------------------------------

source("scripts/00_setup.R")

pacman::p_load(patchwork, hexbin)

dir_processed <- "data/processed"
dir_figs <- "figures"

path_scatter_input <- file.path(dir_processed, "18_refit_no_uparea_scatter_input.rds")
path_allometry_results <- file.path(dir_processed, "18_refit_allometry_results.csv")
path_huc8 <- file.path(dir_processed, "04c_huc08_boundaries.gpkg")

output_png <- file.path(dir_figs, "24_figure1_cumulative_flux_allometry.png")
output_pdf <- file.path(dir_figs, "24_figure1_cumulative_flux_allometry.pdf")

# Easy-to-edit figure settings
hex_bins <- 75
figure_width <- 10
figure_height <- 10
map_crs <- 5070


# 1. Read completed allometry outputs ------------------------------------------

# Each row is a screened, routed reach-month value from the no-area refit.
flux_data <- readRDS(path_scatter_input) %>%
  mutate(metabolism = factor(metabolism, levels = c("GPP", "ER")))

# One slope per HUC8: median of the 12 completed monthly fits.
huc8_slopes <- read_csv(path_allometry_results, show_col_types = FALSE) %>%
  mutate(huc = str_pad(as.character(huc), width = 8, pad = "0")) %>%
  filter(model == "without_uparea", metabolism %in% c("GPP", "ER"), is.finite(slope)) %>%
  group_by(huc, metabolism) %>%
  summarise(annual_slope = median(slope), .groups = "drop") %>%
  mutate(metabolism = factor(metabolism, levels = c("GPP", "ER")))

huc8 <- read_sf(path_huc8, quiet = TRUE) %>%
  rename(huc = huc8) %>%
  select(huc) %>%
  st_make_valid() %>%
  st_transform(map_crs)

map_data <- huc8 %>%
  inner_join(huc8_slopes, by = "huc")

slope_limits <- range(map_data$annual_slope)


# 2. Top row: cumulative-area relationships ------------------------------------

scatter_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )

p_scatter <- ggplot(
  flux_data,
  aes(x = cumulative_area_km2, y = cumulative_flux_g_o2_d)
) +
  geom_hex(bins = hex_bins) +
  geom_abline(slope = 1, intercept = 0, 
            linetype = "dashed") + 
  scale_x_log10(labels = scales::label_number_auto()) +
  scale_y_log10(labels = scales::label_number_auto()) +
  scale_fill_viridis_c(
    option = "D",
    trans = "log10",
    name = "Reach-month count") +
  facet_wrap(
    ~metabolism,
    nrow = 1,
    scales = "free_y",
    labeller = as_labeller(c(GPP = "GPP", ER = "ER"))) +
  labs(
    x = expression("Cumulative upstream area (km"^2*")"),
    y = expression("Cumulative metabolism (g O"[2]*" d"^-1*")")) +
  scatter_theme


# 3. Bottom row: annual HUC8 slopes --------------------------------------------

p_map <- ggplot(map_data) +
  geom_sf(aes(fill = annual_slope), colour = NA) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 1,
    limits = slope_limits,
    oob = scales::squish,
    name = expression("Annual HUC8 slope")
  ) +
  facet_wrap(
    ~metabolism,
    nrow = 1,
    labeller = as_labeller(c(GPP = "GPP HUC8 scaling", ER = "ER HUC8 scaling"))
  ) +
  coord_sf(datum = NA) +
  theme_void(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right"
  )


# 4. Combine and save ----------------------------------------------------------

figure_1 <- p_scatter / p_map +
  plot_layout(heights = c(1, 0.75), guides = "collect")

ggsave(output_png, figure_1,
       width = figure_width, height = figure_height, dpi = 320, bg = "white")
ggsave(output_pdf, figure_1,
       width = figure_width, height = figure_height, bg = "white")
