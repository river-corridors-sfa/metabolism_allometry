## 30_allometric_slope_vs_upstream_area.R
##
## Compare each HUC8's annual allometric scaling slope against the largest
## cumulative upstream area represented in that HUC8.
##
## Run from the project root:
## Rscript scripts/30_allometric_slope_vs_upstream_area.R

source("scripts/00_setup.R")

pacman::p_load(hexbin)

dir_processed <- "data/processed"
dir_figs <- "figures"

path_allometry <- file.path(dir_processed, "18_refit_allometry_results.csv")
path_huc8 <- file.path(dir_processed, "04c_huc08_summary.csv")

output_png <- file.path(dir_figs, "30_allometric_slope_vs_upstream_area.png")
output_pdf <- file.path(dir_figs, "30_allometric_slope_vs_upstream_area.pdf")

# Annual slope is the median of the twelve completed monthly HUC8 fits,
# matching the annual-slope definition used in Figure 1.
annual_slopes <- read_csv(path_allometry, show_col_types = FALSE) %>%
  mutate(huc = str_pad(as.character(huc), width = 8, pad = "0")) %>%
  filter(
    model == "without_uparea",
    metabolism %in% c("GPP", "ER"),
    is.finite(slope)
  ) %>%
  group_by(huc, metabolism) %>%
  summarise(annual_slope = median(slope), .groups = "drop")

huc8_area <- read_csv(path_huc8, show_col_types = FALSE) %>%
  transmute(
    huc = str_pad(as.character(huc8), width = 8, pad = "0"),
    max_uparea_km2
  )

plot_data <- annual_slopes %>%
  inner_join(huc8_area, by = "huc") %>%
  filter(
    is.finite(max_uparea_km2), max_uparea_km2 > 0,
    is.finite(annual_slope)
  ) %>%
  mutate(metabolism = factor(metabolism, levels = c("GPP", "ER")))

p <- ggplot(plot_data, aes(x = max_uparea_km2, y = annual_slope)) +
  geom_hex(bins = 45) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  scale_x_log10(labels = scales::label_number_auto()) +
  scale_fill_viridis_c(
    option = "D",
    trans = "log10",
    name = "HUC8 count"
  ) +
  facet_wrap(
    ~metabolism,
    nrow = 1,
    labeller = as_labeller(c(GPP = "GPP", ER = "ER"))
  ) +
  labs(
    x = expression("Maximum cumulative upstream area in HUC8 (km"^2*")"),
    y = expression("Annual HUC8 allometric scaling slope (" * beta * ")")
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave(output_png, p, width = 10, height = 5.5, dpi = 320, bg = "white")
ggsave(output_pdf, p, width = 10, height = 5.5, bg = "white")
