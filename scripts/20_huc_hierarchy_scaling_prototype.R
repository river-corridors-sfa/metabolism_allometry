#!/usr/bin/env Rscript
# Prototype: one parent HUC4, with the same reach-scale allometry refit in its
# nested HUC6, HUC8, HUC10, and HUC12 units.
suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(nhdplusTools); library(broom)
})

huc4_id <- "1802"                 # best local HUC8 coverage in a first screen
month_id <- 7                      # July; retain GPP and ER as separate panels
min_reaches <- 5                  # exploratory lower-HUC fit threshold
min_log10_area_range <- 1
out <- "data/processed/huc_hierarchy_prototype"
cache <- "data/cache/wbd_hierarchy"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

huc8_all <- st_read("data/processed/04c_huc08_boundaries.gpkg", quiet = TRUE) |>
  mutate(huc8 = as.character(huc8), huc4 = str_sub(huc8, 1, 4))
huc4_geometry <- huc8_all |>
  filter(huc4 == huc4_id) |>
  summarise(geometry = st_union(geom))

get_cached_huc <- function(level) {
  level_code <- str_pad(level, 2, pad = "0")
  path <- file.path(cache, paste0("huc", level_code, "_", huc4_id, ".rds"))
  if (!file.exists(path)) {
    x <- get_huc(AOI = huc4_geometry, type = paste0("huc", level_code))
    saveRDS(x, path)
  }
  readRDS(path) |>
    mutate(huc_id = as.character(.data[[paste0("huc", level)]]), level = paste0("HUC", level)) |>
    select(level, huc_id, name, states, geometry)
}

huc_units <- bind_rows(map(c(6, 10, 12), get_cached_huc)) |>
  bind_rows(huc8_all |>
    transmute(level = "HUC8", huc_id = huc8, name = huc8_name, states, geometry = geom)) |>
  filter(st_intersects(geometry, huc4_geometry, sparse = FALSE)[, 1])

# Assign reaches by flowline midpoint; this preserves the HUC8 source assignment
# for a check and avoids attributing a long boundary-crossing flowline by a
# vertex. The allometry inputs retain their project-standard cumulative flux and
# cumulative contributing area.
rivers <- st_read("data/processed/merit_data/03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp", quiet = TRUE) |>
  mutate(comid = as.character(comid), huc = as.character(huc)) |>
  filter(str_sub(huc, 1, 4) == huc4_id)
reach_points <- st_line_sample(st_transform(rivers, 5070), sample = .5) |>
  st_as_sf() |>
  st_transform(4326) |>
  mutate(comid = rivers$comid) |>
  select(comid)

assignments <- map_dfr(unique(huc_units$level), \(lev) {
  units <- filter(huc_units, level == lev)
  st_join(reach_points, units |> select(huc_id), join = st_within, left = FALSE) |>
    st_drop_geometry() |>
    distinct(comid, huc_id) |>
    mutate(level = lev)
})

inputs <- readRDS("data/processed/07_huc08_allometry_input.rds") |>
  filter(month == month_id, cumulative_area_km2 > 0, cumulative_flux_g_o2_d > 0,
         pct_area_missing < 10,
         case_when(metabolism_metric == "GPP" ~ pct_gpp_metabolism_missing < 10,
                   metabolism_metric == "ER" ~ pct_er_metabolism_missing < 10,
                   TRUE ~ FALSE)) |>
  transmute(comid = downstream_comid, metabolism_metric, cumulative_area_km2,
            cumulative_flux_g_o2_d)

fit_slope <- function(d) {
  n <- nrow(d)
  area_range <- diff(range(log10(d$cumulative_area_km2)))
  if (n < min_reaches || area_range < min_log10_area_range) {
    return(tibble(n_reaches = n, log10_area_range = area_range, slope = NA_real_,
                  slope_se = NA_real_, r_squared = NA_real_))
  }
  fit <- lm(log10(cumulative_flux_g_o2_d) ~ log10(cumulative_area_km2), data = d)
  tibble(n_reaches = n, log10_area_range = area_range,
         slope = coef(fit)[2], slope_se = summary(fit)$coefficients[2, 2],
         r_squared = summary(fit)$r.squared)
}

child_fits <- inputs |>
  inner_join(assignments, by = "comid", relationship = "many-to-many") |>
  group_by(level, huc_id, metabolism_metric) |>
  group_modify(~ fit_slope(.x)) |>
  ungroup() |>
  mutate(level = factor(level, levels = c("HUC6", "HUC8", "HUC10", "HUC12")))

parent_fit <- inputs |>
  filter(comid %in% unique(assignments$comid)) |>
  group_by(metabolism_metric) |>
  group_modify(~ fit_slope(.x)) |>
  ungroup() |>
  mutate(level = factor("HUC4", levels = c("HUC4", "HUC6", "HUC8", "HUC10", "HUC12")),
         huc_id = huc4_id)

all_fits <- bind_rows(
  parent_fit |> mutate(level = "HUC4"),
  child_fits |> mutate(level = as.character(level))
) |>
  mutate(level = factor(level, levels = c("HUC4", "HUC6", "HUC8", "HUC10", "HUC12")))
write_csv(all_fits, file.path(out, "huc_hierarchy_slopes.csv"))

summary <- all_fits |>
  group_by(metabolism_metric, level) |>
  summarise(n_units_total = n(), n_valid_fits = sum(!is.na(slope)),
            median_slope = median(slope, na.rm = TRUE), .groups = "drop")
write_csv(summary, file.path(out, "huc_hierarchy_summary.csv"))
write_csv(assignments, file.path(out, "reach_huc_assignments.csv"))

# For comparing many HUC2s, robust quantile ribbons convey both the central
# tendency and changing spread without allowing a few extreme fits to dominate.
gpp_ribbon <- all_fits |>
  filter(metabolism_metric == "GPP", !is.na(slope)) |>
  group_by(level) |>
  summarise(n_units = n(), median_slope = median(slope),
            q10 = quantile(slope, .10), q25 = quantile(slope, .25),
            q75 = quantile(slope, .75), q90 = quantile(slope, .90),
            .groups = "drop") |>
  mutate(level_index = match(as.character(level), c("HUC4", "HUC6", "HUC8", "HUC10", "HUC12")))
write_csv(gpp_ribbon, file.path(out, "gpp_hierarchy_quantiles.csv"))

p_ribbon <- ggplot(gpp_ribbon, aes(level_index, median_slope)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), fill = "#77AADD", alpha = .20) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#77AADD", alpha = .55) +
  geom_line(linewidth = .8, colour = "#225577") +
  geom_point(size = 2.2, colour = "#225577") +
  geom_point(data = filter(gpp_ribbon, level == "HUC4"), shape = 21, size = 3.4,
             fill = "#EE7733", colour = "black", stroke = .35) +
  scale_x_continuous(breaks = 1:5, labels = c("HUC4", "HUC6", "HUC8", "HUC10", "HUC12")) +
  labs(title = paste0("GPP scaling across nested HUC levels: HUC4 ", huc4_id),
       subtitle = "Line = median slope; dark ribbon = interquartile range; light ribbon = 10th–90th percentiles",
       x = NULL, y = "Slope of log10(cumulative GPP) ~ log10(drainage area)") +
  theme_bw() + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave("figures/20_huc_hierarchy_gpp_ribbon.png", p_ribbon, width = 8, height = 4.5, dpi = 240)

plot_data <- filter(all_fits, !is.na(slope))
p <- ggplot() +
  geom_boxplot(data = filter(plot_data, level != "HUC4"),
    aes(level, slope), width = .58, outlier.alpha = .25, fill = "#77AADD") +
  geom_point(data = filter(plot_data, level != "HUC4"),
    aes(level, slope), position = position_jitter(width = .10), alpha = .32, size = .9) +
  geom_point(data = filter(plot_data, level == "HUC4"), aes(level, slope),
    shape = 21, fill = "#EE7733", colour = "black", stroke = .4, size = 3) +
  scale_x_discrete(limits = c("HUC4", "HUC6", "HUC8", "HUC10", "HUC12")) +
  facet_wrap(~ metabolism_metric, nrow = 1) +
  labs(title = paste0("Hierarchical allometry prototype: HUC4 ", huc4_id),
       subtitle = "July cumulative reach flux ~ cumulative drainage area; HUC4 is one point, nested units are distributions",
       x = NULL, y = "Slope of log10(cumulative flux) ~ log10(drainage area)") +
  theme_bw() + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave("figures/20_huc_hierarchy_scaling_prototype.png", p, width = 10, height = 5, dpi = 240)

message("HUC4 ", huc4_id, " prototype written to ", out)
