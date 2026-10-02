#!/usr/bin/env Rscript
# National extension of the hierarchical GPP-scaling prototype. Each HUC2 is
# fitted independently; child units inherit only the reaches spatially assigned
# to that HUC2's WBD HUC12 polygons.
suppressPackageStartupMessages({ library(tidyverse); library(sf); library(nhdplusTools) })

min_reaches <- 5
min_log10_area_range <- 1
month_id <- 7
out <- "data/processed/huc2_hierarchy_gpp"
cache <- "data/cache/wbd_hierarchy/huc12_by_huc2"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

huc8 <- st_read("data/processed/04c_huc08_boundaries.gpkg", quiet = TRUE) |>
  mutate(huc8 = as.character(huc8), huc2 = str_sub(huc8, 1, 2))
huc2_ids <- huc8 |> st_drop_geometry() |> distinct(huc2) |> arrange(huc2) |> pull(huc2)
requested_huc2 <- commandArgs(trailingOnly = TRUE)
if (length(requested_huc2)) huc2_ids <- intersect(huc2_ids, requested_huc2)

inputs <- readRDS("data/processed/07_huc08_allometry_input.rds") |>
  filter(month == month_id, metabolism_metric == "GPP", cumulative_area_km2 > 0,
         cumulative_flux_g_o2_d > 0, pct_area_missing < 10,
         pct_gpp_metabolism_missing < 10,
         str_sub(huc, 1, 2) %in% huc2_ids) |>
  transmute(comid = downstream_comid, huc2 = str_sub(huc, 1, 2),
            cumulative_area_km2, cumulative_flux_g_o2_d)

rivers <- st_read("data/processed/merit_data/03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp", quiet = TRUE) |>
  mutate(comid = as.character(comid), huc2 = str_sub(as.character(huc), 1, 2))

fit_slope <- function(d) {
  n <- nrow(d); area_range <- diff(range(log10(d$cumulative_area_km2)))
  if (n < min_reaches || area_range < min_log10_area_range)
    return(tibble(n_reaches = n, log10_area_range = area_range, slope = NA_real_, r_squared = NA_real_))
  fit <- lm(log10(cumulative_flux_g_o2_d) ~ log10(cumulative_area_km2), data = d)
  tibble(n_reaches = n, log10_area_range = area_range,
         slope = unname(coef(fit)[2]), r_squared = summary(fit)$r.squared)
}

get_huc12 <- function(huc2_id) {
  path <- file.path(cache, paste0("huc12_", huc2_id, ".rds"))
  if (!file.exists(path)) {
    aoi <- huc8 |> filter(huc2 == huc2_id) |> st_make_valid() |>
      st_transform(5070) |> summarise(geometry = st_union(geom)) |>
      st_simplify(dTolerance = 1000) |> st_transform(4326)
    saveRDS(get_huc(AOI = aoi, type = "huc12"), path)
  }
  readRDS(path) |> st_make_valid() |> transmute(huc12 = as.character(huc12), geometry)
}

assign_one_huc2 <- function(huc2_id) {
  river_part <- filter(rivers, huc2 == huc2_id)
  points <- st_line_sample(st_transform(river_part, 5070), sample = .5) |>
    st_as_sf() |> st_transform(4326) |> mutate(comid = river_part$comid) |> select(comid)
  sf_use_s2(FALSE)
  on.exit(sf_use_s2(TRUE), add = TRUE)
  st_join(points, get_huc12(huc2_id), join = st_within, left = FALSE) |>
    st_drop_geometry() |> distinct(comid, huc12) |>
    transmute(comid, huc2 = huc2_id, huc6 = str_sub(huc12, 1, 6),
              huc8 = str_sub(huc12, 1, 8), huc10 = str_sub(huc12, 1, 10), huc12)
}

# Caches avoid repeating the slow WBD spatial join after an interrupted run.
assignments <- map_dfr(huc2_ids, \(id) {
  path <- file.path(out, paste0("assignments_huc2_", id, ".rds"))
  if (!file.exists(path)) saveRDS(assign_one_huc2(id), path)
  readRDS(path)
})
saveRDS(assignments, file.path(out, "huc2_reach_huc12_assignments.rds"))
write_csv(assignments, file.path(out, "huc2_reach_huc12_assignments.csv"))

child_long <- inputs |>
  inner_join(assignments, by = c("comid", "huc2"), relationship = "many-to-many") |>
  pivot_longer(c(huc6, huc8, huc10, huc12), names_to = "level", values_to = "huc_id") |>
  mutate(level = str_to_upper(level)) |>
  group_by(huc2, level, huc_id) |>
  group_modify(~ fit_slope(.x)) |>
  ungroup()

parent <- inputs |>
  group_by(huc2) |>
  group_modify(~ fit_slope(.x)) |>
  ungroup() |>
  transmute(huc2, level = "HUC2", huc_id = huc2, n_reaches, log10_area_range, slope, r_squared)

fits <- bind_rows(parent, child_long) |>
  mutate(level = factor(level, levels = c("HUC2", "HUC6", "HUC8", "HUC10", "HUC12"))) |>
  arrange(huc2, level, huc_id)
saveRDS(fits, file.path(out, "huc2_hierarchy_gpp_fits.rds"))
write_csv(fits, file.path(out, "huc2_hierarchy_gpp_fits.csv"))

quantiles <- fits |>
  filter(!is.na(slope)) |>
  group_by(huc2, level) |>
  summarise(n_units = n(), median_slope = median(slope), q25 = quantile(slope, .25),
            q75 = quantile(slope, .75), .groups = "drop") |>
  mutate(level_index = match(as.character(level), c("HUC2", "HUC6", "HUC8", "HUC10", "HUC12")))
saveRDS(quantiles, file.path(out, "huc2_hierarchy_gpp_quantiles.rds"))
write_csv(quantiles, file.path(out, "huc2_hierarchy_gpp_quantiles.csv"))

# Requested overlay: all parent points are identical; each HUC2's center line
# is distinguishable, while a single light ribbon encodes its IQR.
overall_means <- quantiles |>
  group_by(level_index) |>
  summarise(mean_slope = mean(median_slope), .groups = "drop")

p <- ggplot(quantiles, aes(level_index, median_slope, group = huc2)) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = huc2), alpha = .10, colour = NA) +
  geom_line(aes(colour = huc2), linewidth = .65, alpha = .85) +
  geom_point(data = filter(quantiles, level == "HUC2"), colour = "#333333", size = 1.8) +
  geom_smooth(data = overall_means, aes(level_index, mean_slope), inherit.aes = FALSE,
              method = "lm", se = FALSE, colour = "black", linetype = "dashed",
              linewidth = .85) +
  scale_x_continuous(breaks = 1:5, labels = c("HUC2", "HUC6", "HUC8", "HUC10", "HUC12")) +
  scale_colour_viridis_d(end = .9, guide = guide_legend(ncol = 3)) +
  scale_fill_viridis_d(end = .9, guide = "none") +
  labs(title = "GPP scaling across nested HUC levels", subtitle =
    "Lines = HUC2 median slopes; light ribbons = HUC2 interquartile ranges; black dash = LM through across-HUC2 means",
    x = NULL, y = "Slope of log10(cumulative GPP) ~ log10(drainage area)", colour = "HUC2") +
  theme_bw() + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
    legend.position = "bottom")
ggsave("figures/21_huc2_hierarchy_gpp_overlay.png", p, width = 11, height = 6.5, dpi = 240)

message("Wrote hierarchy results for HUC2 ", paste(huc2_ids, collapse = ", "), " to ", out)
