# 05_examine_first_huc8_allometry.R

rm(list = ls())

source("scripts/00_setup.R")


# 1. Settings ------------------------------------------------------------------

target_huc8 <- "13020202"

dir.create("figures/05_huc8_allometry", recursive = TRUE, showWarnings = FALSE)


# 2. Read inputs ---------------------------------------------------------------

huc8 <- nhdplusTools::get_huc(id = target_huc8, type = "huc08")

ggplot() + 
  geom_sf(data = huc8) + 
  theme_map() 

flowlines_sf <- read_sf("data/processed/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp") |>
  clean_names() |>
  mutate(
    comid = as.character(comid),
    nxt_dw = as.character(nxt_dw))

catchments <- read_sf("data/dryad_data/MERIT_Hydro_CONUS_catchments/MERIT_Hydro_CONUS_catchments.shp") |>
  clean_names() |>
  st_drop_geometry() |>
  mutate(comid = as.character(comid)) |>
  select(comid, unitarea)

metabolism <- readRDS("data/processed/03_maavara_monthly_metabolism.rds") |>
  mutate(comid = as.character(comid))

upstream_lookup <- readRDS("data/processed/04_upstream_lookup.rds") |>
  mutate(
    downstream_comid = as.character(downstream_comid),
    upstream_comid = as.character(upstream_comid)
  )

flowlines <- flowlines_sf |>
  st_drop_geometry() |>
  left_join(catchments, by = "comid")


# 3. Identify MERIT reaches intersecting HUC8 ----------------------------------

huc8 <- st_transform(huc8, st_crs(flowlines_sf))

huc8_comids <- flowlines_sf |>
  filter(lengths(st_intersects(geometry, huc8)) > 0) |>
  st_drop_geometry() |>
  distinct(comid)

x <- flowlines_sf |>
  filter(lengths(st_intersects(geometry, huc8)) > 0)

ggplot() + 
  geom_sf(data = huc8, fill = "forestgreen", alpha = 0.1) +
  geom_sf(data = x, color = "lightblue") + 
  theme_map() 

huc8_flowlines <- flowlines |>
  semi_join(huc8_comids, by = "comid")

huc8_outlets <- huc8_flowlines |>
  filter(!nxt_dw %in% huc8_flowlines$comid) |>
  arrange(desc(uparea))


focal_outlet <- huc8_outlets |>
  slice(1) |>
  pull(comid)

focal_comids <- upstream_lookup |>
  filter(downstream_comid == focal_outlet) |>
  distinct(comid = upstream_comid)

huc8_outlets |>
  select(comid, nxt_dw, uparea, order, lngthkm) |>
  print(n = 20)


# 4. Area QC -------------------------------------------------------------------

area_check <- upstream_lookup |>
  left_join(
    flowlines |> select(upstream_comid = comid, unitarea),
    by = "upstream_comid") |>
  group_by(downstream_comid) |>
  summarise(
    sum_unitarea = sum(unitarea, na.rm = TRUE),
    .groups = "drop") |>
  left_join(
    flowlines |> select(downstream_comid = comid, uparea),
    by = "downstream_comid") |>
  mutate(
    area_ratio = sum_unitarea / uparea,
    area_qc = if_else(area_ratio < 0.9 | area_ratio > 1.1, "fail", "pass"))

# saveRDS(area_check, "data/processed/05_area_accumulation_qc.rds")


# 5. Assemble cumulative metabolism -------------------------------------------

local_metabolism <- metabolism |>
  semi_join(focal_comids, by = "comid")

x1 <- x %>% 
  left_join(local_metabolism %>% filter(month == "Jan"), by = "comid")

p1 <- ggplot() + 
  geom_sf(data = huc8, fill = "forestgreen", alpha = 0.1) +
  geom_sf(data = x1, aes(color = nep)) + 
  scale_color_viridis_c() + 
  theme_map() + 
  theme(legend.position = "bottom")

cumulative_metabolism <- upstream_lookup |>
  semi_join(focal_comids, by = c("downstream_comid" = "comid")) |>
  semi_join(focal_comids, by = c("upstream_comid" = "comid")) |>
  left_join(
    local_metabolism |> select(upstream_comid = comid, month, gpp, er, nep),
    by = "upstream_comid"
  ) |>
  group_by(downstream_comid, month) |>
  summarise(
    cumulative_gpp = sum(gpp, na.rm = TRUE),
    cumulative_er = sum(er, na.rm = TRUE),
    cumulative_nep = sum(nep, na.rm = TRUE),
    n_upstream_reaches = n_distinct(upstream_comid),
    n_upstream_with_gpp = n_distinct(upstream_comid[!is.na(gpp)]),
    n_upstream_with_er = n_distinct(upstream_comid[!is.na(er)]),
    .groups = "drop"
  ) |>
  left_join(
    flowlines |> select(downstream_comid = comid, uparea, order, lngthkm, unitarea),
    by = "downstream_comid"
  ) |>
  left_join(area_check, by = "downstream_comid") |>
  mutate(in_target_huc8 = downstream_comid %in% huc8_comids$comid)

x2 <- x %>% 
  left_join(cumulative_metabolism %>% 
              filter(month == "Jan") %>% 
              dplyr::select(downstream_comid, cumulative_nep), 
            by = c("comid" = "downstream_comid"))


p1 <- ggplot() + 
  geom_sf(data = huc8, fill = "forestgreen", alpha = 0.1) +
  geom_sf(data = x1, aes(color = nep)) + 
  geom_sf_text_repel(data = x1, aes(label = round(nep, 1)), size = 2) +
  scale_color_viridis_c() + 
  theme_map() + 
  theme(legend.position = "bottom")

p2 <- ggplot() + 
  geom_sf(data = huc8, fill = "forestgreen", alpha = 0.1) +
  geom_sf(data = x2, aes(color = cumulative_nep)) + 
  geom_sf_text_repel(data = x2, 
                     max.overlaps = 15,
                     aes(label = round(cumulative_nep, 1)), 
                     size = 2) +
  scale_color_viridis_c() + 
  theme_map() + 
  theme(legend.position = "bottom")

plot_grid(p1, p2, nrow = 1)
ggsave("figures/05_huc8_allometry/01_comparison.png", width = 8, height = 5)


#saveRDS(cumulative_metabolism, "data/processed/05_first_huc8_cumulative_metabolism.rds")


# 6. QC summaries --------------------------------------------------------------

network_qc <- cumulative_metabolism |>
  distinct(downstream_comid, uparea, order, area_ratio, area_qc, in_target_huc8) |>
  summarise(
    target_huc8 = target_huc8,
    focal_outlet = focal_outlet,
    n_huc8_reaches = nrow(huc8_flowlines),
    n_huc8_outlets = nrow(huc8_outlets),
    n_focal_reaches = n(),
    pct_focal_reaches_in_huc8 = mean(in_target_huc8) * 100,
    max_order = max(order, na.rm = TRUE),
    min_uparea = min(uparea, na.rm = TRUE),
    median_uparea = median(uparea, na.rm = TRUE),
    max_uparea = max(uparea, na.rm = TRUE),
    pct_area_qc_fail = mean(area_qc == "fail", na.rm = TRUE) * 100,
    median_area_ratio = median(area_ratio, na.rm = TRUE)
  )

metabolism_qc <- cumulative_metabolism |>
  group_by(month) |>
  summarise(
    n_reaches = n(),
    n_positive_gpp = sum(cumulative_gpp > 0, na.rm = TRUE),
    n_positive_er = sum(cumulative_er > 0, na.rm = TRUE),
    n_positive_nep = sum(cumulative_nep > 0, na.rm = TRUE),
    min_nep = min(cumulative_nep, na.rm = TRUE),
    median_nep = median(cumulative_nep, na.rm = TRUE),
    max_nep = max(cumulative_nep, na.rm = TRUE),
    min_gpp_coverage = min(n_upstream_with_gpp / n_upstream_reaches, na.rm = TRUE),
    min_er_coverage = min(n_upstream_with_er / n_upstream_reaches, na.rm = TRUE),
    .groups = "drop"
  )

network_qc
metabolism_qc

write_csv(huc8_outlets, "data/processed/05_first_huc8_outlet_candidates.csv")
write_csv(network_qc, "data/processed/05_first_huc8_network_qc.csv")
write_csv(metabolism_qc, "data/processed/05_first_huc8_metabolism_qc.csv")


# 7. Allometric fits -----------------------------------------------------------

allometry_data <- cumulative_metabolism |>
  filter(
    area_qc == "pass",
    uparea > 0,
    cumulative_gpp > 0,
    cumulative_er > 0
  ) |>
  select(downstream_comid, month, uparea, cumulative_gpp, cumulative_er) |>
  pivot_longer(
    c(cumulative_gpp, cumulative_er),
    names_to = "metabolism_var",
    values_to = "cumulative_metabolism"
  )

allometry_fits <- allometry_data |>
  group_by(month, metabolism_var) |>
  summarise(
    n = n(),
    fit = list(lm(log10(cumulative_metabolism) ~ log10(uparea))),
    .groups = "drop"
  ) |>
  mutate(
    intercept = map_dbl(fit, \(x) coef(x)[1]),
    slope = map_dbl(fit, \(x) coef(x)[2]),
    r2 = map_dbl(fit, \(x) summary(x)$r.squared)
  ) |>
  select(-fit)

saveRDS(allometry_fits, "data/processed/05_first_huc8_allometry_fits.rds")
write_csv(allometry_fits, "data/processed/05_first_huc8_allometry_fits.csv")

allometry_fits


# 8. Figures -------------------------------------------------------------------

flowlines_sf |>
  semi_join(huc8_comids, by = "comid") |>
  ggplot() +
  geom_sf(data = huc8, fill = NA, linewidth = 0.7) +
  geom_sf(linewidth = 0.25, alpha = 0.7) +
  geom_sf(
    data = flowlines_sf |> filter(comid == focal_outlet),
    linewidth = 1.2
  ) +
  theme_bw()

ggsave("figures/05_huc8_allometry/05_huc8_flowline_selection.png", width = 7, height = 6)


cumulative_metabolism |>
  distinct(downstream_comid, uparea, n_upstream_reaches, area_qc) |>
  ggplot(aes(uparea, n_upstream_reaches, color = area_qc)) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_d(end = 0.8) +
  labs(x = "Upstream area", y = "Number of upstream reaches", color = "Area QC") +
  theme_bw()

ggsave("figures/05_huc8_allometry/05_network_size_check.png", width = 7, height = 5)


cumulative_metabolism |>
  distinct(downstream_comid, uparea, area_ratio, area_qc) |>
  ggplot(aes(uparea, area_ratio, color = area_qc)) +
  geom_hline(yintercept = 1, linewidth = 0.4) +
  geom_point(alpha = 0.5, size = 0.7) +
  scale_x_log10() +
  scale_color_viridis_d(end = 0.8) +
  coord_cartesian(ylim = c(0, 1.2)) +
  labs(x = "Upstream area", y = "sum(unitarea upstream) / uparea", color = "Area QC") +
  theme_bw()

ggsave("figures/05_huc8_allometry/05_area_qc_check.png", width = 7, height = 5)


allometry_data |>
  ggplot(aes(uparea, cumulative_metabolism)) +
  geom_point(alpha = 0.35, size = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(metabolism_var ~ month, scales = "free_y") +
  labs(x = "Upstream area", y = "Cumulative metabolism") +
  theme_bw()

ggsave("figures/05_huc8_allometry/05_loglog_gpp_er_by_month.png", width = 12, height = 6)


cumulative_metabolism |>
  filter(area_qc == "pass", uparea > 0) |>
  ggplot(aes(uparea, cumulative_nep)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 0.6) +
  scale_x_log10() +
  facet_wrap(~ month, scales = "free_y") +
  labs(x = "Upstream area", y = "Cumulative NEP") +
  theme_bw()

ggsave("figures/05_huc8_allometry/05_nep_signed_by_month.png", width = 12, height = 6)


allometry_fits |>
  ggplot(aes(month, slope, group = metabolism_var, color = metabolism_var)) +
  geom_hline(yintercept = 1, linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line() +
  scale_color_viridis_d(end = 0.8) +
  labs(x = NULL, y = "Log-log allometric slope", color = NULL) +
  theme_bw()

ggsave("figures/05_huc8_allometry/05_monthly_slopes.png", width = 7, height = 5)


allometry_fits |>
  ggplot(aes(month, r2, group = metabolism_var, color = metabolism_var)) +
  geom_hline(yintercept = 0.8, linewidth = 0.4) +
  geom_point(size = 2) +
  geom_line() +
  scale_color_viridis_d(end = 0.8) +
  labs(x = NULL, y = expression(R^2), color = NULL) +
  theme_bw()

ggsave("figures/05_huc8_allometry/05_monthly_r2.png", width = 7, height = 5)