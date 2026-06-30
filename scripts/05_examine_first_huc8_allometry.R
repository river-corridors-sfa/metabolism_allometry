# 05_examine_first_huc8_allometry.R

rm(list = ls())

source("scripts/00_setup.R")


# 1. Settings ------------------------------------------------------------------

target_huc8 <- "13020202"

dir.create("figures/05_huc8_allometry", recursive = TRUE, showWarnings = FALSE)


# 2. Read inputs ---------------------------------------------------------------

## Gives us the HUC boundary
huc8 <- nhdplusTools::get_huc(id = target_huc8, type = "huc08")

## Gives us river linestrings and some metrics
flowlines_sf <- read_sf("data/processed/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp") |>
  clean_names() |>
  mutate(
    comid = as.character(comid),
    nxt_dw = as.character(nxt_dw))

## Gives us area for each catchment
catchments <- read_sf("data/dryad_data/MERIT_Hydro_CONUS_catchments/MERIT_Hydro_CONUS_catchments.shp") |>
  clean_names() |>
  st_drop_geometry() |>
  mutate(comid = as.character(comid)) |>
  select(comid, unitarea)

## Metabolism estimates
metabolism <- readRDS("data/processed/03_maavara_monthly_metabolism.rds") |>
  mutate(comid = as.character(comid))

## Pull in width and length
comid_characteristics <- read_csv("data/dryad_data/scaling_dataset_Qmonth_06082024.csv") %>% 
  clean_names() %>% 
  dplyr::select(comid, unitarea, lengthkm, merit_width, order, merit_qmean) %>% 
  mutate(area_m2 = lengthkm * 1000 *  merit_width, 
         comid = as.character(comid))

## Calculate scaled metabolism - NEP in g O2/m2/day, area (length * width) in m2
## gives g O2/day
metabolism_scaled <- left_join(metabolism, 
                               comid_characteristics, 
                               by = "comid") %>% 
  mutate(nep_gO2_d = nep * area_m2)

## Lookup table to route reaches
upstream_lookup <- readRDS("data/processed/04_upstream_lookup.rds") |>
  mutate(
    downstream_comid = as.character(downstream_comid),
    upstream_comid = as.character(upstream_comid)
  )

##
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

local_metabolism <- metabolism_scaled |>
  semi_join(focal_comids, by = "comid")

x1 <- x %>% 
  left_join(local_metabolism %>% filter(month == "Jan"), by = "comid")

p1 <- ggplot() + 
  geom_sf(data = huc8, fill = "forestgreen", alpha = 0.1) +
  geom_sf(data = x1, aes(color = nep_gO2_d)) + 
  scale_color_viridis_c() + 
  theme_map() + 
  theme(legend.position = "bottom")

cumulative_metabolism <- upstream_lookup |>
  semi_join(focal_comids, by = c("downstream_comid" = "comid")) |>
  semi_join(focal_comids, by = c("upstream_comid" = "comid")) |>
  left_join(
    local_metabolism |> select(upstream_comid = comid, month, gpp, er, nep_gO2_d),
    by = "upstream_comid"
  ) |>
  group_by(downstream_comid, month) |>
  summarise(
    cumulative_gpp = sum(gpp, na.rm = TRUE),
    cumulative_er = sum(er, na.rm = TRUE),
    cumulative_nep_g02_d = sum(nep_gO2_d, na.rm = TRUE),
    n_upstream_reaches = n_distinct(upstream_comid),
    n_upstream_with_gpp = n_distinct(upstream_comid[!is.na(gpp)]),
    n_upstream_with_er = n_distinct(upstream_comid[!is.na(er)]),
    .groups = "drop"
  ) |>
  left_join(
    flowlines |> select(downstream_comid = comid, uparea, order, lngthkm, unitarea),
    by = "downstream_comid"
  ) |>
  left_join(
    area_check |> select(downstream_comid, sum_unitarea, area_ratio, area_qc),
    by = "downstream_comid"
  ) |>
  mutate(in_target_huc8 = downstream_comid %in% huc8_comids$comid)

x2 <- x %>% 
  left_join(cumulative_metabolism %>% 
              filter(month == "Jan") %>% 
              dplyr::select(downstream_comid, cumulative_nep_g02_d), 
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
  geom_sf(data = x2, aes(color = cumulative_nep_g02_d)) + 
  geom_sf_text_repel(data = x2, 
                     max.overlaps = 15,
                     aes(label = round(cumulative_nep_g02_d, 1)), 
                     size = 2) +
  scale_color_viridis_c() + 
  theme_map() + 
  theme(legend.position = "bottom")

plot_grid(p1, p2, nrow = 1)
ggsave("figures/05_huc8_allometry/01_comparison.png", width = 8, height = 5)


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
    n_positive_nep = sum(cumulative_nep_g02_d > 0, na.rm = TRUE),
    min_nep = min(cumulative_nep_g02_d, na.rm = TRUE),
    median_nep = median(cumulative_nep_g02_d, na.rm = TRUE),
    max_nep = max(cumulative_nep_g02_d, na.rm = TRUE),
    min_gpp_coverage = min(n_upstream_with_gpp / n_upstream_reaches, na.rm = TRUE),
    min_er_coverage = min(n_upstream_with_er / n_upstream_reaches, na.rm = TRUE),
    .groups = "drop"
  )


# 7. Allometric fits -----------------------------------------------------------

allometry_data <- cumulative_metabolism |>
  filter(
    area_qc == "pass",
    uparea > 0,
    cumulative_gpp > 0,
    cumulative_er > 0
  ) |>
  select(downstream_comid, month, uparea, cumulative_nep_g02_d) |>
  pivot_longer(
    c(cumulative_nep_g02_d),
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

allometry_fits




