# 04_build_merit_flow_network.R

rm(list = ls())

source("scripts/00_setup.R")


# 1. Read flowline attributes --------------------------------------------------

flowlines <- read_sf("data/processed/riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS.shp") |>
  clean_names() |>
  st_drop_geometry() |>
  mutate(comid = as.character(comid),
    nxt_dw = as.character(nxt_dw))


# 2. Build downstream edge table -----------------------------------------------

edges_downstream <- flowlines |>
  transmute(
    from_comid = comid,
    to_comid = nxt_dw
  ) |>
  filter(
    !is.na(to_comid),
    to_comid != "0",
    to_comid %in% flowlines$comid
  )

saveRDS(edges_downstream, "data/processed/04_edges_downstream.rds")


# 3. Build upstream lookup table ------------------------------------------------

upstream_lookup <- flowlines |>
  transmute(
    downstream_comid = comid,
    upstream_comid = comid
  )

new_upstream <- upstream_lookup
i <- 1

# This loop takes a very long time, potentially about an hour
while (nrow(new_upstream) > 0) {
  
  message("Network step ", i, ": ", nrow(new_upstream), " new upstream links")
  
  new_upstream <- new_upstream |>
    left_join(
      edges_downstream,
      by = c("upstream_comid" = "to_comid")
    ) |>
    filter(!is.na(from_comid)) |>
    transmute(
      downstream_comid,
      upstream_comid = from_comid
    ) |>
    anti_join(upstream_lookup, by = c("downstream_comid", "upstream_comid"))
  
  upstream_lookup <- bind_rows(upstream_lookup, new_upstream)
  
  i <- i + 1
}

saveRDS(upstream_lookup, "data/processed/04_upstream_lookup.rds")


# 4. Basic network checks -------------------------------------------------------

## n_edges = 172,629 and n_from_comids = 172,629 means every routed reach has 
## exactly one downstream neighbor in the retained CONUS network. 
## n_to_comids = 86,466 means many reaches share downstream targets, consistent 
## with tributary merging. The difference between 175,116 flowlines and 172,629 
## edges means 2,487 reaches do not route to another retained reach, matching 
## the outlet count below.
edges_downstream |>
  summarise(
    n_edges = n(),
    n_from_comids = n_distinct(from_comid),
    n_to_comids = n_distinct(to_comid))

## Headwaters/outlets check: n_headwaters = 88,650 means ~50.6% of reaches have 
##no upstream reach inside the retained network. That seems high but not 
##necessarily wrong for MERIT: many first-order reaches plus clipped coastal/boundary 
##fragments. n_outlets = 2,487 means these reaches have no downstream reach 
##inside your retained dataset. That is also plausible for CONUS because you 
##have many terminal coastal outlets and edge-of-domain exits. The key point is 
##that n_outlets exactly equals n_flowlines - n_edges, so this check is 
##internally consistent.
flowlines |>
  summarise(
    n_flowlines = n(),
    n_headwaters = sum(!comid %in% edges_downstream$to_comid),
    n_outlets = sum(!comid %in% edges_downstream$from_comid)
  )

## Upstream lookup check: n_downstream_comids = 175,116 confirms every flowline 
## got an upstream membership set, including itself. min = 1 is correct because 
## each reach includes itself. median = 1 means at least half of all reaches 
## have no upstream retained reaches, which agrees with the large headwater count. 
## max = 63,993 means the largest downstream reach accumulates ~64k upstream reaches; 
## that is plausible for a major river outlet and indicates the recursion is 
## working beyond immediate neighbors.
upstream_lookup |>
  count(downstream_comid, name = "n_upstream_reaches") |>
  summarise(
    n_downstream_comids = n(),
    min_upstream_reaches = min(n_upstream_reaches),
    median_upstream_reaches = median(n_upstream_reaches),
    max_upstream_reaches = max(n_upstream_reaches))








# 5. Compare accumulated unitarea to MERIT uparea ------------------------------

catchments <- read_sf("data/dryad_data/MERIT_Hydro_CONUS_catchments/MERIT_Hydro_CONUS_catchments.shp") |>
  clean_names() |>
  st_drop_geometry() |>
  mutate(comid = as.character(comid)) |>
  select(comid, unitarea)

flowlines <- flowlines |>
  left_join(catchments, by = "comid")

area_check <- upstream_lookup |>
  left_join(
    flowlines |> select(upstream_comid = comid, unitarea),
    by = "upstream_comid"
  ) |>
  group_by(downstream_comid) |>
  summarise(
    sum_unitarea = sum(unitarea, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    flowlines |> select(downstream_comid = comid, uparea),
    by = "downstream_comid"
  ) |>
  mutate(area_ratio = sum_unitarea / uparea)

#median_area_ratio = 1, p05 = 1, p95 = 1: For at least 90% of COMIDs, 
# sum(unitarea of all upstream reaches) exactly reproduces uparea. That validates 
# both the nxt_dw routing and the recursive upstream_lookup. 
## n_bad_ratio = 1,824 out of 175,116 - Only ~1.04% are outside your tolerance
area_check |>
  summarise(
    n = n(),
    median_area_ratio = median(area_ratio, na.rm = TRUE),
    p05_area_ratio = quantile(area_ratio, 0.05, na.rm = TRUE),
    p95_area_ratio = quantile(area_ratio, 0.95, na.rm = TRUE),
    n_bad_ratio = sum(area_ratio < 0.9 | area_ratio > 1.1, na.rm = TRUE))


## This is still mostly good, but the failures are not random rounding noise. 
## They are cases where your routed lookup under-accumulates area relative to uparea.
## Likely cause: clipped-domain / boundary inflow - These are probably reaches 
## whose uparea includes drainage area outside your retained CONUS flowline set, 
## or whose upstream tributaries were removed during trimming. Your upstream_lookup 
## can only sum retained COMIDs; uparea may still reflect the full MERIT upstream drainage.
area_check |>
  filter(area_ratio < 0.9 | area_ratio > 1.1) |>
  summarise(
    n = n(),
    min_ratio = min(area_ratio, na.rm = TRUE),
    median_ratio = median(area_ratio, na.rm = TRUE),
    max_ratio = max(area_ratio, na.rm = TRUE),
    median_uparea = median(uparea, na.rm = TRUE),
    median_sum_unitarea = median(sum_unitarea, na.rm = TRUE)
  )


# 6. Check that upstream area increases downstream -----------------------------

## Should be 0s (except n_edges) - good news!
edges_downstream |>
  left_join(
    flowlines |> select(from_comid = comid, uparea_from = uparea),
    by = "from_comid"
  ) |>
  left_join(
    flowlines |> select(to_comid = comid, uparea_to = uparea),
    by = "to_comid"
  ) |>
  summarise(
    n_edges = n(),
    n_area_decreases = sum(uparea_to < uparea_from, na.rm = TRUE),
    pct_area_decreases = mean(uparea_to < uparea_from, na.rm = TRUE) * 100
  )


