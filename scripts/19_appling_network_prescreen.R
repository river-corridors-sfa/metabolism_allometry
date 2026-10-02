#!/usr/bin/env Rscript
# Provisional Appling network screen.  This is deliberately not a substitute
# for NHDPlus/NLDI verification; see the report it writes for the limitation.

suppressPackageStartupMessages({
  library(tidyverse); library(janitor); library(igraph); library(lubridate)
  library(sf); library(cowplot)
})

dir.create("data/processed/appling_network_screen", recursive = TRUE, showWarnings = FALSE)
out <- "data/processed/appling_network_screen"
input <- "data/raw/dryad_data/appling_data_06202024.csv"
edges_file <- "data/processed/merit_data/05_edges_downstream.rds"
stopifnot(file.exists(input), file.exists(edges_file))

# Appling's released daily table carries central estimates, 95% intervals,
# split R-hat, and n_eff.  `usable` is a transparent day-level screen, not the
# unavailable release-level High/Medium/Low confidence label.  `log_usable`
# additionally requires positive GPP and negative ER (ER is signed negative).
daily <- read_csv(input, show_col_types = FALSE) |>
  clean_names() |>
  transmute(
    site_id = site_name, date = mdy(date), gpp, er,
    gpp_rhat, gpp_n_eff, er_rhat, er_n_eff,
    longitude, latitude, merit_comid = as.character(comid),
    drainage_area_km2 = uparea
  ) |>
  mutate(
    usable = if_all(c(gpp, er, gpp_rhat, gpp_n_eff, er_rhat, er_n_eff), is.finite) &
      gpp_rhat <= 1.1 & er_rhat <= 1.1,
    log_usable = usable & gpp > 0 & er < 0,
    er_magnitude = -er
  )

checks <- daily |>
  summarise(
    input_rows = n(), sites = n_distinct(site_id),
    duplicate_site_dates = sum(duplicated(paste(site_id, date))),
    missing_dates = sum(is.na(date)),
    finite_gpp_er = sum(is.finite(gpp) & is.finite(er)),
    usable_day_rows = sum(usable), log_usable_day_rows = sum(log_usable),
    positive_er_rows = sum(er > 0, na.rm = TRUE), negative_gpp_rows = sum(gpp < 0, na.rm = TRUE)
  )
write_csv(checks, file.path(out, "input_checks.csv"))
stopifnot(checks$duplicate_site_dates == 0, checks$missing_dates == 0)

sites <- daily |>
  distinct(site_id, merit_comid, longitude, latitude, drainage_area_km2) |>
  arrange(site_id)
stopifnot(nrow(sites) == n_distinct(daily$site_id))

# This is only a pre-screen: MERIT COMIDs identify the reaches and `nxt_dw`
# supplies a directed MERIT routing graph.  It is not an NHDPlus COMID join.
edges <- readRDS(edges_file) |>
  transmute(from = as.character(from_comid), to = as.character(to_comid)) |>
  distinct()
g <- graph_from_data_frame(edges, directed = TRUE)
member <- components(as_undirected(g, mode = "collapse"))$membership
sites <- sites |>
  mutate(merit_component = unname(member[merit_comid]),
         routing_prescreen = !is.na(merit_component))
write_csv(sites, file.path(out, "gauge_metadata_prescreen.csv"))

# Directed reachability is calculated only between the observed sites. It
# distinguishes a pair on a directed path from tributaries that only share a
# downstream confluence.
site_vertices <- sites |>
  filter(routing_prescreen) |>
  pull(merit_comid) |>
  unique()
d <- distances(g, v = site_vertices, to = site_vertices, mode = "out")
relationships <- as_tibble(as.table(d), .name_repair = "minimal") |>
  setNames(c("upstream_merit_comid", "downstream_merit_comid", "routing_distance")) |>
  filter(is.finite(routing_distance), routing_distance > 0) |>
  left_join(sites |> select(site_id, upstream_merit_comid = merit_comid),
            by = "upstream_merit_comid", relationship = "many-to-many") |>
  rename(upstream_site_id = site_id) |>
  left_join(sites |> select(site_id, downstream_merit_comid = merit_comid),
            by = "downstream_merit_comid", relationship = "many-to-many") |>
  rename(downstream_site_id = site_id) |>
  transmute(upstream_site_id, downstream_site_id, upstream_merit_comid,
            downstream_merit_comid, routing_distance,
            evidence = "MERIT nxt_dw routing pre-screen; NHDPlus/NLDI verification pending")
write_csv(relationships, file.path(out, "merit_prescreen_relationships.csv"))

# For every calendar day, retain the complete simultaneously log-usable set
# within each weakly connected component.  Collapsing identical sets gives all
# observed candidate subsets, rather than requiring every gauge in a component.
sets <- daily |>
  filter(log_usable) |>
  inner_join(sites |> select(site_id, merit_component, routing_prescreen), by = "site_id") |>
  filter(routing_prescreen) |>
  group_by(merit_component, date) |>
  summarise(site_set = str_c(sort(unique(site_id)), collapse = " | "),
            n_gauges = n_distinct(site_id), .groups = "drop") |>
  group_by(merit_component, site_set, n_gauges) |>
  arrange(date, .by_group = TRUE) |>
  mutate(run = cumsum(c(TRUE, diff(date) != 1))) |>
  summarise(shared_log_usable_days = n(),
            first_shared_date = min(date), last_shared_date = max(date),
            longest_uninterrupted_days = max(as.integer(table(run))), .groups = "drop") |>
  separate_rows(site_set, sep = " \\| ") |>
  left_join(sites |> select(site_id, drainage_area_km2), by = c("site_set" = "site_id")) |>
  group_by(merit_component, n_gauges, shared_log_usable_days,
           first_shared_date, last_shared_date, longest_uninterrupted_days) |>
  summarise(gauge_ids = str_c(sort(site_set), collapse = " | "),
            min_drainage_area_km2 = min(drainage_area_km2),
            max_drainage_area_km2 = max(drainage_area_km2),
            drainage_area_ratio = max_drainage_area_km2 / min_drainage_area_km2,
            .groups = "drop") |>
  arrange(desc(n_gauges), desc(longest_uninterrupted_days), desc(drainage_area_ratio))
write_csv(sets, file.path(out, "candidate_subsets_merit_prescreen.csv"))

# A component summary guards against interpreting overlapping subsets as three
# independent networks. No row is labelled a verified winner in this script.
components <- sites |>
  filter(routing_prescreen) |>
  count(merit_component, name = "component_gauges") |>
  left_join(sets |> group_by(merit_component) |> slice_max(n_gauges, n = 1,
    with_ties = FALSE) |> ungroup(), by = "merit_component") |>
  arrange(desc(n_gauges), desc(longest_uninterrupted_days), desc(drainage_area_ratio))
write_csv(components, file.path(out, "component_rankings_merit_prescreen.csv"))

unresolved <- sites |>
  filter(!routing_prescreen) |>
  transmute(site_id, merit_comid, reason = "MERIT COMID absent from local nxt_dw graph; needs NHDPlus/NLDI match")
write_csv(unresolved, file.path(out, "unresolved_sites.csv"))

# Top three are distinct routing components, so overlapping subsets are never
# presented as separate networks. Their gauge-level detail is retained exactly
# as used in the simultaneous-date calculation.
top_three <- components |>
  filter(component_gauges >= 2) |>
  slice_head(n = 3) |>
  mutate(network_rank = row_number())
write_csv(top_three, file.path(out, "top_three_merit_networks.csv"))

top_gauges <- top_three |>
  select(network_rank, merit_component, gauge_ids) |>
  separate_rows(gauge_ids, sep = " \\| ") |>
  rename(site_id = gauge_ids) |>
  left_join(sites, by = "site_id") |>
  left_join(daily |>
    group_by(site_id) |>
    summarise(first_daily_date = min(date), last_daily_date = max(date),
              usable_days = sum(usable), log_usable_days = sum(log_usable), .groups = "drop"),
    by = "site_id") |>
  arrange(network_rank, drainage_area_km2)
write_csv(top_gauges, file.path(out, "top_three_gauge_details.csv"))

# Map actual local-MERIT paths, not straight-line connectors. For each observed
# set, retain only nearest observed upstream/downstream links; their shortest
# directed graph paths are then extracted from the local flowline layer.
direct_links <- function(ids) {
  ids <- sites |> filter(site_id %in% ids) |> pull(merit_comid) |> unique()
  rel <- relationships |>
    filter(upstream_merit_comid %in% ids, downstream_merit_comid %in% ids) |>
    distinct(upstream_merit_comid, downstream_merit_comid)
  rel |> filter(!map2_lgl(upstream_merit_comid, downstream_merit_comid, \(u, v) {
    any(ids != u & ids != v &
      is.finite(d[u, ids]) & d[u, ids] > 0 &
      is.finite(d[ids, v]) & d[ids, v] > 0)
  }))
}
path_ids <- top_three |>
  select(network_rank, gauge_ids) |>
  mutate(gauge_ids = str_split(gauge_ids, " \\| "),
         links = map(gauge_ids, direct_links),
         path_comids = map(links, \(z) unique(unlist(map2(z$upstream_merit_comid,
           z$downstream_merit_comid, \(u, v) names(shortest_paths(g, u, to = v)$vpath[[1]]))))))

riv_file <- "data/processed/merit_data/03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp"
if (file.exists(riv_file)) {
  riv_all <- st_read(riv_file, quiet = TRUE) |>
    clean_names() |>
    mutate(comid = as.character(comid), huc8 = str_sub(huc, 1, 8))
  river_paths <- riv_all |>
    inner_join(path_ids |> select(network_rank, path_comids) |> unnest(path_comids) |>
      rename(comid = path_comids), by = "comid")
  p <- map(1:3, \(rank) ggplot() +
    geom_sf(data = filter(river_paths, network_rank == rank), colour = "#2878B5", linewidth = .25) +
    geom_sf(data = st_as_sf(filter(top_gauges, network_rank == rank),
      coords = c("longitude", "latitude"), crs = 4326), colour = "#D55E00", size = 1.4) +
    labs(title = paste("Candidate", rank),
      subtitle = "Blue: local MERIT routed paths; orange: gauges") + theme_void())
  p <- plot_grid(plotlist = p, nrow = 1) |>
    ggdraw() + draw_label("Top local-MERIT routing candidates — not NHDPlus verification",
      x = .5, y = .99, hjust = .5, vjust = 1, fontface = "bold")
  ggsave("figures/19_appling_networks_merit_prescreen.png", p, width = 11, height = 4, dpi = 220)

  # HUC8 context maps use all local flowlines in candidate HUC8s. Stream order
  # controls both alpha and width; traced candidate paths are overlaid in blue.
  huc_file <- "data/processed/04c_huc08_boundaries.gpkg"
  if (file.exists(huc_file)) {
    candidate_hucs <- top_gauges |>
      mutate(merit_comid = as.character(merit_comid)) |>
      left_join(st_drop_geometry(riv_all) |> select(comid, huc8), by = c("merit_comid" = "comid")) |>
      distinct(network_rank, huc8) |>
      filter(!is.na(huc8))
    all_huc8 <- st_read(huc_file, quiet = TRUE) |>
      mutate(huc8 = as.character(huc8)) |>
      mutate(huc4 = str_sub(huc8, 1, 4))
    huc8 <- all_huc8 |>
      inner_join(candidate_hucs, by = "huc8")
    candidate_huc4 <- candidate_hucs |>
      mutate(huc4 = str_sub(huc8, 1, 4)) |>
      distinct(network_rank, huc4)
    huc4 <- all_huc8 |>
      group_by(huc4) |>
      summarise(states = str_c(sort(unique(states)), collapse = ", "), .groups = "drop") |>
      inner_join(candidate_huc4, by = "huc4")
    huc_maps <- map(1:3, \(rank) {
      boundary <- filter(huc8, network_rank == rank)
      context <- st_filter(riv_all, boundary) |>
        mutate(order = pmax(order, 1))
      ggplot() +
        geom_sf(data = boundary, fill = "#F7F7F7", colour = "#777777", linewidth = .25) +
        geom_sf(data = filter(huc4, network_rank == rank), fill = NA, colour = "#222222",
          linewidth = .6, linetype = "dashed") +
        geom_sf(data = context, aes(alpha = order, linewidth = order), colour = "#5C7893") +
        scale_alpha_continuous(range = c(.08, .70), guide = "none") +
        scale_linewidth_continuous(range = c(.08, .75), guide = "none") +
        geom_sf(data = filter(river_paths, network_rank == rank), colour = "#0072B2", linewidth = .5) +
        geom_sf(data = st_as_sf(filter(top_gauges, network_rank == rank),
          coords = c("longitude", "latitude"), crs = 4326), colour = "#D55E00", size = 1.4) +
        labs(title = paste("Candidate", rank), subtitle = "Gray = HUC8; dashed black = HUC4; darker/wider = higher stream order") +
        theme_void()
    })
    huc_plot <- plot_grid(plotlist = huc_maps, nrow = 1) |>
      ggdraw() + draw_label("Candidate HUC8 maps: local MERIT stream order and routed paths",
        x = .5, y = .99, hjust = .5, vjust = 1, fontface = "bold")
    ggsave("figures/19_appling_candidate_huc8_stream_order.png", huc_plot,
      width = 13, height = 4.5, dpi = 220)
    write_csv(candidate_hucs, file.path(out, "top_three_candidate_huc8s.csv"))
    write_csv(st_drop_geometry(huc4) |> arrange(network_rank, huc4),
      file.path(out, "top_three_candidate_huc4s.csv"))

    # Requested focused view: Candidate 2 without its southern-most gauge and
    # that gauge's otherwise single-gauge HUC8. Keep this derivative selection
    # separate from the ranked candidate table.
    focus_all <- top_gauges |>
      filter(network_rank == 2) |>
      mutate(merit_comid = as.character(merit_comid)) |>
      left_join(st_drop_geometry(riv_all) |> select(comid, huc8), by = c("merit_comid" = "comid"))
    excluded <- focus_all |> slice_min(latitude, n = 1, with_ties = FALSE)
    excluded_huc8 <- excluded$huc8
    focus <- focus_all |> filter(site_id != excluded$site_id, huc8 != excluded_huc8)
    focus_boundary <- all_huc8 |> filter(huc8 %in% unique(focus$huc8))
    focus_context <- st_filter(riv_all, focus_boundary) |> mutate(order = pmax(order, 1))
    focus_links <- direct_links(focus$site_id)
    focus_path_comids <- unique(unlist(map2(focus_links$upstream_merit_comid,
      focus_links$downstream_merit_comid, \(u, v) names(shortest_paths(g, u, to = v)$vpath[[1]]))))
    focus_paths <- filter(riv_all, comid %in% focus_path_comids)
    focus_xlim <- range(focus$longitude) + c(-.25, .25)
    focus_ylim <- range(focus$latitude) + c(-.20, .20)
    p_focus <- ggplot() +
      geom_sf(data = focus_boundary, fill = "#F7F7F7", colour = "#777777", linewidth = .3) +
      geom_sf(data = focus_context, aes(alpha = order, linewidth = order), colour = "#5C7893") +
      scale_alpha_continuous(range = c(.08, .70), guide = "none") +
      scale_linewidth_continuous(range = c(.08, .75), guide = "none") +
      geom_sf(data = focus_paths, colour = "#0072B2", linewidth = .65) +
      geom_sf(data = st_as_sf(focus, coords = c("longitude", "latitude"), crs = 4326),
        colour = "#D55E00", size = 2) +
      coord_sf(xlim = focus_xlim, ylim = focus_ylim, expand = TRUE) +
      labs(title = "Candidate 2: focused HUC8 context", subtitle = paste0(
        "Excluded southern-most site ", excluded$site_id, " and HUC8 ", excluded_huc8,
        "; darker/wider = higher stream order")) + theme_void()
    ggsave("figures/19_candidate2_focused_huc8_stream_order.png", p_focus,
      width = 7, height = 6, dpi = 240)
    write_csv(bind_rows(
      excluded |> transmute(selection = "excluded", site_id, huc8, latitude, longitude),
      focus |> transmute(selection = "retained", site_id, huc8, latitude, longitude)),
      file.path(out, "candidate2_focused_gauge_selection.csv"))
  }
}

message("Wrote MERIT routing pre-screen to ", out)
