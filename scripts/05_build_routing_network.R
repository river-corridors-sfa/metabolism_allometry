## This script builds the routing tables. Since it's a pretty straightforward 
## (but heavy) step, we'll do it in a regular script. 

# 0. Setup ---------------------------------------------------------------------

source("scripts/00_setup.R")

tic("run 05") #2250s, 37.5 minutes

pacman::p_load(knitr)

dir_merit <- "data/processed/merit_data"

## Input paths
path_cat <- file.path(dir_merit, "03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp")
path_riv <- file.path(dir_merit, "03_riv_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.shp")

## Build output paths
path_edges_downstream <- file.path(dir_merit, "05_edges_downstream.rds")
path_upstream_lookup <- file.path(dir_merit, "05_upstream_lookup.rds")


# 1. Read in data --------------------------------------------------------------

tic("Read filtered MERIT data") #11s

merit_cat <- read_sf(path_cat) %>%
  clean_names() %>%
  mutate(comid = as.character(comid))

merit_riv <- read_sf(path_riv) %>%
  clean_names() %>%
  st_drop_geometry() %>%
  mutate(
    comid = as.character(comid),
    nxt_dw = as.character(nxt_dw))

toc()

## Safety checks
stopifnot("comid" %in% names(merit_riv))
stopifnot("nxt_dw" %in% names(merit_riv))


# 2. Build downstream lookup ---------------------------------------------------

tic("Build downstream edge table")

edges_downstream <- merit_riv %>%
  transmute(
    from_comid = comid,
    to_comid = nxt_dw
  ) %>%
  filter(
    !is.na(to_comid),
    to_comid != "0",
    to_comid %in% merit_riv$comid
  )

base::saveRDS(edges_downstream, path_edges_downstream)
toc()


# 3. Build upstream lookup -----------------------------------------------------

tic("Build full upstream lookup")

upstream_lookup <- merit_riv %>%
  transmute(
    downstream_comid = comid,
    upstream_comid = comid
  )

new_upstream <- upstream_lookup
network_step <- 1

tic("Build upstream lookup")

while (nrow(new_upstream) > 0) {
  
  message("Network step ", network_step, ": ", nrow(new_upstream), " new upstream links")
  
  new_upstream <- new_upstream %>%
    left_join(edges_downstream, by = c("upstream_comid" = "to_comid")) %>%
    filter(!is.na(from_comid)) %>%
    transmute(downstream_comid,
              upstream_comid = from_comid
    ) %>%
    anti_join(upstream_lookup, by = c("downstream_comid", "upstream_comid"))
  
  upstream_lookup <- bind_rows(upstream_lookup, new_upstream)
  
  network_step <- network_step + 1
}

base::saveRDS(upstream_lookup, path_upstream_lookup)

toc()


# 4. Basic routing checks ------------------------------------------------------

network_summary <- tibble(
  n_flowlines = nrow(merit_riv),
  n_edges = nrow(edges_downstream),
  n_from_comids = n_distinct(edges_downstream$from_comid),
  n_to_comids = n_distinct(edges_downstream$to_comid),
  n_headwaters = sum(!merit_riv$comid %in% edges_downstream$to_comid),
  n_outlets = sum(!merit_riv$comid %in% edges_downstream$from_comid),
  n_upstream_lookup_rows = nrow(upstream_lookup),
  n_downstream_comids = n_distinct(upstream_lookup$downstream_comid)
)

upstream_summary <- upstream_lookup %>%
  count(downstream_comid, name = "n_upstream_reaches") %>%
  summarise(
    min_upstream_reaches = min(n_upstream_reaches),
    median_upstream_reaches = median(n_upstream_reaches),
    max_upstream_reaches = max(n_upstream_reaches),
    .groups = "drop"
  )

network_summary %>%
  kable()

upstream_summary %>%
  kable()

toc()



