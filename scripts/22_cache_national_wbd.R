#!/usr/bin/env Rscript
# Cache the official USGS national WBD GeoPackage locally.  `curl --continue-at`
# makes this safe to re-run after an interrupted transfer.
suppressPackageStartupMessages(library(fs))

url <- paste0("https://prd-tnm.s3.amazonaws.com/StagedProducts/",
              "Hydrography/WBD/National/GPKG/WBD_National_GPKG.zip")
cache_dir <- "data/cache/wbd_bulk"
zip_path <- file.path(cache_dir, "WBD_National_GPKG.zip")
extract_dir <- file.path(cache_dir, "WBD_National_GPKG")
dir_create(cache_dir, recurse = TRUE)

if (!file_exists(zip_path) || file_size(zip_path) < 5e9) {
  status <- system2("curl", c("--fail", "--location", "--retry", "5",
                               "--continue-at", "-", "--output", zip_path, url))
  if (status != 0) stop("Official WBD download did not complete; rerun this script to resume.")
}

if (!dir_exists(extract_dir)) {
  dir_create(extract_dir)
  unzip(zip_path, exdir = extract_dir)
}

gpkg <- dir_ls(extract_dir, glob = "*.gpkg", recurse = TRUE)
if (length(gpkg) != 1) stop("Expected exactly one GeoPackage after extraction.")
message("Cached official WBD GeoPackage: ", gpkg)
