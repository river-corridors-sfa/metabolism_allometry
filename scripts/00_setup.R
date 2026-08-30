# 00_setup.R
# Monthly metabolism scaling: proof-of-concept setup

## Start with a clean environment
rm(list = ls())


# 1. Packages ---------------------------------------------------------------

## Load packages
require(pacman)

pacman::p_load(tidyverse,
               sf,
               here,
               ggthemes,
               nhdplusTools,
               ggsflabel,
               rnaturalearth,
               cowplot,
               tictoc,
               janitor)


# 3. other standards -----------------------------------------------------------

common_crs <- 4326   # WGS84

