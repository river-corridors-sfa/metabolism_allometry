# 00_setup.R
# Monthly metabolism scaling: proof-of-concept setup

## Start with a clean environment
rm(list = ls())


# 1. Packages ---------------------------------------------------------------

## Load packages
require(pacman)

pacman::p_load(tidyverse, janitor, # keep it tidy
               sf, ggthemes, ggsflabel, # spatial things
               here,
               nhdplusTools, # USGS data
               rnaturalearth, # map spatial features
               cowplot, # combine plots
               tictoc) # time things


# 3. other standards -----------------------------------------------------------

common_crs <- 4326   # WGS84

