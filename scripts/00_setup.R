# 00_setup.R
# Monthly metabolism scaling: proof-of-concept setup

# 1. Packages ---------------------------------------------------------------

require(pacman)

pacman::p_load(tidyverse,
               sf,
               here,
               ggthemes,
               janitor)


# 2. Paths ------------------------------------------------------------------

dir_dryad <- file.path("data/dryad_data")


# 3. other standards -----------------------------------------------------------

common_crs <- 4326   # WGS84
