# 00_setup.R
# Monthly metabolism scaling: proof-of-concept setup

# 1. Packages ---------------------------------------------------------------

require(pacman)

pacman::p_load(tidyverse,
               sf,
               here,
               janitor)


# 2. Paths ------------------------------------------------------------------

dir_project <- here::here()

dir_data <- file.path(dir_project, "data")
dir_raw <- file.path(dir_data, "raw")
dir_interim <- file.path(dir_data, "interim")
dir_processed <- file.path(dir_data, "processed")

dir_figures <- file.path(dir_project, "figures")
dir_outputs <- file.path(dir_project, "outputs")

dir.create(dir_raw, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_interim, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_outputs, recursive = TRUE, showWarnings = FALSE)


# 3. Input files ------------------------------------------------------------

path_huc8_raw <- file.path(dir_raw, "huc8_conus")
path_nep_raw <- file.path(dir_raw, "monthly_nep_flowlines")

path_huc8_clean <- file.path(dir_processed, "huc8_conus.gpkg")
path_flowlines_clean <- file.path(dir_processed, "monthly_nep_flowlines.rds")


# 4. Proof-of-concept settings ---------------------------------------------

target_huc8 <- NA_character_

target_crs <- 5070   # NAD83 / Conus Albers
month_levels <- month.name


# 5. Plot settings ----------------------------------------------------------

theme_set(
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
)

scale_fill_nep <- function(...) {
  scale_fill_viridis_c(option = "magma", ...)
}

scale_color_nep <- function(...) {
  scale_color_viridis_c(option = "magma", ...)
}