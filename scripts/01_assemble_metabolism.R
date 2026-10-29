## This script assembles the raw Maavara monthly metabolism prediction files into
## one long-form dataset with one row per COMID-month. We keep this step separate
## from spatial processing so later scripts can use a consistent metabolism table
## without repeatedly reading and reshaping the raw monthly CSVs.


# 0. Setup ---------------------------------------------------------------------

source("scripts/00_setup.R")

tic("run 01") # 7s

# 0. Settings ------------------------------------------------------------------

## Paths
dir_er_raw <- "data/processed/17_maavara_refits/without_uparea/ER"
dir_gpp_raw <- "data/processed/17_maavara_refits/without_uparea/GPP"
dir_processed <- "data/processed"
path_baseline <- file.path(dir_processed, "01_assembled_maavara_monthly_metabolism.rds")
path_metabolism <- file.path(dir_processed, "01_assembled_no_uparea_monthly_metabolism.rds")
path_log <- file.path(dir_processed, "01_assembled_no_uparea_monthly_metabolism_log.csv")

dir.create(dir_processed, recursive = TRUE, showWarnings = FALSE)

## map month numbers to month abbreviations
month_lookup <- tibble(
  month_name = month.abb,
  month = 1:12
)


# 1. Identify input files -------------------------------------------------------

## Create lists of files to read
er_files <- list.files(dir_er_raw, pattern = "\\.csv$", full.names = TRUE)
gpp_files <- list.files(dir_gpp_raw, pattern = "\\.csv$", full.names = TRUE)

## Create a tibble for mapping months to ER monthly files
er_file_tbl <- tibble(path = er_files) %>%
  mutate(month_name = str_extract(basename(path), paste(month.abb, collapse = "|"))) %>%
  left_join(month_lookup, by = "month_name") %>%
  arrange(month)

## Create a tibble for mapping months to GPP monthly files
gpp_file_tbl <- tibble(path = gpp_files) %>%
  mutate(month_name = str_extract(basename(path), paste(month.abb, collapse = "|"))) %>%
  left_join(month_lookup, by = "month_name") %>%
  arrange(month)

## Safeties in case months are missing (they're not)
stopifnot(nrow(er_file_tbl) == 12, nrow(gpp_file_tbl) == 12,
          identical(sort(er_file_tbl$month), 1:12),
          identical(sort(gpp_file_tbl$month), 1:12))
stopifnot(!any(is.na(er_file_tbl$month)), !any(is.na(gpp_file_tbl$month)))


# 2. Read and standardize GPP and ER ----------------------------------------------------

## Function for reading each file
read_metabolism_file <- function(path, metabolism_col) {
  month_name <- str_extract(basename(path), paste(month.abb, collapse = "|"))
  process <- if (grepl("^er_", metabolism_col)) "ER" else "GPP"
  expected <- paste0(process, "_pred_meanQ_", month_name)
  raw <- read_csv(path, show_col_types = FALSE, col_types = cols(.default = col_double(), COMID = col_character()))
  if (!identical(names(raw), c("COMID", expected))) stop("Unexpected columns in ", path)
  if (anyNA(raw$COMID) || anyDuplicated(raw$COMID)) stop("Invalid COMID keys in ", path)
  tibble(comid = raw$COMID, month = match(month_name, month.abb),
         !!metabolism_col := raw[[expected]])
}


## Create ER dataset
er_monthly <- er_files %>%
  map(read_metabolism_file, metabolism_col = "er_gO2_m2_d") %>%
  bind_rows()

## Create GPP dataset
gpp_monthly <- gpp_files %>%
  map(read_metabolism_file, metabolism_col = "gpp_gO2_m2_d") %>%
  bind_rows()


# 4. Assemble metabolism --------------------------------------------------------

## Bind together
metabolism <- full_join(gpp_monthly, er_monthly, by = c("comid", "month")) %>%
  arrange(comid, month)
baseline <- readRDS(path_baseline) %>%
  mutate(comid = as.character(comid)) %>% select(comid, month)
if (anyDuplicated(baseline) || anyDuplicated(metabolism[c("comid", "month")]) ||
    nrow(anti_join(baseline, metabolism, by = c("comid", "month"))) ||
    nrow(anti_join(metabolism, baseline, by = c("comid", "month")))) {
  stop("No-uparea and original COMID-month coverage differ; investigate before ingestion")
}
if (any(!metabolism$month %in% 1:12)) stop("Invalid month")

## Diagnostics
metabolism_log <- tibble(
  dataset = c("ER", "GPP", "assembled"),
  n_files = c(nrow(er_file_tbl), nrow(gpp_file_tbl), NA_integer_),
  n_rows = c(nrow(er_monthly), nrow(gpp_monthly), nrow(metabolism)),
  n_comids = c(n_distinct(er_monthly$comid), n_distinct(gpp_monthly$comid), n_distinct(metabolism$comid)),
  n_comid_months = c(n_distinct(paste(er_monthly$comid, er_monthly$month)),
                     n_distinct(paste(gpp_monthly$comid, gpp_monthly$month)),
                     n_distinct(paste(metabolism$comid, metabolism$month))),
  n_missing_gpp = c(NA_integer_, sum(is.na(gpp_monthly$gpp_gO2_m2_d)), sum(is.na(metabolism$gpp_gO2_m2_d))),
  n_missing_er = c(sum(is.na(er_monthly$er_gO2_m2_d)), NA_integer_, sum(is.na(metabolism$er_gO2_m2_d))))


# 5. Checks --------------------------------------------------------------------

metabolism_log

metabolism %>%
  count(month) %>%
  mutate(month_name = month.abb[month]) %>%
  select(month, month_name, n)

metabolism %>%
  summarize(
    n_rows = n(),
    n_comids = n_distinct(comid),
    n_months = n_distinct(month),
    n_duplicate_comid_months = n() - n_distinct(paste(comid, month)),
    pct_missing_gpp = mean(is.na(gpp_gO2_m2_d)) * 100,
    pct_missing_er = mean(is.na(er_gO2_m2_d)) * 100)


# 6. Write outputs -------------------------------------------------------------

saveRDS(metabolism, path_metabolism)
write_csv(metabolism_log, path_log)

toc()
