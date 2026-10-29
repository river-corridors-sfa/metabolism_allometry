# 08_huc8_fit_stability.R
# Run from the project root after 07. Tests deletion sensitivity, not confidence.
# Sections end before confluences and at HUC8 / retained-network boundaries.
# Cumulative flux and area NEVER change when a section leaves the regression.
# Input from 07 must already exclude negative local rates and reach-months with
# >=10% missing upstream area or metabolism coverage.
# Network dependence: Ver Hoef & Peterson (2010), doi:10.1198/jasa.2009.ap08248;
# Kojevnikov (2021), doi:10.48550/arXiv.2101.12312. This is an adaptation.

local({
  library(tidyverse)
  library(igraph)
  library(progress)

  # 1. Settings ---------------------------------------------------------------
  test_mode <- FALSE
  test_hucs <- NULL                  # Or c("17090001", "17090002", "17090003")
  max_blocks <- 400L                 # Inf = omit every section once
  seed <- 20261002L
  min_reaches <- 20L                 # Same eligibility rules as 07
  min_area_range_log10 <- 1
  missing_area_threshold <- 10       # Must match 07's strict <10% filters
  missing_metabolism_threshold <- 10
  r2_cutoff <- 0.8
  slope_bounds <- c(0.95, 1.05)

  # 400 is a runtime choice, NOT a published river-network cutoff:
  # 1.96 * sqrt(0.25 / 400) = 0.049, approximate worst-case sampling error
  # for a section-level change frequency with 400 valid sampled deletions.
  # This does not bound missed extreme slopes or provide a slope CI.
  dir_processed <- "data/processed"
  dir_merit <- file.path(dir_processed, "merit_data")
  path_output <- file.path(dir_processed, if (test_mode) {
    "08_huc08_fit_stability_test.csv"
  } else {
    "08_huc08_fit_stability.csv"
  })

  # 2. Read only the needed attributes ---------------------------------------
  # RDS must be read in full once; narrow immediately. A fresh R session avoids
  # retaining 07's large spatial tables and expanded upstream lookup alongside it.
  message("Reading cumulative input from 07...")
  input_07 <- readRDS(file.path(dir_processed, "07_huc08_allometry_input.rds"))
  required_07 <- c("huc", "downstream_comid", "month", "metabolism_metric",
                   "cumulative_area_km2", "cumulative_flux_g_o2_d",
                   "pct_area_missing", "pct_metabolism_missing")
  missing_07 <- setdiff(required_07, names(input_07))
  if (length(missing_07)) {
    stop("07 input lacks the new coverage field(s): ",
         paste(missing_07, collapse = ", "), ". Rerun 07 before 08.")
  }
  if (any(!is.finite(input_07$pct_area_missing) |
          input_07$pct_area_missing >= missing_area_threshold) ||
      any(!is.finite(input_07$pct_metabolism_missing) |
          input_07$pct_metabolism_missing >= missing_metabolism_threshold) ||
      any(!is.finite(input_07$cumulative_area_km2) |
          input_07$cumulative_area_km2 <= 0) ||
      any(!is.finite(input_07$cumulative_flux_g_o2_d) |
          input_07$cumulative_flux_g_o2_d <= 0)) {
    stop("07 input violates the positive-flux or <10% coverage rules.")
  }
  input <- input_07 %>%
    transmute(
      huc = str_pad(as.character(huc), 8, pad = "0"),
      comid = as.character(downstream_comid), month, metabolism_metric,
      x = log10(cumulative_area_km2), y = log10(cumulative_flux_g_o2_d),
      pct_metabolism_missing
    )
  rm(input_07)
  stopifnot(all(is.finite(input$x)), all(is.finite(input$y)))

  if (test_mode) {
    if (is.null(test_hucs)) {
      # First 10 HUC8 IDs, selected deterministically.
      test_hucs <- input %>% distinct(huc) %>% arrange(huc) %>%
        slice_head(n = 10) %>% pull(huc)
    }
    input <- input %>% filter(huc %in% str_pad(test_hucs, 8, pad = "0"))
  }
  stopifnot(nrow(input) > 0)
  hucs <- sort(unique(input$huc))
  message("Selected ", length(hucs), " HUC8s; ", max_blocks, " section cap.")

  # 07 takes HUC from the catchment layer; 04c takes it from flowlines.
  # Boundary reaches can therefore have different HUC labels in the two files.
  # Use 07's labels for sections so they match the HUC8s being refit.
  key <- input %>% distinct(comid, huc)
  stopifnot(!anyDuplicated(key$comid))
  retained_comids <- readRDS(file.path(dir_processed, "04c_comid_huc08_key.rds")) %>%
    transmute(comid = as.character(comid))
  stopifnot(!anyDuplicated(retained_comids$comid),
            nrow(anti_join(key, retained_comids, by = "comid")) == 0)
  rm(retained_comids)

  # Count incoming links before cropping: external tributaries still identify
  # confluences. Cutting incoming confluence links starts the downstream section
  # at the confluence reach; no reach appears in two sections.
  message("Defining sections from retained routing...")
  edges <- readRDS(file.path(dir_merit, "05_edges_downstream.rds")) %>%
    transmute(from_comid = as.character(from_comid),
              to_comid = as.character(to_comid))
  stopifnot(!anyDuplicated(edges$from_comid))
  edges <- edges %>%
    add_count(to_comid, name = "n_incoming") %>%
    filter(n_incoming == 1) %>%
    inner_join(key, by = c("from_comid" = "comid")) %>%
    inner_join(key, by = c("to_comid" = "comid"), suffix = c("", "_to")) %>%
    filter(huc == huc_to) %>% select(from_comid, to_comid)
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE,
                                          vertices = key %>% select(comid))
  stopifnot(igraph::is_dag(graph))
  membership <- igraph::components(graph, mode = "weak")$membership
  blocks <- tibble(comid = names(membership), block = unname(membership))
  input <- input %>% left_join(blocks, by = "comid") %>% select(-comid)
  rm(key, edges, graph, membership, blocks)
  invisible(gc())

  # 3. Fit and summarize deletions -------------------------------------------
  classify <- function(slope, r2) {
    case_when(
      !is.finite(slope) | !is.finite(r2) ~ NA_character_,
      r2 < r2_cutoff ~ "low_r2",
      slope < slope_bounds[1] ~ "sublinear",
      slope <= slope_bounds[2] ~ "near_linear",
      TRUE ~ "superlinear"
    )
  }
  direction <- function(slope) {
    case_when(slope < 1 ~ "below_1", slope > 1 ~ "above_1",
              is.finite(slope) ~ "at_1", TRUE ~ NA_character_)
  }

  summarize_fit <- function(data, selected_blocks) {
    n_total <- nrow(data)
    n_unique <- n_distinct(data$x)
    area_range <- diff(range(data$x))
    # Center first to avoid cancellation in sum-of-squares arithmetic.
    data <- data %>% mutate(dx = x - mean(x), dy = y - mean(y))
    section_stats <- data %>% group_by(block) %>% summarise(
      n = n(), sx = sum(dx), sy = sum(dy),
      sxx = sum(dx^2), sxy = sum(dx * dy), syy = sum(dy^2),
      xmin = min(x), xmax = max(x), .groups = "drop"
    )
    totals <- section_stats %>% summarise(across(c(n, sx:syy), sum))
    estimate <- function(n, sx, sy, sxx, sxy, syy) {
      xx <- sxx - sx^2 / n
      xy <- sxy - sx * sy / n
      yy <- syy - sy^2 / n
      tibble(slope = if_else(xx > 0, xy / xx, NA_real_),
             r2 = if_else(xx > 0 & yy > 0, pmin(1, pmax(0, xy^2 / (xx * yy))), NA_real_))
    }
    baseline <- with(totals, estimate(n, sx, sy, sxx, sxy, syy))
    baseline_ok <- n_total >= min_reaches && n_unique >= min_reaches &&
      area_range >= min_area_range_log10 && all(is.finite(unlist(baseline)))

    # Count x values that disappear entirely when a section is removed.
    lost_values <- data %>% distinct(block, x) %>% add_count(x) %>%
      filter(n == 1) %>% count(block, name = "n_unique_lost")
    # The second section minimum/maximum replaces the first when it is omitted.
    lower <- section_stats %>% arrange(xmin) %>% slice_head(n = 2)
    upper <- section_stats %>% arrange(desc(xmax)) %>% slice_head(n = 2)
    refits <- section_stats %>% filter(block %in% selected_blocks) %>%
      left_join(lost_values, by = "block") %>% mutate(
        n_removed = n,
        n_unique_remaining = n_unique - replace_na(n_unique_lost, 0L),
        area_range_remaining =
          if_else(block == upper$block[1], upper$xmax[2], upper$xmax[1]) -
          if_else(block == lower$block[1], lower$xmin[2], lower$xmin[1]),
        across(c(n, sx:syy), ~ totals[[cur_column()]] - .x)
      )
    coefficients <- with(refits, estimate(n, sx, sy, sxx, sxy, syy))
    refits <- bind_cols(refits, coefficients) %>% mutate(
      valid = baseline_ok & n >= min_reaches &
        n_unique_remaining >= min_reaches &
        area_range_remaining >= min_area_range_log10 &
        is.finite(slope) & is.finite(r2)
    )
    valid_refits <- refits %>% filter(valid)
    has_refits <- nrow(valid_refits) > 0
    base_class <- classify(baseline$slope, baseline$r2)
    class_changes <- if (has_refits) {
      mean(classify(valid_refits$slope, valid_refits$r2) != base_class)
    } else NA_real_
    direction_changes <- if (has_refits) {
      mean(direction(valid_refits$slope) != direction(baseline$slope))
    } else NA_real_
    tibble(
      n_reaches = n_total, n_unique_areas = n_unique,
      log10_area_range = area_range,
      median_pct_metabolism_missing = median(data$pct_metabolism_missing),
      max_pct_metabolism_missing = max(data$pct_metabolism_missing),
      baseline_slope = baseline$slope, baseline_r2 = baseline$r2,
      baseline_intercept = mean(data$y) - baseline$slope * mean(data$x),
      baseline_class = base_class,
      baseline_direction = direction(baseline$slope),
      n_blocks = nrow(section_stats), n_tested = nrow(refits),
      n_valid = nrow(valid_refits), n_invalid = sum(!refits$valid),
      n_untested = nrow(section_stats) - nrow(refits),
      slope_min = if (has_refits) min(valid_refits$slope) else NA_real_,
      slope_max = if (has_refits) max(valid_refits$slope) else NA_real_,
      max_abs_slope_change = if (has_refits) {
        max(abs(valid_refits$slope - baseline$slope))
      } else NA_real_,
      r2_min = if (has_refits) min(valid_refits$r2) else NA_real_,
      r2_max = if (has_refits) max(valid_refits$r2) else NA_real_,
      max_removed_fraction = if (nrow(refits)) max(refits$n_removed) / n_total else NA_real_,
      class_change_fraction = class_changes,
      direction_change_fraction = direction_changes,
      status = case_when(
        !baseline_ok ~ "insufficient_baseline",
        !has_refits ~ "insufficient_refits",
        class_changes > 0 | direction_changes > 0 ~ "class_changes_detected",
        TRUE ~ "no_class_changes_detected"
      )
    )
  }

  # 4. Process one HUC8 at a time; write only its small summary ----------------
  # An integer row index avoids repeatedly scanning the full CONUS table.
  rows_by_huc <- split(seq_len(nrow(input)), input$huc)
  n_fits <- input %>% distinct(huc, month, metabolism_metric) %>% nrow()
  bar <- progress::progress_bar$new(
    format = "[:bar] :percent | :current/:total fits | elapsed :elapsed | ETA :eta",
    total = n_fits, clear = FALSE, show_after = 0
  )
  first_write <- TRUE
  walk(hucs, function(huc_i) {
    data <- input[rows_by_huc[[huc_i]], ]
    # HUC-specific seeds make the test and full runs use identical sections.
    set.seed(as.integer(huc_i) + seed)
    selected_blocks <- sort(unique(data$block))
    if (length(selected_blocks) > max_blocks) {
      selected_blocks <- sample(selected_blocks, max_blocks, replace = FALSE)
    }
    result <- data %>% group_by(huc, month, metabolism_metric) %>%
      group_modify(~ {
        answer <- summarize_fit(.x, selected_blocks)
        bar$tick()
        answer
      }) %>% ungroup()
    write_csv(result, path_output, append = !first_write,
              col_names = first_write)
    first_write <<- FALSE
  })
  message("Saved: ", path_output)

  # Optional test-mode figure (disabled for the full CONUS run).
  # if (test_mode) {
  #   result <- read_csv(path_output, col_types = cols(huc = col_character()))
  #   print(result %>% count(status))
  #   print(ggplot(result, aes(factor(month), baseline_slope)) +
  #     geom_hline(yintercept = c(slope_bounds, 1), linetype = "dashed", color = "grey60") +
  #     geom_linerange(aes(ymin = slope_min, ymax = slope_max), na.rm = TRUE) +
  #     geom_point(aes(color = status)) +
  #     facet_grid(metabolism_metric ~ huc) + theme_bw() +
  #     labs(x = "Month", y = "Slope", color = "Deletion result",
  #          subtitle = "Points: full fits; ranges: valid section deletions (not confidence intervals)"))
  # }
  invisible(NULL)
})
