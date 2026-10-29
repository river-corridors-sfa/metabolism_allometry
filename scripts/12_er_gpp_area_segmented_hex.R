# 12_er_gpp_area_segmented_hex.R
# Run from the project root after 07. Fit one log-log breakpoint separately
# for GPP and ER using the full 07 reach-month input (not the sampled RF
# predictions from 18), saving a two-panel viridis hex figure.

local({
  library(tidyverse)
  library(hexbin)
  library(segmented)

  path_input <- "data/processed/07_no_uparea_huc08_allometry_input.rds"
  path_figure <- "figures/figure_3_no_uparea_er_gpp_area_segmented_hex.png"
  path_summary <- "figures/12_no_uparea_er_gpp_area_segmented_fit.csv"
  
  hex_bins <- 65
  
  if (!file.exists(path_input)) stop("Missing input: ", path_input, ". Run 07 first.")
  dir.create("figures", recursive = TRUE, showWarnings = FALSE)

  input <- readRDS(path_input)
  required <- c("downstream_comid", "month", "metabolism_metric",
                "cumulative_area_km2", "cumulative_flux_g_o2_d",
                "pct_area_missing", "pct_metabolism_missing")
  if (!all(required %in% names(input))) {
    stop("07 input lacks required columns: ",
         paste(setdiff(required, names(input)), collapse = ", "))
  }
  if (any(!is.finite(input$pct_area_missing) |
          input$pct_area_missing >= 10) ||
      any(!is.finite(input$pct_metabolism_missing) |
          input$pct_metabolism_missing >= 10)) {
    stop("07 input violates the <10% area or metabolism coverage rules.")
  }
  plot_data <- input %>%
    transmute(comid = as.character(downstream_comid), month,
              metabolism = metabolism_metric,
              cumulative_area_km2, cumulative_flux_g_o2_d,
              log_area = log10(cumulative_area_km2),
              log_flux = log10(cumulative_flux_g_o2_d))
  if (any(!plot_data$metabolism %in% c("GPP", "ER")) ||
      any(!is.finite(plot_data$log_area)) ||
      any(!is.finite(plot_data$log_flux)) ||
      anyDuplicated(plot_data[c("comid", "month", "metabolism")])) {
    stop("07 input has invalid or duplicate reach-month-process rows.")
  }
  rm(input)
  message("Using ", nrow(plot_data), " reach-month-process observations from 07.")

  summaries <- list()
  predictions <- list()
  breakpoints <- list()
  for (metric in c("GPP", "ER")) {
    points <- plot_data %>% filter(metabolism == metric)
    if (nrow(points) < 30 || n_distinct(points$log_area) < 10) {
      stop("Too few observations for the ", metric, " segmented fit.")
    }
    linear_fit <- lm(log_flux ~ log_area, data = points)
    segmented_fit <- segmented(linear_fit, seg.Z = ~ log_area,
                               psi = list(log_area = median(points$log_area)))
    breakpoint_km2 <- 10^unname(segmented_fit$psi[1, "Est."])
    if (!is.finite(breakpoint_km2)) stop(metric, " fit returned an invalid breakpoint.")
    slopes <- slope(segmented_fit)$log_area[, "Est."]
    prediction <- tibble(log_area = seq(min(points$log_area), max(points$log_area),
                                        length.out = 400))
    prediction$cumulative_flux_g_o2_d <-
      10^as.numeric(predict(segmented_fit, newdata = prediction))
    prediction$cumulative_area_km2 <- 10^prediction$log_area
    prediction$metabolism <- metric
    predictions[[metric]] <- prediction
    breakpoints[[metric]] <- tibble(
      metabolism = metric,
      breakpoint_km2 = breakpoint_km2,
      label_y = 10^quantile(points$log_flux, 0.9),
      label = paste0(scales::comma(round(breakpoint_km2)), " km²")
    )

    summaries[[metric]] <- tibble(
      input_dataset = basename(path_input),
      metabolism = metric, n_reach_months = nrow(points),
      n_reaches = n_distinct(points$comid), breakpoint_km2 = breakpoint_km2,
      slope_below_breakpoint = slopes[1], slope_above_breakpoint = slopes[2],
      linear_aic = AIC(linear_fit), segmented_aic = AIC(segmented_fit)
    )

  }
  breakpoints <- bind_rows(breakpoints)
  figure <- ggplot(plot_data, aes(x = cumulative_area_km2,
                                  y = cumulative_flux_g_o2_d)) +
    geom_hex(bins = hex_bins) +
    geom_vline(data = breakpoints, aes(xintercept = breakpoint_km2),
               inherit.aes = FALSE, linetype = "dashed",
               colour = "grey75", linewidth = 0.8) +
    geom_line(data = bind_rows(predictions),
              aes(x = cumulative_area_km2, y = cumulative_flux_g_o2_d),
              inherit.aes = FALSE, colour = "white", linewidth = 1.5) +
    geom_line(data = bind_rows(predictions),
              aes(x = cumulative_area_km2, y = cumulative_flux_g_o2_d),
              inherit.aes = FALSE, colour = "#253746", linewidth = 0.8) +
    geom_text(data = breakpoints,
              aes(x = breakpoint_km2, y = label_y, label = label),
              inherit.aes = FALSE, hjust = -0.1, colour = "grey25", size = 3.5) +
    facet_wrap(~ metabolism, nrow = 1, scales = "free_y") +
    scale_x_log10(labels = scales::label_number_auto()) +
    scale_y_log10(labels = scales::label_number_auto()) +
    scale_fill_viridis_c(option = "D", trans = "log10",
                         name = "Reach-month count") +
    labs(x = expression("Cumulative upstream area (km"^2*")"),
         y = expression("Cumulative flux (g O"[2]*" d"^-1*")"),
         caption = "Full 07 input: negative local rates set to NA; <10% missing upstream area and metabolism coverage.") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
  ggsave(path_figure, figure, width = 12, height = 5.5, dpi = 300, bg = "white")
  message("Saved: ", path_figure)
  write_csv(bind_rows(summaries), path_summary)
  message("Saved: ", path_summary)

  # 4. Supplemental figures: one breakpoint per Level I ecoregion ------------
  path_key <- "data/processed/04c_no_uparea_comid_huc08_key.rds"
  if (!file.exists(path_key)) stop("Missing reach/ecoregion lookup: ", path_key)
  reach_regions <- readRDS(path_key) %>%
    dplyr::select(comid, l1_ecoregion) %>%
    distinct(comid, .keep_all = TRUE)
  for (metric in c("GPP", "ER")) {
  path_supplement <- file.path("figures", paste0(
    "figure_s2_no_uparea_", tolower(metric), "_ecoregion_segmented_hex.png"))
  region_data <- plot_data %>%
    filter(metabolism == metric) %>%
    inner_join(reach_regions, by = "comid", relationship = "many-to-one") %>%
    filter(!is.na(l1_ecoregion), nzchar(l1_ecoregion))
  if (nrow(region_data) == 0) stop("No ", metric, " observations have an ecoregion assignment.")

  region_predictions <- list()
  region_breakpoints <- list()
  fitted_regions <- character()
  for (region in sort(unique(region_data$l1_ecoregion))) {
    points <- region_data %>% filter(l1_ecoregion == region)
    if (nrow(points) < 30 || n_distinct(points$log_area) < 10) {
      warning("Skipping ", metric, " / ", region,
              ": too few observations for a segmented fit.")
      next
    }
    region_linear <- lm(log_flux ~ log_area, data = points)
    region_fit <- tryCatch(
      segmented(region_linear, seg.Z = ~ log_area,
                psi = list(log_area = median(points$log_area))),
      error = function(e) {
        warning("Skipping ", metric, " / ", region, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(region_fit)) next
    region_break <- 10^unname(region_fit$psi[1, "Est."])
    if (!is.finite(region_break)) {
      warning("Skipping ", metric, " / ", region, ": invalid breakpoint.")
      next
    }
    region_grid <- tibble(log_area = seq(min(points$log_area),
                                         max(points$log_area), length.out = 300))
    region_grid$cumulative_flux_g_o2_d <-
      10^as.numeric(predict(region_fit, newdata = region_grid))
    region_grid$cumulative_area_km2 <- 10^region_grid$log_area
    region_grid$l1_ecoregion <- region
    region_predictions[[region]] <- region_grid
    region_breakpoints[[region]] <- tibble(
      l1_ecoregion = region,
      breakpoint_km2 = region_break,
      label_y = 10^quantile(points$log_flux, 0.9),
      label = paste0(scales::comma(round(region_break)), " km²")
    )
    fitted_regions <- c(fitted_regions, region)
  }
  if (length(fitted_regions) == 0) stop("No ", metric, " ecoregion segmented fits succeeded.")
  region_data <- region_data %>% filter(l1_ecoregion %in% fitted_regions)

  supplement <- ggplot(region_data,
                       aes(x = cumulative_area_km2, y = cumulative_flux_g_o2_d)) +
    geom_hex(bins = hex_bins) +
    geom_vline(data = bind_rows(region_breakpoints),
               aes(xintercept = breakpoint_km2), inherit.aes = FALSE,
               linetype = "dashed", colour = "grey75", linewidth = 0.7) +
    geom_line(data = bind_rows(region_predictions),
              aes(x = cumulative_area_km2, y = cumulative_flux_g_o2_d),
              inherit.aes = FALSE, colour = "white", linewidth = 1.3) +
    geom_line(data = bind_rows(region_predictions),
              aes(x = cumulative_area_km2, y = cumulative_flux_g_o2_d),
              inherit.aes = FALSE, colour = "#253746", linewidth = 0.7) +
    geom_text(data = bind_rows(region_breakpoints),
              aes(x = breakpoint_km2, y = label_y, label = label),
              inherit.aes = FALSE, hjust = -0.1, colour = "grey25", size = 3) +
    facet_wrap(~ l1_ecoregion, scales = "free_y", ncol = 3) +
    scale_x_log10(labels = scales::label_number_auto()) +
    scale_y_log10(labels = scales::label_number_auto()) +
    scale_fill_viridis_c(option = "D", trans = "log10",
                         name = "Reach-month count") +
    labs(x = expression("Cumulative upstream area (km"^2*")"),
         y = paste0("Cumulative ", metric, " (g O2 d-1)"),
         caption = "Full 07 input: negative local rates set to NA; <10% missing upstream area and metabolism coverage.") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
  ggsave(path_supplement, supplement, width = 13,
         height = 3.6 * ceiling(length(fitted_regions) / 3),
         dpi = 300, bg = "white")
  message("Saved: ", path_supplement)
  }
})
