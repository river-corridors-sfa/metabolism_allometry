#!/usr/bin/env Rscript
# Monthly scale-dependence diagnostic for GPP, using completed local WBD
# assignments for HUC2s 01--12. Each regional value is a median of child-unit
# slopes, so each HUC2 receives equal weight in the continental summaries.
suppressPackageStartupMessages(library(tidyverse))

min_reaches <- 5
min_log10_area_range <- 1
out <- "data/processed/huc2_hierarchy_gpp"
dir.create("figures", showWarnings = FALSE)

assignments <- readRDS(file.path(out, "huc2_reach_huc12_assignments.rds"))
huc2_ids <- sort(unique(assignments$huc2))

inputs <- readRDS("data/processed/07_huc08_allometry_input.rds") |>
  filter(metabolism_metric == "GPP", cumulative_area_km2 > 0,
         cumulative_flux_g_o2_d > 0, pct_area_missing < 10,
         pct_gpp_metabolism_missing < 10, str_sub(huc, 1, 2) %in% huc2_ids) |>
  transmute(comid = downstream_comid, month, huc2 = str_sub(huc, 1, 2),
            log_area = log10(cumulative_area_km2),
            log_flux = log10(cumulative_flux_g_o2_d))

slope_stats <- function(d, groups) {
  d |>
    group_by(across(all_of(groups))) |>
    summarise(n_reaches = n(), min_area = min(log_area), max_area = max(log_area),
              mean_area = mean(log_area), mean_flux = mean(log_flux),
              mean_area_flux = mean(log_area * log_flux),
              mean_area_sq = mean(log_area^2), mean_flux_sq = mean(log_flux^2),
              .groups = "drop") |>
    mutate(log10_area_range = max_area - min_area,
           ss_area = mean_area_sq - mean_area^2,
           ss_flux = mean_flux_sq - mean_flux^2,
           cov_area_flux = mean_area_flux - mean_area * mean_flux,
           slope = if_else(n_reaches >= min_reaches & log10_area_range >= min_log10_area_range & ss_area > 0,
                           cov_area_flux / ss_area, NA_real_),
           r_squared = if_else(!is.na(slope) & ss_flux > 0,
                               cov_area_flux^2 / (ss_area * ss_flux), NA_real_)) |>
    select(all_of(groups), n_reaches, log10_area_range, slope, r_squared)
}

parent_fits <- slope_stats(inputs, c("huc2", "month")) |>
  mutate(level = "HUC2", huc_id = huc2) |>
  select(huc2, month, level, huc_id, everything())

child_fits <- inputs |>
  inner_join(assignments, by = c("comid", "huc2"), relationship = "many-to-many") |>
  pivot_longer(c(huc6, huc8, huc10, huc12), names_to = "level", values_to = "huc_id") |>
  mutate(level = str_to_upper(level)) |>
  slope_stats(c("huc2", "month", "level", "huc_id"))

fits <- bind_rows(parent_fits, child_fits) |>
  mutate(level = factor(level, levels = c("HUC2", "HUC6", "HUC8", "HUC10", "HUC12"))) |>
  arrange(huc2, month, level, huc_id)

regional <- fits |>
  filter(!is.na(slope)) |>
  group_by(huc2, month, level) |>
  summarise(n_units = n(), median_slope = median(slope), .groups = "drop")

tiles <- regional |>
  group_by(month, level) |>
  summarise(n_huc2 = n(), mean_slope = mean(median_slope),
            sd_across_huc2 = sd(median_slope),
            iqr_across_huc2 = IQR(median_slope), .groups = "drop") |>
  mutate(month = factor(month, levels = 1:12, labels = month.abb),
         level = factor(level, levels = c("HUC2", "HUC6", "HUC8", "HUC10", "HUC12")))

saveRDS(fits, file.path(out, "huc2_month_level_gpp_fits.rds"))
write_csv(fits, file.path(out, "huc2_month_level_gpp_fits.csv"))
saveRDS(regional, file.path(out, "huc2_month_level_gpp_regional_medians.rds"))
write_csv(regional, file.path(out, "huc2_month_level_gpp_regional_medians.csv"))
saveRDS(tiles, file.path(out, "huc2_month_level_gpp_tiles.rds"))
write_csv(tiles, file.path(out, "huc2_month_level_gpp_tiles.csv"))

base <- theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        plot.title = element_text(face = "bold"))

p_mean <- ggplot(tiles, aes(level, month, fill = mean_slope)) +
  geom_tile(colour = "white", linewidth = .35) +
  geom_text(aes(label = sprintf("%.2f", mean_slope)), size = 3.1) +
  scale_fill_viridis_c(option = "C", name = "Mean\nexponent") +
  labs(title = "Mean regional GPP scaling exponent by month and nested HUC level",
       subtitle = "Unweighted mean of HUC2-specific median slopes; HUC2s 01--12") + base

p_sd <- ggplot(tiles, aes(level, month, fill = sd_across_huc2)) +
  geom_tile(colour = "white", linewidth = .35) +
  geom_text(aes(label = sprintf("%.2f", sd_across_huc2)), size = 3.1) +
  scale_fill_viridis_c(option = "magma", name = "SD across\nHUC2s") +
  labs(title = "Between-region variability in GPP scaling exponent",
       subtitle = "SD of HUC2-specific median slopes; each tile has 12 HUC2s") + base

ggsave("figures/23_huc2_month_level_gpp_mean_tiles.png", p_mean, width = 8.4, height = 6.4, dpi = 240)
ggsave("figures/23_huc2_month_level_gpp_sd_tiles.png", p_sd, width = 8.4, height = 6.4, dpi = 240)
message("Wrote monthly HUC-level GPP tile summaries for HUC2s ", paste(huc2_ids, collapse = ", "))
