rm(list = ls())
source("scripts/00_setup.R")

# Local reach production, using the CONUS-filtered Maavara no-uparea refit
# used by the current main figures. Rates are g O2 m^-2 d^-1. ER is a
# positive respiration magnitude; negative predictions are invalid.
path_rates <- "data/processed/03_filtered_no_uparea_monthly_metabolism.rds"
path_attrs <- "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv"
path_figure <- "figures/33_where_river_metabolism_occurs.png"
stopifnot(file.exists(path_rates), file.exists(path_attrs))

rates <- readRDS(path_rates) %>%
  transmute(comid = as.character(comid), month = as.integer(month),
            gpp = gpp_g_o2_m2_d, er = er_g_o2_m2_d)
stopifnot(!anyNA(rates[c("comid", "month")]),
          all(rates$month %in% 1:12),
          !anyDuplicated(rates[c("comid", "month")]))

# lengthkm is km; merit_width is m; uparea is drainage area in km2.
# The product is LOCAL river water surface area in m2, never catchment area.
attrs <- readr::read_csv(path_attrs,
                         col_select = c(COMID, lengthkm, merit_width, uparea),
                         col_types = cols(COMID = col_character(),
                                          .default = col_double()),
                         show_col_types = FALSE) %>%
  transmute(comid = COMID, drainage_km2 = uparea,
            river_area_m2 = lengthkm * 1000 * merit_width)
stopifnot(!anyNA(attrs$comid), !anyDuplicated(attrs$comid))

monthly <- rates %>%
  left_join(attrs, by = "comid") %>%
  mutate(days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[month],
         valid = is.finite(gpp) & gpp >= 0 & is.finite(er) & er >= 0 &
           is.finite(drainage_km2) & drainage_km2 > 0 &
           is.finite(river_area_m2) & river_area_m2 > 0)
coverage <- monthly %>%
  group_by(comid) %>%
  summarise(n_months = n_distinct(month), n_valid = sum(valid), .groups = "drop")
complete_ids <- coverage %>% filter(n_months == 12, n_valid == 12) %>% pull(comid)
message("Full-year complete reaches: ", length(complete_ids), " / ", nrow(coverage),
        "; excluded incomplete reaches: ", nrow(coverage) - length(complete_ids),
        "; invalid GPP reach-months: ", sum(!is.finite(monthly$gpp) | monthly$gpp < 0),
        "; invalid ER reach-months: ", sum(!is.finite(monthly$er) | monthly$er < 0))
if (!length(complete_ids)) stop("No reaches have complete valid annual coverage")

annual <- monthly %>%
  filter(comid %in% complete_ids) %>%
  group_by(comid, drainage_km2) %>%
  summarise(river_area_m2 = sum(river_area_m2 * days) / sum(days),
            gpp_g_o2_y = sum(gpp * river_area_m2 * days),
            er_g_o2_y = sum(er * river_area_m2 * days), .groups = "drop") %>%
  arrange(drainage_km2, comid)
stopifnot(nrow(annual) == length(complete_ids),
          all(is.finite(annual$river_area_m2)),
          sum(annual$gpp_g_o2_y) > 0, sum(annual$er_g_o2_y) > 0)

curve <- annual %>%
  mutate(area_pct = 100 * cumsum(river_area_m2) / sum(river_area_m2),
         gpp_pct = 100 * cumsum(gpp_g_o2_y) / sum(gpp_g_o2_y),
         er_pct = 100 * cumsum(er_g_o2_y) / sum(er_g_o2_y))
stopifnot(all(abs(tail(curve[c("area_pct", "gpp_pct", "er_pct")], 1) - 100) < 1e-8))

# Split the reach crossing 80% proportionally so the highlighted interval
# represents exactly the final 20% of surface area, including its boundary.
boundary <- function(y) approx(c(0, curve$area_pct), c(0, y), xout = 80)$y
gpp_final <- 100 - boundary(curve$gpp_pct)
er_final <- 100 - boundary(curve$er_pct)
message(sprintf("Final 20%% of river surface area: GPP %.2f%%; ER %.2f%%",
                gpp_final, er_final))

plot_data <- bind_rows(
  tibble(area_pct = c(0, curve$area_pct), metabolism_pct = c(0, curve$gpp_pct),
         process = "GPP"),
  tibble(area_pct = c(0, curve$area_pct), metabolism_pct = c(0, curve$er_pct),
         process = "ER"))
labels <- tibble(process = c("GPP", "ER"),
                 metabolism_pct = c(64, 58),
                 text = sprintf(c("GPP: %.1f%%", "ER: %.1f%%"),
                                c(gpp_final, er_final)))

figure <- ggplot(plot_data, aes(area_pct, metabolism_pct, colour = process)) +
  annotate("rect", xmin = 80, xmax = 100, ymin = -Inf, ymax = Inf,
           fill = "grey94", alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60",
              linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 1.05) +
  geom_vline(xintercept = 80, colour = "grey65", linewidth = 0.4) +
  geom_label(data = labels,
             aes(x = 82, y = metabolism_pct, label = text),
             hjust = 0, fill = "white", linewidth = 0,
             size = 3.7, show.legend = FALSE) +
  scale_colour_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Cumulative river surface area (%)",
       y = "Cumulative annual metabolism (%)", colour = NULL,
       caption = stringr::str_wrap(paste0("CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea predictions (n = ",
                        format(nrow(annual), big.mark = ","), "). Reaches ordered by upstream drainage area; ",
                        "shading marks the final 20% of river surface area."), width = 100)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", legend.justification = "left",
        axis.line = element_line(colour = "#253746"),
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(12, 20, 12, 12))
ggsave(path_figure, figure, width = 8.2, height = 6.4, dpi = 300, bg = "white")
