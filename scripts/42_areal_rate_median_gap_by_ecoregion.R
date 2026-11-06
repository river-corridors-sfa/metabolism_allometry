rm(list = ls())
source("scripts/00_setup.R")

# Two heaps give the exact median of every prefix in O(n log n) time.
Rcpp::cppFunction('
#include <queue>
#include <vector>
#include <functional>
Rcpp::NumericVector running_median_cpp(Rcpp::NumericVector x) {
  std::priority_queue<double> lower;
  std::priority_queue<double, std::vector<double>, std::greater<double>> upper;
  Rcpp::NumericVector out(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    double value = x[i];
    if (lower.empty() || value <= lower.top()) lower.push(value);
    else upper.push(value);
    if (lower.size() > upper.size() + 1) {
      upper.push(lower.top()); lower.pop();
    } else if (upper.size() > lower.size()) {
      lower.push(upper.top()); upper.pop();
    }
    out[i] = lower.size() == upper.size()
      ? (lower.top() + upper.top()) / 2.0 : lower.top();
  }
  return out;
}')
stopifnot(isTRUE(all.equal(running_median_cpp(c(5, 1, 9, 2)),
                           c(5, 3, 5, 3.5))))

# CONUS-filtered Maavara no-uparea predictions; compare annual mean areal
# rates directly, without multiplying rates by river surface area.
rates <- readRDS("data/processed/03_filtered_no_uparea_monthly_metabolism.rds") %>%
  transmute(comid = as.character(comid), month = as.integer(month),
            gpp = gpp_g_o2_m2_d, er = er_g_o2_m2_d)
stopifnot(!anyNA(rates[c("comid", "month")]),
          all(rates$month %in% 1:12), !anyDuplicated(rates[c("comid", "month")]))
attrs <- readr::read_csv(
  "data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",
  col_select = c(COMID, lengthkm, merit_width, uparea),
  col_types = cols(COMID = col_character(), .default = col_double()),
  show_col_types = FALSE) %>%
  transmute(comid = COMID, drainage_km2 = uparea,
            river_area_m2 = lengthkm * 1000 * merit_width)
regions <- readRDS("data/processed/04c_no_uparea_comid_huc08_key.rds") %>%
  transmute(comid = as.character(comid), ecoregion = l1_ecoregion)
stopifnot(!anyDuplicated(attrs$comid), !anyDuplicated(regions$comid))

monthly <- rates %>% left_join(attrs, by = "comid") %>%
  left_join(regions, by = "comid") %>%
  mutate(days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[month],
         valid = is.finite(gpp) & gpp >= 0 & is.finite(er) & er >= 0 &
           is.finite(drainage_km2) & drainage_km2 > 0 &
           is.finite(river_area_m2) & river_area_m2 > 0 & !is.na(ecoregion))
complete_ids <- monthly %>% group_by(comid) %>%
  summarise(complete = n_distinct(month) == 12 && all(valid), .groups = "drop") %>%
  filter(complete) %>% pull(comid)
message("Complete GPP and ER reaches: ", length(complete_ids), " / ",
        n_distinct(monthly$comid), "; excluded: ",
        n_distinct(monthly$comid) - length(complete_ids))
if (!length(complete_ids)) stop("No complete reaches")

annual <- monthly %>% filter(comid %in% complete_ids) %>%
  group_by(comid, ecoregion, drainage_km2) %>%
  summarise(river_area_m2 = sum(river_area_m2 * days) / sum(days),
            GPP = sum(gpp * days) / sum(days),
            ER = sum(er * days) / sum(days), .groups = "drop")
stopifnot(nrow(annual) == length(complete_ids))

# X retains cumulative river surface area (km2). Y is the exact running
# unweighted median annual mean rate across reaches seen so far, in
# g O2 m^-2 d^-1. No rate is multiplied by river area.
ordered <- annual %>% group_by(ecoregion) %>%
  arrange(drainage_km2, comid, .by_group = TRUE) %>%
  mutate(x = cumsum(river_area_m2) / 1e6,
         gpp_y = running_median_cpp(GPP),
         er_y = running_median_cpp(ER)) %>%
  ungroup()
stopifnot(all(is.finite(ordered$x)), all(is.finite(ordered$gpp_y)),
          all(is.finite(ordered$er_y)), all(ordered$x > 0))
points <- ordered %>% select(ecoregion, x, gpp_y, er_y) %>%
  arrange(ecoregion, x)

# Insert exact linear crossings before building fill runs. This prevents a
# colored ribbon from spilling across a point where GPP and ER switch rank.
make_bands <- function(part) {
  gap <- part$gpp_y - part$er_y
  crossing <- which(head(gap, -1) * tail(gap, -1) < 0)
  if (length(crossing)) {
    fraction <- -gap[crossing] / (gap[crossing + 1] - gap[crossing])
    extra <- tibble(
      ecoregion = part$ecoregion[1],
      x = part$x[crossing] + fraction * (part$x[crossing + 1] - part$x[crossing]),
      gpp_y = part$gpp_y[crossing] + fraction *
        (part$gpp_y[crossing + 1] - part$gpp_y[crossing]),
      er_y = part$er_y[crossing] + fraction *
        (part$er_y[crossing + 1] - part$er_y[crossing]))
    part <- bind_rows(part, extra) %>% arrange(x)
  }
  gap <- part$gpp_y - part$er_y
  dominant <- if_else(head(gap, -1) + tail(gap, -1) >= 0, "GPP", "ER")
  run <- cumsum(c(1L, as.integer(diff(match(dominant, c("GPP", "ER"))) != 0)))
  map_dfr(unique(run), function(id) {
    interval <- which(run == id)
    part[min(interval):(max(interval) + 1), ] %>%
      mutate(higher = dominant[min(interval)],
             band = paste(ecoregion, id, sep = "_"))
  })
}
bands <- points %>% group_split(ecoregion) %>% map_dfr(make_bands)
stopifnot(all(is.finite(bands$x)), all(bands$higher %in% c("GPP", "ER")))

line_data <- points %>%
  pivot_longer(c(gpp_y, er_y), names_to = "process", values_to = "y") %>%
  mutate(process = recode(process, gpp_y = "GPP", er_y = "ER"))
figure <- ggplot() +
  geom_ribbon(data = bands,
              aes(x = x, ymin = pmin(gpp_y, er_y), ymax = pmax(gpp_y, er_y),
                  fill = higher, group = band), alpha = 0.18) +
  geom_line(data = line_data, aes(x, y, colour = process, group = process),
            linewidth = 0.85) +
  facet_wrap(vars(ecoregion), ncol = 3, scales = "free") +
  scale_colour_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_fill_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_x_continuous(labels = scales::label_number(big.mark = ","),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(labels = scales::label_number(big.mark = ","),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Cumulative river surface area (km²)",
       y = "Running median annual areal rate (g O₂ m⁻² d⁻¹)",
       colour = "Process", fill = "Higher running median",
       caption = stringr::str_wrap(paste0(
         "CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea GPP and ER predictions (n = ",
         format(nrow(annual), big.mark = ","), "). Reaches ordered by upstream drainage area within each Level I ecoregion. ",
         "Y is the unweighted median of annual mean reach rates added so far; no rate is multiplied by river area. Panel axes vary."), 120)) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F2F4F5", colour = NA),
        strip.text = element_text(colour = "#253746", face = "bold", size = 9),
        panel.spacing = unit(0.9, "lines"),
        legend.position = "top", legend.box = "horizontal",
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(10, 10, 10, 10))
ggsave("figures/42_areal_rate_median_gap_by_ecoregion.png", figure,
       width = 12, height = 8.5, dpi = 300, bg = "white")
