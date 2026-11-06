rm(list = ls())
source("scripts/00_setup.R")

# Complete-year, matched local reach predictions. ER is a positive magnitude.
# Stream order comes from the existing no-uparea reach key. Each order contributes
# one mean annual areal rate, so this is an order-level rate index, not river flux.
rates <- readRDS("data/processed/03_filtered_no_uparea_monthly_metabolism.rds") %>%
  transmute(comid = as.character(comid), month = as.integer(month),
            gpp = gpp_g_o2_m2_d, er = er_g_o2_m2_d)
key <- readRDS("data/processed/04c_no_uparea_comid_huc08_key.rds") %>%
  transmute(comid = as.character(comid), ecoregion = l1_ecoregion,
            stream_order = as.integer(stream_order))
stopifnot(!anyDuplicated(rates[c("comid", "month")]),
          !anyDuplicated(key$comid), all(rates$month %in% 1:12))

monthly <- rates %>% left_join(key, by = "comid") %>%
  mutate(days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[month],
         valid = is.finite(gpp) & gpp >= 0 & is.finite(er) & er >= 0 &
           !is.na(ecoregion) & !is.na(stream_order) & stream_order >= 1)
complete_ids <- monthly %>% group_by(comid) %>%
  summarise(complete = n_distinct(month) == 12 && all(valid), .groups = "drop") %>%
  filter(complete) %>% pull(comid)
if (!length(complete_ids)) stop("No complete matched reaches")
message("Complete matched reaches: ", length(complete_ids), " / ",
        n_distinct(monthly$comid))

annual <- monthly %>% filter(comid %in% complete_ids) %>%
  group_by(comid, ecoregion, stream_order) %>%
  summarise(gpp_rate = weighted.mean(gpp, days),
            er_rate = weighted.mean(er, days), .groups = "drop")
stopifnot(nrow(annual) == length(complete_ids))

# Each stream order has equal weight in the cumulative index. This preserves
# the previous rate-based comparison without letting abundant low-order reaches
# dominate solely because there are more of them.
by_order <- annual %>% group_by(ecoregion, stream_order) %>%
  summarise(mean_gpp = mean(gpp_rate), mean_er = mean(er_rate),
            n_reaches = n(), .groups = "drop") %>%
  group_by(ecoregion) %>% arrange(stream_order, .by_group = TRUE) %>%
  mutate(total_er = sum(mean_er),
         gpp_y = 100 * cumsum(mean_gpp) / total_er,
         er_y = 100 * cumsum(mean_er) / total_er) %>% ungroup()
stopifnot(all(is.finite(by_order$total_er) & by_order$total_er > 0),
          all(is.finite(by_order$gpp_y)), all(is.finite(by_order$er_y)))
ends <- by_order %>% group_by(ecoregion) %>% slice_tail(n = 1) %>% ungroup()
stopifnot(all(abs(ends$er_y - 100) < 1e-8))
message("Final cumulative rates as % of summed order-mean ER:")
print(ends %>% select(ecoregion, stream_order, gpp_y, er_y), n = Inf)

points <- bind_rows(
  by_order %>% transmute(ecoregion, x = stream_order, gpp_y, er_y),
  by_order %>% distinct(ecoregion) %>% mutate(x = 0L, gpp_y = 0, er_y = 0)) %>%
  arrange(ecoregion, x)

# Insert any crossings so the light gap fill changes process at the crossing.
make_ribbons <- function(part) {
  gap <- part$gpp_y - part$er_y
  crossing <- which(head(gap, -1) * tail(gap, -1) < 0)
  if (length(crossing)) {
    fraction <- -gap[crossing] / (gap[crossing + 1] - gap[crossing])
    part <- bind_rows(part, tibble(
      ecoregion = part$ecoregion[1],
      x = part$x[crossing] + fraction * (part$x[crossing + 1] - part$x[crossing]),
      gpp_y = part$gpp_y[crossing] + fraction *
        (part$gpp_y[crossing + 1] - part$gpp_y[crossing]),
      er_y = part$er_y[crossing] + fraction *
        (part$er_y[crossing + 1] - part$er_y[crossing]))) %>% arrange(x)
  }
  gap <- part$gpp_y - part$er_y
  higher <- if_else(head(gap, -1) + tail(gap, -1) >= 0, "GPP", "ER")
  runs <- cumsum(c(1L, as.integer(diff(match(higher, c("GPP", "ER"))) != 0)))
  map_dfr(unique(runs), function(id) {
    interval <- which(runs == id)
    part[min(interval):(max(interval) + 1), ] %>%
      mutate(higher = higher[min(interval)],
             ribbon = paste(ecoregion, id, sep = "_"))
  })
}
ribbons <- points %>% group_split(ecoregion) %>% map_dfr(make_ribbons)
lines <- points %>%
  pivot_longer(c(gpp_y, er_y), names_to = "process", values_to = "y") %>%
  mutate(process = recode(process, gpp_y = "GPP", er_y = "ER"))
upper <- max(105, ceiling(max(by_order$gpp_y) / 10) * 10)
figure <- ggplot() +
  geom_ribbon(data = ribbons,
              aes(x = x, ymin = pmin(gpp_y, er_y), ymax = pmax(gpp_y, er_y),
                  fill = higher, group = ribbon), alpha = 0.16) +
  geom_line(data = lines, aes(x, y, colour = process, group = process),
            linewidth = 0.85) +
  geom_point(data = lines %>% filter(x > 0),
             aes(x, y, colour = process), size = 1.25) +
  facet_wrap(vars(ecoregion), ncol = 3) +
  scale_colour_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_fill_manual(values = c(GPP = "#009E73", ER = "#0072B2")) +
  scale_x_continuous(breaks = 1:8, expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = c(0, 50, 100), expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(0, 8), ylim = c(0, upper)) +
  labs(x = "Stream order", y = "Cumulative annual areal rate (% of summed order-mean ER)",
       colour = "Process", fill = "Higher cumulative total",
       caption = stringr::str_wrap(paste0(
         "CONUS-filtered MERIT reaches with complete 12-month Maavara no-uparea GPP and ER predictions (n = ",
         format(nrow(annual), big.mark = ","), "). Within each Level I ecoregion, each stream order contributes its mean annual reach rate; ",
         "both curves are divided by the sum of ER order means. Orders are equally weighted regardless of reach count or river area. ",
         "ER ends at 100%; these curves summarize rates, not total river metabolism."), 120)) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "#F2F4F5", colour = NA),
        strip.text = element_text(colour = "#253746", face = "bold", size = 9),
        panel.spacing = unit(0.9, "lines"), legend.position = "top",
        legend.box = "horizontal",
        plot.caption = element_text(hjust = 0, colour = "grey35", size = 9),
        plot.margin = margin(10, 10, 10, 10))
ggsave("figures/48_cumulative_gpp_er_by_stream_order_ecoregion.png", figure,
       width = 12, height = 8.5, dpi = 300, bg = "white")
