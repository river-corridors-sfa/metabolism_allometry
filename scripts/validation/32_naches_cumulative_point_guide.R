# Companion to Step 2 of 32_naches_normalization_exploration.Rmd.
# Reads saved rates and routing; no modeling or routing is rerun.
suppressPackageStartupMessages({library(data.table); library(sf); library(ggplot2); library(patchwork)})
huc <- "17030002"
chosen <- data.table(comid=c("78018875","78018496","78018213"),
                     stage=c("Small","Medium","Large"), letter=c("S","M","L"))
fit <- as.data.table(readRDS("data/processed/07_huc08_allometry_input.rds"))
fit <- unique(fit[huc == "17030002" & month == 7L & metabolism_metric == "GPP",
                  .(comid=downstream_comid, area_km2=cumulative_area_km2)])
links <- as.data.table(readRDS("data/processed/merit_data/05_upstream_lookup.rds"))
links <- links[downstream_comid %in% fit$comid,
               .(comid=downstream_comid, source_comid=upstream_comid)]
rates <- as.data.table(readRDS("data/processed/01_assembled_maavara_monthly_metabolism.rds"))
rates <- rates[month == 7L & comid %in% fit$comid,
               .(source_comid=comid, rate=as.numeric(gpp_gO2_m2_d))]
rates[!is.finite(rate) | rate < 0, rate := NA_real_]
stopifnot(!anyDuplicated(rates$source_comid))
# Match the Rmd's Step 2: sum raw upstream areal rates, without area weighting.
v <- merge(links, rates, by="source_comid", all.x=TRUE)
v <- v[, .(gpp_sum=if(all(is.na(rate))) NA_real_ else sum(rate,na.rm=TRUE), n_upstream=.N), by=comid]
v <- merge(v, fit, by="comid", all.y=TRUE)
v <- v[is.finite(gpp_sum) & gpp_sum > 0 & is.finite(area_km2) & area_km2 > 0]
stopifnot(nrow(v)==51L)
v[, `:=`(x=log10(area_km2*1e6), y=log10(gpp_sum))]
marked <- merge(v, chosen, by="comid")
stopifnot(nrow(marked)==3L, setequal(marked$n_upstream, c(1L,3L,51L)))
read_huc <- function(prefix) {
  base <- paste0(prefix,"_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered")
  st_transform(st_read(file.path("data/processed/merit_data",paste0(base,".shp")),
                       query=sprintf('SELECT * FROM "%s" WHERE huc = \'%s\'',base,huc),
                       quiet=TRUE),5070)
}
cat <- read_huc("03_cat"); riv <- read_huc("03_riv")
stopifnot(nrow(cat)==51L,nrow(riv)==51L)
cat <- merge(cat,v[,.(comid,y)],by="comid",all.x=TRUE)
riv <- merge(riv,v[,.(comid,y)],by="comid",all.x=TRUE)
cols <- c("#F1C40F","#F39C12","#E74C3C","#A42A75","#40358B","#163B73")
color_scale <- scale_color_gradientn(colors=cols,limits=range(v$y),
                                    name=expression(log[10]*" cumulative GPP"))
map_panel <- function(i) {
  row <- chosen[i]; ids <- links[comid==row$comid,source_comid]
  point <- marked[comid==row$comid]
  ggplot() +
    geom_sf(data=cat,fill="#ECECEC",color="white",linewidth=.18) +
    geom_sf(data=cat[cat$comid %in% ids,],fill="#C8DCE8",color="white",linewidth=.2) +
    geom_sf(data=riv,color="#999999",linewidth=.3) +
    geom_sf(data=riv[riv$comid %in% ids,],aes(color=y),linewidth=1.15) +
    color_scale + coord_sf(datum=NA,expand=FALSE) +
    labs(title=paste0(row$letter,"  ",row$stage," upstream set"),
         subtitle=sprintf("%d reach%s | %.1f km² upstream\nSum of raw GPP rates = %.2f g O2/m²/day",
                          point$n_upstream,ifelse(point$n_upstream==1,"","es"),
                          point$area_km2,point$gpp_sum)) +
    theme_void(base_size=11) +
    theme(plot.title=element_text(face="bold",size=13),
          plot.subtitle=element_text(size=9.5,lineheight=1.1),
          legend.position="none",plot.margin=margin(5,5,5,5))
}
scatter <- ggplot(v,aes(x,y)) +
  geom_hline(yintercept=0,linetype="dashed",color="#A35C43") +
  geom_smooth(method="lm",formula=y~x,se=FALSE,color="#19354A",linewidth=.8) +
  geom_point(aes(color=y),size=3) +
  geom_point(data=marked,shape=21,fill=NA,color="black",stroke=1.2,size=6) +
  geom_text(data=marked,aes(label=letter),color="black",fontface="bold",size=3.5,nudge_y=.12) +
  color_scale +
  labs(x=expression(log[10]*" cumulative upstream watershed area (m"^2*")"),
       y="log10 cumulative GPP rate (g O2/m2/day)") +
  theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),legend.position="right")
figure <- ((map_panel(1)+map_panel(2)+map_panel(3))/scatter) +
  plot_layout(heights=c(1.05,1)) +
  plot_annotation(caption="Naches HUC8 17030002, July. Recreates Step 2 of 32_naches_normalization_exploration.Rmd from saved inputs.")
dir.create("figures",showWarnings=FALSE)
fwrite(marked[,.(comid,stage,n_upstream,area_km2,gpp_sum)],"figures/32_naches_cumulative_point_values.csv")
ggsave("figures/32_naches_cumulative_point_guide.png",figure,width=15,height=10,dpi=220,bg="white")
ggsave("figures/32_naches_cumulative_point_guide.pdf",figure,width=15,height=10,device=cairo_pdf,bg="white")
print(marked)
