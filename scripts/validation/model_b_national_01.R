#!/usr/bin/env Rscript
# National one-off Model B assessment. Run from repository root.
suppressPackageStartupMessages(library(data.table))
out <- "reports/no_uparea_migration/model_b_national_assessment"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
read_input <- function(path, metric=NULL) {
  x <- as.data.table(readRDS(path))
  x[, `:=`(huc=sprintf("%08d",as.integer(huc)),downstream_comid=as.character(downstream_comid))]
  if (!is.null(metric)) x <- x[metabolism_metric==metric]
  x[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2,cumulative_flux_g_o2_d,pct_area_missing,pct_metabolism_missing)]
}
original <- read_input("data/processed/07_huc08_allometry_input.rds")
modified_gpp <- read_input("data/processed/07_no_uparea_huc08_allometry_input.rds","GPP")
stopifnot(!anyDuplicated(original,by=c("downstream_comid","month","metabolism_metric")))
catdb <- as.data.table(foreign::read.dbf("data/processed/merit_data/03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.dbf",as.is=TRUE))
catdb[, COMID:=as.character(comid)]
stopifnot(!anyDuplicated(catdb$COMID))
focal <- catdb[,.(downstream_comid=COMID,huc=sprintf("%08d",as.integer(huc)))]
stopifnot(!anyDuplicated(focal$downstream_comid))
scaling <- fread("data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",select=c("COMID","lengthkm","merit_width","uparea"))
scaling[,COMID:=as.character(COMID)]
attr <- merge(catdb[,.(COMID,unitarea)],scaling,by="COMID",all.x=TRUE,sort=FALSE)
attr[,stream_area_m2:=lengthkm*1000*merit_width]
setnames(attr,"COMID","upstream_comid")
rm(catdb,scaling); gc()
lk <- as.data.table(readRDS("data/processed/merit_data/05_upstream_lookup.rds"))
lk[,`:=`(downstream_comid=as.character(downstream_comid),upstream_comid=as.character(upstream_comid))]
lk <- lk[downstream_comid %chin% focal$downstream_comid]
stopifnot(!anyDuplicated(lk,by=c("downstream_comid","upstream_comid")))
lk <- merge(lk,attr[,.(upstream_comid,unitarea)],by="upstream_comid",all.x=TRUE,sort=FALSE)
cat("Route pairs without catchment area:",sum(is.na(lk$unitarea)),"\n")
area <- lk[,.(cumulative_area_km2=sum(unitarea,na.rm=TRUE)),by=downstream_comid]
area <- merge(area,focal,by="downstream_comid",all.x=TRUE,sort=FALSE)
area <- merge(area,attr[,.(downstream_comid=upstream_comid,uparea)],by="downstream_comid",all.x=TRUE,sort=FALSE)
area[,pct_area_missing:=pmax(0,1-cumulative_area_km2/uparea)*100]
rm(focal); gc()
replay <- list()
for (m in seq_along(months)) {
  cat("Model B ER accumulation",months[m],"\n"); flush.console()
  p <- file.path("reports/no_uparea_migration/er_qc_sensitivity/conus/near_zero",paste0("scaling_ER_pred_",months[m],".csv"))
  if (!file.exists(p)) stop("Model B ER input missing: ",p)
  rates <- fread(p)
  setnames(rates,c(names(rates)[1:2]),c("upstream_comid","rate"))
  rates[,upstream_comid:=as.character(upstream_comid)]
  rates <- merge(rates,attr[,.(upstream_comid,stream_area_m2)],by="upstream_comid",all=FALSE,sort=FALSE)
  rates[rate<0,rate:=NA_real_]
  rates[,flux:=rate*stream_area_m2]
  z <- merge(lk,rates[,.(upstream_comid,flux)],by="upstream_comid",all.x=TRUE,sort=FALSE)
  agg <- z[,.(cumulative_flux_g_o2_d=if (all(is.na(flux))) NA_real_ else sum(flux,na.rm=TRUE),covered_area=sum(unitarea[!is.na(flux)],na.rm=TRUE)),by=downstream_comid]
  agg <- merge(agg,area,by="downstream_comid",all.x=TRUE,sort=FALSE)
  agg[,`:=`(month=m,metabolism_metric="ER",pct_metabolism_missing=(1-covered_area/cumulative_area_km2)*100)]
  agg <- agg[!is.na(huc)&is.finite(cumulative_area_km2)&cumulative_area_km2>0 & is.finite(cumulative_flux_g_o2_d)&cumulative_flux_g_o2_d>0 & is.finite(pct_area_missing)&pct_area_missing<10 & is.finite(pct_metabolism_missing)&pct_metabolism_missing<10,
             .(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2,cumulative_flux_g_o2_d,pct_area_missing,pct_metabolism_missing)]
  saveRDS(agg,file.path(out,sprintf("model_b_er_screened_%02d.rds",m)))
  replay[[m]] <- data.table(month=m,screened_er=nrow(agg))
  rm(rates,z,agg); gc(verbose=FALSE)
}
fwrite(rbindlist(replay),file.path(out,"accumulation_counts.csv"))
modified_er <- rbindlist(lapply(seq_along(months),function(m) readRDS(file.path(out,sprintf("model_b_er_screened_%02d.rds",m)))))
saveRDS(modified_er,file.path(out,"model_b_er_screened.rds"))
modified <- rbind(modified_gpp,modified_er)
stopifnot(!anyDuplicated(modified,by=c("downstream_comid","month","metabolism_metric")))
fwrite(modified,file.path(out,"modified_screened_reaches.csv"))
fwrite(original,file.path(out,"original_screened_reaches.csv"))
cat("Saved national screened data\n")
