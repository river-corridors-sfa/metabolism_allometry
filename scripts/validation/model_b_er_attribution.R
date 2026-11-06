#!/usr/bin/env Rscript
# Three-version ER audit. Run from repository root; no production writes.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2)
out <- "reports/no_uparea_migration/model_b_er_attribution"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
paths <- list(
  original="data/raw/dryad_data/scaling_ER_pred",
  unfiltered="data/processed/17_maavara_refits/without_uparea/ER",
  model_b="reports/no_uparea_migration/er_qc_sensitivity/conus/near_zero")
read_pred <- function(version,mo) {
  p <- file.path(paths[[version]],paste0("scaling_ER_pred_",mo,".csv"))
  if (!file.exists(p)) stop("Missing ",version," ER prediction: ",p)
  x <- fread(p)
  stopifnot(ncol(x)==2L,!anyDuplicated(x[[1]]))
  setnames(x,c("COMID","rate"))
  x[,COMID:=as.character(COMID)]
  x
}
sc <- fread("data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",select=c("COMID","lengthkm","merit_width"))
sc[,COMID:=as.character(COMID)]
sc[,stream_area_m2:=lengthkm*1000*merit_width]
stopifnot(!anyDuplicated(sc$COMID))
lk <- as.data.table(readRDS("data/processed/merit_data/05_upstream_lookup.rds"))
lk[,`:=`(downstream_comid=as.character(downstream_comid),upstream_comid=as.character(upstream_comid))]
stopifnot(!anyDuplicated(lk,by=c("downstream_comid","upstream_comid")))
catdb <- as.data.table(foreign::read.dbf("data/processed/merit_data/03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.dbf",as.is=TRUE))
catdb[,COMID:=as.character(comid)]
attr <- merge(catdb[,.(COMID,unitarea)],sc[,.(COMID,stream_area_m2)],by="COMID",all=TRUE)
attr[,good_geometry:=is.finite(unitarea)&unitarea>0&is.finite(stream_area_m2)&stream_area_m2>0]
bad_geometry <- attr[good_geometry==FALSE,COMID]
rm(catdb,sc,attr); gc()
flags <- vector("list",12)
rate_counts <- vector("list",12)
for (i in seq_along(months)) {
  mo <- months[i]; cat("Complete-network audit",mo,"\n"); flush.console()
  r <- lapply(names(paths),read_pred,mo=mo)
  names(r) <- names(paths)
  common_ids <- Reduce(intersect,lapply(r,`[[`,"COMID"))
  good <- Reduce(intersect,lapply(r,function(x) x[is.finite(rate)&rate>0,COMID]))
  good <- setdiff(good,bad_geometry)
  bad_ids <- setdiff(unique(lk$upstream_comid),good)
  bad_focal <- unique(lk[upstream_comid %chin% bad_ids,downstream_comid])
  good_focal <- setdiff(unique(lk$downstream_comid),bad_focal)
  flags[[i]] <- data.table(downstream_comid=good_focal,month=i)
  rate_counts[[i]] <- rbindlist(lapply(names(r),function(v) data.table(month=i,version=v,n=nrow(r[[v]]),n_missing=sum(!is.finite(r[[v]]$rate)),n_nonpositive=sum(is.finite(r[[v]]$rate)&r[[v]]$rate<=0),n_positive=sum(is.finite(r[[v]]$rate)&r[[v]]$rate>0),n_common_ids=length(common_ids),n_complete_focal_networks=length(good_focal))))
}
complete <- rbindlist(flags)
fwrite(rbindlist(rate_counts),file.path(out,"local_rate_coverage.csv"))
fwrite(complete,file.path(out,"complete_network_focals.csv"))
rm(lk,flags,rate_counts); gc()
screened <- list(
  original=as.data.table(readRDS("data/processed/07_huc08_allometry_input.rds")),
  unfiltered=as.data.table(readRDS("data/processed/07_no_uparea_huc08_allometry_input.rds")),
  model_b=as.data.table(readRDS("reports/no_uparea_migration/model_b_national_assessment/model_b_er_screened.rds")))
for (v in names(screened)) {
  z <- screened[[v]][metabolism_metric=="ER"]
  z[,`:=`(huc=sprintf("%08d",as.integer(huc)),downstream_comid=as.character(downstream_comid))]
  screened[[v]] <- z[,.(huc,downstream_comid,month,cumulative_area_km2,cumulative_flux_g_o2_d,pct_area_missing,pct_metabolism_missing)]
  stopifnot(!anyDuplicated(screened[[v]],by=c("downstream_comid","month")))
}
base <- Reduce(function(a,b) merge(a,b,by=c("huc","downstream_comid","month"),all=FALSE,suffixes=c("","_other")),
  lapply(names(screened),function(v) {z <- copy(screened[[v]]);setnames(z,setdiff(names(z),c("huc","downstream_comid","month")),paste0(setdiff(names(z),c("huc","downstream_comid","month")),"_",v));z}))
base <- merge(base,complete,by=c("downstream_comid","month"),all=FALSE)
stopifnot(max(abs(base$cumulative_area_km2_original-base$cumulative_area_km2_unfiltered))<1e-5,
          max(abs(base$cumulative_area_km2_original-base$cumulative_area_km2_model_b))<1e-5,
          max(base$pct_metabolism_missing_original,base$pct_metabolism_missing_unfiltered,base$pct_metabolism_missing_model_b)<1e-8)
fwrite(base[,.(huc,downstream_comid,month,area=cumulative_area_km2_original,original_flux=cumulative_flux_g_o2_d_original,unfiltered_flux=cumulative_flux_g_o2_d_unfiltered,model_b_flux=cumulative_flux_g_o2_d_model_b)],file.path(out,"fixed_network_reaches.csv"))
fit <- function(a,y) {
  n <- length(a); nu <- uniqueN(a); span <- log10(max(a)/min(a))
  if (n<20L||nu<20L||span<1) return(list(n_reaches=n,n_unique_areas=nu,area_span=span,slope=NA_real_,intercept=NA_real_,r2=NA_real_))
  x <- log10(a); z <- log10(y); xx <- x-mean(x); zz <- z-mean(z); slope <- sum(xx*zz)/sum(xx^2)
  list(n_reaches=n,n_unique_areas=nu,area_span=span,slope=slope,intercept=mean(z)-slope*mean(x),r2=sum(xx*zz)^2/sum(xx^2)/sum(zz^2))
}
fixed <- rbindlist(lapply(names(paths),function(v) {
  z <- base[,fit(cumulative_area_km2_original,get(paste0("cumulative_flux_g_o2_d_",v))),by=.(huc,month)]
  z[,version:=v]
  z
}))
fwrite(fixed,file.path(out,"fixed_network_fits.csv"))
own <- rbindlist(lapply(names(paths),function(v) {
  z <- screened[[v]][,fit(cumulative_area_km2,cumulative_flux_g_o2_d),by=.(huc,month)]
  z[,version:=v]
  z
}))
fwrite(own,file.path(out,"own_coverage_fits.csv"))
saved <- list(original="data/processed/07_huc08_allometry_results.csv",unfiltered="data/processed/07_no_uparea_huc08_allometry_results.csv",model_b="reports/no_uparea_migration/model_b_national_assessment/four_design_fits.csv")
replay <- rbindlist(lapply(names(saved),function(v) {
  z <- fread(saved[[v]],colClasses=list(character="huc"))[metabolism_metric=="ER"]
  if (v=="model_b") z <- z[design=="modified_all" & is.finite(slope)]
  y <- merge(own[version==v & is.finite(slope)],z,by=c("huc","month"),suffixes=c("_replay","_saved"))
  data.table(version=v,eligible=nrow(y),saved=nrow(z),max_slope_error=max(abs(y$slope_replay-y$slope_saved)))
}))
stopifnot(all(replay$eligible==replay$saved),all(replay$max_slope_error<1e-8))
fwrite(replay,file.path(out,"fit_replay_validation.csv"))
cat("Three-version fixed-network and own-coverage fits complete\n")
