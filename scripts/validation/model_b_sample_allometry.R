#!/usr/bin/env Rscript
# Recompute routed cumulative metabolism only for six selected HUC8s.
suppressPackageStartupMessages(library(data.table))
out <- "reports/no_uparea_migration/model_b_comparability"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# Stratified, deterministic selection: two HUC8s per longitude band at the
# 25th and 75th percentiles of maximum upstream area among modest-size HUC8s
# with all 24 original monthly process fits.
s <- fread("data/processed/04c_huc08_summary.csv",colClasses=list(character="huc8"))
f <- fread("data/processed/07_huc08_allometry_results.csv",colClasses=list(character="huc"))
count <- f[,.(fits=.N),by=huc]
s <- merge(s,count,by.x="huc8",by.y="huc")
s <- s[n_flowlines>=40 & n_flowlines<=150 & fits==24 & log10_uparea_range>=1.5]
s[, region:=fifelse(centroid_longitude< -105,"West",fifelse(centroid_longitude< -90,"Central","East"))]
chosen <- rbindlist(lapply(c("West","Central","East"),function(region_i) {
  z <- s[region==region_i][order(max_uparea_km2)]
  z[round(c(.25,.75)*(nrow(z)-1))+1]
}))
stopifnot(nrow(chosen)==6L,!anyDuplicated(chosen$huc8))
fwrite(chosen[,.(huc8,huc8_name,states,region,centroid_longitude,n_flowlines,max_uparea_km2)],file.path(out,"sample_huc8s.csv"))
cat("Sample HUC8s:",paste(chosen$huc8,collapse=", "),"\n")

cat_dbf <- as.data.table(foreign::read.dbf("data/processed/merit_data/03_cat_pfaf_7_MERIT_Hydro_v07_Basins_v01_CONUS_filtered.dbf",as.is=TRUE))
cat_dbf[,`:=`(COMID=as.character(comid),huc8=sprintf("%08s",as.character(huc)))]
cat_dbf[, huc8:=gsub(" ","0",huc8,fixed=TRUE)]
stopifnot(!anyDuplicated(cat_dbf$COMID))
focal <- cat_dbf[huc8 %chin% chosen$huc8,COMID]
stopifnot(length(focal)>0)
rm(s,f,count)

cat("Reading saved upstream lookup\n")
lk <- as.data.table(readRDS("data/processed/merit_data/05_upstream_lookup.rds"))
lk[,`:=`(downstream_comid=as.character(downstream_comid),upstream_comid=as.character(upstream_comid))]
lk <- lk[downstream_comid %chin% focal]
stopifnot(nrow(lk)>0,!anyDuplicated(lk,by=c("downstream_comid","upstream_comid")))
up_ids <- unique(lk$upstream_comid)
attr <- cat_dbf[,.(COMID,huc8,unitarea)]
setnames(attr,"COMID","upstream_comid")
lk <- merge(lk,attr[,.(upstream_comid,unitarea)],by="upstream_comid",all.x=TRUE,sort=FALSE)
stopifnot(!anyNA(lk$unitarea))
cum_area <- lk[,.(cumulative_area_km2=sum(unitarea),n_upstream_reaches=uniqueN(upstream_comid)),by=downstream_comid]
rm(cat_dbf)

scaling <- fread("data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv",select=c("COMID","lengthkm","merit_width","uparea"))
scaling[,COMID:=as.character(COMID)]
attr <- merge(attr,scaling,by.x="upstream_comid",by.y="COMID",all.x=TRUE,sort=FALSE)
attr[,stream_area_m2:=lengthkm*1000*merit_width]
stopifnot(!anyNA(attr[upstream_comid %chin% up_ids]$stream_area_m2))
focal_attr <- attr[upstream_comid %chin% focal,.(downstream_comid=upstream_comid,huc=huc8,uparea)]
cum_area <- merge(cum_area,focal_attr,by="downstream_comid",all.x=TRUE,sort=FALSE)
cum_area[,pct_area_missing:=pmax(0,1-cumulative_area_km2/uparea)*100]
stopifnot(!anyNA(cum_area$huc))
rm(scaling)

read_rate <- function(path,metric) {
  z <- fread(path)
  setnames(z,names(z)[2],metric)
  z[,COMID:=as.character(COMID)]
  z[COMID %chin% up_ids]
}

b_all <- vector("list",12)
for (i in seq_along(months)) {
  mo <- months[i]
  cat("Sample cumulative",mo,"\n")
  g <- read_rate(file.path("data/processed/17_maavara_refits/without_uparea/GPP",paste0("scaling_GPP_pred_",mo,".csv")),"gpp")
  e <- read_rate(file.path("reports/no_uparea_migration/er_qc_sensitivity/conus/near_zero",paste0("scaling_ER_pred_",mo,".csv")),"er")
  rates <- merge(g,e,by="COMID",all=TRUE,sort=FALSE)
  setnames(rates,"COMID","upstream_comid")
  rates <- merge(rates,attr[,.(upstream_comid,stream_area_m2)],by="upstream_comid",all.x=TRUE,sort=FALSE)
  rates[gpp<0,gpp:=NA_real_]; rates[er<0,er:=NA_real_]
  rates[,`:=`(gpp_flux=gpp*stream_area_m2,er_flux=er*stream_area_m2)]
  z <- merge(lk,rates[,.(upstream_comid,gpp_flux,er_flux)],by="upstream_comid",all.x=TRUE,sort=FALSE)
  agg <- z[,.(gpp_flux=if (all(is.na(gpp_flux))) NA_real_ else sum(gpp_flux,na.rm=TRUE),
               er_flux=if (all(is.na(er_flux))) NA_real_ else sum(er_flux,na.rm=TRUE),
               gpp_covered_area=sum(unitarea[!is.na(gpp_flux)]),
               er_covered_area=sum(unitarea[!is.na(er_flux)])),by=downstream_comid]
  agg <- merge(agg,cum_area,by="downstream_comid",all.x=TRUE,sort=FALSE)
  agg[,month:=i]
  b_all[[i]] <- rbind(
    agg[,.(huc,downstream_comid,month,metabolism_metric="GPP",cumulative_area_km2,cumulative_flux_g_o2_d=gpp_flux,pct_area_missing,pct_metabolism_missing=(1-gpp_covered_area/cumulative_area_km2)*100)],
    agg[,.(huc,downstream_comid,month,metabolism_metric="ER",cumulative_area_km2,cumulative_flux_g_o2_d=er_flux,pct_area_missing,pct_metabolism_missing=(1-er_covered_area/cumulative_area_km2)*100)])
}
b <- rbindlist(b_all)
b <- b[is.finite(cumulative_area_km2)&cumulative_area_km2>0 &
         is.finite(cumulative_flux_g_o2_d)&cumulative_flux_g_o2_d>0 &
         is.finite(pct_area_missing)&pct_area_missing<10 &
         is.finite(pct_metabolism_missing)&pct_metabolism_missing<10]
stopifnot(!anyDuplicated(b,by=c("huc","downstream_comid","month","metabolism_metric")))

cat("Reading saved original and six-feature 07 inputs for validation\n")
original <- readRDS("data/processed/07_huc08_allometry_input.rds")
original <- as.data.table(original[original$huc %in% chosen$huc8,])
six <- readRDS("data/processed/07_no_uparea_huc08_allometry_input.rds")
six <- as.data.table(six[six$huc %in% chosen$huc8 & six$metabolism_metric=="GPP",])
original[,downstream_comid:=as.character(downstream_comid)]
six[,downstream_comid:=as.character(downstream_comid)]
common_gpp <- merge(b[metabolism_metric=="GPP"],six,by=c("huc","downstream_comid","month","metabolism_metric"),suffixes=c("_recomputed","_saved"))
validation <- data.table(recomputed_gpp_rows=nrow(b[metabolism_metric=="GPP"]),saved_gpp_rows=nrow(six),matched=nrow(common_gpp),
                         max_area_diff=max(abs(common_gpp$cumulative_area_km2_recomputed-common_gpp$cumulative_area_km2_saved)),
                         max_flux_relative_diff=max(abs(common_gpp$cumulative_flux_g_o2_d_recomputed-common_gpp$cumulative_flux_g_o2_d_saved)/common_gpp$cumulative_flux_g_o2_d_saved))
fwrite(validation,file.path(out,"sample_gpp_replay_validation.csv"))
if (validation$recomputed_gpp_rows!=validation$saved_gpp_rows || validation$matched!=validation$saved_gpp_rows || validation$max_area_diff>1e-5 || validation$max_flux_relative_diff>1e-8) stop("Sample GPP replay does not match saved 07 input")
rm(six,common_gpp)

key_cols <- c("huc","downstream_comid","month","metabolism_metric")
common <- merge(original,b,by=key_cols,suffixes=c("_original","_B"))
stopifnot(max(abs(common$cumulative_area_km2_original-common$cumulative_area_km2_B))<1e-5)
paired <- common[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2=cumulative_area_km2_original,
                   original_flux=cumulative_flux_g_o2_d_original,B_flux=cumulative_flux_g_o2_d_B,
                   ratio_B_original=cumulative_flux_g_o2_d_B/cumulative_flux_g_o2_d_original)]
fwrite(paired,file.path(out,"sample_cumulative_paired.csv"))

fits_input <- rbindlist(list(
  original[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2,cumulative_flux_g_o2_d,design="original_all")],
  b[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2,cumulative_flux_g_o2_d,design="B_all")],
  common[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2=cumulative_area_km2_original,cumulative_flux_g_o2_d=cumulative_flux_g_o2_d_original,design="original_common")],
  common[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2=cumulative_area_km2_B,cumulative_flux_g_o2_d=cumulative_flux_g_o2_d_B,design="B_common")]))
fwrite(fits_input,file.path(out,"sample_allometry_fit_input.csv"))
fit_one <- function(z) {
  n <- nrow(z); narea <- uniqueN(z$cumulative_area_km2)
  span <- log10(max(z$cumulative_area_km2)/min(z$cumulative_area_km2))
  if (n<20 || narea<20 || span<1) return(list(n_reaches=n,n_areas=narea,log10_area_range=span,median_flux=median(z$cumulative_flux_g_o2_d),slope=NA_real_,r2=NA_real_,simple_class="insufficient"))
  model <- lm(log10(cumulative_flux_g_o2_d)~log10(cumulative_area_km2),data=z)
  slope <- unname(coef(model)[2]); r2 <- summary(model)$r.squared
  class <- if (r2<.8) "low_r2" else if (slope<.95) "sublinear" else if (slope<=1.05) "near_linear" else "superlinear"
  list(n_reaches=n,n_areas=narea,log10_area_range=span,median_flux=median(z$cumulative_flux_g_o2_d),slope=slope,r2=r2,simple_class=class)
}
fits <- fits_input[,fit_one(.SD),by=.(huc,month,metabolism_metric,design)]
fwrite(fits,file.path(out,"sample_allometry_fits.csv"))
wide <- dcast(fits,huc+month+metabolism_metric~design,value.var=c("n_reaches","median_flux","slope","r2","simple_class"))
wide[,`:=`(controlled_slope_delta=slope_B_common-slope_original_common,
           coverage_slope_delta=slope_B_all-slope_B_common,
           total_slope_delta=slope_B_all-slope_original_all,
           controlled_class_change=simple_class_B_common!=simple_class_original_common,
           expanded_class_change=simple_class_B_all!=simple_class_original_all,
           B_only_reaches=n_reaches_B_all-n_reaches_B_common,
           original_only_reaches=n_reaches_original_all-n_reaches_original_common)]
fwrite(wide,file.path(out,"sample_allometry_comparison.csv"))
cat("Saved sample allometry comparisons\n")
