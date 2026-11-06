#!/usr/bin/env Rscript
# Run after model_b_national_01.R; reuses saved original and modified GPP outputs.
suppressPackageStartupMessages(library(data.table))
out <- "reports/no_uparea_migration/model_b_national_assessment"
keys <- c("huc","downstream_comid","month","metabolism_metric")
o <- fread(file.path(out,"original_screened_reaches.csv"),colClasses=list(character=c("huc","downstream_comid")))
b <- fread(file.path(out,"modified_screened_reaches.csv"),colClasses=list(character=c("huc","downstream_comid")))
e <- b[metabolism_metric=="ER"]
saveRDS(e,file.path(out,"model_b_er_screened.rds"))
six <- as.data.table(readRDS("data/processed/07_no_uparea_huc08_allometry_input.rds"))
six <- six[six$metabolism_metric=="GPP",]
six[,`:=`(huc=sprintf("%08d",as.integer(huc)),downstream_comid=as.character(downstream_comid))]
check <- merge(b[metabolism_metric=="GPP"],six,by=keys,suffixes=c("_here","_saved"))
stopifnot(nrow(check)==nrow(six),nrow(check)==nrow(b[metabolism_metric=="GPP"]),
          max(abs(check$cumulative_flux_g_o2_d_here-check$cumulative_flux_g_o2_d_saved)/check$cumulative_flux_g_o2_d_saved)<1e-10)
sample <- fread("reports/no_uparea_migration/model_b_comparability/sample_cumulative_paired.csv",colClasses=list(character=c("huc","downstream_comid")))
sample <- sample[metabolism_metric=="ER"]
validation <- merge(e,sample,by=keys)
stopifnot(nrow(validation)==nrow(sample),
          max(abs(validation$cumulative_flux_g_o2_d-validation$B_flux)/sample$B_flux)<1e-8)
fwrite(data.table(sample_er_rows=nrow(sample),matched=nrow(validation),max_relative_flux_error=max(abs(validation$cumulative_flux_g_o2_d-validation$B_flux)/sample$B_flux)),file.path(out,"six_huc8_validation.csv"))
common <- merge(o,b,by=keys,suffixes=c("_original","_modified"))
stopifnot(max(abs(common$cumulative_area_km2_original-common$cumulative_area_km2_modified))<1e-5)
common[,`:=`(log_area=log10(cumulative_area_km2_original),log_flux_ratio=log10(cumulative_flux_g_o2_d_modified/cumulative_flux_g_o2_d_original))]
fwrite(common[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2=cumulative_area_km2_original,original_flux=cumulative_flux_g_o2_d_original,modified_flux=cumulative_flux_g_o2_d_modified,log_flux_ratio)],file.path(out,"common_reach_flux.csv"))
fit <- function(z) {
  n <- nrow(z); a <- z$cumulative_area_km2; y <- z$cumulative_flux_g_o2_d
  span <- log10(max(a)/min(a)); nu <- uniqueN(a)
  if (n<20 || nu<20 || span<1) return(list(n_reaches=n,n_unique_areas=nu,log10_area_range=span,median_flux=median(y),slope=NA_real_,intercept=NA_real_,r2=NA_real_,median_area=median(a)))
  x <- log10(a); y <- log10(y); xx <- x-mean(x); yy <- y-mean(y)
  slope <- sum(xx*yy)/sum(xx^2)
  list(n_reaches=n,n_unique_areas=nu,log10_area_range=span,median_flux=median(z$cumulative_flux_g_o2_d),slope=slope,intercept=mean(y)-slope*mean(x),r2=sum(xx*yy)^2/sum(xx^2)/sum(yy^2),median_area=median(a))
}
z <- rbindlist(list(o[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2,cumulative_flux_g_o2_d,design="original_all")],
                    b[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2,cumulative_flux_g_o2_d,design="modified_all")],
                    common[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2=cumulative_area_km2_original,cumulative_flux_g_o2_d=cumulative_flux_g_o2_d_original,design="original_common")],
                    common[,.(huc,downstream_comid,month,metabolism_metric,cumulative_area_km2=cumulative_area_km2_modified,cumulative_flux_g_o2_d=cumulative_flux_g_o2_d_modified,design="modified_common")]))
fits <- z[,fit(.SD),by=.(huc,month,metabolism_metric,design)]
fwrite(fits,file.path(out,"four_design_fits.csv"))
saved <- fread("data/processed/07_huc08_allometry_results.csv",colClasses=list(character="huc"))
replay <- merge(fits[design=="original_all" & is.finite(slope)],saved,by=c("huc","month","metabolism_metric"),suffixes=c("_replay","_saved"))
stopifnot(nrow(replay)==nrow(saved),max(abs(replay$slope_replay-replay$slope_saved))<1e-8)
fwrite(data.table(matched=nrow(replay),max_slope_error=max(abs(replay$slope_replay-replay$slope_saved))),file.path(out,"original_fit_validation.csv"))
wide <- dcast(fits,huc+month+metabolism_metric~design,value.var=c("n_reaches","median_flux","slope","intercept","r2","median_area"))
wide[,`:=`(total_slope_delta=slope_modified_all-slope_original_all,prediction_slope_delta=slope_modified_common-slope_original_common,
           coverage_slope_delta=slope_modified_all-slope_modified_common,original_coverage_slope_delta=slope_original_common-slope_original_all,
           total_intercept_delta=intercept_modified_all-intercept_original_all,prediction_intercept_delta=intercept_modified_common-intercept_original_common,
           total_r2_delta=r2_modified_all-r2_original_all,prediction_r2_delta=r2_modified_common-r2_original_common,
           original_only_reaches=n_reaches_original_all-n_reaches_original_common,modified_only_reaches=n_reaches_modified_all-n_reaches_modified_common,
           median_flux_ratio=median_flux_modified_common/median_flux_original_common)]
fwrite(wide,file.path(out,"paired_fit_comparison.csv"))
trend <- common[,{
  x <- log_area; y <- log_flux_ratio
  list(n=.N,log_ratio_area_slope=if (.N>=20 && sd(x)>0) cov(x,y)/var(x) else NA_real_,median_log_ratio=median(y))
},by=.(huc,month,metabolism_metric)]
fwrite(trend,file.path(out,"common_flux_ratio_area_trend.csv"))
cat("National paired fits and validation complete\n")
