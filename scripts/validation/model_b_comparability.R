#!/usr/bin/env Rscript
# Read-only paired coverage and prediction comparison for Maavara vs six-feature GPP + Model B ER.
suppressPackageStartupMessages({library(data.table)})

root <- normalizePath(".")
out <- file.path(root, "reports/no_uparea_migration/model_b_comparability")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
raw <- file.path(root,"data/raw/dryad_data")
refit <- file.path(root,"data/processed/17_maavara_refits/without_uparea")
qcer <- file.path(root,"reports/no_uparea_migration/er_qc_sensitivity/conus/near_zero")
months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

scaling <- fread(file.path(raw,"scaling_dataset_Qmonth_06082024.csv"),select=c("COMID","uparea","order"))
stopifnot(nrow(scaling)==175116L,!anyDuplicated(scaling$COMID))
scaling[, COMID:=as.character(COMID)]
key <- as.data.table(readRDS("data/processed/04c_comid_huc08_key.rds"))[,.(COMID=as.character(comid),huc8=as.character(huc8))]
stopifnot(!anyDuplicated(key$COMID))
geo <- fread("data/processed/04c_huc08_summary.csv",colClasses=list(character="huc8"))[,.(huc8,longitude=centroid_longitude)]
scaling <- merge(scaling,key,by="COMID",all.x=TRUE,sort=FALSE)
scaling <- merge(scaling,geo,by="huc8",all.x=TRUE,sort=FALSE)
stopifnot(nrow(scaling)==175116L)
scaling[, region:=fifelse(is.na(longitude),"Unassigned HUC8",fifelse(longitude< -105,"West",fifelse(longitude< -90,"Central","East")))]
scaling[, area_bin:=cut(uparea,c(-Inf,50,100,1000,10000,Inf),labels=c("<=50","50-100","100-1000","1000-10000",">10000"))]

read_pred <- function(path,value) {
  x <- fread(path)
  stopifnot(nrow(x)==175116L,!anyDuplicated(x$COMID))
  setnames(x,names(x)[2],value)
  x[, COMID:=as.character(COMID)]
  x
}
monthly <- vector("list",12)
for (i in seq_along(months)) {
  mo <- months[i]
  cat("Paired coverage",mo,"\n")
  files <- list(
    read_pred(file.path(raw,"scaling_GPP_pred",paste0("scaling_GPP_pred_",mo,".csv")),"gpp_original"),
    read_pred(file.path(raw,"scaling_ER_pred",paste0("scaling_ER_pred_",mo,".csv")),"er_original"),
    read_pred(file.path(refit,"GPP",paste0("scaling_GPP_pred_",mo,".csv")),"gpp_modified"),
    read_pred(file.path(qcer,paste0("scaling_ER_pred_",mo,".csv")),"er_modified")
  )
  d <- Reduce(function(a,b) merge(a,b,by="COMID",all=TRUE,sort=FALSE),files)
  d <- merge(d,scaling,by="COMID",all.x=TRUE,sort=FALSE)
  stopifnot(nrow(d)==175116L,!anyNA(d$uparea))
  d[, month:=i]
  d[, original_both_est:=is.finite(gpp_original)&is.finite(er_original)]
  d[, modified_both_est:=is.finite(gpp_modified)&is.finite(er_modified)]
  d[, original_both_positive:=original_both_est & gpp_original>0 & er_original>0]
  d[, modified_both_positive:=modified_both_est & gpp_modified>0 & er_modified>0]
  for (metric in c("gpp","er")) {
    a <- paste0(metric,"_original"); b <- paste0(metric,"_modified")
    d[, (paste0(metric,"_paired_positive")):=is.finite(get(a))&get(a)>0&is.finite(get(b))&get(b)>0]
  }
  monthly[[i]] <- d
}
d <- rbindlist(monthly)
stopifnot(nrow(d)==175116L*12L,!anyDuplicated(d,by=c("COMID","month")))
loss <- d[original_both_positive & !modified_both_positive]
loss_reason <- loss[,.(n=.N,unique_comids=uniqueN(COMID)),by=.(
  modified_gpp_status=fifelse(!is.finite(gpp_modified),"missing",fifelse(gpp_modified<=0,"nonpositive","positive")),
  modified_er_status=fifelse(!is.finite(er_modified),"missing",fifelse(er_modified<=0,"nonpositive","positive")))]
fwrite(loss_reason,file.path(out,"positive_coverage_loss_reasons.csv"))
gain <- d[!original_both_positive & modified_both_positive]
gain_reason <- gain[,.(n=.N,unique_comids=uniqueN(COMID)),by=.(
  original_gpp_status=fifelse(!is.finite(gpp_original),"missing",fifelse(gpp_original<=0,"nonpositive","positive")),
  original_er_status=fifelse(!is.finite(er_original),"missing",fifelse(er_original<=0,"nonpositive","positive")))]
fwrite(gain_reason,file.path(out,"positive_coverage_gain_reasons.csv"))

coverage <- function(x, original, modified, label) {
  data.table(definition=label,total=nrow(x),original=sum(x[[original]]),modified=sum(x[[modified]]),
             overlap=sum(x[[original]]&x[[modified]]),
             original_only=sum(x[[original]]&!x[[modified]]),
             modified_only=sum(!x[[original]]&x[[modified]]))
}
coverage_all <- rbind(coverage(d,"original_both_est","modified_both_est","both_nonmissing"),
                      coverage(d,"original_both_positive","modified_both_positive","both_positive"))
coverage_all[, retained_pct:=100*overlap/original]
coverage_all[, modified_pct_of_total:=100*modified/total]
fwrite(coverage_all,file.path(out,"coverage_total.csv"))

coverage_strata <- rbindlist(lapply(c("month","area_bin","region"),function(g) {
  rbindlist(lapply(split(d,by=g,keep.by=TRUE),function(s) {
    z <- rbind(coverage(s,"original_both_est","modified_both_est","both_nonmissing"),
               coverage(s,"original_both_positive","modified_both_positive","both_positive"))
    z[,`:=`(grouping=g,level=as.character(s[[g]][1]))]
    z
  }))
}))
coverage_strata[, retained_pct:=fifelse(original>0,100*overlap/original,NA_real_)]
fwrite(coverage_strata,file.path(out,"coverage_by_stratum.csv"))

reach <- d[,.(original_est_months=sum(original_both_est),modified_est_months=sum(modified_both_est),
              original_positive_months=sum(original_both_positive),modified_positive_months=sum(modified_both_positive)),by=.(COMID,uparea,order,area_bin,region,huc8)]
stopifnot(nrow(reach)==175116L)
reach_summary <- rbindlist(lapply(c("est","positive"),function(s) {
  a <- paste0("original_",s,"_months"); b <- paste0("modified_",s,"_months")
  rbindlist(lapply(c(1L,12L),function(threshold) {
    data.table(definition=s,months_required=threshold,total=nrow(reach),
               original=sum(reach[[a]]>=threshold),modified=sum(reach[[b]]>=threshold),
               overlap=sum(reach[[a]]>=threshold & reach[[b]]>=threshold),
               original_only=sum(reach[[a]]>=threshold & reach[[b]]<threshold),
               modified_only=sum(reach[[a]]<threshold & reach[[b]]>=threshold))
  }))
}))
reach_summary[,retained_pct:=100*overlap/original]
fwrite(reach_summary,file.path(out,"coverage_unique_comids.csv"))
fwrite(reach,file.path(out,"reach_month_counts.csv"))

paired_summary <- function(x,metric,g="overall") {
  a <- x[[paste0(metric,"_original")]]; b <- x[[paste0(metric,"_modified")]]
  valid <- is.finite(a)&is.finite(b)&a>0&b>0
  a <- a[valid]; b <- b[valid]
  if (!length(a)) return(data.table(metric=metric,grouping=g,n=0L,correlation=NA_real_,median_ratio=NA_real_,mean_bias=NA_real_,median_relative_pct=NA_real_,median_original=NA_real_,median_modified=NA_real_))
  data.table(metric=metric,grouping=g,n=length(a),correlation=cor(a,b),
             median_ratio=median(b/a),mean_bias=mean(b-a),
             median_relative_pct=100*median((b-a)/a),
             median_original=median(a),median_modified=median(b))
}
paired_overall <- rbind(paired_summary(d,"gpp"),paired_summary(d,"er"))
fwrite(paired_overall,file.path(out,"paired_magnitudes_overall.csv"))
paired_strata <- rbindlist(lapply(c("month","area_bin","region"),function(g) {
  rbindlist(lapply(split(d,by=g,keep.by=TRUE),function(s) {
    z <- rbind(paired_summary(s,"gpp",g),paired_summary(s,"er",g))
    z[,level:=as.character(s[[g]][1])]
    z
  }))
}))
fwrite(paired_strata,file.path(out,"paired_magnitudes_by_stratum.csv"))

# Newly recovered means the original joint positive condition failed; report modified values without a ratio.
new <- d[!original_both_positive & modified_both_positive]
new_summary <- new[,.(reach_months=.N,unique_comids=uniqueN(COMID),
                       median_gpp=median(gpp_modified),median_er=median(er_modified)),by=.(area_bin,region)]
fwrite(new_summary,file.path(out,"newly_recovered_magnitudes.csv"))

# Descriptive trend: one slope per month plus pooled slope; no independent-row inference.
er_pair <- d[er_paired_positive & is.finite(uparea) & uparea>0,
             .(month,log_area=log10(uparea),log_ratio=log(er_modified/er_original))]
trend <- rbind(er_pair[,.(slope=coef(lm(log_ratio~log_area))[2],n=.N),by=month],
               er_pair[,.(month=0L,slope=coef(lm(log_ratio~log_area))[2],n=.N)])
fwrite(trend,file.path(out,"er_log_ratio_area_trend.csv"))

cat("Saved coverage and paired comparison summaries to",out,"\n")
