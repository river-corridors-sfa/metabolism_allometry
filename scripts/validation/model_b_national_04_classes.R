#!/usr/bin/env Rscript
# Apply the saved 10/11 descriptive stability and annual rules to four designs.
suppressPackageStartupMessages(library(data.table))
out <- "reports/no_uparea_migration/model_b_national_assessment"
orig <- fread("data/processed/08_huc08_fit_stability.csv",colClasses=list(character="huc"))
gpp <- fread("data/processed/08_no_uparea_huc08_fit_stability.csv",colClasses=list(character="huc"))[metabolism_metric=="GPP"]
er <- fread(file.path(out,"model_b_er_stability.csv"),colClasses=list(character="huc"))
orig[,design:="original"]
gpp[,design:="modified"]
er[,design:="modified"]
s <- rbindlist(list(orig,gpp,er),fill=TRUE)
stopifnot(!anyDuplicated(s,by=c("huc","month","metabolism_metric","design")))
required <- c("n_reaches","n_unique_areas","log10_area_range","n_valid","baseline_slope","baseline_r2","slope_min","slope_max","r2_min","median_pct_metabolism_missing","max_pct_metabolism_missing")
s[,supported:=rowSums(is.finite(as.matrix(.SD)))==length(required) & n_reaches>=20 & n_unique_areas>=20 & log10_area_range>=1 & n_valid>0,.SDcols=required]
s[,monthly_class:=fifelse(!supported,"Insufficient data",fifelse(pmin(baseline_r2,r2_min)<.8,"Uncertain",fifelse(pmax(baseline_slope,slope_max)<.95,"Sub-linear",fifelse(pmin(baseline_slope,slope_min)>=.95 & pmax(baseline_slope,slope_max)<=1.05,"Linear",fifelse(pmin(baseline_slope,slope_min)>1.05,"Super-linear","Uncertain")))))]
monthly <- s[,.(huc,month,metabolism_metric,design,monthly_class,baseline_slope,baseline_r2,n_reaches)]
fwrite(monthly,file.path(out,"monthly_classes.csv"))
allhuc <- fread("data/processed/04c_huc08_summary.csv",colClasses=list(character="huc8"))[,.(huc=huc8)]
grid <- CJ(huc=allhuc$huc,month=1:12,metabolism_metric=c("GPP","ER"),design=c("original","modified"),unique=TRUE)
grid <- merge(grid,monthly,by=c("huc","month","metabolism_metric","design"),all.x=TRUE)
annual <- grid[,.(months_present=sum(!is.na(monthly_class)),annual_class=if (sum(!is.na(monthly_class))<12) "Insufficient data" else if (uniqueN(monthly_class)>1) "Seasonally variable" else monthly_class[1]),by=.(huc,metabolism_metric,design)]
fwrite(annual,file.path(out,"annual_classes.csv"))
mp <- dcast(grid,huc+month+metabolism_metric~design,value.var="monthly_class")
mp[is.na(original),original:="Insufficient data"]
mp[is.na(modified),modified:="Insufficient data"]
fwrite(mp[,.(n=.N),by=.(metabolism_metric,original,modified)][order(metabolism_metric,original,modified)],file.path(out,"monthly_class_transitions.csv"))
ap <- dcast(annual,huc+metabolism_metric~design,value.var="annual_class")
fwrite(ap[,.(n=.N),by=.(metabolism_metric,original,modified)][order(metabolism_metric,original,modified)],file.path(out,"annual_class_transitions.csv"))
saved <- fread("data/processed/10_no_uparea_monthly_classes.csv",colClasses=list(character="huc"))[metabolism_metric=="GPP"]
chk <- merge(monthly[design=="modified" & metabolism_metric=="GPP"],saved,by=c("huc","month","metabolism_metric"),suffixes=c("_here","_saved"))
stopifnot(nrow(chk)==nrow(saved),all(chk$monthly_class_here==chk$monthly_class_saved))
fwrite(data.table(modified_gpp_class_rows=nrow(chk),mismatches=0L),file.path(out,"gpp_class_validation.csv"))
sample <- fread("reports/no_uparea_migration/model_b_comparability/sample_monthly_final_classes.csv",colClasses=list(character="huc"))
check_er <- merge(monthly[design=="modified" & metabolism_metric=="ER"],sample[design=="B_all" & metabolism_metric=="ER"],by=c("huc","month","metabolism_metric"),suffixes=c("_national","_sample"))
stopifnot(nrow(check_er)==72L,all(check_er$monthly_class_national==check_er$monthly_class_sample),max(abs(check_er$baseline_slope_national-check_er$baseline_slope_sample))<1e-8)
fwrite(data.table(sample_er_monthly_classes=nrow(check_er),mismatches=0L,max_slope_error=max(abs(check_er$baseline_slope_national-check_er$baseline_slope_sample))),file.path(out,"six_huc8_class_validation.csv"))
cat("National class transitions complete\n")
