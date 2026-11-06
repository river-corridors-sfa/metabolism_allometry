#!/usr/bin/env Rscript
# Reconcile saved 08/10/11 ER classes with baseline fit and eligibility gates.
suppressPackageStartupMessages(library(data.table))
out <- "reports/no_uparea_migration/model_b_er_attribution"
own <- fread(file.path(out,"own_coverage_fits.csv"),colClasses=list(character="huc"))
fixed <- fread(file.path(out,"fixed_network_fits.csv"),colClasses=list(character="huc"))
legacy <- fread("reports/no_uparea_migration/step2_scientific_results/monthly_classes_from_08_using_10_11_rules.csv",colClasses=list(character="huc"))[metabolism_metric=="ER"]
legacy[,version:=fifelse(variant=="original","original","unfiltered")]
legacy <- legacy[,.(huc,month,version,monthly_class)]
b <- fread("reports/no_uparea_migration/model_b_national_assessment/monthly_classes.csv",colClasses=list(character="huc"))[metabolism_metric=="ER" & design=="modified",.(huc,month,version="model_b",monthly_class)]
cls <- rbind(legacy,b)
stopifnot(!anyDuplicated(cls,by=c("huc","month","version")))
hucs <- fread("data/processed/04c_huc08_summary.csv",colClasses=list(character="huc8"))$huc8
grid <- CJ(huc=hucs,month=1:12,version=c("original","unfiltered","model_b"),unique=TRUE)
grid <- merge(grid,own[,.(huc,month,version,n_reaches,slope,r2)],by=c("huc","month","version"),all.x=TRUE)
grid <- merge(grid,cls,by=c("huc","month","version"),all.x=TRUE)
grid[,`:=`(eligible=is.finite(slope),baseline_super=is.finite(slope)&slope>1.05&is.finite(r2)&r2>=.8,stable_super=monthly_class=="Super-linear")]
grid[is.na(stable_super),stable_super:=FALSE]
fwrite(grid,file.path(out,"monthly_gates.csv"))
annual <- grid[,.(eligible_months=sum(eligible),baseline_super_months=sum(baseline_super),stable_super_months=sum(stable_super)),by=.(huc,version)]
annual[,gate:=fifelse(eligible_months<12,"eligibility",fifelse(baseline_super_months<12,"slope_or_r2",fifelse(stable_super_months<12,"section_stability","superlinear")))]
a_legacy <- fread("reports/no_uparea_migration/step2_scientific_results/annual_classes_from_08_using_10_11_rules.csv",colClasses=list(character="huc"))[metabolism_metric=="ER"]
a_legacy[,version:=fifelse(variant=="original","original","unfiltered")]
a_b <- fread("reports/no_uparea_migration/model_b_national_assessment/annual_classes.csv",colClasses=list(character="huc"))[metabolism_metric=="ER" & design=="modified",.(huc,version="model_b",annual_class)]
ac <- rbind(a_legacy[,.(huc,version,annual_class=scaling_class)],a_b)
annual <- merge(annual,ac,by=c("huc","version"),all.x=TRUE)
stopifnot(nrow(annual)==length(hucs)*3L,all((annual$gate=="superlinear")== (annual$annual_class=="Super-linear")))
fwrite(annual,file.path(out,"annual_gates_and_classes.csv"))
counts <- annual[,.(huc8s=.N,annual_super=sum(annual_class=="Super-linear"),annual_super_pct=100*mean(annual_class=="Super-linear"),all_12_eligible=sum(eligible_months==12),all_12_baseline_super=sum(baseline_super_months==12),all_12_stable_super=sum(stable_super_months==12)),by=version]
fwrite(counts,file.path(out,"annual_reconciliation.csv"))
fwrite(annual[,.(huc8s=.N),by=.(version,gate)],file.path(out,"annual_gate_counts.csv"))
wide <- dcast(annual,huc~version,value.var=c("annual_class","gate","eligible_months","baseline_super_months","stable_super_months"))
fwrite(wide,file.path(out,"annual_huc8_three_way.csv"))
fixed_summary <- fixed[,.(groups=.N,eligible=sum(is.finite(slope)),median_slope=median(slope,na.rm=TRUE),baseline_super=sum(is.finite(slope)&slope>1.05&is.finite(r2)&r2>=.8)),by=version]
fwrite(fixed_summary,file.path(out,"fixed_network_summary.csv"))
cat("Class and eligibility reconciliation complete\n")
