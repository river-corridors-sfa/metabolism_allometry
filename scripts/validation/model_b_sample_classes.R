#!/usr/bin/env Rscript
# Sample-only reproduction of 08 section deletions and 10/11 monthly/annual rules.
suppressPackageStartupMessages({library(data.table);library(igraph)})
out <- "reports/no_uparea_migration/model_b_comparability"
x <- fread(file.path(out,"sample_allometry_fit_input.csv"),colClasses=list(character=c("huc","downstream_comid")))
edges <- as.data.table(readRDS("data/processed/merit_data/05_edges_downstream.rds"))
edges[,`:=`(from_comid=as.character(from_comid),to_comid=as.character(to_comid))]
stopifnot(!anyDuplicated(edges$from_comid))
incoming <- edges[,.(n_incoming=.N),by=to_comid]
edges <- merge(edges,incoming,by="to_comid")
edges <- edges[n_incoming==1,.(from_comid,to_comid)]

fit_stats <- function(z) {
  if (nrow(z)<20 || uniqueN(z$cumulative_area_km2)<20 ||
      log10(max(z$cumulative_area_km2)/min(z$cumulative_area_km2))<1) return(NULL)
  fit <- lm(log10(cumulative_flux_g_o2_d)~log10(cumulative_area_km2),data=z)
  c(slope=unname(coef(fit)[2]),r2=summary(fit)$r.squared)
}
classify <- function(base,deletions,n,narea,span) {
  if (is.null(base) || !length(deletions) || n<20 || narea<20 || span<1) return("Insufficient data")
  slopes <- c(base["slope"],vapply(deletions,`[[`,numeric(1),"slope"))
  r2s <- c(base["r2"],vapply(deletions,`[[`,numeric(1),"r2"))
  if (min(r2s)<.8) return("Uncertain")
  if (max(slopes)<.95) return("Sub-linear")
  if (min(slopes)>=.95 && max(slopes)<=1.05) return("Linear")
  if (min(slopes)>1.05) return("Super-linear")
  "Uncertain"
}

rows <- list(); idx <- 0L
for (design_i in c("original_all","B_all","original_common","B_common")) {
  data <- x[design==design_i]
  key <- unique(data[,.(comid=downstream_comid,huc)])
  stopifnot(!anyDuplicated(key$comid))
  e <- merge(edges,key,by.x="from_comid",by.y="comid")
  setnames(e,"huc","huc_from")
  e <- merge(e,key,by.x="to_comid",by.y="comid")
  e <- e[huc_from==huc,.(from_comid,to_comid)]
  g <- graph_from_data_frame(e,directed=TRUE,vertices=key$comid)
  stopifnot(is_dag(g))
  membership <- components(g,mode="weak")$membership
  blocks <- data.table(downstream_comid=names(membership),block=as.integer(membership))
  data <- merge(data,blocks,by="downstream_comid",all.x=TRUE)
  stopifnot(!anyNA(data$block))
  groups <- split(data,by=c("huc","month","metabolism_metric"))
  for (z in groups) {
    base <- fit_stats(z)
    del <- lapply(unique(z$block),function(b) fit_stats(z[block!=b]))
    del <- Filter(Negate(is.null),del)
    n <- nrow(z); narea <- uniqueN(z$cumulative_area_km2)
    span <- log10(max(z$cumulative_area_km2)/min(z$cumulative_area_km2))
    idx <- idx+1L
    rows[[idx]] <- data.table(design=design_i,huc=z$huc[1],month=z$month[1],metabolism_metric=z$metabolism_metric[1],
                              n_reaches=n,n_unique_areas=narea,log10_area_range=span,
                              n_blocks=uniqueN(z$block),n_valid=length(del),
                              baseline_slope=if (is.null(base)) NA_real_ else base["slope"],
                              baseline_r2=if (is.null(base)) NA_real_ else base["r2"],
                              slope_min=if (length(del)) min(vapply(del,`[[`,numeric(1),"slope")) else NA_real_,
                              slope_max=if (length(del)) max(vapply(del,`[[`,numeric(1),"slope")) else NA_real_,
                              r2_min=if (length(del)) min(vapply(del,`[[`,numeric(1),"r2")) else NA_real_,
                              monthly_class=classify(base,del,n,narea,span))
  }
  cat("Classified",design_i,"\n")
}
monthly <- rbindlist(rows)
fwrite(monthly,file.path(out,"sample_monthly_final_classes.csv"))
annual <- monthly[,.(months=.N,
                     annual_class=if (.N<12) "Insufficient data" else if (uniqueN(monthly_class)>1) "Seasonally variable" else monthly_class[1]),
                  by=.(huc,metabolism_metric,design)]
fwrite(annual,file.path(out,"sample_annual_final_classes.csv"))

# Check the new sample replay against the actual saved 08/10 results for the
# existing six-feature GPP data. This validates the section logic independently.
saved <- fread("data/processed/10_no_uparea_monthly_classes.csv",colClasses=list(character="huc"))
check <- merge(monthly[design=="B_all" & metabolism_metric=="GPP",.(huc,month,metabolism_metric,monthly_class)],
               saved[,.(huc,month,metabolism_metric,monthly_class_saved=monthly_class)],
               by=c("huc","month","metabolism_metric"))
validation <- data.table(n=nrow(check),mismatches=sum(check$monthly_class!=check$monthly_class_saved))
fwrite(validation,file.path(out,"sample_class_replay_validation.csv"))
if (validation$n!=72L || validation$mismatches!=0L) stop("Sample class replay did not match saved 10 GPP classes")
cat("Saved sample final monthly and annual classes\n")
