## 31_bootstrap_huc8_scaling.R
##
## PLACEHOLDER: bootstrap HUC8 scaling behavior before constructing CONUS maps.
##
## Workflow position
##   07_calculate_huc8_allometric_scaling.Rmd builds cumulative reach values and
##   fits the original model. Script 18 repeats the fit for paired RF variants.
##   Run this script after the chosen cumulative reach table is finalized and
##   before a new HUC8 CONUS map script classifies or displays scaling behavior.
##   Do not bootstrap the already aggregated slopes in results CSV files: sample
##   reaches and refit log10(cumulative flux) ~ log10(cumulative area).
##
## Reference: Regier et al. (2025), JGR: Biogeosciences,
## https://doi.org/10.1029/2024JG008344, section 2.3. They sampled 50 reaches
## spanning stream orders 1,000 times within each basin x cumulative-HEF decile,
## summarized slope/intercept/R2, and used slope percentile intervals plus
## mean R2 >= 0.8 to label scaling behavior.
##
## Decisions to settle during the code audit (the paper's HEF bins do not map
## directly to our HUC8 x month x GPP/ER analysis):
##   1. Use the canonical script 07 input or the no-uparea RF reach-level input
##      from script 18; make the map and bootstrap use the same model variant.
##   2. Define the resampling unit and scheme for connected, nested reaches.
##      Check stream-order coverage and the effect of spatial dependence.
##   3. Set minimum eligible reaches/area range, sample size for small HUC8s,
##      number of iterations, fixed seed, and treatment of failed fits.
##   4. Define annual behavior: summarize monthly bootstrap fits first, then
##      classify HUC8s, or use another explicit annual estimand.
##   5. Decide whether the paper's mean R2 cutoff of 0.8 is appropriate here.
##
## Intended inputs (not yet wired):
##   data/processed/07_huc08_allometry_input.rds
##   or reach-level cumulative data from the selected script 18 RF variant.
## Intended outputs (not yet created):
##   data/processed/31_huc8_bootstrap_replicates.rds
##   data/processed/31_huc8_bootstrap_summary.csv
##   Columns should include HUC8, month, metabolism, model variant, eligible n,
##   successful replicate n, slope/intercept/R2 summaries, slope 2.5/97.5%
##   limits, and classification (uncertain/sublinear/linear/superlinear).

stop("Placeholder only: resolve the bootstrap design during the audit before running.")
