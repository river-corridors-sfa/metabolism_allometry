# PaperBot 3000: no-uparea metabolism migration

Status: **CODE IMPLEMENTED — NOT RUN** (2026-10-29). Implementation base revision: `23c32f385a32eae8226b0a768905faed01b249e0`. No pipeline scripts or validation were executed during implementation.

## Decision and provenance

The primary HUC8 allometry will use the six-feature GPP and ER refits that exclude `uparea` as a random-forest predictor, reducing a circular link between modeled metabolism and drainage-area scaling. `uparea` remains in the network/coverage checks and the allometric analysis retains cumulative catchment area as its predictor. Original downloaded Maavara predictions and existing processed baseline outputs remain untouched for comparison.

Script `17_refit_maavara_with_without_uparea.py` removes exactly `uparea` in `feature_names`; the other retained features include upstream land-cover proportions. Its `run_manifest.json` identifies the September 22, 2026 full-CONUS, all-12-month run, with sources `data/raw/dryad_data/appling_data_06202024.csv` and `data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv`. Production inputs are the 24 CSVs under `data/processed/17_maavara_refits/without_uparea/{GPP,ER}/scaling_{GPP,ER}_pred_{Jan..Dec}.csv`. Each has `COMID` and `{process}_pred_meanQ_{month}`. These are predictions for every scaling-data COMID, **not held-out test predictions**. Grain is one COMID per nominal mean-flow month, with areal rate in g O₂ m⁻² d⁻¹. ER is positive respiration magnitude; script 17 negates the negative StreamMetabolizer training response before fitting. Five source rows per month have incomplete covariates and are retained as missing predictions. Negative model predictions remain negative until the unchanged script 07 recodes them to NA.

Original inputs: `data/raw/dryad_data/scaling_{GPP,ER}_pred/scaling_{GPP,ER}_pred_{Jan..Dec}.csv`. Original assembly: `data/processed/01_assembled_maavara_monthly_metabolism.rds`. Other baseline artifacts: `03_filtered_maavara_monthly_metabolism.rds`, `04c_comid_huc08_key.rds`, `04c_huc08_summary.csv`, `07_huc08_allometry_{input.rds,results.csv}`, `07_negative_rate_recode_summary.csv`, `08_huc08_fit_stability.csv`, and original 10–12 figures in `figures/`. Their exact original run commands, external-layer snapshots, and file hashes were not recorded by the old pipeline. Preserve these files; do not rerun the original scripts from the base revision over them.

## Change and dependency ledger

| Step | Action | Input and output decision |
|---|---|---|
| 00 | Reuse | Package setup, no prediction dependency. |
| 01 | Edited; run | Explicit no-uparea 17 CSVs → `01_assembled_no_uparea_monthly_metabolism.rds` and log. Strict month, columns, unique keys, and original key coverage checks. |
| 02 | Reuse | MERIT crop and geometry do not depend on rates. |
| 03 | Edited reference; **do not rerender for migration** | Full Rmd now points to no-uparea input/output, but its spatial recomputation and shapefile writes are unnecessary. Run `03_refresh_no_uparea_metabolism.R`, which joins the existing approved catchments and writes `03_filtered_no_uparea_monthly_metabolism.rds`; it asserts exact baseline filtered key coverage. Existing 03 spatial layers and coverage remain authoritative. |
| 04 | Edited reference; optional diagnostic only | Single-HUC diagnostic is not production input. No rerun required for the primary pipeline. |
| 04b | Reuse | Ecoregion lookup depends on HUC geometry, not rates. |
| 05 | Reuse | Routing edges and upstream lookup depend on 03 spatial layers, not rates. |
| 04c | Edited; run **after 05 exists** | No-uparea filtered rates → variant-specific COMID key and HUC8 summary. Cached `04c_huc08_boundaries.gpkg` is reused when present; do not refresh it for this migration. |
| 06 | Reuse | Routing diagnostic, no production output required. |
| 07 | Edited; render | Reads no-uparea 01 independently of 03; writes variant-specific screened input, fits, negative-rate summary, and ordered screen attrition. Same network, thresholds, recode, and fit criteria. |
| 08 | Edited; run | No-uparea 07 input and 04c key → `08_no_uparea_huc08_fit_stability.csv`. |
| 09 | None | No production script. |
| 10–11 | Edited; run | Read no-uparea 08 output; write separately named map/bar figures. |
| 12 | Edited; run | Read no-uparea 07 and 04c; write separately named segmented figures and summary. |

The current 07 already uses retained routed reaches and the existing <10% area and metabolism-coverage screens. Its negative-rate convention, ≥20 reaches, ≥20 distinct areas, one-decade range, OLS fit, and simple classes are unchanged. Script 08 retains its section deletion design and classification rules. Script 12 retains its segmented model. No model-selection switch was added.

## User-run order

Run from repository root in an R environment with the packages used by the scripts. Stop at the first error; check that each named output exists and is newer than its input before continuing. Use fresh R sessions for the large 07 and 08 jobs. Rmd 07 expects `scripts/` as its working directory; `rmarkdown::render` supplies this default. Its HTML name is variant-specific so baseline HTML is preserved.

```bash
Rscript scripts/01_assemble_metabolism.R
Rscript scripts/03_refresh_no_uparea_metabolism.R
Rscript scripts/04c_build_huc08_index.R
Rscript -e 'rmarkdown::render("scripts/07_calculate_huc8_allometric_scaling.Rmd", output_file = "07_no_uparea_huc8_allometric_scaling.html")'
Rscript scripts/08_huc8_fit_stability.R
Rscript scripts/10_figure1_scaling_maps.R
Rscript scripts/11_ecoregion_scaling_bars.R
Rscript scripts/12_er_gpp_area_segmented_hex.R
Rscript scripts/validation/validate_no_uparea_migration.R
```

01 succeeds only if all 24 CSVs match their expected schema and the original COMID-month key universe exactly. 03 refresh succeeds only if approved-catchment filtering reproduces original 03 keys. 04c should create its two `04c_no_uparea_*` outputs. 07 should create `07_no_uparea_huc08_allometry_input.rds`, `07_no_uparea_huc08_allometry_results.csv`, `07_no_uparea_negative_rate_recode_summary.csv`, and `07_no_uparea_screen_attrition.csv`. 08 should create its no-uparea CSV. 10–12 should create `figure_1_no_uparea*`, `figure_2_no_uparea*`, `figure_3_no_uparea*`, `figure_s2_no_uparea*`, and `12_no_uparea_er_gpp_area_segmented_fit.csv`. Final validation writes CSV ledgers, `summary.txt`, and two paired plots to `reports/no_uparea_migration/`. Existing original outputs are never target paths in this run order.

Inspect `source_coverage.csv` for unmatched keys, `source_rate_counts.csv` for missing/negative/zero/nonfinite changes, `paired_refit_test_skill.csv` for paired test R²/RMSE/MAE changes, `screened_composition.csv` and `no_uparea_screen_attrition.csv` for filtering effects, `paired_allometry_fits.csv` and `classification_transitions.csv` for changed fits/classes, and `segmented_comparison.csv` for breakpoints/slopes. The model-skill table compares **two new refits on the same seeded 70/30 split** (147,108 test rows per process in saved metrics). It does not establish test skill relative to the published original forest; the latter has only full-network prediction comparisons. The saved metric CSV gives identical test counts but not held-out row IDs or predictions, so exact record-level identity is inferred from script 17's shared split and source, rather than independently verifiable from saved test artifacts.

## Completion and validation ledger

| Item | Status | Observed / acceptance work |
|---|---|---|
| 01, 03 refresh, 04c, 07, 08, 10, 11, 12 | USER RUN REQUIRED | No new outputs generated by this migration yet. |
| 00, 02, 04b, 05, 06 | VALIDATED FOR REUSE BY DEPENDENCY | Their existing outputs are independent of rate values; provenance caveats below remain. |
| Predictive skill | USER RUN REQUIRED | Saved paired refit test metrics exist; inspect emitted deltas. No arbitrary deterioration tolerance. |
| Source substitution and coverage | USER RUN REQUIRED | Exact key equality is a hard ingestion gate; inspect rate distributions and negative/missing changes. |
| Attrition and sample composition | USER RUN REQUIRED | Compare baseline/new 01, 03, screened 07, eligible fits, and 08 outputs. Explain all differences before acceptance. Original 07 did not save per-rule exclusions, so its full stepwise rule ledger cannot be reconstructed exactly from final artifacts. |
| Scientific results | USER RUN REQUIRED | Inspect paired slopes, R², sample sizes, classes, changed HUC8s, and segmented summaries. |
| Scientific acceptance | NOT STARTED | Requires user-run outputs and inspection; no result is claimed here. |

Expected source coverage is 175,116 COMIDs × 12 months in each process, based on the earlier audit. The original 01 contains 2,101,392 COMID-months; original filtered 03 contains 1,968,876. Original 07 has 3,815,080 screened reach-month-process rows and 44,450 eligible fits; 08 has 45,456 rows, including ineligible baselines. These are comparison baselines, not expected no-uparea counts. Differences in negative rates can change coverage and fit eligibility; every observed difference needs a quantified cause. The original 03→07 exact per-rule attrition is unavailable because the original 07 saved only its final screened input.

## Issue register (separate from migration)

Unresolved from `PB_pipeline_00_12_audit.md`: route-pair uniqueness and cycle/self-inclusion integrity of the 05 closure; 117-COMID difference between original filtered metabolism and 04c flowline key plus boundary HUC assignment differences; scientific interpretation of negative-rate recoding and coverage thresholds; dependence among nested cumulative reach observations and limits of class/breakpoint inference; incomplete historical manifests and cached external layer provenance. This migration neither changes nor resolves those issues.
