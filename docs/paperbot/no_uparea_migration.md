# PaperBot 3000: no-uparea metabolism migration

Status: **CODE IMPLEMENTED — NOT RUN** (2026-10-29). Implementation base revision: `23c32f385a32eae8226b0a768905faed01b249e0`. No pipeline scripts or validation were executed during implementation.

## Decision and provenance

The primary HUC8 allometry will use the six-feature GPP and ER refits that exclude `uparea` as a random-forest predictor, reducing a circular link between modeled metabolism and drainage-area scaling. `uparea` remains in the network/coverage checks and the allometric analysis retains cumulative catchment area as its predictor. Original downloaded Maavara predictions and existing processed baseline outputs remain untouched for comparison.

Script `01a_refit_maavara_with_without_uparea.py` removes exactly `uparea` in `feature_names`; the other retained features include upstream land-cover proportions. Its `run_manifest.json` identifies the September 22, 2026 full-CONUS, all-12-month run, with sources `data/raw/dryad_data/appling_data_06202024.csv` and `data/raw/dryad_data/scaling_dataset_Qmonth_06082024.csv`. Production inputs are the 24 CSVs under `data/processed/17_maavara_refits/without_uparea/{GPP,ER}/scaling_{GPP,ER}_pred_{Jan..Dec}.csv`. Each has `COMID` and `{process}_pred_meanQ_{month}`. These are predictions for every scaling-data COMID, **not held-out test predictions**. Grain is one COMID per nominal mean-flow month, with areal rate in g O₂ m⁻² d⁻¹. ER is positive respiration magnitude; script 17 negates the negative StreamMetabolizer training response before fitting. Five source rows per month have incomplete covariates and are retained as missing predictions. Negative model predictions remain negative until the unchanged script 07 recodes them to NA.

Original inputs: `data/raw/dryad_data/scaling_{GPP,ER}_pred/scaling_{GPP,ER}_pred_{Jan..Dec}.csv`. Original assembly: `data/processed/01_assembled_maavara_monthly_metabolism.rds`. Other baseline artifacts: `03_filtered_maavara_monthly_metabolism.rds`, `04c_comid_huc08_key.rds`, `04c_huc08_summary.csv`, `07_huc08_allometry_{input.rds,results.csv}`, `07_negative_rate_recode_summary.csv`, `08_huc08_fit_stability.csv`, and original 10–12 figures in `figures/`. Their exact original run commands, external-layer snapshots, and file hashes were not recorded by the old pipeline. Preserve these files; do not rerun the original scripts from the base revision over them.

## Change and dependency ledger

| Step | Action | Input and output decision |
|---|---|---|
| 00 | Reuse | Package setup, no prediction dependency. |
| 01a | Reuse existing outputs | Python refit reads raw Dryad training, scaling, and original prediction files. Its saved full-network predictions are the input to 01b; no refit is needed for this migration run. |
| 01b | Edited; run | Explicit no-uparea 01a CSVs → `01_assembled_no_uparea_monthly_metabolism.rds` and log. Strict month, columns, unique keys, and process-specific key coverage checks against the 24 original **raw** Dryad prediction CSVs. It no longer depends on the old assembled RDS. |
| 02 | Reuse | MERIT crop and geometry do not depend on rates. |
| 03 | Edited; render after 01 | The full Rmd performs its existing spatial filtering, then joins no-uparea 01 rates and writes `03_filtered_no_uparea_monthly_metabolism.rds`. It writes a separate no-uparea spatial-coverage checkpoint. If approved spatial layers already exist, it checks their COMID/HUC and routing keys and preserves them; a spatial-key change stops the run for review. It asserts exact original filtered metabolism key coverage. |
| 04 | Edited reference; optional diagnostic only | Single-HUC diagnostic is not production input. No rerun required for the primary pipeline. |
| 04b | Reuse | Ecoregion lookup depends on HUC geometry, not rates. |
| 04c | Edited; run after 04b | No-uparea filtered rates and 03 spatial coverage → variant-specific COMID key and HUC8 summary. Its retained-edge indicator is computed directly from filtered flowlines using 05's rule, so 04c no longer depends on a later script. Cached `04c_huc08_boundaries.gpkg` is reused when present; do not refresh it for this migration. |
| 05 | Reuse | Routing edges and upstream lookup depend on 03 spatial layers, not rates. In a fresh full build, run 05 after 04c and before 07. |
| 06 | Reuse | Routing diagnostic, no production output required. |
| 07 | Edited; render | Reads no-uparea 01 independently of 03; writes variant-specific screened input, fits, negative-rate summary, and ordered screen attrition. Same network, thresholds, recode, and fit criteria. |
| 08 | Edited; run | No-uparea 07 input and 04c key → `08_no_uparea_huc08_fit_stability.csv`. |
| 09 | None | No production script. |
| 10–11 | Edited; run | Read no-uparea 08 output; write separately named map/bar figures. |
| 12 | Edited; run | Read no-uparea 07 and 04c; write separately named segmented figures and summary. |

The current 07 already uses retained routed reaches and the existing <10% area and metabolism-coverage screens. Its negative-rate convention, ≥20 reaches, ≥20 distinct areas, one-decade range, OLS fit, and simple classes are unchanged. Script 08 retains its section deletion design and classification rules. Script 12 retains its segmented model. No model-selection switch was added.

Readability update: the top-level `local({ ... })` wrappers were removed from 08, 10, 11, and 12. This changes variable scope when sourced interactively but not the `Rscript` calculations or output paths. The separate `x_` exploratory scripts were outside the requested 00–12 production sequence.

## User-run order

Run from repository root in an R environment with the packages used by the scripts. Stop at the first error; check that each named output exists and is newer than its input before continuing. Use fresh R sessions for the large 07 and 08 jobs. Rmd 07 expects `scripts/` as its working directory; `rmarkdown::render` supplies this default. Its HTML name is variant-specific so baseline HTML is preserved.

```bash
Rscript scripts/01b_assemble_metabolism.R
Rscript -e 'rmarkdown::render("scripts/03_filter_HUC8s.Rmd", output_file = "03_no_uparea_filter_HUC8s.html")'
Rscript scripts/04c_build_huc08_index.R
Rscript -e 'rmarkdown::render("scripts/07_calculate_huc8_allometric_scaling.Rmd", output_file = "07_no_uparea_huc8_allometric_scaling.html")'
Rscript scripts/08_huc8_fit_stability.R
Rscript scripts/10_figure1_scaling_maps.R
Rscript scripts/11_ecoregion_scaling_bars.R
Rscript scripts/12_er_gpp_area_segmented_hex.R
Rscript scripts/validation/validate_no_uparea_migration.R
```

01 succeeds only if all 24 refit CSVs match their expected schema and both processes' COMID-month keys match the original raw Dryad CSVs exactly. 03 succeeds only if spatial keys agree with existing approved layers and the filtered metabolism keys match the original 03 output; it creates `03_filtered_no_uparea_monthly_metabolism.rds`, `03_no_uparea_huc_spatial_coverage.rds`, and variant-specific HTML. 04c should create its two `04c_no_uparea_*` outputs. 07 should create `07_no_uparea_huc08_allometry_input.rds`, `07_no_uparea_huc08_allometry_results.csv`, `07_no_uparea_negative_rate_recode_summary.csv`, and `07_no_uparea_screen_attrition.csv`. 08 should create its no-uparea CSV. 10–12 should create `figure_1_no_uparea*`, `figure_2_no_uparea*`, `figure_3_no_uparea*`, `figure_s2_no_uparea*`, and `12_no_uparea_er_gpp_area_segmented_fit.csv`. Final validation writes CSV ledgers, `summary.txt`, and two paired plots to `reports/no_uparea_migration/`. Existing original outputs are never target paths in this run order.

Inspect `source_coverage.csv` for unmatched keys, `source_rate_counts.csv` for missing/negative/zero/nonfinite changes, `paired_refit_test_skill.csv` for paired test R²/RMSE/MAE changes, `screened_composition.csv` and `no_uparea_screen_attrition.csv` for filtering effects, `paired_allometry_fits.csv` and `classification_transitions.csv` for changed fits/classes, and `segmented_comparison.csv` for breakpoints/slopes. The model-skill table compares **two new refits on the same seeded 70/30 split** (147,108 test rows per process in saved metrics). It does not establish test skill relative to the published original forest; the latter has only full-network prediction comparisons. The saved metric CSV gives identical test counts but not held-out row IDs or predictions, so exact record-level identity is inferred from script 17's shared split and source, rather than independently verifiable from saved test artifacts.

## Completion and validation ledger

| Item | Status | Observed / acceptance work |
|---|---|---|
| 01b | OUTPUT GENERATED — VALIDATION PENDING | User ran the assembly on 2026-10-29, before the raw-file coverage-check revision; its log reports 2,101,392 COMID-month rows, 175,116 COMIDs, and 60 missing values per process. Full validation remains pending. |
| 03 | OUTPUT GENERATED — VALIDATION PENDING | User ran it on 2026-10-29; no-uparea filtered metabolism and spatial-coverage checkpoint exist. Full validation remains pending. |
| 04c | OUTPUT GENERATED — VALIDATION PENDING | User ran it on 2026-10-29; new key has 163,956 unique COMIDs, matching baseline count, and summary has 1,894 HUC8 rows, matching baseline count. Full validation remains pending. |
| 07 | OUTPUT GENERATED — VALIDATION PENDING | User ran it on 2026-10-29; screened input, fits, negative-rate log, and ordered attrition ledger exist. The ledger ends at 3,610,917 screened reach-month-process rows versus 3,815,080 in the original. The new fit CSV has 43,778 rows versus 44,450 original fits. Changes require paired validation and explanation before acceptance. |
| 08 | OUTPUT GENERATED — VALIDATION PENDING | User ran it on 2026-10-29; the no-uparea CSV has 45,456 data rows, the same count as the original 08 file. Classification and slope comparisons remain pending. |
| 10, 11 | GENERATED AND CLASS TABLES VERIFIED | Both scripts reran on 2026-10-29 after adding saved monthly/annual class tables; tables agree exactly by key and class. Figure content was not visually reviewed in this step. |
| 12 | GENERATED — VALIDATION PENDING | No-uparea segmented summary and figures exist; full reconciliation is reserved for later review steps. |
| 00, 02, 04b, 05, 06 | VALIDATED FOR REUSE BY DEPENDENCY | Their existing outputs are independent of rate values; provenance caveats below remain. |
| Predictive skill | USER RUN REQUIRED | Saved paired refit test metrics exist; inspect emitted deltas. No arbitrary deterioration tolerance. |
| Source substitution and coverage | USER RUN REQUIRED | Exact key equality is a hard ingestion gate; inspect rate distributions and negative/missing changes. |
| Attrition and sample composition | USER RUN REQUIRED | Compare baseline/new 01, 03, screened 07, eligible fits, and 08 outputs. Explain all differences before acceptance. Original 07 did not save per-rule exclusions, so its full stepwise rule ledger cannot be reconstructed exactly from final artifacts. |
| Scientific results | USER RUN REQUIRED | Inspect paired slopes, R², sample sizes, classes, changed HUC8s, and segmented summaries. |
| Scientific acceptance | NOT STARTED | Requires user-run outputs and inspection; no result is claimed here. |

Expected source coverage is 175,116 COMIDs × 12 months in each process, based on the earlier audit. The original 01 contains 2,101,392 COMID-months; original filtered 03 contains 1,968,876. Original 07 has 3,815,080 screened reach-month-process rows and 44,450 eligible fits; 08 has 45,456 rows, including ineligible baselines. These are comparison baselines, not expected no-uparea counts. Differences in negative rates can change coverage and fit eligibility; every observed difference needs a quantified cause. The original 03→07 exact per-rule attrition is unavailable because the original 07 saved only its final screened input.

## Issue register (separate from migration)

Unresolved from `PB_pipeline_00_12_audit.md`: route-pair uniqueness and cycle/self-inclusion integrity of the 05 closure; 117-COMID difference between original filtered metabolism and 04c flowline key plus boundary HUC assignment differences; scientific interpretation of negative-rate recoding and coverage thresholds; dependence among nested cumulative reach observations and limits of class/breakpoint inference; incomplete historical manifests and cached external layer provenance. This migration neither changes nor resolves those issues.

## Post-migration review — step 1 (2026-10-29)

**Status: CONFIRMED for sample-size accounting; scientific acceptance remains UNRESOLVED.** Initial saved-output audit: `scripts/validation/step1_sample_size_audit.R` and `reports/no_uparea_migration/step1_sample_size/`. The authorized follow-up replayed only 07's accumulation and screens, using saved 05 routing and both saved 01 assemblies: `scripts/validation/step1_filter_attribution.R` and `reports/no_uparea_migration/step1_filter_attribution/`. It did not rerun model fitting, routing construction, or the full pipeline. A separate inexpensive fit-eligibility replay is in `scripts/validation/step1_fit_eligibility.R`.

| Check | GPP | ER | Evidence |
|---|---:|---:|---|
| Original → no-uparea negative source rates, all 01 COMID-months | 3,300 → 8,850 | 71,968 → 152,517 | `rate_counts.csv` |
| Original → no-uparea negative rates after 03 catchment filter | 3,062 → 8,021 | 68,911 → 146,288 | `rate_counts.csv` |
| New missing source values | 60 | 60 | `source_rate_transitions.csv`; 24 per process after 03 |
| Original → no-uparea 07 screened rows | 1,954,228 → 1,951,894 (−2,334; −0.12%) | 1,860,852 → 1,659,023 (−201,829; −10.85%) | `screened_counts.csv` |
| Original-only / no-uparea-only screened keys | 3,941 / 1,607 | 232,091 / 30,262 | `screened_key_changes.csv` |
| Reaches lost / gained completely across all months | 0 / 0 | 5,843 / 1,158 | `reach_changes.csv` |
| Original → no-uparea eligible fits | 22,448 → 22,447 (−1) | 22,002 → 21,331 (−671) | `fit_counts_by_process.csv` |

Both variants retain exactly 2,101,392 source COMID-month keys and 1,968,876 post-03 keys, with 175,116 and 164,073 COMIDs respectively (`stage_counts.csv`). The no-uparea 07 screen ledger starts with 3,933,744 process rows and removes 384 missing HUC, 75,397 nonpositive/missing flux, 20,767 area-coverage, and 226,279 metabolism-coverage rows, leaving 3,610,917; no rows were removed for missing month, nonpositive area, or nonfinite metabolism coverage (`no_uparea_screen_rules.csv`). These are *sequential new-run exclusions*, not differences from original exclusions.

The follow-up replay reproduced **all 24** saved 07 retained-key sets exactly (12 months × two prediction variants); `replay_verification.csv` has zero missing or extra keys throughout. Its no-uparea sequential ledger also matches the saved `07_no_uparea_screen_attrition.csv` exactly. Thus every screened-row difference has an assigned first failing rule under 07's existing filter order (`changed_screened_rows.csv`, `changed_key_rule_transitions.csv`):

| Process and direction | Newly fails / newly passes positive cumulative flux | Newly fails / newly passes <10% metabolism-area coverage | Total changed rows |
|---|---:|---:|---:|
| GPP original-only | 1,460 | 2,481 | 3,941 |
| GPP no-uparea-only | 317 | 1,290 | 1,607 |
| ER original-only | 72,181 | 159,910 | 232,091 |
| ER no-uparea-only | 885 | 29,377 | 30,262 |

The original 07 replay excludes 2,300 rows at positive flux and 95,201 at metabolism coverage; the no-uparea replay excludes 75,397 and 226,279 respectively. Missing HUC, area, and month screening do not explain variant differences; the area-coverage count changes by 12 only because the preceding positive-flux rule changes which rows reach it. Among original-only ER screened rows, 93,884 have a newly negative **focal** ER prediction and six newly missing focal predictions (`screened_change_focal_rate_status.csv`); other losses can arise from changes anywhere in their upstream network. The large ER shift is therefore mechanically consistent with 07's negative-to-NA rule and cumulative coverage screen, but whether this treatment of negative modeled ER is scientifically appropriate remains a decision for the research team.

The separate replay of 07's three fit-eligibility predicates reproduces both saved fit-key sets exactly (`changed_fit_eligibility.csv`, `fit_eligibility_transitions.csv`). Of 888 ER fits lost, 392 fall below 20 screened reaches and 496 fall below one decade of area range; of 217 ER fits gained, 32 newly reach 20 screened reaches and 185 newly span one decade. GPP loses three fits (two reach count, one area range) and gains two (one each). No changed fit is attributed first to the 20-distinct-area rule. These transitions explain the net 671 ER and one GPP fit losses. The accounting verifies the implementation outcome; it does not establish an acceptable scientific effect size or resolve the dependence among routed cumulative observations.

**2026-10-30 exploratory interpolation sensitivity (scientific decision UNRESOLVED).** A read-only scenario replaced 33,269 eligible negative no-uparea reach-month rates using positive direct upstream and downstream neighbors, averaging multiple upstream rates first. Replaying 07 screens added 39,572 rows (38,988 ER, 584 GPP), closing 19.38% of the net screened-row gap; 29,549 of these rows were specifically original-only losses. No production input or threshold changed. Exact assumptions, counts, and limitations are in `reports/no_uparea_migration/step1_neighbor_imputation_sensitivity/README.md`; code is `scripts/validation/step1_neighbor_imputation_sensitivity.R`. This does not evaluate altered slopes or classifications and is not a recommendation to adopt interpolation.

**2026-10-30 ER source/median follow-up.** The archived ER notebook and refit code train on all 490,358 rows after flipping raw ER signs; 43,010 raw responses were positive, yielding negative training targets under the intended positive-magnitude convention. Of those, 34,112 exceed +0.5 and 25,596 have intervals entirely above zero. The full-network median remains positive and changes only modestly: +3.613 original → +3.622 with-uparea refit → +3.748 no-uparea refit. The large positive-to-negative median shift applies only to the selected sign-flipping rows. Full provenance, quality-field checks, predictor associations, and caveats are in `reports/no_uparea_migration/step1_er_label_and_shift/README.md`. Treatment of wrong-sign source responses is a scientific decision; no model or analytical method was changed.

## Post-migration review — step 2 (2026-10-29)

**Status: CONFIRMED for saved comparisons; FIXED for missing 10/11 classification tables; scientific acceptance remains UNRESOLVED.** Read-only comparison code: `scripts/validation/step2_scientific_results.R`. Detailed outputs: `reports/no_uparea_migration/step2_scientific_results/`. Scripts 10 and 11 were rerun only to save their monthly and annual class tables; no model, routing, 07, or 08 run occurred.

The paired refits use the same training source, response, 70/30 split, and random seed in `01a_refit_maavara_with_without_uparea.py`; both saved test rows have `n = 147,108` per process. Exact held-out row IDs/predictions were not saved, so record identity is supported by source code rather than independently checked IDs. These are comparable **with-uparea versus without-uparea refits**, not held-out skill estimates versus the published original model. Saved test metrics (`paired_refit_test_skill.csv`): GPP R² 0.827608 → 0.825142 (−0.002466), RMSE 1.931705 → 1.945472 (+0.71%), MAE 0.903738 → 0.913654 (+1.10%); ER R² 0.777697 → 0.775994 (−0.001703), RMSE 2.605172 → 2.615133 (+0.38%), MAE 1.507041 → 1.514744 (+0.51%). These declines are numerically small but no acceptance tolerance was imposed.

The 07 fit comparison has 22,445 matched GPP and 21,114 matched ER HUC8-month fits; GPP has 3 original-only and 2 new-only fits, ER 888 original-only and 217 new-only (`fit_composition.csv`). Among matched fits, median paired exponent change is +0.00499 GPP and +0.03360 ER; median paired R² change is −0.00905 GPP and −0.01372 ER (`matched_fit_summary.csv`). Exact screened reach-ID comparison shows changed samples in 1,420 matched GPP fits and 14,649 matched ER fits. Even among same-sample fits, median exponent changes are +0.00502 GPP and +0.01883 ER; among changed-sample fits, +0.00492 GPP and +0.04657 ER (`fit_changes_by_sample_composition.csv`). This separates value changes from changed fit support descriptively, not causally.

Script 10 and script 11 now save `data/processed/{10,11}_no_uparea_{monthly,annual}_classes.csv` without changing classification rules. Their generated tables have 45,456 monthly rows and 3,788 annual HUC8-process rows each, with **zero key or class mismatches** against each other and against the audit's reconstruction from 08 (`saved_10_11_class_verification.csv`). Under the final monthly rule, 4,448 ER and 681 GPP HUC8-month classes change. Under the annual rule, 696 ER and 185 GPP HUC8 classes change (`monthly_class_pairs.csv`, `annual_class_pairs.csv`, `changed_annual_huc8s.csv`). Annual ER Super-linear falls 1,350 → 810 while Seasonally variable rises 483 → 990; annual GPP Super-linear falls 1,724 → 1,685 while Seasonally variable rises 129 → 158 (`annual_class_counts.csv`). The dominant ER transition is 574 HUC8s from Super-linear to Seasonally variable (`annual_class_transitions.csv`). These are substantive descriptive shifts for scientific review; the existing rule is heuristic and is not a confidence interval.
