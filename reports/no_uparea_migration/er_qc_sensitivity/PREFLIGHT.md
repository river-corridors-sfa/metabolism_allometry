# ER training-label sensitivity: preflight (2026-10-30)

**Status:** APPROVED AND COMPLETED on 2026-10-30. Four ER forests were fitted and the QC variants' CONUS predictions generated with `.venv-maavara/bin/python scripts/validation/er_qc_sensitivity.py`. Results and figures are in `ER_QC_comparison.Rmd` and its rendered HTML. Production files were not changed.

## Fixed comparison design

The source is `data/raw/dryad_data/appling_data_06202024.csv`. Reproduce the existing row-wise 70/30 `train_test_split` with `random_state=24` before applying any label treatment. ER magnitude is `-ER` in source units. Evaluate every locally fitted forest on the same 134,336 held-out rows with raw `ER <= 0`. Do not alter or filter the held-out set in model-specific ways. Retain the existing RF settings: 500 trees, `max_features=5`, `max_depth=30`, `max_samples=0.4`, `random_state=24`, `oob_score=True`; use the pinned Maavara package environment.

| Fit | Predictors | Training treatment | Training rows | Training sites |
|---|---|---|---:|---:|
| Reproduced seven-predictor baseline | Existing seven, including `uparea` | Existing unfiltered `-ER` | 343,250 | 355 |
| Unfiltered six-predictor baseline | Existing six, excluding `uparea` | Existing unfiltered `-ER` | 343,250 | 355 |
| A: strict exclusion | Same six | Keep raw `ER <= 0`; target `-ER` | 313,012 | 353 |
| B: near-zero treatment | Same six | Keep raw `ER <= +0.5`; set `(0,+0.5]` target to zero, otherwise `-ER` | 319,263 | 354 |

The full source has 43,010 raw-positive ER rows. The original seeded split places 30,238 in training (23,987 above +0.5 and 6,251 in `(0,+0.5]`) and 12,772 in testing. Of the latter, 134,336 rows meet the common physically valid held-out rule. Thus A and B exclude **training** rows only; the original test membership is preserved.

The Maavara published seven-predictor model is a separate benchmark. The archived Dryad `random_forest_CONUS_ER.ipynb` final seven-feature cell reports test Pearson `r = 0.8823235` and `r² = 0.7784948` on its original unfiltered test set; it does not report RMSE, MAE, or bias in that cell. This is Pearson correlation squared, **not** `sklearn.metrics.r2_score`. The paper's separately stated performance, if any, must be quoted with its definition and protocol. Its saved CONUS predictions support a paired network-distribution comparison, but cannot be scored on the local valid-only held-out rows unless individual held-out predictions or a fitted model are available. Do not present its reported r² as a paired test result.

## Why four fits are required

`data/processed/17_maavara_refits/model_metrics.csv` contains aggregate metrics for the two existing local ER baselines, and the 12 monthly CONUS CSVs contain their network predictions. The run did not save the fitted forests or row-level held-out predictions. Those metrics include the 12,772 raw-positive held-out ER responses and therefore cannot be recalculated on the 134,336-row valid-only set. A fair common-test comparison requires refitting the two baselines plus A and B. The existing baseline files remain immutable; all new work belongs under this report directory or a separate scratch output directory.

## Required outputs after approval

1. Save held-out predictions for all four local forests keyed to original source-row index, and aggregate R², RMSE, MAE, bias, and negative/zero prediction counts on the common valid set. Report uncertainty or meaningful paired differences without using a post hoc acceptance threshold.
2. Stratify the same held-out rows by fixed, pooled stream-size and discharge bins, month/season, and geographic region. Report bin counts, observed response support, R² when defined, RMSE, MAE, bias, and prediction sign rates. Show the paired skill change of A and B relative to the unfiltered six-predictor forest.
3. Record retained sites, training-row counts, and predictor ranges/quantiles for each fit, with full-source/held-out comparisons. Note site overlap between training and testing caused by the original row-wise split; this is **not** a site-held-out test.
4. Predict all 12 months on the unchanged CONUS scaling dataset using the existing feature construction. Keep five incomplete covariate rows per month as NA. Compare the four local forests and the downloaded Maavara predictions on identical COMID-month keys: distributions, negative/zero rates, and usable coverage (`finite & >0`) overall and by upstream-area bins, including headwaters. Count gains/losses per COMID-month, not just totals.
5. Save a concise R Markdown report with figures, source links, exact run metadata, and a conditional recommendation based on skill, physical plausibility, and coverage. Do not change scripts 01–12 or any production output.

## Current evidence and limits

The existing unfiltered local ER test metrics on **all** 147,108 test rows are R² 0.777697, RMSE 2.605172, MAE 1.507041, bias 0.006479 for the seven-predictor fit and R² 0.775994, RMSE 2.615133, MAE 1.514744, bias 0.009041 for the six-predictor fit. These are context only, not common-valid-test benchmarks. The [wrong-sign training audit](../step1_er_training_distribution/ER_wrong_sign_training_audit.Rmd) shows positive raw ER across sites, flow, stream sizes, seasons, and regions; QC changes training support. The final scientific decision awaits the comparison.
