# Step 1 sensitivity: interpolate negative no-uparea predictions

Status: **SCENARIO TESTED; scientific choice unresolved** (2026-10-30). Code: `scripts/validation/step1_neighbor_imputation_sensitivity.R`. Production source data, 07 outputs, thresholds, and fits were not changed.

For each month and process, a negative no-uparea areal rate is eligible only when it has a direct downstream reach and one or more direct upstream reaches in saved `05_edges_downstream.rds`, and **every** direct neighbor's original no-uparea rate is positive. At confluences, the direct upstream rates are averaged into one upstream mean; the replacement is `(upstream mean + downstream rate) / 2`. This is a single simultaneous pass: imputed values cannot serve as neighbors for other imputations. Original missing values are not imputed. The calculation averages areal rates, not fluxes, and uses equal weight for each direct upstream reach.

The script replayed the existing 07 cumulative flux/coverage calculation using saved `05_upstream_lookup.rds` and applied the unchanged 07 screens. Every saved no-uparea screened key was retained in the scenario; the additional keys are in `recovered_screened_rows.csv`.

| Process | Negative no-uparea source reach-month rates | Eligible for interpolation | Saved no-uparea screened rows | Scenario screened rows | Additional rows |
|---|---:|---:|---:|---:|---:|
| GPP | 8,850 | 3,437 | 1,951,894 | 1,952,478 | 584 |
| ER | 152,517 | 29,832 | 1,659,023 | 1,698,011 | 38,988 |
| Total | 161,367 | 33,269 | 3,610,917 | 3,650,489 | 39,572 |

The scenario closes **39,572 / 204,163 = 19.38%** of the original-to-no-uparea *net* screened-row gap. Intersecting recovered keys with the earlier `changed_screened_rows.csv` shows that 29,549 (29,116 ER, 433 GPP) were actually original-only losses; all had failed the metabolism-coverage screen, not the positive-flux screen. The other 10,023 additional rows (9,872 ER, 151 GPP) were absent from both original and no-uparea screens. Thus the scenario directly recovers **29,549 / 236,032 = 12.52%** of original-only lost keys. The scenario total remains 164,591 rows below the original 07 total.

This is only a support/count sensitivity. It does not test whether interpolation is scientifically justified, how altered cumulative values change slopes, R², fit eligibility, or final classes. Using adjacent modeled rates can smooth spatial variation and add dependence to already nested cumulative observations. Replacing negative predictions also changes the treatment of possible model extrapolation or genuinely low process values. The equal-weight confluence rule is one possible definition, not an established hydrologic estimator. A research-team decision is required before any production method change.
