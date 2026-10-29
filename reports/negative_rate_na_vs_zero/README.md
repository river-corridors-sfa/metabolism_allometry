# Negative local rates as missing versus zero: HUC8 scaling sensitivity

This one-off reconstructs cumulative GPP and ER from the raw monthly rate
table, retained MERIT catchments, reach dimensions, and full upstream lookup.
It treats every negative local rate as either **missing (NA)** or **zero**.
Both policies route the same nonnegative local flux; they differ only in
whether a formerly negative reach counts as having metabolism coverage.
Both apply the script 07 requirement of less than 10% missing upstream area
and less than 10% missing process-specific metabolism coverage. Watershed
area and routing are identical under both policies.

The missing-policy reconstruction matches every one of the 3,815,080
observations in the regenerated 07 input. The maximum relative cumulative
flux discrepancy is 2.33 × 10⁻¹⁴, attributable to summation order. Its
44,450 supported monthly line fits match the saved 07 slopes and R² to
floating-point precision.

## Monthly allometric slopes

These are HUC8 × month × process `log10(cumulative flux) ~
log10(cumulative upstream area)` slopes using the script 07 requirements:
at least 20 reaches, 20 unique areas, and one decade of area range. The
range is the minimum to maximum across supported monthly fits.

| Process | Negative-rate policy | Eligible reach-months | Supported monthly fits | HUC8s with ≥1 fit | Slope range | Median slope |
|---|---|---:|---:|---:|---:|---:|
| ER | Missing | 1,860,852 | 22,002 | 1,869 | 0.949–2.518 | 1.451 |
| ER | Zero | 1,954,354 | 22,452 | 1,871 | 0.871–2.100 | 1.424 |
| GPP | Missing | 1,954,228 | 22,448 | 1,871 | 0.991–2.239 | 1.475 |
| GPP | Zero | 1,955,927 | 22,452 | 1,871 | 0.991–2.139 | 1.474 |

The zero policy admits 93,502 additional ER and 1,699 additional GPP
reach-month observations. Among HUC8–month groups fitted under both policies,
2,824 ER and 44 GPP slopes differ by more than 0.05. Zero also adds 450
supported monthly ER fits and four GPP fits. The flux values of observations
shared by both policies do not change; these differences arise from which
observations pass the coverage filter.

## Annual HUC8 category assignments

The one-off also repeats script 08's section-deletion stability calculation
for the zero-policy input. It then applies the **same** monthly and annual
classification rules as scripts 10 and 11. The canonical missing-policy 08
results are reused; the zero-policy stability results are saved separately.

| Process | Policy | Insufficient data | Uncertain | Seasonally variable | Super-linear |
|---|---|---:|---:|---:|---:|
| ER | Missing | 26 | 35 | 483 | 1,350 |
| ER | Zero | 24 | 23 | 143 | 1,704 |
| GPP | Missing | 24 | 17 | 129 | 1,724 |
| GPP | Zero | 24 | 18 | 120 | 1,732 |

**375 of 1,894 ER HUC8s (19.8%)** change annual category, compared with
**11 GPP HUC8s (0.58%)**. The largest single transition is 351 ER HUC8s
from *Seasonally variable* under missing rates to *Super-linear* under zero
rates. Two additional HUC8s receive at least one supported monthly ER line
fit under the zero policy; GPP HUC8 fit eligibility is unchanged. The annual
category shifts are primarily a consequence of retaining observations with
originally negative local ER predictions, not a change to cumulative flux
values for observations shared by both policies.

## Queryable outputs

- `reach_month_policy_observations.rds`: each retained catchment × month ×
  process, cumulative area and flux, missing coverage under both policies,
  and `eligible_na`/`eligible_zero` flags. Read with `readRDS()`.
- `na_reconstruction_check.csv`: row-count and numerical agreement with the
  canonical 07 missing-policy input.
- `monthly_huc8_fits.csv`: all HUC8 × month × process fits under both policies,
  including eligibility status, slope, R², and sample size. The policy labels
  are `missing` and `zero`.
- `monthly_fit_comparison.csv`: paired monthly fits with reach-count and slope
  changes (written by the final summarizer).
- `slope_summary.csv`, `paired_slope_change_summary.csv`, and
  `huc8_fit_eligibility.csv`: compact numerical summaries.
- `zero_huc08_allometry_input.rds`: zero-policy input with the same schema
  used by script 08, retained for reproduction.
- `zero_huc08_fit_stability.csv`: zero-policy section-deletion results from
  the same algorithm as script 08.
- `monthly_stability_by_policy.csv`, `monthly_category_by_policy.csv`,
  `annual_huc8_categories_by_policy.csv`, `huc8_assignment_comparison.csv`,
  `huc8_category_transitions.csv`, and `huc8_category_totals.csv`: the exact
  Figure 1–2 category assignments and transitions under each policy.

Reproduce in order from the repository root:

```sh
Rscript scripts/x_28_negative_rate_na_vs_zero_scaling.R
Rscript scripts/x_28b_zero_huc8_fit_stability.R
Rscript scripts/x_28c_summarize_na_vs_zero_assignments.R
```

The stability run is computationally the expensive step. It saves to this
report directory and never replaces the canonical 07 or 08 outputs.
When reading CSVs, import `huc` as character to preserve leading zeros,
for example `readr::read_csv(path, col_types = readr::cols(huc =
readr::col_character()))`.

This report compares the HUC8 **single-line** allometric slopes used by
07/08 and the associated Figure 1–2 category assignments. It does not
compare the pooled two-segment fits from script 12.
