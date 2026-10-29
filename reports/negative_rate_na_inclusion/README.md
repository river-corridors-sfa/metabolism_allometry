# HUC8 inclusion if negative local metabolism rates are missing

`scripts/x_26_negative_rates_to_na_huc8_inclusion.R` subtracts every negative
local GPP or ER flux from the corresponding saved upstream sums, leaving
nonnegative local contributions and cumulative watershed area unchanged. It
then refits affected HUC8 × month × process groups with the settings and
eligibility rules in `scripts/x_19_monthly_huc8_segmented_full_reaches.R`.
Unaffected groups retain their current status. The input universe is the
already saved `07_huc08_allometry_input.rds`, which contains only reaches with
positive original cumulative flux; the table therefore measures the effect on
the *current analysis universe* rather than rebuilding every originally
discarded reach from raw input.

| Fit definition | Currently included HUC8s | Currently included HUC8s excluded after NA | Excluded as % of currently included | HUC8s passing after NA | % of 1,894 evaluated HUC8s passing after NA |
|---|---:|---:|---:|---:|---:|
| At least one GPP or ER monthly fit | 1,807 | 0 | 0.000% | 1,807 | 95.41% |
| At least one GPP monthly fit | 1,797 | 0 | 0.000% | 1,797 | 94.88% |
| At least one ER monthly fit | 1,802 | 2 | 0.111% | 1,801 | 95.09% |

The two HUC8s losing their last ER fit are **10020006** and **11030009**.
One HUC8 (**03040102**) gains an ER fit, explaining why the net ER count
decreases by one rather than two. No HUC8 loses all successful fits.

These small inclusion changes are a direct consequence of the **current**
filters: positive cumulative flux, less than 10% missing *watershed area*, at
least 30 reaches and unique areas, at least one decade of area range, and at
least five reaches/three unique areas on each side of the fitted breakpoint.
There is **no threshold on missing metabolism coverage**. A reach with an NA
local rate contributes no local flux to upstream sums, but it can remain a
focal observation if other upstream reaches provide positive cumulative flux.
No observation was lost from the current positive-flux universe; the 156
formerly successful monthly fits that fail after refitting do so because the
new breakpoint leaves too few observations or unique areas on one side.

The alternative policy “exclude any focal reach with an NA local rate” or
“require complete/near-complete upstream metabolism coverage” would produce a
different and potentially much larger exclusion count. Those are additional
rules, not part of the present HUC8 filters.
