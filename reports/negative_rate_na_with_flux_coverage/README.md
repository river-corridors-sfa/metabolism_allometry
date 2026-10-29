# HUC8 inclusion with negative rates missing and <10% missing metabolism coverage

This sensitivity run applies two rules together: negative local GPP and ER
rates become NA before cumulative flux is calculated, and each focal
reach-month must have **less than 10%** missing metabolism coverage in its
cumulative upstream area. The cutoff is the same strict `<10%` comparison
already used for missing watershed area. Coverage is metric-specific: a
negative GPP rate affects GPP coverage, and a negative ER rate affects ER
coverage. The missing reach's upstream *area* remains in the denominator;
its local metabolism flux is omitted from the numerator.

| Fit definition | Currently included HUC8s | Currently included HUC8s excluded | Excluded as % of currently included | HUC8s passing new rules | % of 1,894 evaluated HUC8s passing |
|---|---:|---:|---:|---:|---:|
| At least one GPP or ER monthly fit | 1,807 | 0 | 0.000% | 1,807 | 95.41% |
| At least one GPP monthly fit | 1,797 | 0 | 0.000% | 1,798 | 94.93% |
| At least one ER monthly fit | 1,802 | 9 | 0.499% | 1,794 | 94.72% |

The nine HUC8s losing their last ER fit are `10020006`, `11030009`,
`12090206`, `13040211`, `13070005`, `17070307`, `17110004`, `18040007`,
and `18060008`. One other HUC8 (`03040102`) gains an ER fit, so net ER
inclusion falls by eight. GPP gains one HUC8 (`12100407`) and loses none.
These gains can occur despite fewer observations because changing inputs
and refitting can move a breakpoint so that both sides newly satisfy the
minimum side-support rule. No HUC8 loses all successful fits.

The coverage filter removes 77,548 of 3,892,374 focal reach-month-process
rows in the saved positive-cumulative-flux input. These rows include some
with missing metabolism coverage before negative rates are recoded. This
analysis uses the same existing saved input universe as the prior NA-only
comparison; it does not reconstruct originally discarded nonpositive
cumulative-flux observations. The HUC8 counts therefore address the effect
on currently evaluated HUC8s.

Run `Rscript scripts/x_27_negative_rates_to_na_with_flux_coverage.R` from
the project root to reproduce `summary.csv`, `huc8_inclusion_comparison.csv`,
and `group_status_comparison.csv`.
