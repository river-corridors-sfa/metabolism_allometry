# Appling observed-network screen — local-data ranking

This repository contains Maavara et al.'s `appling_data_06202024.csv` (342 MB
on disk), not the three authoritative Appling release inputs needed for the
requested final answer: **1. Site data** (including the site-level confidence
flag), **2. Spatial data**, and **8. Metabolism estimates and predictors**.
The local file has 490,358 unique site-date rows for 355 sites and retains the
Appling daily GPP, ER, credible intervals, R-hat, and n_eff fields. Dryad's
README explicitly says the Appling variables retain their original names and
units, but that statement does not replace a file-level checksum or join to the
USGS release.

I attempted to fetch Maavara's two published notebooks from their public Dryad
file-stream URLs on 2026-09-23; Dryad returned HTTP 403. Thus the notebook's
actual filtering cannot be asserted. The analysis below consequently uses an
explicit daily rule, not an inferred notebook rule or absent site-confidence
label.

## What was run

`scripts/19_appling_network_prescreen.R` is reproducible tidyverse R code and
writes compact audit CSVs under `data/processed/appling_network_screen/`. It:

1. validates duplicate site-dates, site-to-MERIT joins, date parsing, and the
   sign distribution of ER;
2. defines a *provisional daily screen* as finite GPP and ER with both daily
   R-hat values <= 1.1; its separate log screen also requires GPP > 0 and ER <
   0, and uses `-ER` only as a positive respiration magnitude;
3. identifies every simultaneous subset within a directed local MERIT
   `nxt_dw` component, retaining total common days, longest uninterrupted run,
   drainage-area range, and ratio; and
4. writes directed upstream/downstream relationships from that local graph.

The threshold is transparent but not a substitute for Appling's published
site-level confidence labels. The Appling paper describes those labels as
depending on model-level convergence, K600 variation, negative GPP, and
positive ER; it reports ER with a negative sign convention. It also explicitly
warns that output users should filter according to their target analysis.

## Local-data candidates (not NHDPlus-verified winners)

The script ranks distinct components by simultaneous log-usable gauge count,
then longest uninterrupted overlap and drainage-area ratio. It writes the
complete result in `top_three_merit_networks.csv` and the included figure uses
blue only for paths traced through local `nxt_dw` routing—not straight-line
schematic links—and orange for the simultaneous gauge set.

`figures/19_appling_candidate_huc8_stream_order.png` adds HUC8 context. Within
each selected HUC8, local MERIT flowlines are drawn with increasing opacity and
width for increasing Strahler order; traced candidate paths are blue and gauges
are orange.

1. Component 333: 75 gauges, one shared usable/loggable day (2016-07-31),
   area ratio 113,290. This is a maximum-count but impractically brief,
   continent-scale branched candidate, not a good first rate~area test.
2. Component 256: 14 gauges, 25 shared days, longest run 2 days, area ratio
   479. It is also brief and branched; it is useful only after directed
   NHDPlus verification and a stable-window subset search.
3. Component 259: 12 gauges, 77 shared days, longest run 4 days, area ratio
   18.5. It is the most credible of the three for an initial local-data
   exploration, but still has a very short continuous overlap.

These are **MERIT-routing candidates, not the answer to the requested
“verified river networks” question**. A MERIT COMID is not an NHDPlus COMID;
the local `nxt_dw` graph supports reproducible directed-path screening but
does not supply the requested independent NHDPlus/NLDI or authoritative
upstream-polygon confirmation.

The required missing inputs/services are:

* USGS item 1 site table and item 7 diagnostics (to retain only eligible
  site-model records and their published confidence labels);
* USGS item 2a points and 2b upstream catchment polygons (to reconcile gauge
  identity and independently check nesting); and
* NHDPlus routing or NLDI navigation for every candidate gauge (to verify the
  directed paths and identify tributary versus strictly nested structures).

After they are available, replace the provisional screen with the release
confidence filter and verify every CSV relationship against NHDPlus/NLDI.
That will determine whether the branched components or smaller strictly nested
sequences are the appropriate first empirical tests.

## Sources

* [Maavara et al. Dryad dataset](https://doi.org/10.5061/dryad.c2fqz61kx)
* [USGS Appling release](https://doi.org/10.5066/F70864KX)
* [Appling et al. 2018 data descriptor](https://doi.org/10.1038/sdata.2018.292)

## Scientific limitation

GPP and ER here are reach-averaged areal rates (g O2 m-2 d-1). Even a verified
`log(rate) ~ log(drainage area)` comparison is not a direct test of total or
cumulative network-metabolism allometry, and repeated days never add spatial
gauges.
