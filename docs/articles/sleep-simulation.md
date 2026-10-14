# Sleep example: actigraphy SQI read on the PSQI components

This vignette walks through the worked example of the manuscript: a
simulated data set (n = 1000) in which seven PSQI component scores (C1
to C7, 0 to 3) and an actigraphy sleep-quality index (SQI, higher =
better sleep) are driven by one latent sleep-disturbance trait. The
generating SQI curve is the quadratic 72 − 8η − 1.5η² with residual SD
6; it is used only as the truth.

The light steps below run when the vignette is built. The bootstrap and
the simulation study take tens of minutes to hours and are shown, not
run; the full analysis scripts are in
`system.file("scripts", package = "caliber")`.

## Data

The data ship with the package as the data set `psqi_sleep` (see
[`?caliber::psqi_sleep`](https://paoloeusebi.github.io/caliber/reference/psqi_sleep.md)).
[`sleep_data()`](https://paoloeusebi.github.io/caliber/reference/sleep_data.md)
returns it without the true trait, which
[`sleep_eta_true()`](https://paoloeusebi.github.io/caliber/reference/sleep_data.md)
keeps apart for evaluation.

[`library`](https://rdrr.io/r/base/library.html)`(``caliber``)`` ``dat`` ``<-`` `[`sleep_data`](https://paoloeusebi.github.io/caliber/reference/sleep_data.md)`(``)`` ``items`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"C"``, ``1``:``7``)`` ``Y`` ``<-`` `[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(``dat``[``items``]``)``; ``sqi`` ``<-`` ``dat``$``SQI`` `[`summary`](https://rdrr.io/r/base/summary.html)`(``dat``$``PSQI_global``)``; `[`cor`](https://rdrr.io/r/stats/cor.html)`(``dat``$``PSQI_global``, ``sqi``)`` ``#> Min. 1st Qu. Median Mean 3rd Qu. Max. `` ``#> 0.000 4.000 6.000 6.622 9.000 21.000`` ``#> [1] -0.702831`

## Step 1: the GRM on the PSQI alone

`grm`` ``<-`` `[`fit_grm`](https://paoloeusebi.github.io/caliber/reference/fit_grm.md)`(``Y``)`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`cbind`](https://rdrr.io/r/base/cbind.html)`(``a ``=`` ``grm``$``a``, b ``=`` ``-``grm``$``d`` ``/`` ``grm``$``a``)``, ``2``)`` ``#> a `` ``#> [1,] 2.31 -0.84 0.98 2.07`` ``#> [2,] 1.28 -1.00 0.44 1.59`` ``#> [3,] 1.35 0.00 1.07 2.19`` ``#> [4,] 1.92 0.38 1.11 1.80`` ``#> [5,] 1.23 -2.33 0.88 3.16`` ``#> [6,] 0.83 1.44 2.41 2.96`` ``#> [7,] 1.44 -0.95 0.99 2.77`

## Step 2: the monotone SQI curve, MS-Fixed

MS-Fixed holds the items at the PSQI-only estimates and integrates the
latent trait out over p(η \| C, SQI) at every boosting iteration. It is
the recommended primary analysis; MS-Joint (`mode = "joint"`) is the
robustness check and MS-EAP (`mode = "eap"`) the naive comparator.
[`fit_three()`](https://paoloeusebi.github.io/caliber/reference/fit_three.md)
fits all three.

`fit`` ``<-`` `[`fit_ms`](https://paoloeusebi.github.io/caliber/reference/fit_ms.md)`(``grm``, ``Y``, ``sqi``, mode ``=`` ``"fixed"``)`` `[`c`](https://rdrr.io/r/base/c.html)`(``iter ``=`` ``fit``$``iter``, df ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``fit``$``df``, ``2``)``, sigma0 ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``fit``$``act``$``sigma``, ``2``)``,`` `` g1 ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(``fit``$``act``$``log_sigma``[``2``]``, ``3``)``)`` ``#> iter df sigma0 g1 `` ``#> 1046.00 4.46 5.74 -0.06`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`mu_sqi`](https://paoloeusebi.github.io/caliber/reference/monospline_act.md)`(``-``3``:``3``, ``fit``$``act``)``, ``1``)`` ``#> [1] 79.5 79.5 78.3 72.7 62.5 48.7 28.7`

## Step 3: interpretation in SQI units

Each GRM threshold b_jk is placed on the SQI scale at s_jk = μ(b_jk).

`thr`` ``<-`` `[`threshold_locations`](https://paoloeusebi.github.io/caliber/reference/threshold_locations.md)`(``fit``$``a``, ``fit``$``d``, ``fit``$``act``)`` `[`dimnames`](https://rdrr.io/r/base/dimnames.html)`(``thr``$``sqi``)`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``items``, `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"k"``, ``1``:``3``)``)`` `[`round`](https://rdrr.io/r/base/Round.html)`(``thr``$``sqi``, ``1``)`` ``#> k1 k2 k3`` ``#> C1 77.6 62.7 48.0`` ``#> C2 78.3 67.7 54.8`` ``#> C3 72.7 61.8 47.0`` ``#> C4 68.5 61.3 51.6`` ``#> C5 79.5 63.6 18.9`` ``#> C6 57.0 45.6 31.3`` ``#> C7 78.1 62.6 40.2`

The probabilistic crosswalk gives the expected PSQI response pattern at
a given SQI, averaging over the latent values compatible with it.

`cw`` ``<-`` `[`crosswalk`](https://paoloeusebi.github.io/caliber/reference/crosswalk_probabilistic.md)`(`[`c`](https://rdrr.io/r/base/c.html)`(``50``, ``60``, ``70``, ``80``)``, ``fit``$``a``, ``fit``$``d``, ``fit``$``act``, ``items``, total_name ``=`` ``"PSQI"``)`` `[`round`](https://rdrr.io/r/base/Round.html)`(``cw``[`[`c`](https://rdrr.io/r/base/c.html)`(``"SQI"``, ``"E_PSQI"``, ``"Estar_PSQI"``)``]``, ``2``)`` ``#> SQI E_PSQI Estar_PSQI`` ``#> 1 50 13.30 14.39`` ``#> 2 60 9.69 11.30`` ``#> 3 70 6.24 7.04`` ``#> 4 80 4.01 NA`

`th`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``item ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``items``, ``3``)``, k ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``1``:``3``, each ``=`` ``7``)``, sqi_est ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``thr``$``sqi``)``)`` ``pv`` ``<-`` `[`plausible_values`](https://paoloeusebi.github.io/caliber/reference/plausible_values.md)`(``fit``, ``grm``, ``sqi``, seed ``=`` ``20261004``)`` `[`wright_map_sqi`](https://paoloeusebi.github.io/caliber/reference/wright_map_sqi.md)`(`[`mu_sqi`](https://paoloeusebi.github.io/caliber/reference/monospline_act.md)`(``pv``, ``fit``$``act``)``, ``th``, ``sleep_item_labels``)`

![Wright map on the SQI scale,
MS-Fixed.](sleep-simulation_files/figure-html/wright-1.png)

Wright map on the SQI scale, MS-Fixed.

## Intervals

Item parameters, thresholds within the data and the crosswalk get
intervals from draws of the normal approximation of the estimates:

`vc`` ``<-`` `[`ms_vcov`](https://paoloeusebi.github.io/caliber/reference/ms_vcov.md)`(``fit``, ``Y``, ``sqi``)`` ``iv`` ``<-`` `[`ms_intervals`](https://paoloeusebi.github.io/caliber/reference/ms_intervals.md)`(``fit``, ``vc``, ``function``(``a``, ``d``, ``act``)`` ``{`` `` ``s`` ``<-`` `[`threshold_locations`](https://paoloeusebi.github.io/caliber/reference/threshold_locations.md)`(``a``, ``d``, ``act``)``$``sqi`` `` `[`setNames`](https://rdrr.io/r/stats/setNames.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``s``)``, `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"s_"``, `[`threshold_names`](https://paoloeusebi.github.io/caliber/reference/threshold_names.md)`(``items``)``)``)`` ``}``, n_draw ``=`` ``2000``, seed ``=`` ``1``)`

These draws condition on the selected boosting path and miss most of the
sampling variance of the curve itself (SD of μ̂(0) 0.26 against a true
0.61), so curve bands come from the nonparametric bootstrap, which
refits every step:

`boot`` ``<-`` `[`ms_bootstrap`](https://paoloeusebi.github.io/caliber/reference/ms_bootstrap.md)`(``Y``, ``sqi``, n_boot ``=`` ``200``, keep_categories ``=`` ``FALSE``, cores ``=`` ``4``)`` `[`boot_intervals`](https://paoloeusebi.github.io/caliber/reference/boot_intervals.md)`(``boot``, ``"MS-Fixed"``, `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"mu("``, ``-``3``:``3``, ``")"``)``)`

Neither kind of interval has been validated by a coverage simulation.
The top thresholds of C5, C6 and C7 lie beyond the data and should be
reported as beyond data support, not with an interval.

## Simulation study

`inst/scripts/sleep_simulation.R` runs six scenarios × 200 replications
(about 3 hours on 4 cores) and `sleep_simulation_summary.R` makes the
tables and figures. One replication:

`d`` ``<-`` `[`simulate_scenario`](https://paoloeusebi.github.io/caliber/reference/simulate_scenario.md)`(``s ``=`` ``4``, r ``=`` ``1``)`` ``# gentle logistic truth`` ``f`` ``<-`` `[`fit_three`](https://paoloeusebi.github.io/caliber/reference/fit_three.md)`(``d``$``Y``, ``d``$``sqi``)`` ``tr`` ``<-`` `[`scenario_truth`](https://paoloeusebi.github.io/caliber/reference/scenario_truth.md)`(``4``, `[`seq`](https://rdrr.io/r/base/seq.html)`(``15``, ``95``, by ``=`` ``0.5``)``)`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(`[`mean`](https://rdrr.io/r/base/mean.html)`(``(`[`threshold_locations`](https://paoloeusebi.github.io/caliber/reference/threshold_locations.md)`(``f``$``` `MS-Fixed` ```$``a``, ``f``$``` `MS-Fixed` ```$``d``, ``f``$``` `MS-Fixed` ```$``act``)``$``sqi`` ``-`` ``tr``$``thr``)``^``2``)``)`
