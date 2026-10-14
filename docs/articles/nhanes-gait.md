# NHANES example: gait speed read on seven mobility items

The real-data example of the manuscript: gait speed from the timed 20-ft
walk (m/s) in NHANES 2001–2002 participants aged 60+ (n = 1,426), read
on seven lower-body mobility items of the physical functioning
questionnaire (0 = no difficulty, 1 = some, 2 = much, 3 = unable). The
latent trait is mobility difficulty; gait speed falls as it rises, so
the non-increasing curve applies as it is. The analysis data ship with
the package as the data set `pf_gait` (see
[`?caliber::pf_gait`](https://paoloeusebi.github.io/caliber/reference/pf_gait.md)),
which
[`gait_data()`](https://paoloeusebi.github.io/caliber/reference/gait_data.md)
returns together with the item matrix; they are built from the
public-use files by `inst/scripts/gait_data.R`.

The full analysis (three methods, 200 bootstrap resamples, sensitivity
variants, figures) is `inst/scripts/gait_main.R`; item dimensionality
checks are `inst/scripts/gait_dimensionality.R` (needs lavaan).

## Data

[`library`](https://rdrr.io/r/base/library.html)`(``caliber``)`` ``g`` ``<-`` `[`gait_data`](https://paoloeusebi.github.io/caliber/reference/gait_data.md)`(``)`` ``items`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"C"``, ``1``:``7``)`` `[`c`](https://rdrr.io/r/base/c.html)`(``n ``=`` ``g``$``n``, mean_speed ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`mean`](https://rdrr.io/r/base/mean.html)`(``g``$``speed``)``, ``2``)``, sd_speed ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`sd`](https://rdrr.io/r/stats/sd.html)`(``g``$``speed``)``, ``2``)``)`` ``#> n mean_speed sd_speed `` ``#> 1426.00 0.92 0.26`` `[`colSums`](https://rdrr.io/r/base/colSums.html)`(`[`is.na`](https://rdrr.io/r/base/NA.html)`(``g``$``Y``)``)`` ``# C1, C2 not asked of people who need equipment to walk`` ``#> C1 C2 C3 C4 C5 C6 C7 `` ``#> 105 109 0 2 0 1 14`

Missing item responses are allowed: they contribute nothing to a
person’s likelihood. `gait_data(exclude_aids = TRUE)` drops the 161
walking-aid users, and `gait_data(skip_unable = TRUE)` scores the
skipped C1 and C2 as unable.

## MS-Fixed

`grm`` ``<-`` `[`fit_grm`](https://paoloeusebi.github.io/caliber/reference/fit_grm.md)`(``g``$``Y``)`` ``fit`` ``<-`` `[`fit_ms`](https://paoloeusebi.github.io/caliber/reference/fit_ms.md)`(``grm``, ``g``$``Y``, ``g``$``speed``, mode ``=`` ``"fixed"``)`` `[`curve_ends`](https://paoloeusebi.github.io/caliber/reference/curve_ends.md)`(``fit``$``act``, ``g``$``speed``)`` ``#> floor ceiling share_below_floor share_above_ceiling r2_approx`` ``#> 1 0.4866731 1.044125 0.06451613 0.326087 0.3296524`

The curve resolves only gait speeds between its floor (about 0.49 m/s)
and ceiling (about 1.04 m/s); about a third of observed speeds are above
the ceiling.

## Thresholds in gait-speed units

`thr`` ``<-`` `[`threshold_locations`](https://paoloeusebi.github.io/caliber/reference/threshold_locations.md)`(``fit``$``a``, ``fit``$``d``, ``fit``$``act``)`` ``flags`` ``<-`` `[`threshold_flags`](https://paoloeusebi.github.io/caliber/reference/threshold_flags.md)`(``g``$``Y``, ``thr``$``sqi``, ``fit``$``act``)`` ``tab`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``item ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``gait_item_labels``, ``3``)``, k ``=`` `[`rep`](https://rdrr.io/r/base/rep.html)`(``1``:``3``, each ``=`` ``7``)``,`` `` speed ``=`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``thr``$``sqi``)``, ``2``)``, ``flags``[`[`c`](https://rdrr.io/r/base/c.html)`(``"n_at_or_above"``, ``"flag_sparse"``, ``"at_floor"``)``]``)`` ``tab``$``report`` ``<-`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``tab``$``at_floor``, `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"below about %.2f m/s"``, `[`mu_sqi`](https://paoloeusebi.github.io/caliber/reference/monospline_act.md)`(``4``, ``fit``$``act``)``)``,`` `` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"%.2f"``, ``tab``$``speed``)``)`` ``tab``[`[`order`](https://rdrr.io/r/base/order.html)`(``tab``$``item``, ``tab``$``k``)``, `[`c`](https://rdrr.io/r/base/c.html)`(``"item"``, ``"k"``, ``"n_at_or_above"``, ``"report"``, ``"flag_sparse"``)``]`` ``#> item k n_at_or_above report`` ``#> 1 C1 Walking a quarter mile 1 377 0.88`` ``#> 8 C1 Walking a quarter mile 2 159 0.71`` ``#> 15 C1 Walking a quarter mile 3 82 0.61`` ``#> 2 C2 Walking up 10 steps 1 289 0.82`` ``#> 9 C2 Walking up 10 steps 2 106 0.66`` ``#> 16 C2 Walking up 10 steps 3 38 0.52`` ``#> 3 C3 Walking between rooms 1 116 0.61`` ``#> 10 C3 Walking between rooms 2 23 below about 0.49 m/s`` ``#> 17 C3 Walking between rooms 3 7 below about 0.49 m/s`` ``#> 4 C4 Standing up from armless chair 1 338 0.79`` ``#> 11 C4 Standing up from armless chair 2 106 0.58`` ``#> 18 C4 Standing up from armless chair 3 38 below about 0.49 m/s`` ``#> 5 C5 Getting in and out of bed 1 217 0.68`` ``#> 12 C5 Getting in and out of bed 2 35 below about 0.49 m/s`` ``#> 19 C5 Getting in and out of bed 3 2 below about 0.49 m/s`` ``#> 6 C6 Stooping, crouching, kneeling 1 704 1.00`` ``#> 13 C6 Stooping, crouching, kneeling 2 295 0.75`` ``#> 20 C6 Stooping, crouching, kneeling 3 134 0.59`` ``#> 7 C7 Standing for long periods 1 575 0.95`` ``#> 14 C7 Standing for long periods 2 295 0.75`` ``#> 21 C7 Standing for long periods 3 168 0.65`` ``#> flag_sparse`` ``#> 1 FALSE`` ``#> 8 FALSE`` ``#> 15 FALSE`` ``#> 2 FALSE`` ``#> 9 FALSE`` ``#> 16 FALSE`` ``#> 3 FALSE`` ``#> 10 TRUE`` ``#> 17 TRUE`` ``#> 4 FALSE`` ``#> 11 FALSE`` ``#> 18 FALSE`` ``#> 5 FALSE`` ``#> 12 FALSE`` ``#> 19 TRUE`` ``#> 6 FALSE`` ``#> 13 FALSE`` ``#> 20 FALSE`` ``#> 7 FALSE`` ``#> 14 FALSE`` ``#> 21 FALSE`

Thresholds within 0.01 m/s of the curve’s floor cannot be told apart and
are reported as “below about 0.5 m/s”.

## Crosswalk and exceedance probabilities

`cw`` ``<-`` `[`crosswalk_probabilistic`](https://paoloeusebi.github.io/caliber/reference/crosswalk_probabilistic.md)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0.4``, ``0.6``, ``0.8``, ``1.0``, ``1.2``)``, ``fit``$``a``, ``fit``$``d``, ``fit``$``act``, ``items``)`` `[`round`](https://rdrr.io/r/base/Round.html)`(`[`cbind`](https://rdrr.io/r/base/cbind.html)`(``speed ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0.4``, ``0.6``, ``0.8``, ``1.0``, ``1.2``)``, ``cw``)``, ``2``)`` ``#> speed C1 C2 C3 C4 C5 C6 C7 total`` ``#> 1 0.4 1.85 1.54 0.59 1.29 0.67 1.96 2.03 9.92`` ``#> 2 0.6 1.20 0.94 0.27 0.76 0.39 1.40 1.40 6.36`` ``#> 3 0.8 0.62 0.45 0.09 0.35 0.18 0.87 0.81 3.37`` ``#> 4 1.0 0.31 0.20 0.02 0.15 0.09 0.55 0.46 1.77`` ``#> 5 1.2 0.17 0.10 0.01 0.07 0.05 0.39 0.29 1.08`` ``ex`` ``<-`` `[`exceedance_probabilistic`](https://paoloeusebi.github.io/caliber/reference/exceedance_probabilistic.md)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0.6``, ``1.0``)``, ``fit``$``a``, ``fit``$``d``, ``fit``$``act``, ``items``)`` `[`round`](https://rdrr.io/r/base/Round.html)`(``ex``[``, `[`paste0`](https://rdrr.io/r/base/paste.html)`(``items``, ``"_1"``)``]``, ``2``)`` ``# P(at least some difficulty)`` ``#> C1_1 C2_1 C3_1 C4_1 C5_1 C6_1 C7_1`` ``#> [1,] 0.62 0.54 0.20 0.50 0.32 0.76 0.69`` ``#> [2,] 0.21 0.15 0.02 0.13 0.08 0.40 0.30`

## Bootstrap intervals

`boot`` ``<-`` `[`ms_bootstrap`](https://paoloeusebi.github.io/caliber/reference/ms_bootstrap.md)`(``g``$``Y``, ``g``$``speed``, n_boot ``=`` ``200``, seed ``=`` ``20261010``,`` `` keep_categories ``=`` ``TRUE``, cores ``=`` ``8``)`` `[`boot_intervals`](https://paoloeusebi.github.io/caliber/reference/boot_intervals.md)`(``boot``, ``"MS-Fixed"``, `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"s_"``, `[`threshold_names`](https://paoloeusebi.github.io/caliber/reference/threshold_names.md)`(``items``)``)``)`

A resample in which an item lost a response category is redrawn (the GRM
then has no threshold for it), so the intervals are conditional on every
category being present. In the manuscript run 29 of 200 resamples needed
a redraw, almost all because of the two “unable” responses for getting
in and out of bed.
