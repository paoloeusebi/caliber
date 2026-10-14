# caliber

Monotone calibration of a digital or performance measure to a multi-item
PRO. One latent trait drives the PRO items (graded response model, items
scored 0 to 3) and the digital measure (non-increasing cubic B-spline
mean, log-linear residual SD). The curve is fitted by monotonic
likelihood boosting (Leitenstorfer and Tutz 2007) in three ways:
MS-Fixed (fixed-item calibration, the recommended primary analysis),
MS-Joint (concurrent calibration) and MS-EAP (EAP plug-in, the naive
comparator).

The curve can instead be fitted with Ramsay’s (1988) monotone regression
splines (`fit_ms(..., curve = "ispline")`): a non-negative combination
of I-splines on a few interior knots at quantiles of the latent N(0, 1),
fitted by maximum likelihood, with the number of knots chosen by AIC.
Boosting stays the default.

The digital measure can be continuous (`family = "gaussian"`, the
default: actigraphy indices, gait speed, distance) or a count
(`family = "poisson"` or `"negbin"`: daily steps, events, deaths). For
counts the spline models log E\[y \| η\], with either curve method; θ of
the negative binomial is estimated by ML. Use Poisson for small counts
(events, deaths) and the negative binomial for large, overdispersed ones
(steps): a Poisson count in the thousands locates the trait more finely
than the quadrature resolves, and
[`fit_ms()`](https://paoloeusebi.github.io/caliber/reference/fit_ms.md)
warns. Normal-approximation intervals
([`ms_vcov()`](https://paoloeusebi.github.io/caliber/reference/ms_vcov.md))
are Gaussian-only; use
[`ms_bootstrap()`](https://paoloeusebi.github.io/caliber/reference/ms_bootstrap.md)
for counts.

## Install (from this folder)

[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``"mirt"``, ``"statmod"``, ``"ggplot2"``, ``"patchwork"``, ``"scales"``, ``"viridisLite"``)``)`` ``roxygen2``::`[`roxygenise`](https://roxygen2.r-lib.org/reference/roxygenize.html)`(``"."``)`` ``# builds man/ (and checks NAMESPACE)`` ``devtools``::`[`install`](https://devtools.r-lib.org/reference/install.html)`(``"."``)`` ``# or: R CMD build . && R CMD INSTALL caliber_0.1.0.tar.gz`

## Main functions

| Step | Functions |
|----|----|
| GRM on the PRO | [`fit_grm()`](https://paoloeusebi.github.io/caliber/reference/fit_grm.md) |
| Curve | [`fit_ms()`](https://paoloeusebi.github.io/caliber/reference/fit_ms.md) (`curve = "boost"` or `"ispline"`), [`fit_three()`](https://paoloeusebi.github.io/caliber/reference/fit_three.md), [`mu_sqi()`](https://paoloeusebi.github.io/caliber/reference/monospline_act.md), [`sigma_sqi()`](https://paoloeusebi.github.io/caliber/reference/monospline_act.md), [`ispline_basis()`](https://paoloeusebi.github.io/caliber/reference/ispline_basis.md), [`ispline_knots()`](https://paoloeusebi.github.io/caliber/reference/ispline_knots.md) |
| Interpretation | [`threshold_locations()`](https://paoloeusebi.github.io/caliber/reference/threshold_locations.md), [`threshold_flags()`](https://paoloeusebi.github.io/caliber/reference/threshold_flags.md), [`crosswalk()`](https://paoloeusebi.github.io/caliber/reference/crosswalk_probabilistic.md), [`crosswalk_probabilistic()`](https://paoloeusebi.github.io/caliber/reference/crosswalk_probabilistic.md), [`exceedance_probabilistic()`](https://paoloeusebi.github.io/caliber/reference/exceedance_probabilistic.md), [`curve_ends()`](https://paoloeusebi.github.io/caliber/reference/curve_ends.md), [`plausible_values()`](https://paoloeusebi.github.io/caliber/reference/plausible_values.md) |
| Intervals | [`ms_vcov()`](https://paoloeusebi.github.io/caliber/reference/ms_vcov.md) + [`ms_intervals()`](https://paoloeusebi.github.io/caliber/reference/ms_intervals.md) (items, thresholds, crosswalk), [`ms_bootstrap()`](https://paoloeusebi.github.io/caliber/reference/ms_bootstrap.md) + [`boot_intervals()`](https://paoloeusebi.github.io/caliber/reference/boot_intervals.md) (curve) |
| Figures | [`wright_map_sqi()`](https://paoloeusebi.github.io/caliber/reference/wright_map_sqi.md), [`plot_expected_scores()`](https://paoloeusebi.github.io/caliber/reference/plot_expected_scores.md) |
| Example data | `psqi_sleep`, `pf_gait` (data sets); [`sleep_data()`](https://paoloeusebi.github.io/caliber/reference/sleep_data.md), [`sleep_eta_true()`](https://paoloeusebi.github.io/caliber/reference/sleep_data.md), [`simulate_sleep_data()`](https://paoloeusebi.github.io/caliber/reference/sleep_data.md), [`gait_data()`](https://paoloeusebi.github.io/caliber/reference/gait_data.md) |

## Data sets

Both are lazy-loaded
([`caliber::psqi_sleep`](https://paoloeusebi.github.io/caliber/reference/psqi_sleep.md),
[`caliber::pf_gait`](https://paoloeusebi.github.io/caliber/reference/pf_gait.md))
and documented in
[`?caliber::psqi_sleep`](https://paoloeusebi.github.io/caliber/reference/psqi_sleep.md)
and
[`?caliber::pf_gait`](https://paoloeusebi.github.io/caliber/reference/pf_gait.md).
`data-raw/datasets.R` rebuilds them from the CSV files in `data-raw/`.

| Data set | Rows × columns | Content |
|----|----|----|
| `psqi_sleep` | 1000 × 11 | Simulated PSQI components C1–C7, PSQI global, actigraphy SQI and the true trait `eta_true`; `simulate_sleep_data(seed = 20261003)` |
| `pf_gait` | 1426 × 15 | NHANES 2001–2002 (MSX_B, PFQ_B, DEMO_B), ages 60+, gait speed and seven mobility items; built by `inst/scripts/gait_data.R`, unweighted |

## Documentation

Rendered vignettes and the function reference are in `docs/` (pkgdown
site), served at <https://paoloeusebi.github.io/caliber/> when GitHub
Pages is enabled for `main` / `docs`.

The digital measure is called `sqi` in the code (from the sleep
example); any measure that falls as the latent trait rises can be used.

## Examples

Vignettes:
[`vignette("sleep-simulation", "caliber")`](https://paoloeusebi.github.io/caliber/articles/sleep-simulation.md),
[`vignette("nhanes-gait", "caliber")`](https://paoloeusebi.github.io/caliber/articles/nhanes-gait.md).
Full analyses (write to `$CALIBER_OUT`, default `output/`):

| Script | What it does | Environment variables |
|----|----|----|
| `sleep_step1.R` | Three approaches on the sleep data, bootstrap | `MS_BOOT` (200), `MS_DRAWS` (2000), `MS_HETERO` (1), `MS_CORES` |
| `sleep_simulation.R` | 6 scenarios × 200 replications | `MS_REPS`, `MS_DRAWS` (0), `MS_HETERO`, `MS_CORES` |
| `sleep_simulation_summary.R` | Tables and figures of the simulation | `MS_HETERO` |
| `gait_data.R` | Builds the gait analysis file from NHANES XPT files | `NHANES_DIR` |
| `gait_dimensionality.R` | Dimensionality and local-dependence checks (needs lavaan) | as gait_main |
| `gait_main.R` | Gait example, bootstrap, figures | `GAIT_BOOT` (200), `GAIT_CORES`, `GAIT_EXCL_AIDS`, `GAIT_SKIP_UNABLE` |

Find them with `system.file("scripts", package = "caliber")`.

## Caveats

Neither the normal-approximation intervals nor the curve bootstrap has
been validated by a coverage simulation. The normal-approximation draws
condition on the selected boosting path and understate the curve’s
variance, so curve bands come from the bootstrap only.

The reported results of the manuscript were produced by the project
scripts of 8 October (sleep) and 11 October (gait), before this package.
With `MS_DRAWS = 0` the package scripts should reproduce them;
`CHECK_NOTES.md` records whether that was confirmed.
