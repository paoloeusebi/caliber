# caliber

Monotone calibration of a digital or performance measure to a multi-item PRO.
One latent trait drives the PRO items (graded response model, items scored 0 to 3) and the
digital measure (non-increasing cubic B-spline mean, log-linear residual SD). The curve is
fitted by monotonic likelihood boosting (Leitenstorfer and Tutz 2007) in three ways:
MS-Fixed (fixed-item calibration, the recommended primary analysis), MS-Joint (concurrent
calibration) and MS-EAP (EAP plug-in, the naive comparator).

The curve can instead be fitted with Ramsay's (1988) monotone regression splines
(`fit_ms(..., curve = "ispline")`): a non-negative combination of I-splines on a few interior
knots at quantiles of the latent N(0, 1), fitted by maximum likelihood, with the number of
knots chosen by AIC. Boosting stays the default.

The digital measure can be continuous (`family = "gaussian"`, the default: actigraphy indices,
gait speed, distance) or a count (`family = "poisson"` or `"negbin"`: daily steps, events,
deaths). For counts the spline models log E[y | η], with either curve method; θ of the negative
binomial is estimated by ML. Use Poisson for small counts (events, deaths) and the negative
binomial for large, overdispersed ones (steps): a Poisson count in the thousands locates the
trait more finely than the quadrature resolves, and `fit_ms()` warns. Normal-approximation
intervals (`ms_vcov()`) are Gaussian-only; use `ms_bootstrap()` for counts.

## Install (from this folder)

```r
install.packages(c("mirt", "statmod", "ggplot2", "patchwork", "scales", "viridisLite"))
roxygen2::roxygenise(".")      # builds man/ (and checks NAMESPACE)
devtools::install(".")         # or: R CMD build . && R CMD INSTALL caliber_0.1.0.tar.gz
```

## Main functions

| Step | Functions |
| --- | --- |
| GRM on the PRO | `fit_grm()` |
| Curve | `fit_ms()` (`curve = "boost"` or `"ispline"`), `fit_three()`, `mu_sqi()`, `sigma_sqi()`, `ispline_basis()`, `ispline_knots()` |
| Interpretation | `threshold_locations()`, `threshold_flags()`, `crosswalk()`, `crosswalk_probabilistic()`, `exceedance_probabilistic()`, `curve_ends()`, `plausible_values()` |
| Intervals | `ms_vcov()` + `ms_intervals()` (items, thresholds, crosswalk), `ms_bootstrap()` + `boot_intervals()` (curve) |
| Figures | `wright_map_sqi()`, `plot_expected_scores()` |
| Example data | `psqi_sleep`, `pf_gait` (data sets); `sleep_data()`, `sleep_eta_true()`, `simulate_sleep_data()`, `gait_data()` |

## Data sets

Both are lazy-loaded (`caliber::psqi_sleep`, `caliber::pf_gait`) and documented in `?caliber::psqi_sleep`
and `?caliber::pf_gait`. `data-raw/datasets.R` rebuilds them from the CSV files in `data-raw/`.

| Data set | Rows × columns | Content |
| --- | --- | --- |
| `psqi_sleep` | 1000 × 11 | Simulated PSQI components C1–C7, PSQI global, actigraphy SQI and the true trait `eta_true`; `simulate_sleep_data(seed = 20261003)` |
| `pf_gait` | 1426 × 15 | NHANES 2001–2002 (MSX_B, PFQ_B, DEMO_B), ages 60+, gait speed and seven mobility items; built by `inst/scripts/gait_data.R`, unweighted |

## Documentation

Rendered vignettes and the function reference are in `docs/` (pkgdown site), served at
<https://paoloeusebi.github.io/caliber/> when GitHub Pages is enabled for `main` / `docs`.

The digital measure is called `sqi` in the code (from the sleep example); any measure that falls
as the latent trait rises can be used.

## Examples

Vignettes: `vignette("sleep-simulation", "caliber")`, `vignette("nhanes-gait", "caliber")`.
Full analyses (write to `$CALIBER_OUT`, default `output/`):

| Script | What it does | Environment variables |
| --- | --- | --- |
| `sleep_step1.R` | Three approaches on the sleep data, bootstrap | `MS_BOOT` (200), `MS_DRAWS` (2000), `MS_HETERO` (1), `MS_CORES` |
| `sleep_simulation.R` | 6 scenarios × 200 replications | `MS_REPS`, `MS_DRAWS` (0), `MS_HETERO`, `MS_CORES` |
| `sleep_simulation_summary.R` | Tables and figures of the simulation | `MS_HETERO` |
| `gait_data.R` | Builds the gait analysis file from NHANES XPT files | `NHANES_DIR` |
| `gait_dimensionality.R` | Dimensionality and local-dependence checks (needs lavaan) | as gait_main |
| `gait_main.R` | Gait example, bootstrap, figures | `GAIT_BOOT` (200), `GAIT_CORES`, `GAIT_EXCL_AIDS`, `GAIT_SKIP_UNABLE` |

Find them with `system.file("scripts", package = "caliber")`.

## Caveats

Neither the normal-approximation intervals nor the curve bootstrap has been validated by a
coverage simulation. The normal-approximation draws condition on the selected boosting path and
understate the curve's variance, so curve bands come from the bootstrap only.

The reported results of the manuscript were produced by the project scripts of 8 October
(sleep) and 11 October (gait), before this package. With `MS_DRAWS = 0` the package scripts
should reproduce them; `CHECK_NOTES.md` records whether that was confirmed.
