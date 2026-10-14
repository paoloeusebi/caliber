# caliber

Monotone calibration of a digital or performance measure to a multi-item PRO.
One latent trait drives the PRO items (graded response model, items scored 0 to 3) and the
digital measure (non-increasing cubic B-spline mean, log-linear residual SD). The curve is
fitted by monotonic likelihood boosting (Leitenstorfer and Tutz 2007) in three ways:
MS-Fixed (fixed-item calibration, the recommended primary analysis), MS-Joint (concurrent
calibration) and MS-EAP (EAP plug-in, the naive comparator).

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
| Curve | `fit_ms()`, `fit_three()`, `mu_sqi()`, `sigma_sqi()` |
| Interpretation | `threshold_locations()`, `threshold_flags()`, `crosswalk()`, `crosswalk_probabilistic()`, `exceedance_probabilistic()`, `curve_ends()`, `plausible_values()` |
| Intervals | `ms_vcov()` + `ms_intervals()` (items, thresholds, crosswalk), `ms_bootstrap()` + `boot_intervals()` (curve) |
| Figures | `wright_map_sqi()`, `plot_expected_scores()` |
| Example data | `sleep_data()`, `simulate_sleep_data()`, `gait_data()` |

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
