# Monte Carlo intervals from the normal approximation

Draws (items, curve coefficients, g) from the normal approximation of
\[ms_vcov()\] and pushes each draw through \`target(a, d, act)\`, which
returns a named numeric vector (thresholds, crosswalk values, ...).
Every curve draw is projected onto the non-increasing cone.

## Usage

``` r
ms_intervals(fit, vc, target, n_draw = 1000, level = 0.95, seed = NULL)
```

## Arguments

- fit:

  a \[fit_ms()\] result

- vc:

  its \[ms_vcov()\]

- target:

  function(a, d, act) returning a named numeric vector

- n_draw:

  number of draws

- level:

  interval level

- seed:

  optional RNG seed

## Value

data.frame(param, est, se, lo, hi)

## Details

Use these intervals for item parameters, latent thresholds, thresholds
on the digital scale within the data and the probabilistic crosswalk. Do
not use them for the curve \\\mu(\eta)\\: they condition on the selected
boosting path and miss most of its sampling variance; use
\[ms_bootstrap()\].
