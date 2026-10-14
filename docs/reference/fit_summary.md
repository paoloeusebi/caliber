# Standard summary of one fit

Item parameters, residual SD, mu(eta) at \`eta_report\`, latent
thresholds b_jk, thresholds on the digital scale s_jk = mu(b_jk) and the
curve on \`curve_grid\` (columns \`curve_1\`, ...), plus the stopping
iteration and df. Used as the default summary of \[ms_bootstrap()\].

## Usage

``` r
fit_summary(
  f,
  items = paste0("C", seq_along(f$a)),
  eta_report = -3:3,
  curve_grid = seq(-3.5, 3.5, by = 0.05)
)
```

## Arguments

- f:

  a \[fit_ms()\] result

- items:

  item names

- eta_report:

  latent values at which mu(eta) is reported

- curve_grid:

  latent grid for curve bands

## Value

named numeric vector
