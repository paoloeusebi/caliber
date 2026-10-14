# Floor, ceiling and approximate R-squared of a fitted curve

The curve is constant outside \\\[-4, 4\]\\, so \\\mu(4)\\ is its floor
and \\\mu(-4)\\ its ceiling: no latent value has a mean outside them.
R-squared is approximate (residual SD at eta = 0 against the SD of y).

## Usage

``` r
curve_ends(act, sqi)
```

## Arguments

- act:

  curve object

- sqi:

  observed digital measure
