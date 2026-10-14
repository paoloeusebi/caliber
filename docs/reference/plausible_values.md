# One plausible value of the latent trait per person

For MS-Joint and MS-Fixed, drawn from \\p(\eta \mid C_i, y_i)\\ on the
quadrature nodes and jittered by a quarter node spacing; for MS-EAP,
from \\p(\eta \mid C_i)\\ with mirt (the digital measure not used).

## Usage

``` r
plausible_values(fit, grm, sqi, seed = NULL)
```

## Arguments

- fit:

  a \[fit_ms()\] result

- grm:

  the \[fit_grm()\] result it used

- sqi:

  the digital measure

- seed:

  RNG seed
