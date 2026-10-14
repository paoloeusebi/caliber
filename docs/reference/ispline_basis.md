# I-spline basis

Monotone splines of Ramsay (1988): each \\I_m\\ rises from 0 to 1 and is
the integral of an M-spline of order \`order\`. Computed as right
cumulative sums of B-splines of order \`order + 1\` on the clamped knot
vector over \\\[-4, 4\]\\; \\\eta\\ is clamped to that range as for
\[bspline_basis()\].

## Usage

``` r
ispline_basis(eta, interior = numeric(0), order = 2)
```

## Arguments

- eta:

  latent values

- interior:

  interior knots (see \[ispline_knots()\])

- order:

  I-spline order k (2 = piecewise quadratic, Ramsay's default)

## Value

length(eta) x (length(interior) + order) matrix

## References

Ramsay JO (1988). Monotone regression splines in action. Statistical
Science 3, 425-461.
