# Cubic B-spline basis of the curve

24 cubic B-splines on 20 equally spaced interior knots over \\\[-4,
4\]\\. \\\eta\\ is clamped to \\\[-4, 4\]\\, so the curve is constant at
its end values outside that range.

## Usage

``` r
bspline_basis(eta)
```

## Arguments

- eta:

  latent values

## Value

length(eta) x 24 matrix
