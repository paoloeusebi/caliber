# Item thresholds on the latent and digital scales

\\b\_{jk} = -d\_{jk} / a_j\\ and \\s\_{jk} = \mu(b\_{jk})\\, the mean of
the digital measure among persons sitting at the threshold.

## Usage

``` r
threshold_locations(a, d, act)
```

## Arguments

- a, d:

  item slopes and intercepts (items x 3)

- act:

  curve object

## Value

list(b, sqi), both items x 3 matrices
