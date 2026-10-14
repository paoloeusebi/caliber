# Exceedance probabilities given the digital measure

\\P(C_j \ge k \mid y = s) = \sum_q P(C_j \ge k \mid \eta_q) p(\eta_q
\mid s)\\, k = 1..3, on a grid (probabilistic, as the crosswalk).

## Usage

``` r
exceedance_probabilistic(grid, a, d, act, items = paste0("C", seq_along(a)))
```

## Arguments

- grid:

  values of the digital measure

- a, d:

  item parameters

- act:

  curve object

- items:

  item names

## Value

grid x (items x 3) matrix, columns item-major within k
