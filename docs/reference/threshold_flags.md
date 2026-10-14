# Flags for thresholds that cannot be read as point estimates

A threshold is flagged as sparse when fewer than \`min_at_or_above\`
persons respond at or above its category, and as at the floor when it
lies within \`floor_tol\` of the curve's floor (such thresholds cannot
be told apart and are reported as "below about" the floor). A category
nobody used is not estimable.

## Usage

``` r
threshold_flags(Y, thr_sqi, act, floor_tol = 0.01, min_at_or_above = 30)
```

## Arguments

- Y:

  item responses

- thr_sqi:

  items x 3 matrix of thresholds on the digital scale

- act:

  curve object

- floor_tol:

  tolerance in units of the digital measure

- min_at_or_above:

  minimum count at or above the category

## Value

data.frame(item, k, n_at_or_above, not_estimable, flag_sparse, at_floor)
