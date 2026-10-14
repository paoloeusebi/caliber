# Expected item scores against the digital measure

Structural and probabilistic crosswalks, estimated (solid) and true
(dashed, columns ending in \`\_true\`), with a rug of the observed
measure.

## Usage

``` r
plot_expected_scores(
  cw,
  sqi,
  items = paste0("C", 1:7),
  x_label = "Actigraphy SQI (digital measure)"
)
```

## Arguments

- cw:

  crosswalk table from \[crosswalk()\] with \`\_true\` columns added

- sqi:

  observed digital measure

- items:

  item names

- x_label:

  x-axis title
