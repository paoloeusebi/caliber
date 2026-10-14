# Person-item targeting plot (Wright map) on the SQI scale, in the style of easyRasch2::RMtargeting(): person histogram (top) and one bar per item split into response-category bands (bottom). The x axis is the digital measure (SQI) instead of logits.

Person-item targeting plot (Wright map) on the SQI scale, in the style
of easyRasch2::RMtargeting(): person histogram (top) and one bar per
item split into response-category bands (bottom). The x axis is the
digital measure (SQI) instead of logits.

## Usage

``` r
wright_map_sqi(
  person_sqi,
  thr,
  item_labels,
  xlim = c(15, 90),
  binwidth = 2,
  category_labels = c("0 (no difficulty)", "1", "2", "3 (severe difficulty)"),
  person_fill = "#0072B2",
  note = NULL,
  x_label = "Actigraphy SQI (digital measure; higher = better sleep)",
  x_breaks = seq(0, 100, by = 10)
)
```

## Arguments

- person_sqi:

  person locations on the digital scale

- thr:

  data.frame with item, k, sqi_est (threshold on the digital scale)

- item_labels:

  labels for the items, in the order of unique(thr\$item)

- xlim:

  x range (widened to include persons and thresholds)

- binwidth:

  histogram bin width

- category_labels:

  labels of the four response categories

- person_fill:

  histogram colour

- note:

  caption text placed under the plot

- x_label:

  x-axis title

- x_breaks:

  x-axis breaks

## Value

a patchwork object
