# Crosswalks from the digital measure to the PRO

\`crosswalk_probabilistic()\` gives \\E\[C_j \mid y = s\]\\ and the
expected total, averaging over \\p(\eta_q \mid s) \propto w_q N(s;
\mu(\eta_q), \sigma(\eta_q)^2)\\ (or the Poisson / negative binomial
probability of \`s\` for a count curve). \`crosswalk()\` adds the
structural crosswalk \\E\[C_j \mid \eta^\*(s)\]\\ (columns \`eta_star\`,
\`Estar\_\*\`; NA outside the curve's range) and prefixes the
probabilistic columns with \`E\_\`.

## Usage

``` r
crosswalk_probabilistic(
  grid,
  a,
  d,
  act,
  items = paste0("C", seq_along(a)),
  total_name = "total"
)

crosswalk(
  grid,
  a,
  d,
  act,
  items = paste0("C", seq_along(a)),
  total_name = "total",
  value_name = "SQI"
)
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

- total_name:

  name of the total-score column

- value_name:

  name of the column holding \`grid\` (\`crosswalk()\` only)
