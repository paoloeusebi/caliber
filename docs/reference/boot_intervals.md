# Percentile intervals from bootstrap replicates

Percentile intervals from bootstrap replicates

## Usage

``` r
boot_intervals(boot, method, params, level = 0.95)
```

## Arguments

- boot:

  output of \[ms_bootstrap()\]

- method:

  one of \[ms_methods\]

- params:

  column names to summarise

- level:

  interval level

## Value

data.frame(param, se, lo, hi)
