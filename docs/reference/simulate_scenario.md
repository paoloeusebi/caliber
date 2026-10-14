# Simulate one replication of a scenario

RNG order: eta, nuisance factor u, items C1..C7, SQI residuals; seed
20261003 + 10000 s + r. SQI is rounded to one decimal.

## Usage

``` r
simulate_scenario(s, r, n = 1000)
```

## Arguments

- s:

  scenario index (1..6, order of \[scenarios\])

- r:

  replication number

- n:

  sample size

## Value

list(Y, sqi, eta)
