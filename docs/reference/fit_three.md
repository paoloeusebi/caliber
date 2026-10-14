# All three methods on one data set

Fits MS-Joint, MS-EAP and MS-Fixed, sharing one PRO-only GRM.

## Usage

``` r
fit_three(Y, sqi, grm = fit_grm(Y), ...)
```

## Arguments

- Y:

  item responses

- sqi:

  the digital measure

- grm:

  output of \[fit_grm()\] (fitted if not given)

- ...:

  passed to \[fit_ms()\]

## Value

named list of three \[fit_ms()\] results
