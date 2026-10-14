# Covariance blocks of a fitted curve and its items

Normal approximation of the estimates of one \[fit_ms()\] result: items
(observed information of the packed parameters), spline coefficients
(linearised boosting operator at the AIC stop, sandwich with squared
residuals) and (g0, g1) (inverse Hessian given the curve). Conditional
on the boosting path; see \[ms_intervals()\].

## Usage

``` r
ms_vcov(fit, Y, sqi)
```

## Arguments

- fit:

  a \[fit_ms()\] result (fitted with \`criterion = "AIC"\`)

- Y, sqi:

  the data it was fitted to

## Value

list(theta, V_items, beta, V_beta, g, V_g, hetero)
