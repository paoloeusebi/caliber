# Covariance blocks of a fitted curve and its items

Normal approximation of the estimates of one \[fit_ms()\] result: items
(observed information of the packed parameters), spline coefficients
(linearised boosting operator at the AIC stop, sandwich with squared
residuals) and (g0, g1) (inverse Hessian given the curve). Conditional
on the boosting path; see \[ms_intervals()\]. For an I-spline fit
(\`curve = "ispline"\`) the curve block is the inverse observed
information of the intercept and the non-zero I-spline coefficients
given items and (g0, g1), with coefficients estimated at zero treated as
fixed (Ramsay 1988), conditional on the selected knots.

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
