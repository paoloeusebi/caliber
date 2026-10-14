# The equation of the digital measure

\`monospline_act()\` builds the curve object used everywhere in the
package: \`list(mu = function(eta), sigma = sigma(0), log_sigma = c(g0,
g1) or NULL, range, beta, type)\`; \`ispline_act()\` builds the same
object for a curve fitted with \`fit_ms(curve = "ispline")\`.
\`mu_sqi()\` and \`sigma_sqi()\` evaluate the mean \\\mu(\eta)\\ and the
residual SD \\\sigma(\eta)\\ (constant, or \\\exp(g_0 + g_1 \eta)\\).

## Usage

``` r
ispline_act(beta, interior, order, sigma, log_sigma = NULL)

monospline_act(beta, sigma, log_sigma = NULL)

mu_sqi(eta, act)

sigma_sqi(eta, act)
```

## Arguments

- beta:

  24 non-increasing spline coefficients

- interior, order:

  \`ispline_act()\`: interior knots and I-spline order of an I-spline
  fit (\`fit_ms(curve = "ispline")\`); \`beta\` are then the
  non-increasing B-spline coefficients of order \`order + 1\` on those
  knots.

- sigma:

  constant residual SD (ignored when \`log_sigma\` is given)

- log_sigma:

  \`c(g0, g1)\` for the log-linear SD, or NULL

- eta:

  latent values

- act:

  a curve object
