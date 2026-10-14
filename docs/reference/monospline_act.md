# The equation of the digital measure

\`monospline_act()\` builds the curve object used everywhere in the
package: \`list(mu = function(eta), sigma = sigma(0), log_sigma = c(g0,
g1) or NULL, range, beta)\`. \`mu_sqi()\` and \`sigma_sqi()\` evaluate
the mean \\\mu(\eta)\\ and the residual SD \\\sigma(\eta)\\ (constant,
or \\\exp(g_0 + g_1 \eta)\\).

## Usage

``` r
monospline_act(beta, sigma, log_sigma = NULL)

mu_sqi(eta, act)

sigma_sqi(eta, act)
```

## Arguments

- beta:

  24 non-increasing spline coefficients

- sigma:

  constant residual SD (ignored when \`log_sigma\` is given)

- log_sigma:

  \`c(g0, g1)\` for the log-linear SD, or NULL

- eta:

  latent values

- act:

  a curve object
