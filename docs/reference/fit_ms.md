# Fit the monotone curve

Monotone B-spline equation of the digital measure by monotonic
likelihood boosting (Leitenstorfer and Tutz 2007, GMBBoost; the default)
or by Ramsay's (1988) monotone regression splines (\`curve =
"ispline"\`), in one of three ways (the boosting steps, tuning values
and stopping rule are the same):

## Usage

``` r
fit_ms(
  grm,
  Y,
  sqi,
  mode = c("fixed", "joint", "eap"),
  lambda = 20,
  shrinkage = 0.1,
  max_iter = 10000,
  patience = 100,
  criterion = "AIC",
  item_every = 5,
  hetero = TRUE,
  curve = c("boost", "ispline"),
  n_knots = "AIC",
  max_knots = 5,
  order = 2,
  family = c("gaussian", "poisson", "negbin")
)
```

## Arguments

- grm:

  output of \[fit_grm()\] on the same items

- Y:

  item responses (0..3, NA allowed), n x J

- sqi:

  the digital measure, length n (non-increasing in eta)

- mode:

  "fixed" (MS-Fixed), "joint" (MS-Joint) or "eap" (MS-EAP)

- lambda:

  ridge penalty of the fitting step (default 20)

- shrinkage:

  step length nu (default 0.1)

- max_iter:

  maximum number of boosting iterations

- patience:

  stop when the AIC has not improved for this many iterations

- criterion:

  "AIC" (default) or "BIC" for the stopping iteration

- item_every:

  MS-Joint: item M-step every this many iterations

- hetero:

  TRUE (default) for log sigma(eta) = g0 + g1 eta, FALSE for a constant
  residual SD (sensitivity analysis)

- curve:

  "boost" (default): monotonic likelihood boosting of the 24-coefficient
  B-spline, as above. "ispline": Ramsay's (1988) monotone regression
  spline, \\\mu(\eta) = \beta_0 - \sum_m a_m I_m(\eta)\\, \\a_m \ge 0\\,
  with I-splines of order \`order\` on \`n_knots\` interior knots at
  quantiles of N(0, 1) (see \[ispline_knots()\]), fitted by maximum
  likelihood under the bounds (EM over the nodes for "fixed" and
  "joint", direct for "eap"). The three ways of handling the latent
  trait are the same. \`lambda\`, \`shrinkage\`, \`patience\` and
  \`item_every\` are then unused, and \`max_iter\` caps the EM
  iterations.

- n_knots:

  \`curve = "ispline"\`: number of interior knots, or "AIC" (default) to
  choose among 0 to \`max_knots\` by \`criterion\`

- max_knots:

  largest number of interior knots tried

- order:

  I-spline order k (default 2, piecewise quadratic I-splines, so the
  curve is a quadratic spline with continuous slope)

- family:

  distribution of the digital measure given the latent trait. "gaussian"
  (default): a continuous measure (actigraphy index, gait speed,
  distance), normal with the mean curve and a constant or log-linear SD.
  "poisson" or "negbin" (negative binomial, Var = mu + mu^2 / theta): a
  count measure (daily steps, events, deaths); the spline is then fitted
  on the log scale, log E\[y \| eta\] = f(eta) non-increasing, by the
  GLM version of the boosting (Leitenstorfer and Tutz 2007: penalised
  Fisher-scoring steps with the Fisher weights of the family,
  deviance-based split choice and AIC stop) or by the I-spline ML fit.
  theta is estimated by ML along the way. \`hetero\` applies to
  "gaussian" only. For count families the approximate hat matrix is
  accumulated on the units (nodes or persons) with the current Fisher
  weights, and \[ms_vcov()\] is not available; use \[ms_bootstrap()\].
  "poisson" suits small counts (deaths, events). A large Poisson count
  (mean in the thousands) locates the latent trait more finely than the
  61 quadrature nodes resolve, which biases the curve; fit_ms() warns
  when this happens. Large counts such as daily steps are overdispersed
  anyway and should use "negbin".

## Value

list(a, d, act, mode, converged, iter, iter_aic, iter_bic, df,
stop_reason, path, item_ll (not for "eap"), R_coef (coefficient operator
at the AIC stop, used by \[ms_vcov()\]), eta_cov (EAP scores, "eap"
only), curve). For \`curve = "ispline"\`, \`iter\` is the number of EM
iterations, \`path\` the knot selection table, \`R_coef\` NULL, and the
result also has n_knots, knots, order and coef (beta0 and the I-spline
coefficients).

## Details

"fixed" MS-Fixed. Items held at the GRM estimates. Each iteration:
E-step (posterior p(eta_q \| C_i, SQI_i) on the nodes), sigma update,
one boosting step on the posterior-weighted node totals. "joint"
MS-Joint. As "fixed", but each iteration also takes an M-step for the
item parameters from the same posterior (expected category counts per
node), so the items and the curve are estimated together and the latent
scale is defined by the PSQI and the SQI jointly. Items start at the GRM
estimates. The item M-step is taken every \`item_every\` iterations
(default 5); skipping it in between keeps every iteration a generalised
EM step, and the items settle long before the AIC stop. "eap" MS-EAP.
The EAP scores from the PSQI-only GRM are treated as the covariate: the
paper's algorithm exactly as published (observed covariate, normal
response).

Residual SD: with hetero = TRUE (default) log sigma(eta) = g0 + g1 eta.
The sigma step then updates (g0, g1) by Newton steps on the
posterior-expected (MS-EAP: observed) log-likelihood given the curve,
units are weighted by 1 / sigma_u^2 in the fitting step, and (g0, g1)
are not counted in tr(H_l) (their number is fixed, so the AIC minimum
does not move). With hetero = FALSE sigma is constant, with its
closed-form update.

Stopping: AIC(l) = Dev_l + 2 tr(H_l), H_l = I - (I - M_l)...(I - M_0)
the approximate hat matrix of the person-level fitted values (paper eq.
11). Dev_l is -2 log L(C, SQI) with eta integrated out ("fixed",
"joint") or -2 log L(SQI \| EAP) ("eap"). For "joint" the 28 item
parameters are not counted in tr(H_l); their number is the same at every
iteration, so the AIC minimum does not move. Search stops when the AIC
has not improved for \`patience\` iterations, no monotone split exists,
or max_iter.

## References

Leitenstorfer F, Tutz G (2007). Generalized monotonic regression based
on B-splines with an application to air pollution data. Biostatistics 8,
654-673.

Ramsay JO (1988). Monotone regression splines in action. Statistical
Science 3, 425-461.
