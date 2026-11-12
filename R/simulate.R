# ---- Data-generating scenarios (misspecification study) -------------------------
# The six scenarios of the earlier study, unchanged, and a seventh in which X
# and Y measure correlated but different traits. They are truths only; no
# analysis model below assumes their shape.

#' True item parameters and simulation scenarios of the sleep example
#'
#' `true_items`: the seven PSQI-component GRM items of the simulated sleep
#' example (slopes a, thresholds b, intercepts d = -a b). `scenarios`: the
#' seven data-generating models of the simulation study, each a list(label,
#' m = true mean curve E[Y | eta], sd = true residual SD function, ld = local
#' dependence, and optionally m_y and rho). In the simulation X is the
#' multi-item scale (items X1..X7, stored as C1..C7) and Y the digital measure.
#' Scenario 7 ("different concepts") has Y driven by its own trait
#' xi = rho eta + sqrt(1 - rho^2) u, correlated rho with the trait eta of X:
#' Y = m_y(xi) + N(0, sd^2), and m is then E[Y | eta], the curve the methods
#' estimate. They are truths only; no analysis model assumes their shape.
#' @format lists
#' @export
true_items <- list(
  a = c(2.2, 1.4, 1.6, 1.8, 1.3, 0.9, 1.5),
  b = cbind(c(-0.8, -0.9, 0.0, 0.4, -2.2, 1.5, -0.9),
            c( 0.9,  0.5, 0.9, 1.1,  0.8, 2.3,  0.9),
            c( 2.2,  1.6, 1.9, 1.8,  3.0, 2.8,  2.5)))
true_items$d <- -true_items$a * true_items$b

# Scenario 7 (different concepts): Y measures a related but different trait,
# xi = rho eta + sqrt(1 - rho^2) u, where the nuisance factor u (drawn in every
# scenario) is Y's own part. The true curve of Y on the trait of X is then the
# average of m_y over xi given eta, computed by mean_given_eta() with
# Gauss-Hermite quadrature over u (exact for a quadratic m_y). With rho = 0.6 it
# is 71.04 - 4.8 eta - 0.54 eta^2, monotone for eta > -4.4, and cor(X total, Y)
# is about -0.42, against -0.71 in the baseline.
y_trait <- function(sc, eta, u)
  if (is.null(sc$rho)) eta else sc$rho * eta + sqrt(1 - sc$rho^2) * u
mean_given_eta <- function(m_y, rho) {
  uq <- statmod::gauss.quad.prob(41, dist = "normal")
  function(e) vapply(e, function(x) sum(uq$weights * m_y(rho * x + sqrt(1 - rho^2) * uq$nodes)), 0)
}

#' @rdname true_items
#' @export
scenarios <- list(
  baseline = list(label = "Baseline (curve turns at eta = -2.67)",
                  m = function(e) 72 - 8 * e - 1.5 * e^2, sd = function(e) 6 + 0 * e, ld = FALSE),
  local_dep = list(label = "Local dependence X4-Y",
                  m = function(e) 72 - 8 * e - 1.5 * e^2, sd = function(e) 4.8 + 0 * e, ld = TRUE),
  hetero = list(label = "Heteroscedastic Y error",
                  m = function(e) 72 - 8 * e - 1.5 * e^2, sd = function(e) 6 * exp(0.4 * e), ld = FALSE),
  gentle_logistic = list(label = "Gentle logistic",
                  m = function(e) 40 + 45 / (1 + exp(1.4 * (e - 0.3))), sd = function(e) 6 + 0 * e, ld = FALSE),
  steep_logistic = list(label = "Steep logistic",
                  m = function(e) 35 + 50 / (1 + exp(3 * (e - 0.8))), sd = function(e) 6 + 0 * e, ld = FALSE),
  plateau = list(label = "Plateau, then decline",
                  m = function(e) 80 - 6 * log(1 + exp(2 * e)), sd = function(e) 6 + 0 * e, ld = FALSE),
  diff_concept = list(label = "Different concepts (trait correlation 0.6)",
                  m = mean_given_eta(function(x) 72 - 8 * x - 1.5 * x^2, 0.6),
                  m_y = function(x) 72 - 8 * x - 1.5 * x^2, rho = 0.6,
                  sd = function(e) 6 + 0 * e, ld = FALSE))
ld_item <- 4; ld_item_load <- 1.2; ld_sqi_load <- -3.6

#' Simulate one replication of a scenario
#'
#' RNG order: eta, nuisance factor u, items C1..C7, SQI residuals; seed
#' 20261003 + 10000 s + r. SQI is rounded to one decimal.
#' @param s scenario index (1..7, order of [scenarios])
#' @param r replication number
#' @param n sample size
#' @return list(Y, sqi, eta)
#' @export
simulate_scenario <- function(s, r, n = 1000) {
  sc <- scenarios[[s]]
  set.seed(20261003 + 10000 * s + r)
  eta <- rnorm(n); u <- rnorm(n)
  Y <- sapply(seq_along(true_items$a), function(j) {
    lin_eta <- if (sc$ld && j == ld_item) eta + ld_item_load * u / true_items$a[j] else eta
    draw_categories(grm_probs(lin_eta, true_items$a[j], true_items$d[j, ]))
  })
  colnames(Y) <- paste0("C", seq_along(true_items$a))
  shift <- if (sc$ld) ld_sqi_load * u else 0
  mean_y <- if (is.null(sc$rho)) sc$m(eta) else sc$m_y(y_trait(sc, eta, u))
  sqi <- round(mean_y + shift + rnorm(n, 0, sc$sd(eta)), 1)
  list(Y = Y, sqi = sqi, eta = eta)
}

#' Truth of a scenario
#'
#' True threshold locations m(b_jk) and the true probabilistic crosswalk with
#' the marginal SQI density on a grid (quadrature: 61 nodes for eta x 21 for u).
#' @param s scenario index
#' @param grid SQI grid
#' @return list(thr, cw, wt)
#' @export
scenario_truth <- function(s, grid) {
  sc <- scenarios[[s]]
  uq <- gauss.quad.prob(21, dist = "normal")
  E <- expand.grid(q = seq_along(nodes), v = seq_along(uq$nodes))
  eta <- nodes[E$q]; u <- uq$nodes[E$v]
  w   <- quad$weights[E$q] * uq$weights[E$v]
  ex  <- sapply(seq_along(true_items$a), function(j) {
    lin <- if (sc$ld && j == ld_item) eta + ld_item_load * u / true_items$a[j] else eta
    grm_expected(lin, true_items$a[j], true_items$d[j, ])
  })
  mean_s <- (if (is.null(sc$rho)) sc$m(eta) else sc$m_y(y_trait(sc, eta, u))) +
    (if (sc$ld) ld_sqi_load * u else 0)
  sd_s   <- sc$sd(eta)
  cw <- t(sapply(grid, function(g) {
    f <- w * dnorm(g, mean_s, sd_s)
    c(colSums(f * ex) / sum(f), dens = sum(f))
  }))
  items <- paste0("C", seq_along(true_items$a))
  cw <- as.data.frame(cw); names(cw)[seq_along(items)] <- items
  cw$PSQI <- rowSums(cw[items])
  list(thr = sc$m(true_items$b), cw = cw, wt = cw$dens / sum(cw$dens))
}


# ---- The single simulated data set of the worked example ------------------------

#' The simulated sleep data set of the worked example
#'
#' `sleep_data()` returns the data set analysed in the manuscript (n = 1000,
#' generated by `simulate_sleep_data(seed = 20261003)`): seven PSQI components
#' C1..C7 (0..3), the PSQI global score and the actigraphy sleep-quality index
#' SQI. `sleep_eta_true()` gives the true latent values (evaluation only).
#' Both are taken from the package data set [psqi_sleep].
#' `simulate_sleep_data()` regenerates it: GRM items of [true_items] and
#' SQI = 72 - 8 eta - 1.5 eta^2 + N(0, 6^2), RNG order eta, items, SQI noise.
#' `sleep_truth()` is that generating model as (a, d, act) for comparisons; its
#' `act$range` is the falling branch of the quadratic.
#' @param n sample size
#' @param seed RNG seed
#' @export
sleep_data <- function() {
  d <- caliber::psqi_sleep
  d[setdiff(names(d), "eta_true")]
}

#' @rdname sleep_data
#' @export
sleep_eta_true <- function()
  caliber::psqi_sleep[c("id", "eta_true")]

#' @rdname sleep_data
#' @export
simulate_sleep_data <- function(n = 1000, seed = 20261003) {
  set.seed(seed)
  eta <- rnorm(n)
  Y <- sapply(seq_along(true_items$a), function(j)
    draw_categories(grm_probs(eta, true_items$a[j], true_items$d[j, ])))
  colnames(Y) <- paste0("C", seq_along(true_items$a))
  sqi <- 72 - 8 * eta - 1.5 * eta^2 + rnorm(n, 0, 6)
  data.frame(id = seq_len(n), Y, PSQI_global = rowSums(Y), SQI = round(sqi, 1),
             eta_true = eta)
}

#' @rdname sleep_data
#' @export
sleep_truth <- function()
  list(a = true_items$a, d = true_items$d,
       act = list(mu = function(e) 72 - 8 * e - 1.5 * e^2, sigma = 6, log_sigma = NULL,
                  range = c(-8 / 3, 6)))

#' PSQI component labels
#' @export
sleep_item_labels <- c("C1 Subjective quality", "C2 Latency", "C3 Duration",
                       "C4 Efficiency", "C5 Disturbances", "C6 Medication",
                       "C7 Daytime dysfunction")
