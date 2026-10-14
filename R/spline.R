# ---- The monotone B-spline equation of the digital measure ---------------------

spline_range  <- c(-4, 4)
spline_degree <- 3
spline_knots  <- spline_range[1] + diff(spline_range) / 21 * (-3:24)
n_basis       <- length(spline_knots) - spline_degree - 1        # 24

#' Cubic B-spline basis of the curve
#'
#' 24 cubic B-splines on 20 equally spaced interior knots over \eqn{[-4, 4]}.
#' \eqn{\eta} is clamped to \eqn{[-4, 4]}, so the curve is constant at its end
#' values outside that range.
#' @param eta latent values
#' @return length(eta) x 24 matrix
#' @export
bspline_basis <- function(eta) {
  x <- pmin(pmax(eta, spline_knots[spline_degree + 1]),
            spline_knots[length(spline_knots) - spline_degree])
  splineDesign(spline_knots, x, ord = spline_degree + 1, outer.ok = TRUE)
}

#' The equation of the digital measure
#'
#' `monospline_act()` builds the curve object used everywhere in the package:
#' `list(mu = function(eta), sigma = sigma(0), log_sigma = c(g0, g1) or NULL,
#' range, beta)`. `mu_sqi()` and `sigma_sqi()` evaluate the mean
#' \eqn{\mu(\eta)} and the residual SD \eqn{\sigma(\eta)} (constant, or
#' \eqn{\exp(g_0 + g_1 \eta)}).
#' @param beta 24 non-increasing spline coefficients
#' @param sigma constant residual SD (ignored when `log_sigma` is given)
#' @param log_sigma `c(g0, g1)` for the log-linear SD, or NULL
#' @param eta latent values
#' @param act a curve object
#' @export
monospline_act <- function(beta, sigma, log_sigma = NULL)
  list(mu = function(eta) drop(bspline_basis(eta) %*% beta),
       sigma = if (is.null(log_sigma)) sigma else exp(log_sigma[1]),
       log_sigma = log_sigma, range = spline_range, beta = beta)

#' @rdname monospline_act
#' @export
mu_sqi <- function(eta, act) act$mu(eta)

#' @rdname monospline_act
#' @export
sigma_sqi <- function(eta, act)
  if (is.null(act$log_sigma)) rep(act$sigma, length(eta)) else
    exp(act$log_sigma[1] + act$log_sigma[2] * eta)

#' Newton steps for (g0, g1) given the curve: maximise
#' sum_u [-N_u (g0 + g1 x_u) - RSS_u / (2 exp(2 (g0 + g1 x_u)))] (concave).
#' @noRd
update_log_sigma <- function(g, x, N, RSS, max_steps = 25, tol = 1e-10) {
  X <- cbind(1, x)
  for (it in 1:max_steps) {
    e2 <- RSS * exp(-2 * drop(X %*% g))
    step <- solve(2 * crossprod(X, X * e2), crossprod(X, e2 - N))
    g <- g + drop(step)
    if (max(abs(step)) < tol) break
  }
  g
}

#' log w_q + log P(C_i | eta_q) + log N(y_i; mu(eta_q), sigma(eta_q)^2), n x Q.
#' @noRd
log_joint_nodes <- function(item_ll, sqi, act) {
  sweep(item_ll, 2, log_wts, "+") +
    dnorm(outer(sqi, mu_sqi(nodes, act), "-"), 0,
          rep(sigma_sqi(nodes, act), each = length(sqi)), log = TRUE)
}
