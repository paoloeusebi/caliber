# ---- Monotone regression splines (Ramsay 1988) -------------------------------------
# The alternative to boosting for the curve: mu(eta) = beta0 - sum_m a_m I_m(eta),
# a_m >= 0, with I-splines I_m of order k (Ramsay 1988, eq. 4-5) on a few interior
# knots at quantiles of N(0, 1) (median, terciles, quartiles, ...; Ramsay's
# preliminary placement), fitted by maximum likelihood under the non-negativity
# constraints. An I-spline of order k is a right cumulative sum of B-splines of
# order k + 1 on the same knots, so the curve is also a B-spline of order k + 1
# with non-increasing coefficients beta_j = beta0 - sum_{m < j} a_m, and every
# function that takes a curve object works unchanged.

#' Interior knots for the I-spline curve
#'
#' Ramsay's preliminary placement: one knot at the median, two at the terciles,
#' three at the quartiles, and so on, here of the latent N(0, 1).
#' @param n_knots number of interior knots
#' @return sorted interior knots on the latent scale
#' @export
ispline_knots <- function(n_knots)
  if (n_knots == 0) numeric(0) else stats::qnorm(seq_len(n_knots) / (n_knots + 1))

ispline_knot_vector <- function(interior, order)
  c(rep(spline_range[1], order + 1), sort(interior), rep(spline_range[2], order + 1))

#' I-spline basis
#'
#' Monotone splines of Ramsay (1988): each \eqn{I_m} rises from 0 to 1 and is
#' the integral of an M-spline of order `order`. Computed as right cumulative
#' sums of B-splines of order `order + 1` on the clamped knot vector over
#' \eqn{[-4, 4]}; \eqn{\eta} is clamped to that range as for [bspline_basis()].
#' @param eta latent values
#' @param interior interior knots (see [ispline_knots()])
#' @param order I-spline order k (2 = piecewise quadratic, Ramsay's default)
#' @return length(eta) x (length(interior) + order) matrix
#' @references Ramsay JO (1988). Monotone regression splines in action.
#'   Statistical Science 3, 425-461.
#' @export
ispline_basis <- function(eta, interior = numeric(0), order = 2) {
  x <- pmin(pmax(eta, spline_range[1]), spline_range[2])
  B <- splineDesign(ispline_knot_vector(interior, order), x, ord = order + 1,
                    outer.ok = TRUE)
  M <- ncol(B)
  I <- t(apply(B[, M:1, drop = FALSE], 1, cumsum))[, M:1, drop = FALSE]
  if (length(x) == 1) I <- matrix(I, 1)
  I[, -1, drop = FALSE]                       # column 1 is the constant 1
}

#' @rdname monospline_act
#' @param interior,order `ispline_act()`: interior knots and I-spline order of
#'   an I-spline fit (`fit_ms(curve = "ispline")`); `beta` are then the
#'   non-increasing B-spline coefficients of order `order + 1` on those knots.
#' @export
ispline_act <- function(beta, interior, order, sigma, log_sigma = NULL) {
  knots <- ispline_knot_vector(interior, order)
  list(mu = function(eta) drop(splineDesign(knots, pmin(pmax(eta, spline_range[1]), spline_range[2]),
                                            ord = order + 1, outer.ok = TRUE) %*% beta),
       sigma = if (is.null(log_sigma)) sigma else exp(log_sigma[1]),
       log_sigma = log_sigma, range = spline_range, beta = beta,
       type = "ispline", interior = interior, order = order)
}

#' A curve object like `act` with new coefficients and residual SD.
#' @noRd
act_with <- function(act, beta, sigma, log_sigma = NULL)
  if (identical(act$type, "ispline")) ispline_act(beta, act$interior, act$order, sigma, log_sigma) else
    monospline_act(beta, sigma, log_sigma)

#' B-spline coefficients from (beta0, a): beta_j = beta0 - sum_{m < j} a_m.
#' @noRd
ispline_beta <- function(beta0, a) beta0 - c(0, cumsum(a))

#' One I-spline fit with fixed knots, on the standardised measure y.
#' Units u: the nodes (MS-Fixed, MS-Joint; posterior-weighted totals N, Sy, Syy,
#' refreshed by EM) or the persons at their EAP scores (MS-EAP; N = 1, one
#' M-step is the ML fit). The M-step maximises the (expected) log-likelihood
#' sum_u [-N_u log s_u - (Syy_u - 2 mu_u Sy_u + mu_u^2 N_u) / (2 s_u^2)]
#' over (beta0, a >= 0, g0, g1) by L-BFGS-B (Ramsay fits with gradient
#' projection under the same bounds).
#' @noRd
ispline_fit_knots <- function(interior, order, mode, grm, Y, y, hetero, tol, max_iter) {
  n  <- length(y)
  x_u <- if (mode == "eap") grm$eap else nodes
  Iu <- ispline_basis(x_u, interior, order); P <- ncol(Iu)
  a <- grm$a; d <- grm$d
  if (mode != "eap") {
    item_ll <- item_loglik_nodes(Y, a, d)
    log_w_item <- sweep(item_ll, 2, log_wts, "+")
    if (mode == "joint") onehot <- item_onehot(Y)
  }
  n_g <- if (hetero) 2 else 1
  unpack <- function(th) list(b0 = th[1], a = th[1 + seq_len(P)],
                              g = if (hetero) th[P + 2:3] else c(th[P + 2], 0))
  negQ <- function(th, N, Sy, Syy) {
    p  <- unpack(th)
    mu <- p$b0 - drop(Iu %*% p$a)
    ls <- p$g[1] + p$g[2] * x_u; v <- exp(2 * ls)
    RSS <- Syy - 2 * mu * Sy + mu^2 * N
    val <- sum(N * ls + RSS / (2 * v))
    r   <- (Sy - mu * N) / v
    dls <- N - RSS / v
    attr(val, "gradient") <- c(-sum(r), drop(crossprod(Iu, r)), sum(dls),
                               if (hetero) sum(x_u * dls))
    val
  }
  th <- c(1, rep(2 / P, P), log(0.8), if (hetero) 0)
  lower <- c(-Inf, rep(0, P), rep(-Inf, n_g))
  loglik_old <- -Inf; iter <- 0; status <- "converged"
  repeat {
    p <- unpack(th)
    mu_u <- p$b0 - drop(Iu %*% p$a); sigma_u <- exp(p$g[1] + p$g[2] * x_u)
    if (mode == "eap") {
      loglik <- sum(dnorm(y, mu_u, sigma_u, log = TRUE))
      N <- rep(1, n); Sy <- y; Syy <- y^2
    } else {
      log_joint <- log_w_item + dnorm(outer(y, mu_u, "-"), 0, rep(sigma_u, each = n), log = TRUE)
      loglik <- sum(row_log_sum_exp(log_joint))
      post <- row_posterior(log_joint)
      N <- colSums(post); Sy <- colSums(post * y); Syy <- colSums(post * y^2)
    }
    if (iter > 0 && (mode == "eap" || loglik - loglik_old < tol * (1 + abs(loglik)))) break
    if (iter >= max_iter) { status <- "max_iter"; break }
    loglik_old <- loglik; iter <- iter + 1
    if (mode == "joint") {
      it <- update_items(post, Y, a, d, onehot); a <- it$a; d <- it$d
      item_ll <- item_loglik_nodes(Y, a, d)
      log_w_item <- sweep(item_ll, 2, log_wts, "+")
    }
    opt <- optim(th, function(t) negQ(t, N, Sy, Syy),
                 function(t) attr(negQ(t, N, Sy, Syy), "gradient"),
                 method = "L-BFGS-B", lower = lower, control = list(factr = 1e3, maxit = 1000))
    th <- opt$par
  }
  p <- unpack(th)
  list(b0 = p$b0, a_coef = p$a, g = p$g, loglik = loglik, iter = iter, status = status,
       a = a, d = d, n_active = sum(p$a > 1e-8), interior = interior)
}

#' The I-spline branch of fit_ms() (curve = "ispline").
#' @noRd
fit_ms_ispline <- function(grm, Y, sqi, mode, n_knots = "AIC", max_knots = 5, order = 2,
                           criterion = "AIC", hetero = TRUE, tol = 1e-9, max_iter = 2000,
                           family = "gaussian") {
  if (family != "gaussian")
    return(fit_ms_ispline_count(grm, Y, sqi, mode, n_knots, max_knots, order, criterion,
                                family, tol, max_iter))
  n_obs  <- length(sqi)
  center <- mean(sqi); scale <- sd(sqi)
  y <- (sqi - center) / scale
  cand <- if (identical(n_knots, "AIC")) 0:max_knots else n_knots
  fits <- lapply(cand, function(k)
    ispline_fit_knots(ispline_knots(k), order, mode, grm, Y, y, hetero, tol, max_iter))
  dev <- sapply(fits, function(f) -2 * (f$loglik - n_obs * log(scale)))
  df  <- sapply(fits, function(f) 1 + f$n_active + if (hetero) 2 else 1)
  path <- data.frame(n_knots = cand, deviance = dev, df = df, aic = dev + 2 * df,
                     bic = dev + log(n_obs) * df,
                     em_iter = sapply(fits, `[[`, "iter"))
  pick <- if (length(cand) == 1) 1 else
    which.min(if (criterion == "BIC") path$bic else path$aic)
  f <- fits[[pick]]

  beta  <- center + scale * ispline_beta(f$b0, f$a_coef)
  g_hat <- c(f$g[1] + log(scale), f$g[2])
  act <- if (hetero) ispline_act(beta, f$interior, order, NA, g_hat) else
    ispline_act(beta, f$interior, order, exp(g_hat[1]))
  list(a = f$a, d = f$d, act = act, mode = mode, curve = "ispline",
       converged = grm$converged && f$status == "converged",
       iter = f$iter, iter_aic = NA_integer_, iter_bic = NA_integer_, df = df[pick],
       n_knots = cand[pick], knots = f$interior, order = order,
       coef = c(beta0 = center + scale * f$b0, setNames(scale * f$a_coef, paste0("a", seq_along(f$a_coef)))),
       stop_reason = f$status, path = path,
       item_ll = if (mode == "eap") NULL else item_loglik_nodes(Y, f$a, f$d),
       R_coef = NULL, eta_cov = if (mode == "eap") grm$eap else NULL, family = "gaussian")
}

#' fit_ms_ispline() for a count family: curve on the log scale.
#' @noRd
fit_ms_ispline_count <- function(grm, Y, y, mode, n_knots, max_knots, order, criterion,
                                 family, tol, max_iter) {
  n_obs <- length(y)
  cand <- if (identical(n_knots, "AIC")) 0:max_knots else n_knots
  fits <- lapply(cand, function(k)
    ispline_fit_knots_count(ispline_knots(k), order, mode, grm, Y, y, family, tol, max_iter))
  dev <- sapply(fits, function(f) -2 * f$loglik)
  df  <- sapply(fits, function(f) 1 + f$n_active + (family == "negbin"))
  path <- data.frame(n_knots = cand, deviance = dev, df = df, aic = dev + 2 * df,
                     bic = dev + log(n_obs) * df, em_iter = sapply(fits, `[[`, "iter"))
  pick <- if (length(cand) == 1) 1 else
    which.min(if (criterion == "BIC") path$bic else path$aic)
  f <- fits[[pick]]
  act <- as_count_act(ispline_act(ispline_beta(f$b0, f$a_coef), f$interior, order, NA),
                      family, f$theta)
  list(a = f$a, d = f$d, act = act, mode = mode, curve = "ispline", family = family,
       converged = grm$converged && f$status == "converged",
       iter = f$iter, iter_aic = NA_integer_, iter_bic = NA_integer_, df = df[pick],
       n_knots = cand[pick], knots = f$interior, order = order,
       coef = c(beta0 = f$b0, setNames(f$a_coef, paste0("a", seq_along(f$a_coef)))),
       stop_reason = f$status, path = path,
       item_ll = if (mode == "eap") NULL else item_loglik_nodes(Y, f$a, f$d),
       R_coef = NULL, eta_cov = if (mode == "eap") grm$eap else NULL)
}

#' Normal approximation of the I-spline coefficients: inverse observed
#' information of (beta0, active a) given items and (g0, g1), with the zero
#' coefficients treated as fixed (Ramsay 1988, sec. 4.1), mapped to the
#' B-spline coefficients.
#' @noRd
ispline_vcov_beta <- function(fit, sqi) {
  act <- fit$act; P <- length(act$beta) - 1
  a_full <- -diff(act$beta); active <- which(a_full > 1e-8 * max(1, abs(act$beta)))
  Tm <- cbind(1, -sapply(active, function(m) as.numeric(seq_len(P + 1) > m)))
  th0 <- c(act$beta[1], a_full[active])
  build <- function(th) act_with(act, drop(Tm %*% th), act$sigma, act$log_sigma)
  negll <- if (fit$mode == "eap") function(th) {
    a2 <- build(th); -sum(dnorm(sqi, mu_sqi(fit$eta_cov, a2), sigma_sqi(fit$eta_cov, a2), log = TRUE))
  } else function(th) -sum(row_log_sum_exp(log_joint_nodes(fit$item_ll, sqi, build(th))))
  V <- solve(optimHess(th0, negll))
  Tm %*% V %*% t(Tm)
}
