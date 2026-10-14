# ---- Monotonic likelihood boosting, three ways ---------------------------------

#' Fitting and selection steps of one boosting iteration (paper eq. 4-9), on
#' "units" u with basis rows Bu (via the cumulative group sums S1, S2), weights
#' N_u and weighted sums Sy_u, Syy_u of the (standardised) SQI. For MS-Joint and
#' MS-Fixed the units are the 61 nodes with posterior-weighted totals; for
#' MS-EAP they are the persons at their EAP scores (N = 1, Sy = y, Syy = y^2).
#' Each unit is weighted by 1 / sigma_u^2, the Fisher weight of a normal mean
#' (sigma_u = residual SD at the unit; with a constant SD this is the paper's
#' update, (B'NB + sigma^2 Lambda)^-1 B' resid).
#' Returns the best split among those that keep beta non-increasing, or NULL.
#' @noRd
select_split <- function(S1, S2, N, Sy, Syy, mu, sigma_u, lambda, shrinkage) {
  M <- ncol(S1) + 1
  v <- sigma_u^2                                     # residual variance per unit
  resid <- (Sy - N * mu) / v
  best  <- NULL
  for (r in 1:(M - 1)) {
    Br  <- cbind(1, S1[, r], S2[, r])
    A   <- solve(crossprod(Br, Br * (N / v)) + lambda * diag(c(0, r, M - r)))
    gam <- shrinkage * drop(A %*% crossprod(Br, resid))
    if (gam[3] >= gam[2]) next                       # would break monotonicity
    new_mu <- mu + drop(Br %*% gam)
    dev    <- sum((Syy - 2 * new_mu * Sy + new_mu^2 * N) / v)
    if (is.null(best) || dev < best$dev)
      best <- list(r = r, gam = gam, Br = Br, A = A, dev = dev)
  }
  best
}

#' Fit the monotone curve
#'
#' Monotone B-spline equation of the digital measure by monotonic likelihood
#' boosting (Leitenstorfer and Tutz 2007, GMBBoost; the default) or by
#' Ramsay's (1988) monotone regression splines (`curve = "ispline"`), in one of
#' three ways (the boosting steps, tuning values and stopping rule are the same):
#'
#'  "fixed"  MS-Fixed. Items held at the GRM estimates. Each iteration: E-step
#'           (posterior p(eta_q | C_i, SQI_i) on the nodes), sigma update, one
#'           boosting step on the posterior-weighted node totals.
#'  "joint"  MS-Joint. As "fixed", but each iteration also takes an M-step for
#'           the item parameters from the same posterior (expected category
#'           counts per node), so the items and the curve are estimated together
#'           and the latent scale is defined by the PSQI and the SQI jointly.
#'           Items start at the GRM estimates. The item M-step is taken every
#'           `item_every` iterations (default 5); skipping it in between keeps
#'           every iteration a generalised EM step, and the items settle long
#'           before the AIC stop.
#'  "eap"    MS-EAP. The EAP scores from the PSQI-only GRM are treated as the
#'           covariate: the paper's algorithm exactly as published (observed
#'           covariate, normal response).
#'
#' Residual SD: with hetero = TRUE (default) log sigma(eta) = g0 + g1 eta. The
#' sigma step then updates (g0, g1) by Newton steps on the posterior-expected
#' (MS-EAP: observed) log-likelihood given the curve, units are weighted by
#' 1 / sigma_u^2 in the fitting step, and (g0, g1) are not counted in tr(H_l)
#' (their number is fixed, so the AIC minimum does not move). With
#' hetero = FALSE sigma is constant, with its closed-form update.
#'
#' Stopping: AIC(l) = Dev_l + 2 tr(H_l), H_l = I - (I - M_l)...(I - M_0) the
#' approximate hat matrix of the person-level fitted values (paper eq. 11).
#' Dev_l is -2 log L(C, SQI) with eta integrated out ("fixed", "joint") or
#' -2 log L(SQI | EAP) ("eap"). For "joint" the 28 item parameters are not
#' counted in tr(H_l); their number is the same at every iteration, so the AIC
#' minimum does not move. Search stops when the AIC has not improved for
#' `patience` iterations, no monotone split exists, or max_iter.
#'
#' @param grm output of [fit_grm()] on the same items
#' @param Y item responses (0..3, NA allowed), n x J
#' @param sqi the digital measure, length n (non-increasing in eta)
#' @param mode "fixed" (MS-Fixed), "joint" (MS-Joint) or "eap" (MS-EAP)
#' @param lambda ridge penalty of the fitting step (default 20)
#' @param shrinkage step length nu (default 0.1)
#' @param max_iter maximum number of boosting iterations
#' @param patience stop when the AIC has not improved for this many iterations
#' @param criterion "AIC" (default) or "BIC" for the stopping iteration
#' @param item_every MS-Joint: item M-step every this many iterations
#' @param hetero TRUE (default) for log sigma(eta) = g0 + g1 eta, FALSE for a
#'   constant residual SD (sensitivity analysis)
#' @param curve "boost" (default): monotonic likelihood boosting of the
#'   24-coefficient B-spline, as above. "ispline": Ramsay's (1988) monotone
#'   regression spline, \eqn{\mu(\eta) = \beta_0 - \sum_m a_m I_m(\eta)},
#'   \eqn{a_m \ge 0}, with I-splines of order `order` on `n_knots` interior
#'   knots at quantiles of N(0, 1) (see [ispline_knots()]), fitted by maximum
#'   likelihood under the bounds (EM over the nodes for "fixed" and "joint",
#'   direct for "eap"). The three ways of handling the latent trait are the
#'   same. `lambda`, `shrinkage`, `patience` and `item_every` are then unused,
#'   and `max_iter` caps the EM iterations.
#' @param family distribution of the digital measure given the latent trait.
#'   "gaussian" (default): a continuous measure (actigraphy index, gait speed,
#'   distance), normal with the mean curve and a constant or log-linear SD.
#'   "poisson" or "negbin" (negative binomial, Var = mu + mu^2 / theta): a
#'   count measure (daily steps, events, deaths); the spline is then fitted on
#'   the log scale, log E[y | eta] = f(eta) non-increasing, by the GLM version
#'   of the boosting (Leitenstorfer and Tutz 2007: penalised Fisher-scoring
#'   steps with the Fisher weights of the family, deviance-based split choice
#'   and AIC stop) or by the I-spline ML fit. theta is estimated by ML along
#'   the way. `hetero` applies to "gaussian" only. For count families the
#'   approximate hat matrix is accumulated on the units (nodes or persons)
#'   with the current Fisher weights, and [ms_vcov()] is not available; use
#'   [ms_bootstrap()]. "poisson" suits small counts (deaths, events). A large
#'   Poisson count (mean in the thousands) locates the latent trait more
#'   finely than the 61 quadrature nodes resolve, which biases the curve;
#'   fit_ms() warns when this happens. Large counts such as daily steps are
#'   overdispersed anyway and should use "negbin".
#' @param n_knots `curve = "ispline"`: number of interior knots, or "AIC"
#'   (default) to choose among 0 to `max_knots` by `criterion`
#' @param max_knots largest number of interior knots tried
#' @param order I-spline order k (default 2, piecewise quadratic I-splines, so
#'   the curve is a quadratic spline with continuous slope)
#' @return list(a, d, act, mode, converged, iter, iter_aic, iter_bic, df,
#'   stop_reason, path, item_ll (not for "eap"), R_coef (coefficient operator
#'   at the AIC stop, used by [ms_vcov()]), eta_cov (EAP scores, "eap" only),
#'   curve). For `curve = "ispline"`, `iter` is the number of EM iterations,
#'   `path` the knot selection table, `R_coef` NULL, and the result also has
#'   n_knots, knots, order and coef (beta0 and the I-spline coefficients).
#' @references Leitenstorfer F, Tutz G (2007). Generalized monotonic regression
#'   based on B-splines with an application to air pollution data.
#'   Biostatistics 8, 654-673.
#'
#'   Ramsay JO (1988). Monotone regression splines in action. Statistical
#'   Science 3, 425-461.
#' @export
fit_ms <- function(grm, Y, sqi, mode = c("fixed", "joint", "eap"),
                   lambda = 20, shrinkage = 0.1, max_iter = 10000,
                   patience = 100, criterion = "AIC", item_every = 5,
                   hetero = TRUE, curve = c("boost", "ispline"),
                   n_knots = "AIC", max_knots = 5, order = 2,
                   family = c("gaussian", "poisson", "negbin")) {
  mode   <- match.arg(mode)
  curve  <- match.arg(curve)
  family <- match.arg(family)
  Y      <- as.matrix(Y)
  if (anyNA(sqi)) stop("The digital measure must have no missing values.")
  if (family != "gaussian") check_count(sqi)
  if (curve == "ispline")
    return(fit_ms_ispline(grm, Y, sqi, mode, n_knots = n_knots, max_knots = max_knots,
                          order = order, criterion = criterion, hetero = hetero,
                          max_iter = min(max_iter, 2000), family = family))
  if (family != "gaussian")
    return(fit_ms_count_boost(grm, Y, sqi, mode, family, lambda, shrinkage, max_iter,
                              patience, criterion, item_every))
  n_obs  <- length(sqi)
  center <- mean(sqi); scale <- sd(sqi)
  y      <- (sqi - center) / scale
  a <- grm$a; d <- grm$d

  Bu <- if (mode == "eap") bspline_basis(grm$eap) else bspline_basis(nodes)
  M  <- ncol(Bu)
  S1 <- t(apply(Bu, 1, cumsum))[, 1:(M - 1)]
  S2 <- rowSums(Bu) - S1

  if (mode != "eap") {
    item_ll    <- item_loglik_nodes(Y, a, d)
    log_w_item <- sweep(item_ll, 2, log_wts, "+")
    G <- diag(n_obs) - 1 / n_obs                       # I - H_0 (intercept only)
  } else {
    # For EAP every M_l = Bu K_l Bu' D_l^-1 has column space in span(Bu), so
    # H_l = Bu R_l with an M x n matrix R_l (exact, no n x n matrix):
    # H_l = H_{l-1} + M_l (I - H_{l-1}). Since 1_n = Bu 1_M, H_0 = 11'/n = Bu R_0.
    R_hat <- matrix(1 / n_obs, M, n_obs)
  }
  # Coefficient operator: beta0 + beta ~= R_coef y, linearised with the
  # posterior weights (and, for MS-Joint, the items) held fixed. Kept at the
  # AIC minimum for the interval estimates (ms_vcov()). For "eap" it is R_hat.
  R_coef <- matrix(1 / n_obs, M, n_obs); R_best <- R_coef
  if (mode == "joint") onehot <- item_onehot(Y)

  x_u   <- if (mode == "eap") grm$eap else nodes      # covariate value of each unit
  beta0 <- 0; beta <- rep(0, M); g <- c(0, 0)          # log sigma = g0 + g1 x, standardised scale
  path <- data.frame(iter = 0:max_iter, deviance = NA_real_, df = NA_real_)
  beta_path <- matrix(NA_real_, max_iter + 1, M + 3)
  item_path <- if (mode == "joint") array(NA_real_, c(max_iter + 1, ncol(Y), 4)) else NULL
  best_aic <- Inf; best_l <- 0; stop_reason <- "max_iter"

  for (l in 0:max_iter) {
    mu_u    <- beta0 + drop(Bu %*% beta)
    sigma_u <- exp(g[1] + g[2] * x_u)

    # E-step (or, for EAP, the residuals at the plug-in scores)
    if (mode == "eap") {
      loglik <- sum(dnorm(y, mu_u, sigma_u, log = TRUE)) - n_obs * log(scale)
      df_l   <- sum(Bu * t(R_hat))                      # tr(Bu R)
      N <- rep(1, n_obs); Sy <- y; Syy <- y^2
    } else {
      log_joint <- log_w_item + dnorm(outer(y, mu_u, "-"), 0,
                                      rep(sigma_u, each = n_obs), log = TRUE)
      loglik    <- sum(row_log_sum_exp(log_joint)) - n_obs * log(scale)
      post      <- row_posterior(log_joint)
      df_l      <- n_obs - sum(diag(G))
      N <- colSums(post); Sy <- colSums(post * y); Syy <- colSums(post * y^2)
    }
    path$deviance[l + 1] <- -2 * loglik
    path$df[l + 1]       <- df_l
    beta_path[l + 1, ]   <- c(beta0, beta, g)
    if (mode == "joint") item_path[l + 1, , ] <- cbind(a, d)

    aic <- path$deviance[l + 1] + 2 * df_l
    if (aic < best_aic) { best_aic <- aic; best_l <- l; R_best <- R_coef }
    if (l - best_l >= patience) { stop_reason <- "patience"; break }
    if (l == max_iter) break

    # residual SD given the current curve
    RSS_u <- Syy - 2 * mu_u * Sy + mu_u^2 * N
    if (hetero) g <- update_log_sigma(g, x_u, N, RSS_u) else g <- c(log(sum(RSS_u) / n_obs) / 2, 0)
    sigma_u <- exp(g[1] + g[2] * x_u)

    # MS-Joint: item M-step from the same posterior
    if (mode == "joint" && l %% item_every == 0) {
      it <- update_items(post, Y, a, d, onehot)
      a <- it$a; d <- it$d
      item_ll    <- item_loglik_nodes(Y, a, d)
      log_w_item <- sweep(item_ll, 2, log_wts, "+")
    }

    best <- select_split(S1, S2, N, Sy, Syy, mu_u, sigma_u, lambda, shrinkage)
    if (is.null(best)) { stop_reason <- "no monotone update"; break }

    beta0 <- beta0 + best$gam[1]
    beta  <- beta + ifelse(seq_len(M) <= best$r, best$gam[2], best$gam[3])
    if (mode == "eap") {
      # M_l = Bu K W', K = nu T A T' (T = grouping, M x 3), W = Bu / sigma_i^2;
      # R <- R + K W' (I - Bu R)
      Tm <- cbind(1, seq_len(M) <= best$r, seq_len(M) > best$r)
      K  <- shrinkage * Tm %*% best$A %*% t(Tm)
      KW <- K %*% t(Bu / sigma_u^2)                     # M x n
      R_hat <- R_hat + KW - (KW %*% Bu) %*% R_hat
      R_coef <- R_hat
    } else {
      # M_l = nu Z A Z_w', Z_w = posterior weights x (grouped basis / sigma_q^2)
      Z  <- post %*% best$Br
      Zw <- post %*% (best$Br / sigma_u^2)
      G  <- G - Z %*% (shrinkage * best$A %*% crossprod(Zw, G))
      # beta0 + beta += K Bu' V^-1 (P' y - N * Bu R y), K = nu T A T'
      Tm <- cbind(1, seq_len(M) <= best$r, seq_len(M) > best$r)
      K  <- shrinkage * Tm %*% best$A %*% t(Tm)
      R_coef <- R_coef + K %*% (t(Bu / sigma_u^2) %*% (t(post) - N * (Bu %*% R_coef)))
    }
  }

  path <- path[!is.na(path$deviance), ]
  path$aic <- path$deviance + 2 * path$df
  path$bic <- path$deviance + log(n_obs) * path$df
  l_aic <- path$iter[which.min(path$aic)]
  l_bic <- path$iter[which.min(path$bic)]
  l_use <- if (criterion == "BIC") l_bic else l_aic

  # Curve at the stop, back on the SQI scale (B-splines sum to one, so beta0
  # folds into the coefficients); items at the stop (MS-Joint)
  bp   <- beta_path[l_use + 1, ]
  beta <- center + scale * (bp[1] + bp[2:(M + 1)])
  if (mode == "joint") { a <- item_path[l_use + 1, , 1]; d <- item_path[l_use + 1, , 2:4] }

  # Residual SD at the stop: ML given the curve (and items); on the SQI scale g0
  # gains log SD(SQI), g1 is unchanged
  g_start <- c(bp[M + 2] + log(scale), bp[M + 3])
  if (mode == "eap") {
    res_eap <- sqi - drop(bspline_basis(grm$eap) %*% beta)
    negll <- function(gg) -sum(dnorm(res_eap, 0, exp(gg[1] + gg[2] * grm$eap), log = TRUE))
    item_ll <- NULL
  } else {
    item_ll <- item_loglik_nodes(Y, a, d)
    negll <- function(gg) -sum(row_log_sum_exp(
      log_joint_nodes(item_ll, sqi, monospline_act(beta, NA, gg))))
  }
  if (hetero) {
    g_hat <- optim(g_start, negll, method = "BFGS", control = list(reltol = 1e-12))$par
    act   <- monospline_act(beta, NA, g_hat)
  } else {
    ls  <- optimize(function(v) negll(c(v, 0)), g_start[1] + c(-1, 1), tol = 1e-8)$minimum
    act <- monospline_act(beta, exp(ls))
  }

  list(a = a, d = d, act = act, mode = mode,
       converged = grm$converged && stop_reason != "max_iter",
       iter = l_use, iter_aic = l_aic, iter_bic = l_bic, df = path$df[l_use + 1],
       stop_reason = stop_reason, path = path, item_ll = item_ll,
       R_coef = if (criterion == "AIC") R_best else NULL,
       eta_cov = if (mode == "eap") grm$eap else NULL, curve = "boost", family = "gaussian")
}

#' All three methods on one data set
#'
#' Fits MS-Joint, MS-EAP and MS-Fixed, sharing one PRO-only GRM.
#' @param Y item responses
#' @param sqi the digital measure
#' @param grm output of [fit_grm()] (fitted if not given)
#' @param ... passed to [fit_ms()]
#' @return named list of three [fit_ms()] results
#' @export
fit_three <- function(Y, sqi, grm = fit_grm(Y), ...) {
  list(`MS-Joint` = fit_ms(grm, Y, sqi, "joint", ...),
       `MS-EAP`   = fit_ms(grm, Y, sqi, "eap", ...),
       `MS-Fixed` = fit_ms(grm, Y, sqi, "fixed", ...))
}

#' Joint deviance of items and digital measure
#'
#' \eqn{-2 \log L(C, y)} with \eqn{\eta} integrated out, for any (a, d, act).
#' @param Y item responses
#' @param sqi the digital measure
#' @param a,d item slopes and intercepts
#' @param act curve object ([monospline_act()])
#' @export
neg2ll_joint <- function(Y, sqi, a, d, act)
  -2 * sum(row_log_sum_exp(log_joint_nodes(item_loglik_nodes(Y, a, d), sqi, act)))


