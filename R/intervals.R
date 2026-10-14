# ---- Interval estimates: normal-approximation draws ------------------------------
# Covariance blocks of the estimates, then Monte Carlo draws from their normal
# approximation, pushed through any target (thresholds, crosswalk):
#   curve     Var(beta) = R diag(e_i^2) R' with R the coefficient operator at the
#             AIC stop (fit_ms) and e_i = y_i - E[mu(eta) | C_i] (MS-Joint,
#             MS-Fixed) or y_i - mu(EAP_i) (MS-EAP). Conditional on the stop;
#             posterior weights and items held fixed.
#   (g0, g1)  inverse Hessian of the log-likelihood given curve and items.
#   items     observed information in the packed parameters (a, d1, log gaps):
#             of C alone for MS-EAP and MS-Fixed (the PRO-only GRM), of
#             (C, y) given the curve for MS-Joint.
# The blocks are treated as independent; for MS-Fixed and MS-EAP the effect of
# the GRM estimates on the fitted curve is not propagated. Each curve draw is
# projected onto the non-increasing cone (isotonic regression).
#
# These draws condition on the boosting path. They are used for item, threshold
# and crosswalk quantities. For the curve mu(eta) itself they capture only
# about 20 to 30% of the sampling variance (the selection of the stopping
# iteration and split points is not seen), so curve bands come from the
# nonparametric bootstrap (ms_bootstrap()). See the package vignettes and
# manuscript supplement S2. Coverage of either kind of interval has not been
# validated by simulation.

#' Gradient of the marginal log-likelihood in the packed item parameters
#' (Fisher identity: the M-step gradient at the posterior).
#' @noRd
item_grad <- function(theta, Y, onehot, sqi = NULL, act = NULL) {
  J <- ncol(Y)
  par <- lapply(1:J, function(j) item_unpack(theta[4 * j - 3:0]))
  a <- sapply(par, `[[`, "a"); d <- t(sapply(par, `[[`, "d"))
  item_ll <- item_loglik_nodes(Y, a, d)
  post <- row_posterior(if (is.null(act)) sweep(item_ll, 2, log_wts, "+") else
                          log_joint_nodes(item_ll, sqi, act))
  unlist(lapply(1:J, function(j)
    item_m_objective(theta[4 * j - 3:0], crossprod(post, onehot[[j]]))$grad))
}

#' Inverse observed information of the packed item parameters (central
#' differences of item_grad).
#' @noRd
item_vcov <- function(theta, Y, sqi = NULL, act = NULL, h = 1e-4) {
  onehot <- item_onehot(Y)
  Jm <- sapply(seq_along(theta), function(k) {
    e <- replace(numeric(length(theta)), k, h)
    (item_grad(theta + e, Y, onehot, sqi, act) - item_grad(theta - e, Y, onehot, sqi, act)) / (2 * h)
  })
  solve(-(Jm + t(Jm)) / 2)
}

#' Covariance blocks of a fitted curve and its items
#'
#' Normal approximation of the estimates of one [fit_ms()] result: items
#' (observed information of the packed parameters), spline coefficients
#' (linearised boosting operator at the AIC stop, sandwich with squared
#' residuals) and (g0, g1) (inverse Hessian given the curve). Conditional on
#' the boosting path; see [ms_intervals()]. For an I-spline fit
#' (`curve = "ispline"`) the curve block is the inverse observed information
#' of the intercept and the non-zero I-spline coefficients given items and
#' (g0, g1), with coefficients estimated at zero treated as fixed (Ramsay
#' 1988), conditional on the selected knots.
#' @param fit a [fit_ms()] result (fitted with `criterion = "AIC"`)
#' @param Y,sqi the data it was fitted to
#' @return list(theta, V_items, beta, V_beta, g, V_g, hetero)
#' @export
ms_vcov <- function(fit, Y, sqi) {
  Y <- as.matrix(Y)
  act <- fit$act; hetero <- !is.null(act$log_sigma)
  if (act_family(act) != "gaussian")
    stop("ms_vcov() is implemented for the Gaussian family; for a count family use ms_bootstrap().")
  if (fit$mode == "eap") {
    e <- sqi - mu_sqi(fit$eta_cov, act)
    negll <- function(gg) -sum(dnorm(e, 0, exp(gg[1] + gg[2] * fit$eta_cov), log = TRUE))
  } else {
    p_c <- row_posterior(sweep(fit$item_ll, 2, log_wts, "+"))           # p(eta | C)
    e   <- sqi - drop(p_c %*% mu_sqi(nodes, act))
    negll <- function(gg) -sum(row_log_sum_exp(
      log_joint_nodes(fit$item_ll, sqi, act_with(act, act$beta, NA, gg))))
  }
  ispline <- identical(act$type, "ispline")
  R <- fit$R_coef
  if (!ispline && is.null(R)) stop("fit_ms() was not run with criterion = \"AIC\".")
  g <- if (hetero) act$log_sigma else c(log(act$sigma), 0)
  V_g <- if (hetero) solve(optimHess(g, negll)) else
    diag(c(1 / optimHess(g[1], function(v) negll(c(v, 0)))[1], 0))
  theta <- unlist(lapply(seq_len(ncol(Y)), function(j) item_pack(fit$a[j], fit$d[j, ])))
  V_items <- if (fit$mode == "joint") item_vcov(theta, Y, sqi, act) else item_vcov(theta, Y)
  V_beta <- if (ispline) ispline_vcov_beta(fit, sqi) else R %*% (e^2 * t(R))
  list(theta = theta, V_items = V_items, beta = act$beta, V_beta = V_beta,
       g = g, V_g = V_g, hetero = hetero)
}

#' n draws from N(mu, V), V positive semi-definite.
#' @noRd
rmvn <- function(n, mu, V) {
  ev <- eigen((V + t(V)) / 2, symmetric = TRUE)
  L  <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(mu))
  sweep(matrix(rnorm(n * length(mu)), n) %*% t(L), 2, mu, "+")
}

#' Monte Carlo intervals from the normal approximation
#'
#' Draws (items, curve coefficients, g) from the normal approximation of
#' [ms_vcov()] and pushes each draw through `target(a, d, act)`, which returns
#' a named numeric vector (thresholds, crosswalk values, ...). Every curve draw
#' is projected onto the non-increasing cone.
#'
#' Use these intervals for item parameters, latent thresholds, thresholds on
#' the digital scale within the data and the probabilistic crosswalk. Do not
#' use them for the curve \eqn{\mu(\eta)}: they condition on the selected
#' boosting path and miss most of its sampling variance; use [ms_bootstrap()].
#' @param fit a [fit_ms()] result
#' @param vc its [ms_vcov()]
#' @param target function(a, d, act) returning a named numeric vector
#' @param n_draw number of draws
#' @param level interval level
#' @param seed optional RNG seed
#' @return data.frame(param, est, se, lo, hi)
#' @export
ms_intervals <- function(fit, vc, target, n_draw = 1000, level = 0.95, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  J  <- length(vc$theta) / 4
  th <- rmvn(n_draw, vc$theta, vc$V_items)
  bt <- rmvn(n_draw, vc$beta, vc$V_beta)
  gg <- rmvn(n_draw, vc$g, vc$V_g)
  draws <- sapply(seq_len(n_draw), function(r) {
    par  <- lapply(1:J, function(j) item_unpack(th[r, 4 * j - 3:0]))
    beta <- -isoreg(-bt[r, ])$yf                       # nearest non-increasing curve
    act  <- if (vc$hetero) act_with(fit$act, beta, NA, gg[r, ]) else act_with(fit$act, beta, exp(gg[r, 1]))
    target(sapply(par, `[[`, "a"), t(sapply(par, `[[`, "d")), act)
  })
  draws <- matrix(draws, ncol = n_draw)
  est <- target(fit$a, fit$d, fit$act)
  q <- (1 - level) / 2
  data.frame(param = names(est), est = unname(est), se = apply(draws, 1, sd),
             lo = apply(draws, 1, quantile, q, na.rm = TRUE),
             hi = apply(draws, 1, quantile, 1 - q, na.rm = TRUE))
}


# ---- Nonparametric bootstrap --------------------------------------------------

#' Standard summary of one fit
#'
#' Item parameters, residual SD, mu(eta) at `eta_report`, latent thresholds
#' b_jk, thresholds on the digital scale s_jk = mu(b_jk) and the curve on
#' `curve_grid` (columns `curve_1`, ...), plus the stopping iteration and df.
#' Used as the default summary of [ms_bootstrap()].
#' @param f a [fit_ms()] result
#' @param items item names
#' @param eta_report latent values at which mu(eta) is reported
#' @param curve_grid latent grid for curve bands
#' @return named numeric vector
#' @export
fit_summary <- function(f, items = paste0("C", seq_along(f$a)), eta_report = -3:3,
                        curve_grid = seq(-3.5, 3.5, by = 0.05)) {
  thr <- threshold_locations(f$a, f$d, f$act)
  tn  <- threshold_names(items)
  c(setNames(f$a, paste0("a_", items)),
    setNames(f$d[, 1], paste0("d1_", items)), setNames(f$d[, 2], paste0("d2_", items)),
    setNames(f$d[, 3], paste0("d3_", items)),
    sigma = f$act$sigma, g1 = if (is.null(f$act$log_sigma)) 0 else f$act$log_sigma[2],
    setNames(mu_sqi(eta_report, f$act), paste0("mu(", eta_report, ")")),
    setNames(c(thr$b), paste0("b_", tn)),
    setNames(c(thr$sqi), paste0("s_", tn)),
    setNames(mu_sqi(curve_grid, f$act), paste0("curve_", seq_along(curve_grid))),
    iter = f$iter, df = f$df)
}

#' Nonparametric bootstrap of the three methods
#'
#' Resamples persons with replacement and refits every step (the PRO-only GRM,
#' the boosting with its own AIC stop, the final residual SD) for each
#' method, then applies `summary_fun` to each fit. With `keep_categories =
#' TRUE` a resample in which an item lost a response category present in the
#' data is redrawn (the GRM then has no threshold for it); the number of
#' redraws and, per threshold, whether the resample was redrawn because that
#' category was lost are returned, so the intervals are conditional on every
#' category being present.
#'
#' The bootstrap captures the variability of the selected boosting path,
#' which the normal approximation misses, but it does not remove smoothing
#' bias, is not guaranteed to be consistent where the monotone curve is flat,
#' and carries no information beyond the data. Its coverage has not been
#' validated by simulation.
#' @param Y item responses
#' @param sqi the digital measure
#' @param n_boot number of resamples
#' @param seed base seed; resample r uses `seed + 100000 + r`
#' @param summary_fun function(fit) returning a named numeric vector
#' @param methods which of [ms_methods] to refit
#' @param keep_categories redraw resamples that lose a category
#' @param cores number of worker processes (1 = sequential)
#' @param ... passed to [fit_ms()]
#' @return data.frame with columns method, rep, redrawn, lost_* and the summaries
#' @export
ms_bootstrap <- function(Y, sqi, n_boot = 200, seed = 20261003,
                         summary_fun = function(f) fit_summary(f, items = colnames(Y)),
                         methods = ms_methods, keep_categories = TRUE,
                         cores = 1, ...) {
  Y <- as.matrix(Y); n <- nrow(Y); J <- ncol(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("C", seq_len(J))
  force(sqi); force(seed); force(summary_fun); force(methods); force(keep_categories)
  tn <- threshold_names(colnames(Y))
  cats_seen <- lapply(1:J, function(j) sort(unique(na.omit(Y[, j]))))
  categories_lost <- function(Yb)
    t(sapply(1:J, function(j) 0:3 %in% cats_seen[[j]] & !(0:3 %in% Yb[, j])))
  dots <- list(...)
  boot_one <- function(r) {
    set.seed(seed + 100000 + r)
    redrawn <- -1
    lost <- setNames(numeric(length(tn)), paste0("lost_", tn))
    repeat {
      idx <- sample.int(n, n, replace = TRUE); redrawn <- redrawn + 1
      if (!keep_categories) break
      gone <- categories_lost(Y[idx, , drop = FALSE])
      if (!any(gone)) break
      lost <- pmax(lost, c(gone[, 2:4]))
    }
    Yb <- Y[idx, , drop = FALSE]; yb <- sqi[idx]
    grm <- fit_grm(Yb)
    do.call(rbind, lapply(methods, function(m) {
      f <- do.call(fit_ms, c(list(grm, Yb, yb, ms_mode[[m]]), dots))
      data.frame(method = m, rep = r, redrawn = redrawn, t(lost),
                 t(summary_fun(f)), check.names = FALSE)
    }))
  }
  safe <- function(r) tryCatch(boot_one(r), error = function(e) conditionMessage(e))
  res <- run_jobs(seq_len(n_boot), safe, cores, export = global_refs(summary_fun))
  ok <- !vapply(res, is.character, logical(1))
  if (!any(ok)) stop("all bootstrap resamples failed; first error: ", res[[1]])
  out <- do.call(rbind, res[ok])
  attr(out, "n_failed") <- sum(!ok)
  out
}

# Names of global-environment objects that fun (and global functions it calls)
# refer to, so PSOCK workers see the same objects as a sequential run.
global_refs <- function(fun) {
  ge <- globalenv(); found <- character(0); todo <- list(fun)
  while (length(todo)) {
    f <- todo[[1]]; todo <- todo[-1]
    if (!is.function(f) || !identical(environment(f), ge)) next
    nm <- setdiff(intersect(all.names(body(f)), ls(ge)), found)
    found <- c(found, nm)
    todo <- c(todo, mget(nm, envir = ge))
  }
  found
}

#' Percentile intervals from bootstrap replicates
#'
#' @param boot output of [ms_bootstrap()]
#' @param method one of [ms_methods]
#' @param params column names to summarise
#' @param level interval level
#' @return data.frame(param, se, lo, hi)
#' @export
boot_intervals <- function(boot, method, params, level = 0.95) {
  bm <- boot[boot$method == method, params, drop = FALSE]
  q <- (1 - level) / 2
  data.frame(param = params, se = apply(bm, 2, sd),
             lo = apply(bm, 2, quantile, q, na.rm = TRUE),
             hi = apply(bm, 2, quantile, 1 - q, na.rm = TRUE), row.names = NULL)
}

#' Apply fun to each job, sequentially or on a PSOCK cluster with caliber loaded.
#' @noRd
run_jobs <- function(jobs, fun, cores = 1, export = character(0)) {
  if (cores <= 1) return(lapply(jobs, fun))
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(caliber)))
  if (length(export)) parallel::clusterExport(cl, export, envir = globalenv())
  parallel::parLapplyLB(cl, jobs, fun)
}
