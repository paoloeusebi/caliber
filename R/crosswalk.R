# ---- Thresholds, crosswalks, exceedance probabilities, persons -------------------

#' Item thresholds on the latent and digital scales
#'
#' \eqn{b_{jk} = -d_{jk} / a_j} and \eqn{s_{jk} = \mu(b_{jk})}, the mean of the
#' digital measure among persons sitting at the threshold.
#' @param a,d item slopes and intercepts (items x 3)
#' @param act curve object
#' @return list(b, sqi), both items x 3 matrices
#' @export
threshold_locations <- function(a, d, act) {
  b <- -d / a
  list(b = b, sqi = matrix(mu_sqi(c(b), act), nrow(b)))
}

#' Latent value with a given mean of the digital measure
#'
#' \eqn{\eta^*(s)} solves \eqn{\mu(\eta) = s} on \eqn{[-4, 4]}; NA outside the
#' curve's range.
#' @param s value of the digital measure
#' @param act curve object
#' @export
structural_eta <- function(s, act) {
  ends <- mu_sqi(act$range, act)
  if (s >= ends[1] || s <= ends[2]) return(NA_real_)
  uniroot(function(e) mu_sqi(e, act) - s, act$range, tol = 1e-10)$root
}

#' Crosswalks from the digital measure to the PRO
#'
#' `crosswalk_probabilistic()` gives \eqn{E[C_j \mid y = s]} and the expected
#' total, averaging over \eqn{p(\eta_q \mid s) \propto w_q N(s; \mu(\eta_q),
#' \sigma(\eta_q)^2)} (or the Poisson / negative binomial probability of `s`
#' for a count curve). `crosswalk()` adds the structural crosswalk
#' \eqn{E[C_j \mid \eta^*(s)]} (columns `eta_star`, `Estar_*`; NA outside the
#' curve's range) and prefixes the probabilistic columns with `E_`.
#' @param grid values of the digital measure
#' @param a,d item parameters
#' @param act curve object
#' @param items item names
#' @param total_name name of the total-score column
#' @param value_name name of the column holding `grid` (`crosswalk()` only)
#' @export
crosswalk_probabilistic <- function(grid, a, d, act, items = paste0("C", seq_along(a)),
                                    total_name = "total") {
  expected_at_nodes <- sapply(seq_along(a), function(j) grm_expected(nodes, a[j], d[j, ]))
  mu_nodes <- mu_sqi(nodes, act); sd_nodes <- sigma_sqi(nodes, act)
  fam <- act_family(act)
  ld  <- if (fam == "gaussian")
    dnorm(outer(grid, mu_nodes, "-"), 0, rep(sd_nodes, each = length(grid)), log = TRUE) else
    count_logdens(grid, mu_nodes, fam, act$theta)
  w   <- exp(sweep(ld, 2, log_wts, "+"))                 # grid x nodes
  out <- (w %*% expected_at_nodes) / rowSums(w)
  colnames(out) <- items
  out <- as.data.frame(out)
  out[[total_name]] <- rowSums(out[items])
  out
}

#' @rdname crosswalk_probabilistic
#' @export
crosswalk <- function(grid, a, d, act, items = paste0("C", seq_along(a)),
                      total_name = "total", value_name = "SQI") {
  structural <- t(sapply(grid, function(s) {
    es <- structural_eta(s, act)
    ex <- if (is.na(es)) rep(NA_real_, length(a)) else
      sapply(seq_along(a), function(j) grm_expected(es, a[j], d[j, ]))
    c(es, ex, sum(ex))
  }))
  colnames(structural) <- c("eta_star", paste0("Estar_", c(items, total_name)))
  prob <- crosswalk_probabilistic(grid, a, d, act, items, total_name)
  names(prob) <- paste0("E_", names(prob))
  out <- data.frame(grid, structural, prob)
  names(out)[1] <- value_name
  out
}

#' Exceedance probabilities given the digital measure
#'
#' \eqn{P(C_j \ge k \mid y = s) = \sum_q P(C_j \ge k \mid \eta_q)
#' p(\eta_q \mid s)}, k = 1..3, on a grid (probabilistic, as the crosswalk).
#' @inheritParams crosswalk_probabilistic
#' @return grid x (items x 3) matrix, columns item-major within k
#' @export
exceedance_probabilistic <- function(grid, a, d, act, items = paste0("C", seq_along(a))) {
  J <- length(a)
  p_nodes  <- do.call(cbind, lapply(1:3, function(k) sapply(1:J, function(j)
    plogis(a[j] * nodes + d[j, k]))))                                  # Q x 3J
  mu_nodes <- mu_sqi(nodes, act); sd_nodes <- sigma_sqi(nodes, act)
  fam <- act_family(act)
  out <- t(sapply(grid, function(s) {
    w <- exp(log_wts + if (fam == "gaussian") dnorm(s, mu_nodes, sd_nodes, log = TRUE) else
      drop(count_logdens(s, mu_nodes, fam, act$theta)))
    colSums(w * p_nodes) / sum(w)
  }))
  colnames(out) <- threshold_names(items)
  out
}

#' Floor, ceiling and approximate R-squared of a fitted curve
#'
#' The curve is constant outside \eqn{[-4, 4]}, so \eqn{\mu(4)} is its floor
#' and \eqn{\mu(-4)} its ceiling: no latent value has a mean outside them.
#' R-squared is approximate (residual SD at eta = 0 against the SD of y).
#' @param act curve object
#' @param sqi observed digital measure
#' @export
curve_ends <- function(act, sqi) {
  fl <- mu_sqi(4, act); ce <- mu_sqi(-4, act)
  data.frame(floor = fl, ceiling = ce, share_below_floor = mean(sqi < fl),
             share_above_ceiling = mean(sqi > ce),
             r2_approx = 1 - (sigma_sqi(0, act) / sd(sqi))^2)
}

#' Flags for thresholds that cannot be read as point estimates
#'
#' A threshold is flagged as sparse when fewer than `min_at_or_above` persons
#' respond at or above its category, and as at the floor when it lies within
#' `floor_tol` of the curve's floor (such thresholds cannot be told apart and
#' are reported as "below about" the floor). A category nobody used is not
#' estimable.
#' @param Y item responses
#' @param thr_sqi items x 3 matrix of thresholds on the digital scale
#' @param act curve object
#' @param floor_tol tolerance in units of the digital measure
#' @param min_at_or_above minimum count at or above the category
#' @return data.frame(item, k, n_at_or_above, not_estimable, flag_sparse, at_floor)
#' @export
threshold_flags <- function(Y, thr_sqi, act, floor_tol = 0.01, min_at_or_above = 30) {
  Y <- as.matrix(Y)
  n_ge <- c(sapply(1:3, function(k) colSums(Y >= k, na.rm = TRUE)))
  fl <- mu_sqi(4, act)
  data.frame(item = rep(colnames(Y), 3), k = rep(1:3, each = ncol(Y)),
             n_at_or_above = n_ge, not_estimable = n_ge == 0,
             flag_sparse = n_ge < min_at_or_above,
             at_floor = c(thr_sqi) - fl < floor_tol)
}

#' One plausible value of the latent trait per person
#'
#' For MS-Joint and MS-Fixed, drawn from \eqn{p(\eta \mid C_i, y_i)} on the
#' quadrature nodes and jittered by a quarter node spacing; for MS-EAP, from
#' \eqn{p(\eta \mid C_i)} with mirt (the digital measure not used).
#' @param fit a [fit_ms()] result
#' @param grm the [fit_grm()] result it used
#' @param sqi the digital measure
#' @param seed RNG seed
#' @export
plausible_values <- function(fit, grm, sqi, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (fit$mode == "eap")
    return(grm$orientation * mirt::fscores(grm$mirt, plausible.draws = 1)[, 1])
  post <- row_posterior(log_joint_nodes(fit$item_ll, sqi, fit$act))
  eta  <- apply(post, 1, function(w) sample(nodes, 1, prob = w))
  eta + rnorm(length(eta), 0, diff(range(nodes)) / length(nodes) / 4)
}
