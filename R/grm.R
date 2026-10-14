# ---- Graded response model building blocks ------------------------------------
# Items are scored 0..3 (four ordered categories). Missing responses (NA) are
# allowed: a missing item contributes nothing to a person's likelihood or to
# the item M-step of MS-Joint; mirt treats them the same way in the GRM.

#' Category probabilities and expected score of one GRM item
#'
#' `grm_probs()` gives \eqn{P(C = 0..3 \mid \eta)} for one item with
#' \eqn{\mathrm{logit}\, P(C \ge k \mid \eta) = a \eta + d_k};
#' `grm_expected()` gives \eqn{E[C \mid \eta]}.
#' @param eta latent values
#' @param a slope
#' @param d intercepts \eqn{d_1 > d_2 > d_3}
#' @return `grm_probs()`: length(eta) x 4 matrix; `grm_expected()`: vector
#' @export
grm_probs <- function(eta, a, d) {
  p_ge_k     <- matrix(sapply(d, function(dk) plogis(a * eta + dk)), ncol = 3)
  p_at_least <- cbind(1, p_ge_k, 0)
  p_at_least[, 1:4, drop = FALSE] - p_at_least[, 2:5, drop = FALSE]
}

#' @rdname grm_probs
#' @export
grm_expected <- function(eta, a, d) drop(grm_probs(eta, a, d) %*% 0:3)

#' One ordinal response (0..3) per row of a probability matrix
#' @param probs n x 4 matrix of category probabilities
#' @return integer vector
#' @export
draw_categories <- function(probs) {
  cum_probs <- t(apply(probs, 1, cumsum))[, 1:3, drop = FALSE]
  rowSums(runif(nrow(probs)) > cum_probs)
}

#' log P(C_i | eta_q) for all persons and nodes (n x Q); missing items skipped.
#' @noRd
item_loglik_nodes <- function(Y, a, d) {
  out <- matrix(0, nrow(Y), length(nodes))
  for (j in seq_len(ncol(Y))) {
    log_probs <- log(pmax(grm_probs(nodes, a[j], d[j, ]), 1e-300))   # Q x 4
    seen <- !is.na(Y[, j])
    out[seen, ] <- out[seen, ] + t(log_probs[, Y[seen, j] + 1, drop = FALSE])
  }
  out
}

#' Graded response model on the PRO items alone
#'
#' Fits a unidimensional logistic GRM with mirt (`itemtype = "graded"`,
#' defaults) and orients the latent trait so that the sum of the slopes is
#' positive (higher \eqn{\eta} = higher item scores). For an item whose top
#' categories nobody used, mirt estimates no threshold for them; the missing
#' intercept is set 20 logits below the previous one, so the category has
#' probability close to zero everywhere, and the threshold should be reported
#' as not estimable.
#'
#' @param Y matrix or data frame of item responses scored 0..3 (NA allowed),
#'   one column per item; column names are used as item names
#' @param SE passed to [mirt::mirt()]
#' @return list(a, d (items x 3), eap, converged, orientation, items, mirt)
#' @export
fit_grm <- function(Y, SE = FALSE) {
  Y <- as.matrix(Y)
  if (any(!is.na(Y) & !(Y %in% 0:3))) stop("Items must be scored 0, 1, 2 or 3 (or NA).")
  items <- colnames(Y)
  if (is.null(items)) items <- colnames(Y) <- paste0("C", seq_len(ncol(Y)))
  fit <- mirt::mirt(as.data.frame(Y), model = 1, itemtype = "graded", SE = SE,
                    verbose = FALSE)
  cf  <- mirt::coef(fit, simplify = TRUE)$items[items, , drop = FALSE]
  a   <- cf[, "a1"]
  eap <- mirt::fscores(fit, method = "EAP")[, "F1"]
  orientation <- if (sum(a) < 0) -1 else 1
  dcols <- paste0("d", 1:3)
  d <- matrix(NA_real_, length(items), 3)
  have <- intersect(dcols, colnames(cf))
  d[, match(have, dcols)] <- cf[, have]
  for (k in 2:3) {
    gone <- is.na(d[, k])
    d[gone, k] <- d[gone, k - 1] - 20
  }
  list(a = unname(orientation * a), d = unname(d), eap = orientation * eap,
       converged = mirt::extract.mirt(fit, "converged"),
       orientation = orientation, items = items, mirt = fit)
}

# Item parameters on an unconstrained scale for the MS-Joint M-step:
# theta_j = (a_j, d_j1, log(d_j1 - d_j2), log(d_j2 - d_j3)); the log gaps keep
# the intercepts ordered.
item_pack   <- function(a, d) c(a, d[1], log(d[1] - d[2]), log(d[2] - d[3]))
item_unpack <- function(th) list(a = th[1],
                                 d = th[2] - c(0, exp(th[3]), exp(th[3]) + exp(th[4])))

#' Expected complete-data log-likelihood of one item, sum_q sum_k R_qk log P_k(eta_q),
#' and its gradient in theta_j. R = Q x 4 posterior-expected category counts.
#' @noRd
item_m_objective <- function(th, R) {
  p  <- item_unpack(th)
  ps <- sapply(p$d, function(dk) plogis(p$a * nodes + dk))       # Q x 3, P(C >= k)
  P  <- pmax(cbind(1, ps) - cbind(ps, 0), 1e-300)                  # Q x 4
  W  <- ps * (1 - ps)                                              # dP*/dlinear
  RP <- R / P
  Wf <- cbind(0, W, 0)                                             # W_0 .. W_4
  g_a <- sum(RP * (Wf[, 1:4] - Wf[, 2:5]) * nodes)
  g_d <- colSums((RP[, 2:4] - RP[, 1:3]) * W)                      # d/d d_1..d_3
  e3 <- exp(th[3]); e4 <- exp(th[4])
  list(value = sum(R * log(P)),
       grad  = c(g_a, sum(g_d), -e3 * (g_d[2] + g_d[3]), -e4 * g_d[3]))
}

#' Category indicators per item (n x 4), zero for a missing response.
#' @noRd
item_onehot <- function(Y)
  lapply(seq_len(ncol(Y)), function(j) { o <- outer(Y[, j], 0:3, "==") * 1; o[is.na(o)] <- 0; o })

#' M-step for all items given the posterior over the nodes (n x Q): each item
#' maximises its expected complete-data log-likelihood (BFGS, started from the
#' current values).
#' @noRd
update_items <- function(post, Y, a, d, onehot) {
  for (j in seq_len(ncol(Y))) {
    R  <- crossprod(post, onehot[[j]])                              # Q x 4
    cache <- NULL                                  # value and gradient in one pass
    obj <- function(t) {
      if (is.null(cache) || !identical(cache$t, t)) cache <<- c(list(t = t), item_m_objective(t, R))
      cache
    }
    th <- optim(item_pack(a[j], d[j, ]), function(t) -obj(t)$value,
                function(t) -obj(t)$grad,
                method = "BFGS", control = list(maxit = 50, reltol = 1e-10))$par
    p <- item_unpack(th); a[j] <- p$a; d[j, ] <- p$d
  }
  list(a = a, d = d)
}
