# ---- Fitting the curve for a count measure ---------------------------------------
# family = "poisson" or "negbin" in fit_ms(). The curve f(eta) = log E[y | eta]
# is the same non-increasing spline as for a continuous measure (boosted cubic
# B-spline, or I-splines), and the three ways of handling eta are the same.

#' Elementwise log p(y_i | log mean f_i) (MS-EAP: one unit per person).
#' @noRd
count_logdens_vec <- function(y, f, family, theta = NULL) {
  if (family == "poisson") return(y * f - exp(f) - lgamma(y + 1))
  lt <- log(theta + exp(f))
  lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) + theta * (log(theta) - lt) + y * (f - lt)
}

#' Marginal log-likelihood and posterior for a count curve on units.
#' @noRd
count_estep <- function(mode, y, f_u, family, theta, log_w_item) {
  if (mode == "eap")
    return(list(loglik = sum(count_logdens_vec(y, f_u, family, theta)),
                post = NULL, N = rep(1, length(y)), Sy = y))
  lj <- log_w_item + count_logdens(y, exp(f_u), family, theta)
  post <- row_posterior(lj)
  list(loglik = sum(row_log_sum_exp(lj)), post = post, N = colSums(post),
       Sy = colSums(post * y))
}

#' Starting value of the NB size: moments, bounded.
#' @noRd
theta_start <- function(y) {
  m <- mean(y); v <- stats::var(y)
  min(max(m^2 / max(v - m, 1e-8 * m^2), 0.05), 1e4)
}

#' Best monotone split of one GLM boosting step (Leitenstorfer and Tutz 2007):
#' penalised Fisher-scoring step on the grouped basis, chosen by deviance.
#' @noRd
select_split_count <- function(S1, S2, N, Sy, f, family, theta, lambda, shrinkage) {
  M <- ncol(S1) + 1
  score <- count_unit_score(f, N, Sy, family, theta)
  W     <- count_unit_weight(f, N, family, theta)
  best  <- NULL
  for (r in 1:(M - 1)) {
    Br  <- cbind(1, S1[, r], S2[, r])
    A   <- solve(crossprod(Br, Br * W) + lambda * diag(c(0, r, M - r)))
    gam <- shrinkage * drop(A %*% crossprod(Br, score))
    if (gam[3] >= gam[2]) next
    new_f <- f + drop(Br %*% gam)
    dev   <- -2 * sum(count_unit_ell(new_f, N, Sy, family, theta))
    if (is.null(best) || dev < best$dev)
      best <- list(r = r, gam = gam, Br = Br, A = A, dev = dev, W = W)
  }
  best
}

#' The boosting branch of fit_ms() for a count family.
#' @noRd
fit_ms_count_boost <- function(grm, Y, y, mode, family, lambda, shrinkage, max_iter,
                               patience, criterion, item_every) {
  n_obs <- length(y)
  a <- grm$a; d <- grm$d
  x_u <- if (mode == "eap") grm$eap else nodes
  Bu <- bspline_basis(x_u); M <- ncol(Bu)
  S1 <- t(apply(Bu, 1, cumsum))[, 1:(M - 1)]
  S2 <- rowSums(Bu) - S1
  log_w_item <- NULL
  if (mode != "eap") {
    item_ll <- item_loglik_nodes(Y, a, d)
    log_w_item <- sweep(item_ll, 2, log_wts, "+")
  }
  if (mode == "joint") onehot <- item_onehot(Y)
  theta <- if (family == "negbin") theta_start(y) else NULL
  beta0 <- log(mean(y)); beta <- rep(0, M)
  # Approximate hat matrix of the fitted log means on the units, H_l = Bu R_l
  # (as for MS-EAP in the Gaussian case), with the Fisher weights of each step
  R <- NULL
  path <- data.frame(iter = 0:max_iter, deviance = NA_real_, df = NA_real_)
  beta_path <- matrix(NA_real_, max_iter + 1, M + 2)
  item_path <- if (mode == "joint") array(NA_real_, c(max_iter + 1, ncol(Y), 4)) else NULL
  best_aic <- Inf; best_l <- 0; stop_reason <- "max_iter"
  for (l in 0:max_iter) {
    f_u <- beta0 + drop(Bu %*% beta)
    es  <- count_estep(mode, y, f_u, family, theta, log_w_item)
    if (is.null(R)) {
      w0 <- count_unit_weight(f_u, es$N, family, theta)
      R  <- matrix(w0 / sum(w0), M, length(w0), byrow = TRUE)
    }
    df_l <- sum(Bu * t(R)) + (family == "negbin")
    path$deviance[l + 1] <- -2 * es$loglik
    path$df[l + 1] <- df_l
    beta_path[l + 1, ] <- c(beta0, beta, if (is.null(theta)) NA else theta)
    if (mode == "joint") item_path[l + 1, , ] <- cbind(a, d)
    aic <- path$deviance[l + 1] + 2 * df_l
    if (aic < best_aic) { best_aic <- aic; best_l <- l }
    if (l - best_l >= patience) { stop_reason <- "patience"; break }
    if (l == max_iter) break

    if (family == "negbin") theta <- update_theta(theta, f_u, es$N, es$Sy, y)
    if (mode == "joint" && l %% item_every == 0) {
      it <- update_items(es$post, Y, a, d, onehot); a <- it$a; d <- it$d
      item_ll <- item_loglik_nodes(Y, a, d)
      log_w_item <- sweep(item_ll, 2, log_wts, "+")
    }
    best <- select_split_count(S1, S2, es$N, es$Sy, f_u, family, theta, lambda, shrinkage)
    if (is.null(best)) { stop_reason <- "no monotone update"; break }
    beta0 <- beta0 + best$gam[1]
    beta  <- beta + ifelse(seq_len(M) <= best$r, best$gam[2], best$gam[3])
    Tm <- cbind(1, seq_len(M) <= best$r, seq_len(M) > best$r)
    K  <- shrinkage * Tm %*% best$A %*% t(Tm)
    KW <- K %*% t(Bu * best$W)
    R  <- R + KW - (KW %*% Bu) %*% R
  }
  path <- path[!is.na(path$deviance), ]
  path$aic <- path$deviance + 2 * path$df
  path$bic <- path$deviance + log(n_obs) * path$df
  l_aic <- path$iter[which.min(path$aic)]
  l_bic <- path$iter[which.min(path$bic)]
  l_use <- if (criterion == "BIC") l_bic else l_aic
  bp   <- beta_path[l_use + 1, ]
  beta <- bp[1] + bp[2:(M + 1)]                     # B-splines sum to one
  if (mode == "joint") { a <- item_path[l_use + 1, , 1]; d <- item_path[l_use + 1, , 2:4] }
  item_ll <- if (mode == "eap") NULL else item_loglik_nodes(Y, a, d)
  theta <- if (family == "negbin") count_final_theta(bp[M + 2], mode, y, Bu, beta, item_ll) else NULL
  list(a = a, d = d, act = as_count_act(monospline_act(beta, NA), family, theta),
       mode = mode, curve = "boost", family = family,
       converged = grm$converged && stop_reason != "max_iter",
       iter = l_use, iter_aic = l_aic, iter_bic = l_bic, df = path$df[l_use + 1],
       stop_reason = stop_reason, path = path, item_ll = item_ll, R_coef = NULL,
       eta_cov = if (mode == "eap") grm$eap else NULL)
}

#' ML of the NB size given the final curve (and items), marginal over eta.
#' @noRd
count_final_theta <- function(theta, mode, y, Bu, beta, item_ll) {
  f_u <- drop(Bu %*% beta)
  negll <- if (mode == "eap") function(lt) -sum(count_logdens_vec(y, f_u, "negbin", exp(lt))) else {
    lw <- sweep(item_ll, 2, log_wts, "+")
    function(lt) -sum(row_log_sum_exp(lw + count_logdens(y, exp(f_u), "negbin", exp(lt))))
  }
  exp(optimize(negll, log(theta) + c(-3, 3), tol = 1e-8)$minimum)
}

#' One I-spline fit with fixed knots for a count measure: f = b0 - Iu a, a >= 0,
#' EM over the nodes (MS-Fixed, MS-Joint) or direct ML (MS-EAP), each M-step by
#' L-BFGS-B in (b0, a, log theta).
#' @noRd
ispline_fit_knots_count <- function(interior, order, mode, grm, Y, y, family, tol, max_iter) {
  n <- length(y)
  x_u <- if (mode == "eap") grm$eap else nodes
  Iu <- ispline_basis(x_u, interior, order); P <- ncol(Iu)
  a <- grm$a; d <- grm$d; log_w_item <- NULL
  if (mode != "eap") {
    item_ll <- item_loglik_nodes(Y, a, d)
    log_w_item <- sweep(item_ll, 2, log_wts, "+")
    if (mode == "joint") onehot <- item_onehot(Y)
  }
  nb <- family == "negbin"
  sum_lgy <- sum(lgamma(y + 1))
  negQ <- function(th, N, Sy) {
    f <- th[1] - drop(Iu %*% th[1 + seq_len(P)])
    t <- if (nb) exp(th[P + 2]) else NULL
    val <- -sum(count_unit_ell(f, N, Sy, family, t))
    sc  <- count_unit_score(f, N, Sy, family, t)
    g   <- c(sum(-sc), -drop(crossprod(Iu, -sc)))
    if (nb) {
      val <- val - (sum(lgamma(y + t)) - n * lgamma(t) + n * t * log(t))
      m <- exp(f)
      dt <- sum(digamma(y + t)) - n * digamma(t) + n * (log(t) + 1) -
        sum(N * log(t + m) + (t * N + Sy) / (t + m))
      g <- c(g, -dt * t)
    }
    attr(val, "gradient") <- g
    val
  }
  th <- c(log(mean(y)) + 0.5, rep(1 / P, P), if (nb) log(theta_start(y)))
  lower <- c(-Inf, rep(0, P), if (nb) log(1e-4))
  loglik_old <- -Inf; iter <- 0; status <- "converged"
  repeat {
    f_u <- th[1] - drop(Iu %*% th[1 + seq_len(P)])
    es <- count_estep(mode, y, f_u, family, if (nb) exp(th[P + 2]), log_w_item)
    loglik <- es$loglik
    if (iter > 0 && (mode == "eap" || loglik - loglik_old < tol * (1 + abs(loglik)))) break
    if (iter >= max_iter) { status <- "max_iter"; break }
    loglik_old <- loglik; iter <- iter + 1
    if (mode == "joint") {
      it <- update_items(es$post, Y, a, d, onehot); a <- it$a; d <- it$d
      item_ll <- item_loglik_nodes(Y, a, d)
      log_w_item <- sweep(item_ll, 2, log_wts, "+")
    }
    opt <- optim(th, function(t) negQ(t, es$N, es$Sy),
                 function(t) attr(negQ(t, es$N, es$Sy), "gradient"),
                 method = "L-BFGS-B", lower = lower, control = list(factr = 1e3, maxit = 1000))
    th <- opt$par
  }
  a_coef <- th[1 + seq_len(P)]
  list(b0 = th[1], a_coef = a_coef, theta = if (nb) exp(th[P + 2]) else NULL,
       loglik = loglik, iter = iter, status = status, a = a, d = d,
       n_active = sum(a_coef > 1e-8), interior = interior)
}
