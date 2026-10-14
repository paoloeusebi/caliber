# ---- Count outcomes: Poisson and negative binomial --------------------------------
# For a count measure (steps, events, deaths) the curve is fitted on the log
# scale: log E[y | eta] = f(eta), f a non-increasing spline, so the mean count
# exp(f(eta)) is non-increasing too. y | eta is Poisson or negative binomial
# (NB2, Var = mu + mu^2 / theta). The Gaussian family (the default) is the
# model used for continuous measures throughout the package.
#
# With eta on the quadrature nodes, the expected log-likelihood needs only the
# posterior-weighted unit totals N_u = sum_i p_iu and Sy_u = sum_i p_iu y_i,
# plus terms that do not depend on the curve:
#   Poisson  sum_u [Sy_u f_u - N_u e^f_u] - sum_i log y_i!
#   NB       sum_u [Sy_u f_u - (theta N_u + Sy_u) log(theta + e^f_u)]
#            + sum_i [lgamma(y_i + theta) - log y_i!] - n lgamma(theta) + n theta log(theta)

count_families <- c("poisson", "negbin")

#' Family of a curve object ("gaussian" unless set).
#' @noRd
act_family <- function(act) if (is.null(act$family)) "gaussian" else act$family

#' log p(y_i | mean mu_q), length(y) x length(mu). Written with lgamma so that
#' non-integer y (a crosswalk grid) gives a smooth value instead of -Inf.
#' @noRd
count_logdens <- function(y, mu, family, theta = NULL) {
  ly <- outer(y, log(mu))
  if (family == "poisson")
    return(sweep(ly, 2, mu, "-") - lgamma(y + 1))
  lt <- log(theta + mu)
  lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
    rep(theta * (log(theta) - lt), each = length(y)) + ly - outer(y, lt)
}

#' Expected log-likelihood per unit (the curve-dependent part), its score and
#' Fisher weight in the log mean f.
#' @noRd
count_unit_ell <- function(f, N, Sy, family, theta = NULL) {
  m <- exp(f)
  if (family == "poisson") return(Sy * f - N * m)
  Sy * f - (theta * N + Sy) * log(theta + m)
}
count_unit_score <- function(f, N, Sy, family, theta = NULL) {
  m <- exp(f)
  if (family == "poisson") Sy - N * m else theta * (Sy - N * m) / (theta + m)
}
count_unit_weight <- function(f, N, family, theta = NULL) {
  m <- exp(f)
  if (family == "poisson") N * m else N * m * theta / (theta + m)
}

#' ML update of the NB size theta given the curve, from the unit totals.
#' @noRd
update_theta <- function(theta, f, N, Sy, y) {
  n <- length(y)
  obj <- function(lt) {
    t <- exp(lt)
    sum(lgamma(y + t)) - n * lgamma(t) + n * t * log(t) - sum((t * N + Sy) * log(t + exp(f)))
  }
  exp(optimize(obj, log(theta) + c(-3, 3), maximum = TRUE, tol = 1e-8)$maximum)
}

#' Turn a spline curve object (fitted on the log scale) into a count curve.
#' @noRd
as_count_act <- function(act, family, theta = NULL) {
  lin <- act$mu
  act$mu <- function(eta) exp(lin(eta))
  act$lin <- lin; act$family <- family; act$theta <- theta
  act$log_sigma <- NULL
  m0 <- exp(lin(0))
  act$sigma <- sqrt(m0 + if (family == "negbin") m0^2 / theta else 0)
  act
}

#' Warn when y pins eta down more finely than the quadrature nodes (spacing
#' about 0.4 near 0): each posterior then sits on one node and the integral
#' over eta, and with it the curve, is badly approximated. Typical of a
#' Poisson fit to large counts (steps).
#' @noRd
warn_node_collapse <- function(fit, y) {
  if (fit$mode == "eap") return(invisible(FALSE))
  lj <- log_joint_nodes(fit$item_ll, y, fit$act)
  pmax <- stats::median(1 / rowSums(exp(lj - apply(lj, 1, max))))
  if (pmax > 0.9) {
    warning(sprintf(paste0(
      "The median posterior puts %.0f%% of its mass on one quadrature node: the ",
      "count measure is too precise for the latent-trait quadrature and the curve ",
      "is unreliable.%s"), 100 * pmax,
      if (fit$family == "poisson") " For large, overdispersed counts use family = \"negbin\"." else ""),
      call. = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

check_count <- function(y) {
  if (any(y < 0) || any(abs(y - round(y)) > 1e-8))
    stop("A count family needs non-negative integer values of the digital measure.")
  invisible(TRUE)
}
