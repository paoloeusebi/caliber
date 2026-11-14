# ---- Equipercentile linking (comparator) -----------------------------------------

#' Equipercentile linking of a continuous measure to a discrete score
#'
#' The standard equipercentile link (Kolen & Brennan, 2014, ch. 2-3) from a
#' continuous digital measure Y to a discrete score X (0..`max_score`, a PRO
#' total or one item), single-group design: a value s of Y is mapped to the
#' score with the same percentile rank. X is presmoothed by a log-linear model
#' that preserves its first `presmooth` moments (Holland & Thayer, 2000;
#' `presmooth = 0` for no smoothing) and continuized by the usual
#' percentile-rank convention (uniform on [x - 0.5, x + 0.5]). Y is continuized
#' by linear interpolation of its empirical distribution function between
#' observed values. Because high Y means a low score in the sleep example,
#' `direction = -1` matches the upper tail of Y to the lower tail of X.
#'
#' The link is a percentile-matching function, not a conditional mean: it is
#' compared with the true crosswalk E[X | Y = s] only to show how far the two
#' differ when the instruments are not perfectly related.
#' @param x integer scores 0..max_score
#' @param y continuous measure
#' @param grid values of Y at which to evaluate the link
#' @param max_score highest possible score
#' @param presmooth number of moments preserved by log-linear presmoothing (0 = none)
#' @param direction -1 if high Y goes with low X, +1 otherwise
#' @return numeric vector: the linked score at each grid value
#' @references Kolen MJ, Brennan RL (2014). Test Equating, Scaling, and Linking, 3rd ed. Springer.
#'   Holland PW, Thayer DT (2000). Univariate and bivariate loglinear models for discrete test
#'   score distributions. J Educ Behav Stat 25:133-183.
#' @export
equipercentile_link <- function(x, y, grid, max_score, presmooth = 4, direction = -1) {
  sc <- 0:max_score
  n  <- tabulate(x + 1, nbins = max_score + 1)
  f  <- if (presmooth > 0) {
    deg <- min(presmooth, sum(n > 0) - 1)
    z   <- (sc - mean(sc)) / sd(sc)
    fit <- stats::glm(n ~ stats::poly(z, deg, raw = TRUE), family = stats::poisson())
    stats::fitted(fit) / sum(stats::fitted(fit))
  } else n / sum(n)
  f  <- pmax(f, 1e-10); f <- f / sum(f)
  Fx <- cumsum(f); Fx[length(Fx)] <- 1; Fx_lo <- c(0, Fx[-length(Fx)])
  ys <- sort(y); m <- length(ys)
  p  <- stats::approx(ys, (seq_len(m) - 0.5) / m, xout = grid, rule = 2, ties = mean)$y
  p  <- ifelse(grid < ys[1], 0, ifelse(grid > ys[m], 1, p))
  if (direction < 0) p <- 1 - p
  vapply(p, function(pp) {                        # inverse of the continuized cdf of X
    k <- min(which(Fx >= pp))
    sc[k] - 0.5 + (pp - Fx_lo[k]) / f[k]
  }, 0)
}
