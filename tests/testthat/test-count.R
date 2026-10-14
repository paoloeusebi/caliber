truth <- function(eta, m0 = 8000) exp(log(m0) - 0.3 * eta - 0.03 * eta^2)  # non-increasing
sim_counts <- function(family, seed = 1, m0 = 8000) {
  d <- sleep_data(); Y <- as.matrix(d[paste0("C", 1:7)])
  eta <- caliber::psqi_sleep$eta_true
  set.seed(seed)
  m <- truth(eta, m0)
  y <- if (family == "poisson") rpois(length(m), m) else rnbinom(length(m), size = 5, mu = m)
  list(Y = Y, y = y, grm = fit_grm(Y))
}

test_that("count log-densities match dpois and dnbinom", {
  y <- c(0, 3, 10); mu <- c(0.5, 4, 12)
  expect_equal(caliber:::count_logdens(y, mu, "poisson"),
               outer(y, mu, function(a, b) dpois(a, b, log = TRUE)))
  expect_equal(caliber:::count_logdens(y, mu, "negbin", 2.5),
               outer(y, mu, function(a, b) dnbinom(a, size = 2.5, mu = b, log = TRUE)))
})

test_that("boosting fits a Poisson count curve (small counts)", {
  s <- sim_counts("poisson", m0 = 3)
  f <- expect_no_warning(fit_ms(s$grm, s$Y, s$y, "fixed", family = "poisson"))
  expect_equal(f$family, "poisson")
  expect_equal(f$act$family, "poisson")
  mu <- mu_sqi(seq(-4, 4, 0.1), f$act)
  expect_true(all(diff(mu) <= 1e-8))
  expect_lt(abs(log(mu_sqi(0, f$act) / truth(0, 3))), 0.15)
  expect_equal(sigma_sqi(0, f$act), sqrt(mu_sqi(0, f$act)))
  thr <- threshold_locations(f$a, f$d, f$act)
  expect_equal(dim(thr$sqi), c(7, 3))
  cw <- crosswalk_probabilistic(c(0, 3, 8), f$a, f$d, f$act)
  expect_true(all(diff(cw$total) < 0))
  expect_error(ms_vcov(f, s$Y, s$y), "Gaussian")
})

test_that("a Poisson fit to large counts warns that the quadrature cannot resolve it", {
  s <- sim_counts("poisson")
  expect_warning(fit_ms(s$grm, s$Y, s$y, "fixed", curve = "ispline", family = "poisson",
                        n_knots = 1), "quadrature")
})

test_that("negative binomial: I-spline and boosting fits recover the curve and the size", {
  s <- sim_counts("negbin")
  f <- expect_no_warning(fit_ms(s$grm, s$Y, s$y, "fixed", curve = "ispline",
                                family = "negbin", n_knots = 1))
  expect_true(f$converged)
  expect_gt(f$act$theta, 3.5); expect_lt(f$act$theta, 7)
  expect_lt(abs(log(mu_sqi(0, f$act) / truth(0))), 0.1)
  expect_true(all(diff(mu_sqi(seq(-4, 4, 0.1), f$act)) <= 1e-8))
  b <- fit_ms(s$grm, s$Y, s$y, "fixed", family = "negbin", max_iter = 300)
  expect_gt(b$act$theta, 3.5); expect_lt(b$act$theta, 7)
  e <- fit_ms(s$grm, s$Y, s$y, "eap", curve = "ispline", family = "poisson", max_knots = 2)
  expect_equal(e$path$n_knots, 0:2)
})

test_that("count families reject non-count data", {
  s <- sim_counts("poisson")
  expect_error(fit_ms(s$grm, s$Y, s$y + 0.5, "fixed", family = "poisson"), "integer")
})
