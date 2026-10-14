test_that("GRM category probabilities sum to one and expected scores rise with eta", {
  p <- grm_probs(seq(-3, 3, 0.5), 1.5, c(1, 0, -1))
  expect_equal(rowSums(p), rep(1, nrow(p)))
  e <- grm_expected(seq(-3, 3, 0.5), 1.5, c(1, 0, -1))
  expect_true(all(diff(e) > 0))
})

test_that("B-spline basis has 24 columns, sums to one and clamps outside [-4, 4]", {
  B <- bspline_basis(c(-5, -4, 0, 4, 5))
  expect_equal(ncol(B), 24)
  expect_equal(rowSums(B), rep(1, 5))
  expect_equal(B[1, ], B[2, ])
  expect_equal(B[4, ], B[5, ])
})

test_that("a non-increasing coefficient vector gives a non-increasing curve", {
  act <- monospline_act(seq(80, 30, length.out = 24), NA, c(log(6), 0.1))
  mu <- mu_sqi(seq(-4, 4, 0.1), act)
  expect_true(all(diff(mu) <= 1e-12))
  expect_equal(sigma_sqi(0, act), 6)
  expect_equal(sigma_sqi(1, act), 6 * exp(0.1))
})

test_that("threshold locations and crosswalks have the expected shape", {
  act <- monospline_act(seq(80, 30, length.out = 24), 6)
  thr <- threshold_locations(true_items$a, true_items$d, act)
  expect_equal(dim(thr$sqi), c(7, 3))
  expect_true(all(thr$sqi[, 1] >= thr$sqi[, 3]))
  cw <- crosswalk(c(40, 60, 80), true_items$a, true_items$d, act, total_name = "PSQI")
  expect_equal(nrow(cw), 3)
  expect_true(all(diff(cw$E_PSQI) < 0))
  ex <- exceedance_probabilistic(c(40, 80), true_items$a, true_items$d, act)
  expect_equal(colnames(ex), threshold_names(paste0("C", 1:7)))
})

test_that("the shipped data sets load", {
  d <- sleep_data()
  expect_equal(nrow(d), 1000)
  expect_true(all(c(paste0("C", 1:7), "SQI") %in% names(d)))
  g <- gait_data()
  expect_equal(g$n, 1426)
  expect_equal(gait_data(exclude_aids = TRUE)$n, 1265)
  expect_equal(gait_data(skip_unable = TRUE)$n_skip_recoded, 103)
})

test_that("the regenerated sleep data match the shipped copy", {
  d <- sleep_data(); s <- simulate_sleep_data()
  expect_equal(d$SQI, s$SQI)
  expect_equal(as.matrix(d[paste0("C", 1:7)]), as.matrix(s[paste0("C", 1:7)]),
               ignore_attr = TRUE)
})
