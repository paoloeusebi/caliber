# =============================================================================
# sleep_simulation.R  --  Simulation study: the three monotone-spline
#                         calibration methods (MS-Joint, MS-EAP, MS-Fixed)
#                         under six data-generating scenarios
#
# Each replication: simulate (n = 1000), fit the PSQI-only GRM once, then the
# three methods. Targets on the SQI scale:
#   1. thresholds s_jk = mu(b_jk) against the true m(b_jk)  (bias, RMSE)
#   2. probabilistic crosswalk E[C_j | SQI = s], s = 15..95 by 0.5, against the
#      truth; error = RMSE over the grid weighted by the true SQI density
# With MS_DRAWS > 0, also 95% normal-approximation intervals for the
# thresholds, mu(eta) and the crosswalk, and their coverage (this coverage
# study was started on 2026-10-11 and stopped after 49 of 1200 replications;
# it has not been completed). The default MS_DRAWS = 0 fits point estimates
# only and reproduces the 8 October run reported in the manuscript (the fits
# use no random numbers, so adding intervals does not change them).
#
# Environment: MS_REPS (default 200), MS_CORES (default all but one),
# MS_HETERO (1 = log-linear residual SD, default; 0 = constant SD, outputs
# end in _constsd), MS_DRAWS (default 0), CALIBER_OUT (default "output"),
# MS_CACHE (cache folder, one .rds per scenario x replication, so an
# interrupted run resumes; delete it after any change to the code).
# Run time: about 3 hours on 4 cores for 1200 replications without intervals.
# Then run sleep_simulation_summary.R.
# =============================================================================

suppressPackageStartupMessages(library(caliber))
items <- paste0("C", 1:7); thr_names <- threshold_names(items); methods <- ms_methods
out <- Sys.getenv("CALIBER_OUT", "output")

n_reps  <- as.integer(Sys.getenv("MS_REPS", "200"))
n_cores <- as.integer(Sys.getenv("MS_CORES", max(1, parallel::detectCores() - 1)))
hetero  <- Sys.getenv("MS_HETERO", "1") == "1"
n_draw  <- as.integer(Sys.getenv("MS_DRAWS", "0"))
suffix  <- paste0(if (hetero) "" else "_constsd", if (n_draw > 0) "_mcint" else "")
cache   <- Sys.getenv("MS_CACHE", file.path(out, "data", paste0("ms_sim_cache", suffix)))
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
cw_grid  <- seq(15, 95, by = 0.5)
eta_grid <- seq(-3, 3, by = 0.1)
truths   <- lapply(seq_along(scenarios), scenario_truth, grid = cw_grid)

run_rep <- function(s, r) {
  file <- file.path(cache, sprintf("s%d_r%03d.rds", s, r))
  if (file.exists(file)) return(readRDS(file))
  t0  <- Sys.time()
  dat <- simulate_scenario(s, r)
  grm <- fit_grm(dat$Y)
  tr  <- truths[[s]]
  out <- lapply(methods, function(m) {
    t1 <- Sys.time()
    f  <- fit_ms(grm, dat$Y, dat$sqi, switch(m, "MS-Joint" = "joint",
                                             "MS-EAP" = "eap", "MS-Fixed" = "fixed"), hetero = hetero)
    secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    thr  <- threshold_locations(f$a, f$d, f$act)$sqi
    cw   <- crosswalk_probabilistic(cw_grid, f$a, f$d, f$act, items, "PSQI")
    cw_err <- sapply(c(items, "PSQI"), function(k) sqrt(sum(tr$wt * (cw[[k]] - tr$cw[[k]])^2)))
    n_thr <- length(thr); n_eta <- length(eta_grid); n_cw <- length(cw_grid)
    if (n_draw > 0) {
      t2 <- Sys.time()
      iv <- ms_intervals(f, ms_vcov(f, dat$Y, dat$sqi), function(a, d, act)
              c(threshold_locations(a, d, act)$sqi, mu_sqi(eta_grid, act),
                unlist(crosswalk_probabilistic(cw_grid, a, d, act, items, "PSQI"))),
              n_draw = n_draw, seed = 20261003 + 10000 * s + r)
      secs_iv <- as.numeric(difftime(Sys.time(), t2, units = "secs"))
    } else {                                   # point estimates only: NA intervals
      iv <- data.frame(se = rep(NA_real_, n_thr + n_eta + 8 * n_cw), lo = NA_real_, hi = NA_real_)
      secs_iv <- 0
    }
    thr_iv <- iv[seq_len(n_thr), ]
    mu_iv  <- iv[n_thr + seq_len(n_eta), ]
    cw_cov <- sapply(seq_along(c(items, "PSQI")), function(k) {
      j  <- n_thr + n_eta + (k - 1) * n_cw + seq_len(n_cw)
      tv <- tr$cw[[c(items, "PSQI")[k]]]
      c(sum(tr$wt * (tv >= iv$lo[j] & tv <= iv$hi[j])), sum(tr$wt * (iv$hi[j] - iv$lo[j])))
    })
    row <- data.frame(scenario = names(scenarios)[s], rep = r, method = m,
                      as.list(setNames(c(thr), paste0("s_", thr_names))),
                      as.list(setNames(cw_err, paste0("cw_", c(items, "PSQI")))),
                      mu_m3 = mu_sqi(-3, f$act), mu_0 = mu_sqi(0, f$act),
                      mu_p3 = mu_sqi(3, f$act), sigma = f$act$sigma,
                      g1 = if (hetero) f$act$log_sigma[2] else 0, a_C4 = f$a[4],
                      iter = f$iter, df = f$df, stop = f$stop_reason,
                      converged = f$converged, secs = secs, secs_iv = secs_iv,
                      as.list(setNames(thr_iv$se, paste0("se_", thr_names))),
                      as.list(setNames(thr_iv$lo, paste0("lo_", thr_names))),
                      as.list(setNames(thr_iv$hi, paste0("hi_", thr_names))),
                      as.list(setNames(cw_cov[1, ], paste0("cwcov_", c(items, "PSQI")))),
                      as.list(setNames(cw_cov[2, ], paste0("cwwid_", c(items, "PSQI")))),
                      check.names = FALSE)
    list(row = row, mu = mu_sqi(eta_grid, f$act), cw = as.matrix(cw),
         mu_lo = mu_iv$lo, mu_hi = mu_iv$hi)
  })
  res <- list(rows = do.call(rbind, lapply(out, `[[`, "row")),
              mu = sapply(out, `[[`, "mu"), cw = lapply(out, `[[`, "cw"),
              mu_lo = sapply(out, `[[`, "mu_lo"), mu_hi = sapply(out, `[[`, "mu_hi"),
              secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  colnames(res$mu) <- colnames(res$mu_lo) <- colnames(res$mu_hi) <- methods
  names(res$cw) <- methods
  saveRDS(res, file)
  res
}

jobs <- expand.grid(r = seq_len(n_reps), s = seq_along(scenarios))
jobs <- jobs[order(jobs$r, jobs$s), ]          # interleave scenarios: early results for all
cat(sprintf("%d scenarios x %d replications on %d cores; cache %s\n",
            length(scenarios), n_reps, n_cores, cache))
t0 <- Sys.time()
job_list <- Map(c, jobs$s, jobs$r)
run_safe <- function(p) tryCatch(run_rep(p[1], p[2]), error = function(e) conditionMessage(e))
if (n_cores > 1) {
  cl <- parallel::makeCluster(n_cores)
  invisible(parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(caliber))))
  parallel::clusterExport(cl, c("run_rep", "cache", "cw_grid", "eta_grid", "truths", "hetero",
                                "n_draw", "items", "thr_names", "methods"))
  res <- parallel::parLapplyLB(cl, job_list, run_safe)
  parallel::stopCluster(cl)
} else res <- lapply(job_list, run_safe)
failed <- sapply(res, is.character)
cat(sprintf("Done in %.1f min; %d of %d replications ran\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(!failed), length(res)))
if (any(failed)) print(table(unlist(res[failed])))

ok <- res[!failed]
saveRDS(list(rows = do.call(rbind, lapply(ok, `[[`, "rows")),
             mu = lapply(ok, `[[`, "mu"), cw = lapply(ok, `[[`, "cw"),
             mu_lo = lapply(ok, `[[`, "mu_lo"), mu_hi = lapply(ok, `[[`, "mu_hi"),
             jobs = jobs[!failed, ], cw_grid = cw_grid, eta_grid = eta_grid,
             truths = truths, n_reps = n_reps, hetero = hetero, n_draw = n_draw),
        file.path(out, "data", paste0("ms_sim_results", suffix, ".rds")))
cat(sprintf("Wrote %s/data/ms_sim_results%s.rds\n", out, suffix))
