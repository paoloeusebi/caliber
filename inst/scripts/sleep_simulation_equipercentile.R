# =============================================================================
# sleep_simulation_equipercentile.R  --  Equipercentile linking as a comparator
#                                        in the simulation study
#
# Same ten scenarios, same replications and seeds as sleep_simulation.R (the
# data sets are identical: simulate_scenario(s, r)), and the same crosswalk
# error: RMSE over Y = 15..95 by 0.5 against the true E[X | Y = s], weighted by
# the true density of Y. Equipercentile linking (equipercentile_link()) maps Y
# to the X total by percentile ranks, with log-linear presmoothing of the X
# total preserving 4 moments (primary) and without smoothing (sensitivity). For
# a like-for-like comparison at item level each item X1..X7 is also linked to Y
# on its own (0..3, presmoothing with 2 moments). Equipercentile linking gives
# no item thresholds, so only the crosswalk is compared.
#
# No item response model is fitted, so this script needs neither mirt nor
# the results of sleep_simulation.R; it takes about 2 minutes on one core.
# Outputs: data/ms_sim_equipercentile.csv (per replication),
#          data/ms_sim_equipercentile_summary.csv. sleep_simulation_summary.R
# adds them to its crosswalk comparison when present.
# Environment: MS_REPS (default 200), CALIBER_OUT (default "output").
# =============================================================================

suppressPackageStartupMessages(library(caliber))
if (utils::packageVersion("caliber") < "0.2.0" || length(scenarios) != 10)
  stop("caliber is out of date (", utils::packageVersion("caliber"), ", ", length(scenarios),
       " scenarios); reinstall it: remotes::install_github(\"paoloeusebi/caliber\")")
items <- paste0("C", 1:7)
out <- Sys.getenv("CALIBER_OUT", "output")
dir.create(file.path(out, "data"), recursive = TRUE, showWarnings = FALSE)
n_reps  <- as.integer(Sys.getenv("MS_REPS", "200"))
cw_grid <- seq(15, 95, by = 0.5)
truths  <- lapply(seq_along(scenarios), scenario_truth, grid = cw_grid)

wrmse <- function(est, true, wt) sqrt(sum(wt * (est - true)^2))

t0 <- Sys.time()
rows <- do.call(rbind, lapply(seq_along(scenarios), function(s) {
  tr <- truths[[s]]
  do.call(rbind, lapply(seq_len(n_reps), function(r) {
    dat   <- simulate_scenario(s, r)
    total <- rowSums(dat$Y)
    ep4   <- equipercentile_link(total, dat$sqi, cw_grid, 21, presmooth = 4)
    ep0   <- equipercentile_link(total, dat$sqi, cw_grid, 21, presmooth = 0)
    it    <- sapply(items, function(k)
      wrmse(equipercentile_link(dat$Y[, k], dat$sqi, cw_grid, 3, presmooth = 2), tr$cw[[k]], tr$wt))
    data.frame(scenario = names(scenarios)[s], rep = r,
               cw_PSQI_eq = wrmse(ep4, tr$cw$PSQI, tr$wt),
               cw_PSQI_eq_raw = wrmse(ep0, tr$cw$PSQI, tr$wt),
               cw_items_eq = mean(it),
               bias_PSQI_eq = sum(tr$wt * (ep4 - tr$cw$PSQI)),
               cor_xy = cor(total, dat$sqi))
  }))
}))
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

mc <- function(v) c(mean(v), sd(v) / sqrt(length(v)))
summ <- do.call(rbind, lapply(split(rows, factor(rows$scenario, names(scenarios))), function(d) {
  a <- mc(d$cw_PSQI_eq); b <- mc(d$cw_PSQI_eq_raw); i <- mc(d$cw_items_eq)
  data.frame(scenario = d$scenario[1], method = "Equipercentile", reps = nrow(d),
             cw_psqi = a[1], cw_psqi_mcse = a[2], cw_psqi_raw = b[1], cw_psqi_raw_mcse = b[2],
             cw_items = i[1], cw_items_mcse = i[2], cw_psqi_bias = mean(d$bias_PSQI_eq),
             cor_xy = mean(d$cor_xy))
}))
write.csv(rows, file.path(out, "data", "ms_sim_equipercentile.csv"), row.names = FALSE)
write.csv(summ, file.path(out, "data", "ms_sim_equipercentile_summary.csv"), row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)
