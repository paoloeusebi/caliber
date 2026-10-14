# =============================================================================
# gait_main.R  --  Real-data example: gait speed (timed 20-ft walk) read on
#                  seven lower-body mobility items (NHANES 2001-2002, ages 60+)
#
#   MS-Fixed  main analysis      (items from the GRM on the PRO alone, then fixed)
#   MS-Joint  robustness check   (items and curve estimated together)
#   MS-EAP    supplement         (curve regressed on the PRO-only EAP scores)
#
# The model code is the caliber package. In its code "sqi" is the digital
# measure; here it is gait speed in m/s. eta = mobility difficulty,
# higher = more difficulty; gait speed is non-increasing in eta. Residual SD
# log-linear in eta (the main model).
#
# Direction of interpretation: from gait speed to the PRO. Outputs are
#   (i)   item thresholds in gait-speed units, s_jk = mu(b_jk)
#   (ii)  the probabilistic crosswalk E[C_j | speed] and E[total | speed]
#   (iii) P(C_j >= k | speed)
# with nonparametric bootstrap intervals (persons resampled, every step
# refitted, including the GRM and the AIC stop). A resample in which an item
# lost a response category is redrawn, because the GRM then has no threshold
# for it (C5 has 2 "unable" responses); the number redrawn is reported, and
# per threshold the share of resamples redrawn because that category was lost.
# The intervals are therefore conditional on every category being present.
#
# Environment: GAIT_BOOT (resamples, default 200), GAIT_CORES,
#   GAIT_EXCL_AIDS=1 for the sensitivity analysis without walking-aid users
#   (outputs end in _noaids),
#   GAIT_SKIP_UNABLE=1 for the sensitivity analysis in which C1 and C2, not
#   asked of people who need special equipment to walk, are scored "unable"
#   for them (outputs end in _skipunable).
#
# Run:  Rscript gait_main.R  (data: gait_data(), shipped with the package)
# Outputs: $CALIBER_OUT/data/gait_*.csv, $CALIBER_OUT/figures/gait_*.png
# Run time: about 20 min for 200 resamples on 15 cores.
# =============================================================================

suppressPackageStartupMessages({ library(caliber); library(ggplot2); library(patchwork) })
items <- paste0("C", 1:7); n_items <- 7; thr_names <- threshold_names(items)
methods <- ms_methods; method_long <- ms_method_labels; method_cols <- ms_method_colours
item_labels <- gait_item_labels; category_labels <- gait_category_labels

n_boot  <- as.integer(Sys.getenv("GAIT_BOOT", "200"))
n_cores <- as.integer(Sys.getenv("GAIT_CORES", max(1, parallel::detectCores() - 1)))
seed    <- 20261010
out     <- Sys.getenv("CALIBER_OUT", "output")
dir.create(file.path(out, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
od <- function(f) file.path(out, "data", f); of <- function(f) file.path(out, "figures", f)
g <- gait_data(exclude_aids = Sys.getenv("GAIT_EXCL_AIDS", "0") == "1",
               skip_unable  = Sys.getenv("GAIT_SKIP_UNABLE", "0") == "1")
Y <- g$Y; speed <- g$speed; n <- g$n; tag <- g$tag
cat(sprintf("n = %d%s%s%s | gait speed mean %.3f, SD %.3f, range %.2f to %.2f m/s\n", n,
            if (g$exclude_aids) " (walking-aid users excluded)" else "",
            if (g$skip_unable) sprintf(" (C1, C2 scored unable for %d people who need equipment to walk)",
                                       g$n_skip_recoded) else "",
            "", mean(speed), sd(speed),
            min(speed), max(speed)))
cat(sprintf("No difficulty on every item answered: %.1f%% (%d of %d); among those who answered all 7: %.1f%%\n",
            100 * mean(rowSums(Y, na.rm = TRUE) == 0), sum(rowSums(Y, na.rm = TRUE) == 0), n,
            100 * mean(rowSums(Y[complete.cases(Y), ]) == 0)))


# ---- Fit ------------------------------------------------------------------------

grm  <- fit_grm(Y)
t0   <- Sys.time()
fits <- fit_three(Y, speed, grm)
cat(sprintf("Three fits in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
for (m in methods) {
  f <- fits[[m]]
  cat(sprintf("%-8s AIC stop at iteration %4d (BIC %4d; %s) | df %.2f | sigma(0) %.4f | g1 %.3f | -2logL(C, speed) %.2f\n",
              m, f$iter, f$iter_bic, f$stop_reason, f$df, f$act$sigma, f$act$log_sigma[2],
              neg2ll_joint(Y, speed, f$a, f$d, f$act)))
}

# The curve is constant outside [-4, 4], so mu(-4) is its ceiling and mu(4) its
# floor: no eta has a mean gait speed outside them. R2 is approximate (residual
# SD at eta = 0 against the SD of gait speed).
curve_ends <- do.call(rbind, lapply(methods, function(m)
  data.frame(method = m, caliber::curve_ends(fits[[m]]$act, speed))))
rownames(curve_ends) <- methods
for (m in methods) with(curve_ends[m, ], cat(sprintf(
  "%-8s curve floor mu(4) %.3f, ceiling mu(-4) %.3f m/s | observed speeds below floor %.1f%%, above ceiling %.1f%% | R2 ~ %.3f\n",
  m, floor, ceiling, 100 * share_below_floor, 100 * share_above_ceiling, r2_approx)))
floor_tol       <- 0.01     # m/s: a threshold this close to the floor is "at the floor"
min_at_or_above <- 30       # persons at or above a category; fewer is "sparse"


# ---- Reported quantities and bootstrap ---------------------------------------------

eta_report <- -1:3                              # mu(eta) reported at these points
curve_grid <- seq(-2, 3.5, by = 0.05)           # for the curve figure
speed_grid <- seq(0.20, 1.60, by = 0.01)        # crosswalk grid, m/s
speed_report <- seq(0.3, 1.4, by = 0.1)         # crosswalk reported with intervals

gsumm <- function(f) {
  thr <- threshold_locations(f$a, f$d, f$act)
  cw  <- crosswalk_probabilistic(speed_report, f$a, f$d, f$act, items)
  c(setNames(f$a, paste0("a_", items)),
    setNames(f$d[, 1], paste0("d1_", items)), setNames(f$d[, 2], paste0("d2_", items)),
    setNames(f$d[, 3], paste0("d3_", items)),
    sigma0 = f$act$sigma, g1 = f$act$log_sigma[2],
    setNames(mu_sqi(eta_report, f$act), paste0("mu(", eta_report, ")")),
    setNames(c(thr$b), paste0("b_", thr_names)),
    setNames(c(thr$sqi), paste0("s_", thr_names)),
    setNames(cw$total, paste0("Etotal(", speed_report, ")")),
    setNames(mu_sqi(curve_grid, f$act), paste0("curve_", seq_along(curve_grid))),
    setNames(crosswalk_probabilistic(speed_grid, f$a, f$d, f$act, items)$total,
             paste0("cwtot_", seq_along(speed_grid))),
    iter = f$iter, df = f$df)
}

# A resample in which an item lost a response category that the data have is
# redrawn (ms_bootstrap(keep_categories = TRUE)); lost_Cj_k = 1 when resample r
# was redrawn at least once because item j had lost category k.
cat(sprintf("\nBootstrap: %d resamples x 3 methods on %d cores ...\n", n_boot, n_cores))
t0 <- Sys.time()
boot <- ms_bootstrap(Y, speed, n_boot, seed, summary_fun = gsumm, keep_categories = TRUE,
                     cores = n_cores)
boot_list <- split(boot, boot$rep)
redrawn <- boot$redrawn[boot$method == "MS-Fixed"]
cat(sprintf("Done in %.1f min; %d of %d resamples usable; %d resamples needed a redraw (%d redraws in all)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            length(boot_list), n_boot, sum(redrawn > 0), sum(redrawn)))
redraw_share <- colMeans(boot[boot$method == "MS-Fixed", paste0("lost_", thr_names)])
cat("Share of resamples redrawn because the item had lost the category:",
    if (any(redraw_share > 0)) paste(sprintf("%s %.1f%%", sub("lost_", "", names(redraw_share)),
                                             100 * redraw_share)[redraw_share > 0], collapse = ", ") else "none", "\n")

report <- grep("^(curve_|cwtot_)|^(iter|df)$", names(gsumm(fits[[1]])), invert = TRUE, value = TRUE)
point  <- lapply(fits, gsumm)
qlo <- function(x) quantile(x, 0.025); qhi <- function(x) quantile(x, 0.975)

estimates <- do.call(rbind, lapply(methods, function(m) {
  bm <- boot[boot$method == m, report]
  data.frame(method = m, param = report, est = unname(point[[m]][report]),
             se = unname(apply(bm, 2, sd)), lo95 = unname(apply(bm, 2, qlo)),
             hi95 = unname(apply(bm, 2, qhi)))
}))
write.csv(estimates, od(sprintf("gait_estimates%s.csv", tag)), row.names = FALSE)


# ---- (i) Thresholds in gait-speed units ---------------------------------------------
# A threshold is flagged when fewer than 30 persons are at or above its category
# (flag_sparse) or when it lies within 0.01 m/s of the curve's floor mu(4)
# (at_floor): thresholds at the floor cannot be told apart and are reported as
# "below about" the floor. The flag of the simulation is kept as well
# (poorly_identified: bootstrap 95% interval wider than 2 sigma(0), the +-1 SD
# noise band of one person's gait speed). boot_redraw_share: share of resamples
# redrawn because the item had lost that category.

thresholds <- do.call(rbind, lapply(methods, function(m) {
  e  <- estimates[estimates$method == m, ]
  gp <- function(p, col) e[[col]][match(paste0(p, "_", thr_names), e$param)]
  data.frame(method = m, item = rep(items, 3), label = rep(item_labels, 3),
             k = rep(1:3, each = n_items),
             n_at_or_above = c(sapply(1:3, function(k) colSums(Y >= k, na.rm = TRUE))),
             b_est = gp("b", "est"), b_se = gp("b", "se"),
             b_lo95 = gp("b", "lo95"), b_hi95 = gp("b", "hi95"),
             sqi_est = gp("s", "est"), sqi_se = gp("s", "se"),
             sqi_lo95 = gp("s", "lo95"), sqi_hi95 = gp("s", "hi95"),
             poorly_identified = gp("s", "hi95") - gp("s", "lo95") > 2 * fits[[m]]$act$sigma,
             flag_sparse = c(sapply(1:3, function(k) colSums(Y >= k, na.rm = TRUE))) < min_at_or_above,
             at_floor = gp("s", "est") - curve_ends[m, "floor"] < floor_tol,
             curve_floor = curve_ends[m, "floor"],
             boot_redraw_share = unname(redraw_share))
}))
# A category nobody used has no threshold (see fit_grm() in caliber):
# it is reported as NA and left out of the comparisons. The Wright maps keep
# the placeholder, which sits at the slowest end of the curve.
thr_out <- thresholds
thr_out[thr_out$n_at_or_above == 0,
        c(grep("^(b|sqi)_", names(thr_out), value = TRUE), "poorly_identified", "at_floor")] <- NA
thr_out <- thr_out[order(thr_out$method, thr_out$item, thr_out$k), ]
names(thr_out) <- sub("^sqi_", "speed_", names(thr_out))
write.csv(thr_out, od(sprintf("gait_thresholds_speed%s.csv", tag)), row.names = FALSE)

cat("\nThresholds in gait-speed units (m/s), MS-Fixed: estimate [95% bootstrap interval]\n")
tf <- thr_out[thr_out$method == "MS-Fixed", ]
flag_text <- function(t) {
  x <- cbind(ifelse(t$at_floor %in% TRUE, "at floor", NA), ifelse(t$flag_sparse, "n < 30", NA),
             ifelse(t$poorly_identified %in% TRUE, "interval > 2 sigma(0)", NA))
  ifelse(t$n_at_or_above == 0, "not estimable", apply(x, 1, function(r) paste(na.omit(r), collapse = "; ")))
}
print(data.frame(item = tf$label, k = tf$k, n_ge_k = tf$n_at_or_above, b = round(tf$b_est, 2),
                 speed = ifelse(tf$at_floor %in% TRUE, sprintf("below ~%.2f m/s", tf$curve_floor),
                                sprintf("%.2f [%.2f, %.2f]", tf$speed_est, tf$speed_lo95, tf$speed_hi95)),
                 flag = flag_text(tf), redrawn = sprintf("%.1f%%", 100 * tf$boot_redraw_share)),
      row.names = FALSE)


# ---- (ii) Crosswalk and (iii) exceedance probabilities ---------------------------------

cws <- lapply(fits, function(f) crosswalk(speed_grid, f$a, f$d, f$act, items, value_name = "gait_speed"))
cw_long <- do.call(rbind, lapply(methods, function(m) {
  bm <- boot[boot$method == m, paste0("cwtot_", seq_along(speed_grid))]
  x  <- cws[[m]]
  data.frame(method = m, x, E_total_lo95 = apply(bm, 2, qlo), E_total_hi95 = apply(bm, 2, qhi))
}))
write.csv(cw_long, od(sprintf("gait_crosswalk%s.csv", tag)), row.names = FALSE)

exceed <- do.call(rbind, lapply(methods, function(m) {
  p <- exceedance_probabilistic(speed_grid, fits[[m]]$a, fits[[m]]$d, fits[[m]]$act, items)
  data.frame(method = m, gait_speed = rep(speed_grid, ncol(p)),
             item = rep(rep(items, 3), each = length(speed_grid)),
             k = rep(rep(1:3, each = n_items), each = length(speed_grid)), prob = c(p))
}))
write.csv(exceed, od(sprintf("gait_exceedance%s.csv", tag)), row.names = FALSE)

cat("\nExpected item scores (0-3) and total (0-21) at selected gait speeds, MS-Fixed (probabilistic)\n")
sel <- cw_long[cw_long$method == "MS-Fixed" & round(cw_long$gait_speed * 100) %in% c(40, 60, 80, 100, 120, 140), ]
print(round(sel[c("gait_speed", paste0("E_", items), "E_total", "E_total_lo95", "E_total_hi95")], 2),
      row.names = FALSE)


# ---- MS-Fixed against MS-Joint (robustness) and MS-EAP (supplement) -----------------

eta_fine <- seq(-3, 3, by = 0.01)
central  <- speed_grid >= quantile(speed, .01) & speed_grid <= quantile(speed, .99)
ref <- "MS-Fixed"
comparison <- do.call(rbind, lapply(methods, function(m) {
  f <- fits[[m]]; e <- estimates[estimates$method == m, ]
  ci <- function(p) sprintf("%.3f to %.3f", e$lo95[e$param == p], e$hi95[e$param == p])
  th <- thresholds[thresholds$method == m, ]; th0 <- thresholds[thresholds$method == ref, ]
  seen <- th$n_at_or_above > 0; th <- th[seen, ]; th0 <- th0[seen, ]
  data.frame(method = m, label = method_long[[m]], iter = f$iter, df = f$df,
             boot_iter_median = median(boot$iter[boot$method == m]),
             sigma0 = f$act$sigma, sigma0_ci = ci("sigma0"),
             g1 = f$act$log_sigma[2], g1_ci = ci("g1"),
             mu_m1 = mu_sqi(-1, f$act), mu_0 = mu_sqi(0, f$act), mu_1 = mu_sqi(1, f$act),
             mu_2 = mu_sqi(2, f$act), mu_3 = mu_sqi(3, f$act),
             neg2ll_C_speed = neg2ll_joint(Y, speed, f$a, f$d, f$act),
             n_thr_poorly_identified = sum(th$poorly_identified),
             maxdiff_mu_vs_fixed = max(abs(mu_sqi(eta_fine, f$act) - mu_sqi(eta_fine, fits[[ref]]$act))),
             maxdiff_thr_speed_vs_fixed = max(abs(th$sqi_est - th0$sqi_est)),
             maxdiff_a_vs_fixed = max(abs(f$a - fits[[ref]]$a)),
             maxdiff_b_vs_fixed = max(abs(th$b_est - th0$b_est)),
             maxdiff_Eitem_vs_fixed = max(abs(as.matrix(cws[[m]][central, paste0("E_", items)]) -
                                              as.matrix(cws[[ref]][central, paste0("E_", items)]))),
             maxdiff_Etotal_vs_fixed = max(abs(cws[[m]]$E_total - cws[[ref]]$E_total)[central]),
             n_thr_sparse = sum(th$flag_sparse), n_thr_at_floor = sum(th$at_floor),
             curve_ends[m, -1], row.names = NULL)
}))
write.csv(comparison, od(sprintf("gait_comparison%s.csv", tag)), row.names = FALSE)
cat("\nComparison of the three methods (differences are against MS-Fixed; crosswalk over the central 98% of speeds)\n")
print(t(data.frame(lapply(comparison[-(1:2)], function(x) if (is.numeric(x)) round(x, 3) else x),
                   row.names = comparison$method)), quote = FALSE)

cat("\nThreshold gait speeds by method (estimate; MS-Fixed 95% interval)\n")
wide <- reshape(thresholds[thresholds$n_at_or_above > 0, c("method", "item", "k", "sqi_est")], idvar = c("item", "k"),
                timevar = "method", direction = "wide")
wide <- wide[order(wide$item, wide$k), ]
print(data.frame(item = wide$item, k = wide$k, Fixed = round(wide[["sqi_est.MS-Fixed"]], 3),
                 Joint = round(wide[["sqi_est.MS-Joint"]], 3), EAP = round(wide[["sqi_est.MS-EAP"]], 3)),
      row.names = FALSE)


# ---- Figures ----------------------------------------------------------------------------

frame_theme <- theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = NA, colour = "grey30"),
        plot.caption = element_text(hjust = 0), legend.position = "bottom")
sample_text <- sprintf("NHANES 2001-2002, ages 60+%s%s, n = %d",
                       if (g$exclude_aids) ", walking-aid users excluded" else "",
                       if (g$skip_unable) ", skipped walking items scored unable" else "", n)

# 1. Gait-speed curves mu(eta) with bootstrap 95% pointwise bands, one panel per method
curves <- do.call(rbind, lapply(methods, function(m) {
  bm <- boot[boot$method == m, paste0("curve_", seq_along(curve_grid))]
  data.frame(method = m, eta = curve_grid, est = mu_sqi(curve_grid, fits[[m]]$act),
             lo = apply(bm, 2, qlo), hi = apply(bm, 2, qhi))
}))
curves$panel <- factor(paste0(curves$method, ": ", method_long[curves$method]),
                       levels = paste0(methods, ": ", method_long[methods]))
write.csv(curves[c("method", "eta", "est", "lo", "hi")], od(sprintf("gait_curves%s.csv", tag)), row.names = FALSE)
thr_rug <- thresholds[c("method", "b_est")]
thr_rug$panel <- factor(paste0(thr_rug$method, ": ", method_long[thr_rug$method]), levels = levels(curves$panel))
p_curves <- ggplot(curves, aes(eta)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = method), alpha = 0.18) +
  geom_line(aes(y = est, colour = method), linewidth = 0.9) +
  geom_rug(data = thr_rug, aes(x = b_est), inherit.aes = FALSE, sides = "b", colour = "grey30") +
  facet_wrap(~ panel) +
  scale_colour_manual(values = method_cols, guide = "none") +
  scale_fill_manual(values = method_cols, guide = "none") +
  labs(x = "Latent mobility difficulty (eta; N(0, 1) metric)", y = "E[gait speed | eta] (m/s)",
       caption = paste(strwrap(sprintf(paste0(
         "%s. Solid: estimated non-increasing B-spline. Shaded: bootstrap 95%% pointwise interval ",
         "(%d resamples, every step refitted). Ticks: the 21 item thresholds b_jk. ",
         "Residual SD log-linear in eta."), sample_text, n_boot), 150), collapse = "\n")) +
  frame_theme
ggsave(of(sprintf("gait_speed_curves%s.png", tag)), p_curves, width = 11, height = 4.8, dpi = 200, bg = "white")

# 2. Wright maps on the gait-speed scale, one per method (wright_map_sqi() of
#    ms_plots.R; its x axis is rescaled from SQI points to m/s afterwards, and
#    the person panel is titled as what it shows: model-implied mean speeds
#    mu(eta), which stay between the curve's floor and ceiling, not observed speeds)
for (m in methods) {
  f <- fits[[m]]
  eta_pv  <- plausible_values(f, grm, speed, seed = seed + 1)
  pv_text <- if (m == "MS-EAP") "one plausible value of eta per person from p(eta | items) (mirt; gait speed not used)" else
    "one plausible value of eta per person from p(eta | items, gait speed)"
  how <- switch(m,
    "MS-Joint" = "GRM items and the monotone gait-speed curve estimated together, eta integrated out",
    "MS-EAP"   = "GRM on the items alone; monotone gait-speed curve regressed on the EAP scores",
    "MS-Fixed" = "GRM on the items alone, items then fixed; monotone gait-speed curve with eta integrated out")
  person_speed <- mu_sqi(eta_pv, f$act)
  th   <- thresholds[thresholds$method == m, ]
  xlim <- range(c(0.3, 1.3, person_speed, th$sqi_est))
  wm <- wright_map_sqi(
    person_sqi = person_speed, thr = th, item_labels = item_labels, xlim = xlim, binwidth = 0.02,
    category_labels = category_labels,
    note = paste0("Note. ", m, " (", tolower(method_long[[m]]), "): ", how, " (", sample_text, "). ",
                  "Persons: ", pv_text, ", mapped to the model-implied mean gait speed mu(eta), ",
                  sprintf(paste0("which lies between the curve's floor (%.2f) and ceiling (%.2f m/s). Observed gait ",
                                 "speeds range from %.2f to %.2f m/s; %.1f%% are above the ceiling and %.1f%% ",
                                 "below the floor.\n"), curve_ends[m, "floor"], curve_ends[m, "ceiling"],
                          min(speed), max(speed), 100 * curve_ends[m, "share_above_ceiling"],
                          100 * curve_ends[m, "share_below_floor"]),
                  "Thresholds: GRM thresholds b_jk placed at mu(b_jk). Each band spans the gait speeds at ",
                  "which that category is the median response. Outer bands fade because they are open-ended."))
  wm <- suppressMessages(wm & scale_x_continuous(limits = xlim, expand = c(0, 0), breaks = seq(0, 2, by = 0.1),
                                                 oob = scales::oob_keep))
  wm[[1]] <- wm[[1]] + labs(title = sprintf("Model-implied mean gait speed mu(eta): mean %.2f, SD %.2f m/s",
                                            mean(person_speed), sd(person_speed)))
  wm[[2]] <- wm[[2]] + labs(x = "Gait speed, timed 20-ft walk (m/s; higher = faster)")
  ggsave(of(sprintf("gait_wright_map_%s%s.png", sub("MS-", "", m), tag)), wm,
         width = 10, height = 7.5, dpi = 200, bg = "white")
}

# 3. Probabilistic crosswalk E[C_j | gait speed], one panel per method
item_cols <- setNames(hcl.colors(n_items, "Dark 3"), item_labels)
cw_plot <- do.call(rbind, lapply(methods, function(m) do.call(rbind, lapply(seq_len(n_items), function(j)
  data.frame(method = m, speed = speed_grid, score = cws[[m]][[paste0("E_", items[j])]], item = item_labels[j])))))
cw_plot$panel <- factor(paste0(cw_plot$method, ": ", method_long[cw_plot$method]),
                        levels = paste0(methods, ": ", method_long[methods]))
cw_plot$item <- factor(cw_plot$item, levels = item_labels)
p_cw <- ggplot(cw_plot, aes(speed, score, colour = item)) +
  geom_line(linewidth = 0.9) +
  geom_rug(data = data.frame(speed = speed), aes(x = speed), inherit.aes = FALSE, sides = "b", alpha = 0.1) +
  facet_wrap(~ panel) +
  scale_colour_manual(values = item_cols, name = NULL) +
  coord_cartesian(xlim = c(0.2, 1.6), ylim = c(0, 3)) +
  labs(x = "Gait speed (m/s)", y = "Expected item score, E[C | gait speed] (0 = no difficulty, 3 = unable)",
       caption = sample_text) +
  frame_theme + guides(colour = guide_legend(nrow = 2))
ggsave(of(sprintf("gait_expected_scores%s.png", tag)), p_cw, width = 12, height = 5.6, dpi = 200, bg = "white")

# 4. Main figure (MS-Fixed): P(C_j >= k | gait speed) per item, and the expected total
ex <- exceed[exceed$method == "MS-Fixed", ]
ex$item  <- factor(item_labels[match(ex$item, items)], levels = item_labels)
ex$level <- factor(c("At least some difficulty", "At least much difficulty", "Unable")[ex$k],
                   levels = c("At least some difficulty", "At least much difficulty", "Unable"))
floor_fixed <- curve_ends["MS-Fixed", "floor"]      # speeds below the curve's floor are shaded
floor_shade <- annotate("rect", xmin = -Inf, xmax = floor_fixed, ymin = -Inf, ymax = Inf,
                        fill = "grey40", alpha = 0.15)
p_ex <- ggplot(ex, aes(gait_speed, prob, colour = level)) +
  floor_shade +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ item, nrow = 2) +
  scale_colour_manual(values = setNames(viridisLite::mako(3, begin = 0.75, end = 0.2), levels(ex$level)), name = NULL) +
  coord_cartesian(xlim = c(0.2, 1.6), ylim = c(0, 1)) +
  labs(x = "Gait speed (m/s)", y = "P(response at or above the level | gait speed)") +
  frame_theme
tot <- cw_long[cw_long$method == "MS-Fixed", ]
p_tot <- ggplot(tot, aes(gait_speed)) +
  floor_shade +
  geom_ribbon(aes(ymin = E_total_lo95, ymax = E_total_hi95), fill = method_cols[["MS-Fixed"]], alpha = 0.2) +
  geom_line(aes(y = E_total), colour = method_cols[["MS-Fixed"]], linewidth = 1) +
  geom_rug(data = data.frame(speed = speed), aes(x = speed), inherit.aes = FALSE, sides = "b", alpha = 0.1) +
  coord_cartesian(xlim = c(0.2, 1.6)) +
  labs(x = "Gait speed (m/s)", y = "Expected total of the 7 items (0-21)",
       caption = paste(strwrap(sprintf(paste0(
         "MS-Fixed, %s. Top: probability of reporting at least each level of difficulty at a given gait ",
         "speed. Bottom: expected total with bootstrap 95%% interval (%d resamples); ticks are observed speeds. ",
         "Shaded: gait speeds below the curve's floor (%.2f m/s; %.1f%% of observed speeds), slower than the ",
         "mean speed at any level of difficulty, where the items no longer separate people."),
         sample_text, n_boot, floor_fixed, 100 * curve_ends["MS-Fixed", "share_below_floor"]), 150),
         collapse = "\n")) +
  frame_theme
ggsave(of(sprintf("gait_interpretation_Fixed%s.png", tag)), p_ex / p_tot + plot_layout(heights = c(3, 2)),
       width = 12, height = 9.5, dpi = 200, bg = "white")

write.csv(boot[, !grepl("^(curve_|cwtot_)", names(boot))], od(sprintf("gait_bootstrap%s.csv", tag)), row.names = FALSE)
cat(sprintf("\nWrote %s/data/gait_*%s.csv and %s/figures/gait_*%s.png\n", out, tag, out, tag))
