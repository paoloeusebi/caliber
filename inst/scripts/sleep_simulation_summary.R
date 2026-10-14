# =============================================================================
# sleep_simulation_summary.R  --  Tables and figures from the simulation
#
# Run after sleep_simulation.R, with the same MS_HETERO, MS_DRAWS and CALIBER_OUT.
# Outputs: data/ms_sim_summary*.csv, data/ms_sim_threshold_summary*.csv,
#          figures/ms_sim_thresholds*.png, ms_sim_sqi_curves*.png,
#          ms_sim_crosswalk_curves*.png, ms_sim_crosswalk_error*.png
# =============================================================================

suppressPackageStartupMessages({ library(caliber); library(ggplot2); library(patchwork) })
items <- paste0("C", 1:7); n_items <- 7; thr_names <- threshold_names(items)
methods <- ms_methods; method_long <- ms_method_labels; method_cols <- ms_method_colours
out <- Sys.getenv("CALIBER_OUT", "output")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

suffix <- paste0(if (Sys.getenv("MS_HETERO", "1") == "1") "" else "_constsd",
                 if (as.integer(Sys.getenv("MS_DRAWS", "0")) > 0) "_mcint" else "")
res  <- readRDS(file.path(out, "data", paste0("ms_sim_results", suffix, ".rds")))
rows <- res$rows
sc_names  <- names(scenarios)
sc_labels <- sapply(scenarios, `[[`, "label")
true_thr  <- lapply(res$truths, function(t) c(t$thr))      # item-major within k
s_cols    <- paste0("s_", thr_names)
k_of      <- rep(1:3, each = n_items)

# Threshold errors per replication
err <- t(sapply(seq_len(nrow(rows)), function(i)
  unlist(rows[i, s_cols]) - true_thr[[match(rows$scenario[i], sc_names)]]))

rmse_mc <- function(mse) {            # pooled RMSE and its Monte Carlo SE (delta method)
  r <- sqrt(mean(mse)); c(r, sd(mse) / (2 * r * sqrt(length(mse))))
}

summary_tab <- do.call(rbind, lapply(sc_names, function(sc) do.call(rbind, lapply(methods, function(m) {
  i <- rows$scenario == sc & rows$method == m
  e <- err[i, , drop = FALSE]
  all <- rmse_mc(rowMeans(e^2)); k3 <- rmse_mc(rowMeans(e[, k_of == 3, drop = FALSE]^2))
  k12 <- rmse_mc(rowMeans(e[, k_of < 3, drop = FALSE]^2))
  cw_items <- rowMeans(as.matrix(rows[i, paste0("cw_", items)]))
  data.frame(scenario = sc, method = m, reps = sum(i),
             converged = mean(rows$converged[i]),
             thr_bias = mean(e), thr_rmse = all[1], thr_rmse_mcse = all[2],
             thr_rmse_k12 = k12[1], thr_rmse_k3 = k3[1],
             cw_items = mean(cw_items), cw_items_mcse = sd(cw_items) / sqrt(sum(i)),
             cw_psqi = mean(rows$cw_PSQI[i]), cw_psqi_mcse = sd(rows$cw_PSQI[i]) / sqrt(sum(i)),
             mu_m3 = mean(rows$mu_m3[i]), mu_0 = mean(rows$mu_0[i]), mu_p3 = mean(rows$mu_p3[i]),
             true_mu_m3 = scenarios[[sc]]$m(-3), true_mu_0 = scenarios[[sc]]$m(0),
             true_mu_p3 = scenarios[[sc]]$m(3),
             sigma0 = mean(rows$sigma[i]), g1 = mean(rows$g1[i]),
             g1_p025 = unname(quantile(rows$g1[i], .025)), g1_p975 = unname(quantile(rows$g1[i], .975)),
             a_C4 = mean(rows$a_C4[i]),
             iter_median = median(rows$iter[i]), df_mean = mean(rows$df[i]),
             secs_mean = mean(rows$secs[i]))
}))))
write.csv(summary_tab, file.path(out, "data", paste0("ms_sim_summary", suffix, ".csv")), row.names = FALSE)

thr_tab <- do.call(rbind, lapply(sc_names, function(sc) do.call(rbind, lapply(methods, function(m) {
  e <- err[rows$scenario == sc & rows$method == m, , drop = FALSE]
  data.frame(scenario = sc, method = m, threshold = thr_names, k = k_of,
             true_sqi = true_thr[[match(sc, sc_names)]], bias = colMeans(e),
             rmse = sqrt(colMeans(e^2)))
}))))
write.csv(thr_tab, file.path(out, "data", paste0("ms_sim_threshold_summary", suffix, ".csv")), row.names = FALSE)

cat("Pooled threshold RMSE (MC SE) / crosswalk wRMSE mean over C1..C7 (PSQI total)\n")
for (sc in sc_names) {
  s <- summary_tab[summary_tab$scenario == sc, ]
  cat(sprintf("%-16s %s\n", sc, paste(sprintf("%s %.2f (%.2f) k3 %.2f | cw %.3f (%.3f)",
      s$method, s$thr_rmse, s$thr_rmse_mcse, s$thr_rmse_k3, s$cw_items, s$cw_psqi), collapse = " || ")))
}

base_theme <- theme_minimal(base_size = 12) +
  theme(panel.border = element_rect(fill = NA, colour = "grey30"),
        legend.position = "bottom", plot.caption = element_text(hjust = 0))
method_shapes <- c("MS-Joint" = 16, "MS-EAP" = 17, "MS-Fixed" = 15)
method_lab <- setNames(paste0(methods, " (", method_long[methods], ")"), methods)
sc_factor <- function(x) factor(unname(sc_labels[x]), levels = rev(sc_labels))

# 1. Threshold RMSE by scenario: all thresholds, lower two, top one
long <- rbind(
  data.frame(summary_tab[c("scenario", "method")], what = "All 21 thresholds", rmse = summary_tab$thr_rmse,
             lo = summary_tab$thr_rmse - 1.96 * summary_tab$thr_rmse_mcse,
             hi = summary_tab$thr_rmse + 1.96 * summary_tab$thr_rmse_mcse),
  data.frame(summary_tab[c("scenario", "method")], what = "Lower two (k = 1, 2)",
             rmse = summary_tab$thr_rmse_k12, lo = NA, hi = NA),
  data.frame(summary_tab[c("scenario", "method")], what = "Top threshold (k = 3)",
             rmse = summary_tab$thr_rmse_k3, lo = NA, hi = NA))
long$what <- factor(long$what, levels = unique(long$what))
long$sc   <- sc_factor(long$scenario)
p1 <- ggplot(long, aes(rmse, sc, colour = method, shape = method)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", na.rm = TRUE,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  facet_wrap(~ what, nrow = 1) +
  scale_colour_manual(values = method_cols, labels = method_lab, breaks = methods, name = NULL) +
  scale_shape_manual(values = method_shapes, labels = method_lab, breaks = methods, name = NULL) +
  expand_limits(x = 0) +
  labs(x = "RMSE of threshold location on the SQI scale (SQI points)", y = NULL,
       title = if (res$hetero) NULL else "Constant residual SD (sensitivity analysis)",
       caption = sprintf(paste0("%d replications per scenario, n = 1000. Thresholds: mu-hat(b-hat_jk) ",
                                "against the true m(b_jk). Bars: 95%% Monte Carlo interval."), res$n_reps)) +
  base_theme
ggsave(file.path(out, "figures", paste0("ms_sim_thresholds", suffix, ".png")), p1, width = 12, height = 5, dpi = 200, bg = "white")

# 2. Mean estimated SQI curve (with the 5-95% range over replications) vs truth
curves <- do.call(rbind, lapply(seq_along(res$mu), function(i)
  data.frame(scenario = res$jobs$s[i], eta = res$eta_grid,
             stack(as.data.frame(res$mu[[i]])))))
names(curves)[names(curves) == "values"] <- "mu"; names(curves)[names(curves) == "ind"] <- "method"
agg <- aggregate(mu ~ scenario + method + eta, curves, function(v)
  c(mean = mean(v), lo = unname(quantile(v, .05)), hi = unname(quantile(v, .95))))
agg <- data.frame(agg[1:3], agg$mu)
agg$sc <- factor(unname(sc_labels[agg$scenario]), levels = sc_labels)
truth_curves <- do.call(rbind, lapply(seq_along(scenarios), function(s)
  data.frame(sc = factor(unname(sc_labels[s]), levels = sc_labels), eta = res$eta_grid,
             mu = scenarios[[s]]$m(res$eta_grid))))
p2 <- ggplot(agg, aes(eta)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = method), alpha = 0.10) +
  geom_line(aes(y = mean, colour = method), linewidth = 0.8) +
  geom_line(data = truth_curves, aes(y = mu), linetype = "dashed", colour = "grey15") +
  facet_wrap(~ sc, ncol = 3) +
  scale_colour_manual(values = method_cols, labels = method_lab, breaks = methods, name = NULL) +
  scale_fill_manual(values = method_cols, labels = method_lab, breaks = methods, name = NULL) +
  coord_cartesian(ylim = c(10, 95)) +
  labs(x = "Latent sleep disturbance (eta)", y = "mu(eta) = E[SQI | eta]",
       caption = "Solid: mean estimate over replications; shaded: 5th to 95th percentile; dashed: truth m(eta).") +
  base_theme
ggsave(file.path(out, "figures", paste0("ms_sim_sqi_curves", suffix, ".png")), p2, width = 12, height = 7.5, dpi = 200, bg = "white")

# 3. Mean estimated E[PSQI total | SQI] vs truth
cwc <- do.call(rbind, lapply(seq_along(res$cw), function(i) do.call(rbind, lapply(methods, function(m)
  data.frame(scenario = res$jobs$s[i], method = m, SQI = res$cw_grid, psqi = res$cw[[i]][[m]][, "PSQI"])))))
cwa <- aggregate(psqi ~ scenario + method + SQI, cwc, mean)
cwa$sc <- factor(unname(sc_labels[cwa$scenario]), levels = sc_labels)
cw_truth <- do.call(rbind, lapply(seq_along(scenarios), function(s)
  data.frame(sc = factor(unname(sc_labels[s]), levels = sc_labels), SQI = res$cw_grid,
             psqi = res$truths[[s]]$cw$PSQI, dens = res$truths[[s]]$wt)))
p3 <- ggplot(cwa, aes(SQI, psqi)) +
  geom_area(data = cw_truth, aes(y = dens / max(dens) * 4), fill = "grey85") +
  geom_line(aes(colour = method), linewidth = 0.8) +
  geom_line(data = cw_truth, linetype = "dashed", colour = "grey15") +
  facet_wrap(~ sc, ncol = 3) +
  scale_colour_manual(values = method_cols, labels = method_lab, breaks = methods, name = NULL) +
  labs(x = "Actigraphy SQI", y = "E[PSQI global | SQI]",
       caption = "Solid: mean estimated probabilistic crosswalk; dashed: truth; grey: true SQI density (scaled).") +
  base_theme
ggsave(file.path(out, "figures", paste0("ms_sim_crosswalk_curves", suffix, ".png")), p3, width = 12, height = 7.5, dpi = 200, bg = "white")

# 4. Crosswalk error by component
cwe <- do.call(rbind, lapply(sc_names, function(sc) do.call(rbind, lapply(methods, function(m) {
  i <- rows$scenario == sc & rows$method == m
  data.frame(scenario = sc, method = m, component = items,
             wrmse = colMeans(rows[i, paste0("cw_", items)]))
}))))
cwe$sc <- factor(unname(sc_labels[cwe$scenario]), levels = sc_labels)
p4 <- ggplot(cwe, aes(component, wrmse, colour = method, shape = method, group = method)) +
  geom_line(linewidth = 0.5, alpha = 0.6) + geom_point(size = 2.4) +
  facet_wrap(~ sc, ncol = 3) +
  scale_colour_manual(values = method_cols, labels = method_lab, breaks = methods, name = NULL) +
  scale_shape_manual(values = method_shapes, labels = method_lab, breaks = methods, name = NULL) +
  expand_limits(y = 0) +
  labs(x = NULL, y = "Weighted RMSE of E[C_j | SQI]",
       caption = "RMSE over SQI 15 to 95, weighted by the true SQI density; mean over replications.") +
  base_theme
ggsave(file.path(out, "figures", paste0("ms_sim_crosswalk_error", suffix, ".png")), p4, width = 12, height = 7, dpi = 200, bg = "white")

cat(sprintf("\nWrote %s/data/ms_sim_*summary%s.csv and %s/figures/ms_sim_*%s.png\n", out, suffix, out, suffix))
