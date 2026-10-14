# =============================================================================
# gait_dimensionality.R  --  Checks on the seven mobility items before the
#                            main fit: 1 versus 2 factors, local independence,
#                            and the fit of the graded response model
#
#   1. Polychoric correlations and their eigenvalues
#   2. GRM (mirt): 1 factor against 2 correlated factors
#      {C1..C5: walking, steps, rooms, chair, bed} and {C6, C7: stooping,
#      standing}, and against an exploratory 2-factor model
#   3. The same comparison as an ordinal CFA (lavaan, WLSMV), for CFI/TLI/RMSEA/SRMR
#   4. Local independence under the 1-factor GRM: Yen's Q3 for every item pair
#   5. Fit of the 1-factor GRM: C2 statistic with RMSEA, SRMSR, CFI, TLI;
#      item fit S-X2; item parameters
#
# Missing responses (C1 and C2 are skipped for walking-aid users) are kept:
# mirt uses all available responses, lavaan pairwise polychorics. The C2 and
# S-X2 statistics need complete responses and are computed on those persons.
#
# Run:  Rscript gait_dimensionality.R   (needs lavaan)
# Outputs: $CALIBER_OUT/data/gait_dim_*.csv
# =============================================================================

suppressPackageStartupMessages({ library(caliber); library(mirt); library(lavaan) })
items <- paste0("C", 1:7); n_items <- 7; item_labels <- gait_item_labels
out <- Sys.getenv("CALIBER_OUT", "output")
dir.create(file.path(out, "data"), recursive = TRUE, showWarnings = FALSE)

g <- gait_data(exclude_aids = Sys.getenv("GAIT_EXCL_AIDS", "0") == "1",
               skip_unable  = Sys.getenv("GAIT_SKIP_UNABLE", "0") == "1")
Y <- g$Y; tag <- g$tag
cat(sprintf("n = %d (%d with all 7 items)\n", g$n, sum(complete.cases(Y))))
Yo <- as.data.frame(lapply(as.data.frame(Y), ordered))

# ---- 1. Polychoric correlations ---------------------------------------------------
R <- unclass(lavCor(Yo, ordered = items, missing = "pairwise"))
ev <- eigen(R)$values
cat("\nPolychoric correlations\n"); print(round(R, 2))
cat(sprintf("Eigenvalues: %s | first / second = %.1f | first explains %.0f%%\n",
            paste(round(ev, 2), collapse = ", "), ev[1] / ev[2], 100 * ev[1] / n_items))
write.csv(round(R, 4), file.path(out, "data", sprintf("gait_dim_polychoric%s.csv", tag)))

# ---- 2. GRM: 1 factor vs 2 factors (mirt) ---------------------------------------
m1  <- mirt(as.data.frame(Y), 1, itemtype = "graded", SE = TRUE, verbose = FALSE)
# The factor correlation is bounded at 0.995: without the bound it can run to 1
# (the two factors are then the same factor) and mirt stops with an error.
spec2 <- mirt.model("F1 = 1-5\nF2 = 6-7\nCOV = F1*F2")
sv2   <- mirt(as.data.frame(Y), spec2, itemtype = "graded", pars = "values")
sv2$ubound[sv2$name == "COV_21"] <- 0.995
m2c <- mirt(as.data.frame(Y), spec2, itemtype = "graded", pars = sv2, verbose = FALSE)
m2e <- mirt(as.data.frame(Y), 2, itemtype = "graded", verbose = FALSE)
phi <- coef(m2c, simplify = TRUE)$cov[2, 1]
ic  <- function(m) c(logLik = extract.mirt(m, "logLik"), npar = extract.mirt(m, "nest"),
                     AIC = extract.mirt(m, "AIC"), BIC = extract.mirt(m, "BIC"))
lrt <- function(m) { x <- 2 * (extract.mirt(m, "logLik") - extract.mirt(m1, "logLik"))
                     df <- extract.mirt(m, "nest") - extract.mirt(m1, "nest")
                     c(LR = x, LR_df = df, LR_p = pchisq(x, df, lower.tail = FALSE)) }
grm_compare <- data.frame(
  model = c("GRM 1 factor", "GRM 2 correlated factors (C1-C5 | C6-C7)", "GRM 2 factors, exploratory"),
  rbind(c(ic(m1), LR = NA, LR_df = NA, LR_p = NA), c(ic(m2c), lrt(m2c)), c(ic(m2e), lrt(m2e))),
  factor_cor = c(NA, phi, NA))
cat("\nGRM, 1 vs 2 factors (LR test against 1 factor)\n")
print(cbind(grm_compare[1], round(grm_compare[-1], 3)), row.names = FALSE)
cat("Exploratory 2-factor loadings (oblimin)\n")
print(round(summary(m2e, rotate = "oblimin", verbose = FALSE)$rotF, 2))

# ---- 3. Ordinal CFA (lavaan, WLSMV) ----------------------------------------------
cfa_fit <- function(model) cfa(model, Yo, ordered = items, estimator = "WLSMV",
                               missing = "pairwise", std.lv = TRUE)
c1 <- cfa_fit("F =~ C1 + C2 + C3 + C4 + C5 + C6 + C7")
c2 <- cfa_fit("F1 =~ C1 + C2 + C3 + C4 + C5\nF2 =~ C6 + C7")
fm <- c("chisq.scaled", "df.scaled", "cfi.scaled", "tli.scaled", "rmsea.scaled",
        "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", "srmr")
cfa_compare <- data.frame(model = c("CFA 1 factor", "CFA 2 correlated factors (C1-C5 | C6-C7)"),
                          rbind(fitMeasures(c1, fm), fitMeasures(c2, fm)),
                          factor_cor = c(NA, lavInspect(c2, "cor.lv")[2, 1]))
cat("\nOrdinal CFA (WLSMV, pairwise polychorics)\n")
print(cbind(cfa_compare[1], round(cfa_compare[-1], 3)), row.names = FALSE)
cat("Standardised loadings, 1 factor:", round(lavInspect(c1, "std")$lambda[, 1], 2), "\n")
mi <- modindices(c1, sort. = TRUE); mi <- mi[mi$op == "~~" & mi$lhs != mi$rhs, ]
cat("Largest residual correlations suggested (modification indices):\n")
print(head(mi[c("lhs", "rhs", "mi", "sepc.all")], 5), row.names = FALSE)

write.csv(grm_compare, file.path(out, "data", sprintf("gait_dim_grm_1v2%s.csv", tag)), row.names = FALSE)
write.csv(cfa_compare, file.path(out, "data", sprintf("gait_dim_cfa_1v2%s.csv", tag)), row.names = FALSE)

# ---- 4. Local independence: Yen's Q3 under the 1-factor GRM ----------------------
q3 <- residuals(m1, type = "Q3", verbose = FALSE)
q3_pairs <- do.call(rbind, lapply(1:(n_items - 1), function(i) do.call(rbind, lapply((i + 1):n_items, function(j)
  data.frame(item1 = items[i], item2 = items[j], Q3 = q3[i, j], polychoric = R[i, j])))))
q3_mean <- mean(q3_pairs$Q3)
q3_pairs$Q3_minus_mean <- q3_pairs$Q3 - q3_mean
q3_pairs$flag <- q3_pairs$Q3_minus_mean > 0.2      # Christensen, Makransky & Horton (2017)
q3_pairs <- q3_pairs[order(-q3_pairs$Q3), ]
cat(sprintf("\nYen's Q3 (mean %.3f); flagged when Q3 - mean > 0.2\n", q3_mean))
print(cbind(q3_pairs[1:2], round(q3_pairs[3:5], 3), flag = q3_pairs$flag), row.names = FALSE)
write.csv(q3_pairs, file.path(out, "data", sprintf("gait_dim_q3%s.csv", tag)), row.names = FALSE)

# ---- 5. Fit of the 1-factor GRM ----------------------------------------------------
# C2 needs complete responses. The GRM is refitted to those persons first:
# dropping the incomplete rows (mostly walking-aid users, the most limited
# people) from the full-sample fit would compare a model for everyone with the
# data of a less limited subgroup.
m1_cc  <- mirt(as.data.frame(Y[complete.cases(Y), ]), 1, itemtype = "graded", verbose = FALSE)
c2stat <- M2(m1_cc, type = "C2")
cat(sprintf("\n1-factor GRM, C2 statistic (persons with all 7 items, n = %d)\n", sum(complete.cases(Y))))
print(round(c2stat, 3), row.names = FALSE)
ifit <- itemfit(m1, fit_stats = "S_X2", na.rm = TRUE)
ifit$p_holm <- p.adjust(ifit$p.S_X2, "holm")
cat("Item fit S-X2\n"); print(cbind(ifit[1], round(ifit[-1], 3)), row.names = FALSE)

cf <- coef(m1, simplify = TRUE, IRTpars = TRUE)$items
se <- do.call(rbind, lapply(coef(m1, IRTpars = TRUE, printSE = TRUE)[items], function(x) x["SE", ]))
grm_items <- data.frame(item = items, label = item_labels, a = cf[, "a"], a_se = se[, "a"],
                        b1 = cf[, "b1"], b1_se = se[, "b1"], b2 = cf[, "b2"], b2_se = se[, "b2"],
                        b3 = cf[, "b3"], b3_se = se[, "b3"])
cat("\n1-factor GRM item parameters (a = slope, b = thresholds on eta)\n")
print(cbind(grm_items[1], round(grm_items[-(1:2)], 2)), row.names = FALSE)
cat(sprintf("Empirical reliability of the EAP scores: %.3f | converged: %s\n",
            empirical_rxx(fscores(m1, method = "EAP", full.scores.SE = TRUE)),
            extract.mirt(m1, "converged")))

write.csv(data.frame(n_complete = sum(complete.cases(Y)), c2stat,
                     empirical_reliability = empirical_rxx(fscores(m1, method = "EAP", full.scores.SE = TRUE))),
          file.path(out, "data", sprintf("gait_dim_grm_fit%s.csv", tag)), row.names = FALSE)
write.csv(ifit, file.path(out, "data", sprintf("gait_dim_itemfit%s.csv", tag)), row.names = FALSE)
write.csv(grm_items, file.path(out, "data", sprintf("gait_dim_grm_items%s.csv", tag)), row.names = FALSE)
cat(sprintf("\nWrote %s/data/gait_dim_*%s.csv\n", out, tag))
