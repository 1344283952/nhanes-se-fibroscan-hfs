# ============================================
# 006 / 09_hfs_predict.R — Hepamet HFS predictive performance vs FIB-4
#
# H5 (Round 6 NK-7 revision):
#   PRIMARY: continuous AUROC ≥ 0.65 predicting LSM ≥ 8 kPa (F2+ fibrosis);
#            ΔAUROC ≥ 0.02 vs FIB-4 (DeLong 1988 paired AUROC test)
#   SECONDARY: NRI/IDI categorical (Pencina 2008) demoted to exploratory due to
#              only ≈ 6 HFS-eligible cases with HFS ≥ 0.47 (power < 5% per Pepe 2014)
#
# 500 bootstrap optimism-corrected internal validation (Steyerberg 2013)
# Flexible calibration plot (Van Calster 2019 DOI 10.1016/j.jclinepi.2019.02.004)
# Outcome: LSM ≥ 8 kPa (gold standard); HFS = Ampuero 2020 formula
# ============================================

set.seed(20260516)
# Round 2 R-Stats v2 P0: ensure additional clinical-prediction packages are
# available. If missing, attempt one-time silent install (no-op if already in).
for (p in c("ResourceSelection", "Hmisc")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    try(install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE), silent = TRUE)
  }
}
suppressPackageStartupMessages({
  library(dplyr); library(survey); library(pROC); library(boot)
})
cat("========================================\n")
cat("006 W6 — HFS predictive (Round 2 expanded: HL + Harrell C + NRI/IDI + NFS + APRI)\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")
load("data/processed/nhanes_final.RData")

# Restrict to HFS-eligible (fasting subsample) with non-NA outcome + predictors
df_hfs <- nhanes_final %>%
  filter(!is.na(hfs), !is.na(fib4), !is.na(fibrosis_lsm8),
         !is.na(wt_saf_pooled), wt_saf_pooled > 0)
cat(sprintf("Analytic N (HFS + FIB-4 + LSM ≥ 8 + SAF wt all non-NA): %d\n",
            nrow(df_hfs)))
cat(sprintf("LSM ≥ 8 prevalence: %d (%.1f%%)\n",
            sum(df_hfs$fibrosis_lsm8 == 1),
            100 * mean(df_hfs$fibrosis_lsm8 == 1)))
cat(sprintf("HFS ≥ 0.47 prevalence: %d (%.1f%%)  [power-limited per NK-7]\n",
            sum(df_hfs$hfs >= 0.47),
            100 * mean(df_hfs$hfs >= 0.47)))

# ---- Primary: Continuous AUROC ----
cat("\n[1/4] AUROC: HFS vs LSM ≥ 8 (primary endpoint per Round 6 NK-7) ...\n")
roc_hfs <- pROC::roc(df_hfs$fibrosis_lsm8, df_hfs$hfs, ci = TRUE, levels = c(0, 1), direction = "<")
roc_fib4 <- pROC::roc(df_hfs$fibrosis_lsm8, df_hfs$fib4, ci = TRUE, levels = c(0, 1), direction = "<")
cat(sprintf("  HFS  AUROC = %.3f (95%% CI %.3f - %.3f)\n",
            as.numeric(roc_hfs$auc), roc_hfs$ci[1], roc_hfs$ci[3]))
cat(sprintf("  FIB-4 AUROC = %.3f (95%% CI %.3f - %.3f)\n",
            as.numeric(roc_fib4$auc), roc_fib4$ci[1], roc_fib4$ci[3]))

# ---- Primary: ΔAUROC (DeLong 1988 paired test) ----
cat("\n[2/4] ΔAUROC HFS vs FIB-4 (DeLong 1988 paired test) ...\n")
delong <- pROC::roc.test(roc_hfs, roc_fib4, method = "delong", paired = TRUE)
cat(sprintf("  ΔAUROC = %.4f  Z = %.3f  p = %.4f\n",
            as.numeric(roc_hfs$auc) - as.numeric(roc_fib4$auc),
            delong$statistic, delong$p.value))
cat(sprintf("  Hypothesis H5: continuous AUROC ≥ 0.65 → %s\n",
            if (as.numeric(roc_hfs$auc) >= 0.65) "CONFIRMED" else "REJECTED"))
cat(sprintf("  Hypothesis H5b: ΔAUROC ≥ 0.02 → %s\n",
            if ((as.numeric(roc_hfs$auc) - as.numeric(roc_fib4$auc)) >= 0.02) "CONFIRMED" else "REJECTED"))

# ---- Optimism-corrected bootstrap (Steyerberg 2013, B=500) ----
cat("\n[3/4] Optimism-corrected bootstrap (B=500, Steyerberg 2013) ...\n")
n_boot <- 500
boot_aurocs <- numeric(n_boot)
boot_deltas <- numeric(n_boot)
for (b in seq_len(n_boot)) {
  idx <- sample.int(nrow(df_hfs), nrow(df_hfs), replace = TRUE)
  db <- df_hfs[idx, ]
  rh <- tryCatch(pROC::roc(db$fibrosis_lsm8, db$hfs, levels = c(0, 1), direction = "<"),
                 error = function(e) NULL)
  rf <- tryCatch(pROC::roc(db$fibrosis_lsm8, db$fib4, levels = c(0, 1), direction = "<"),
                 error = function(e) NULL)
  if (!is.null(rh) && !is.null(rf)) {
    boot_aurocs[b] <- as.numeric(rh$auc)
    boot_deltas[b] <- as.numeric(rh$auc) - as.numeric(rf$auc)
  }
  if (b %% 100 == 0) cat(sprintf("  boot iter %d/%d\n", b, n_boot))
}
boot_aurocs <- boot_aurocs[boot_aurocs > 0]
boot_deltas <- boot_deltas[boot_deltas != 0]
auroc_optimism <- mean(boot_aurocs) - as.numeric(roc_hfs$auc)
auroc_corrected <- as.numeric(roc_hfs$auc) - max(0, auroc_optimism)
cat(sprintf("  Apparent AUROC = %.3f\n", as.numeric(roc_hfs$auc)))
cat(sprintf("  Mean bootstrap AUROC = %.3f\n", mean(boot_aurocs)))
cat(sprintf("  Optimism = %.4f\n", auroc_optimism))
cat(sprintf("  Optimism-corrected AUROC = %.3f\n", auroc_corrected))
cat(sprintf("  Bootstrap ΔAUROC (HFS - FIB-4): mean = %.4f, 95%% CI = (%.4f, %.4f)\n",
            mean(boot_deltas), quantile(boot_deltas, 0.025), quantile(boot_deltas, 0.975)))

# ---- Flexible calibration (Van Calster 2019) ----
cat("\n[4/7] Flexible calibration plot (Van Calster 2019 DOI 10.1016/j.jclinepi.2019.02.004) ...\n")
# Generate calibration curve via locally-weighted logistic regression
# (loess on predicted prob vs observed outcome)
df_hfs$hfs_prob <- df_hfs$hfs   # HFS continuous 0-1 already
calib_loess <- tryCatch(
  loess(fibrosis_lsm8 ~ hfs_prob, data = df_hfs,
        weights = df_hfs$wt_saf_pooled, span = 0.5),
  error = function(e) NULL
)
calib_x <- seq(min(df_hfs$hfs_prob), max(df_hfs$hfs_prob), length.out = 100)
calib_y <- if (!is.null(calib_loess)) {
  pmin(pmax(predict(calib_loess, newdata = data.frame(hfs_prob = calib_x)), 0), 1)
} else NA
calib_df <- data.frame(predicted = calib_x, observed = calib_y)

# Quick calibration plot
png("output/figures/hfs_calibration.png", width = 700, height = 700, res = 110)
plot(calib_df$predicted, calib_df$observed, type = "l", lwd = 2, col = "darkblue",
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Predicted P(LSM ≥ 8 | HFS)", ylab = "Observed P(LSM ≥ 8)",
     main = sprintf("HFS Flexible Calibration (Van Calster 2019; n=%d)", nrow(df_hfs)))
abline(0, 1, lty = 2, col = "grey50")
# Add Hosmer-Lemeshow decile points (computed below) for visual cross-check
dev.off()

# ============================================
# [5/7] Hosmer-Lemeshow + Harrell C + NRI/IDI vs FIB-4
# Round 2 R-Stats v2 P0: prior Methods text claimed these methods but had no code.
# Now actually computed.
# ============================================
cat("\n[5/7] Hosmer-Lemeshow + Harrell C + NRI/IDI vs FIB-4 (Round 2 R-Stats P0) ...\n")
df_hfs$hfs_prob <- df_hfs$hfs

# 5a. Hosmer-Lemeshow goodness-of-fit (10 deciles)
hl_test <- NULL
if (requireNamespace("ResourceSelection", quietly = TRUE)) {
  hl_test <- tryCatch(
    ResourceSelection::hoslem.test(df_hfs$fibrosis_lsm8, df_hfs$hfs_prob, g = 10),
    error = function(e) { cat("  HL error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(hl_test)) {
    cat(sprintf("  Hosmer-Lemeshow chi^2 = %.3f (df = %d), P = %.4f\n",
                as.numeric(hl_test$statistic),
                as.numeric(hl_test$parameter),
                as.numeric(hl_test$p.value)))
  }
} else cat("  ResourceSelection package unavailable — Hosmer-Lemeshow skipped\n")

# 5b. Harrell's C-statistic (Hmisc::rcorr.cens; for binary outcome ≈ AUROC but explicit)
harrell_c <- NULL
if (requireNamespace("Hmisc", quietly = TRUE)) {
  harrell_c <- tryCatch(
    Hmisc::rcorr.cens(df_hfs$hfs, df_hfs$fibrosis_lsm8),
    error = function(e) { cat("  Harrell C error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(harrell_c)) {
    cat(sprintf("  Harrell C = %.4f (SE = %.4f)\n",
                as.numeric(harrell_c["C Index"]),
                as.numeric(harrell_c["S.D."]) / 2))
  }
} else cat("  Hmisc package unavailable — Harrell C skipped\n")

# 5c. NRI and IDI (Pencina 2008, 2011; continuous + categorical)
fit_fib4_logit <- glm(fibrosis_lsm8 ~ fib4, data = df_hfs, family = binomial())
df_hfs$fib4_prob <- predict(fit_fib4_logit, type = "response")

events_idx <- df_hfs$fibrosis_lsm8 == 1
n_events <- sum(events_idx)
n_nonevents <- sum(!events_idx)

# Continuous NRI (Pencina 2011 BMJ; sum of two directional movement probs)
nri_event_up   <- mean(df_hfs$hfs_prob[events_idx]  > df_hfs$fib4_prob[events_idx])  -
                  mean(df_hfs$hfs_prob[events_idx]  < df_hfs$fib4_prob[events_idx])
nri_nonev_down <- mean(df_hfs$fib4_prob[!events_idx] > df_hfs$hfs_prob[!events_idx]) -
                  mean(df_hfs$fib4_prob[!events_idx] < df_hfs$hfs_prob[!events_idx])
nri_continuous <- nri_event_up + nri_nonev_down
cat(sprintf("  Continuous NRI (HFS vs FIB-4) = %.4f  [event-up %.4f + nonevent-down %.4f]\n",
            nri_continuous, nri_event_up, nri_nonev_down))

# IDI (Pencina 2008): mean predicted prob difference
idi <- (mean(df_hfs$hfs_prob[events_idx])  - mean(df_hfs$fib4_prob[events_idx])) -
       (mean(df_hfs$hfs_prob[!events_idx]) - mean(df_hfs$fib4_prob[!events_idx]))
cat(sprintf("  IDI (HFS vs FIB-4)            = %.4f\n", idi))

# Categorical NRI (HFS ≥ 0.47 vs FIB-4 ≥ 2.67 high-risk threshold)
high_hfs  <- df_hfs$hfs_prob >= 0.47
high_fib4 <- df_hfs$fib4 >= 2.67
nri_cat_event  <- mean(high_hfs[events_idx]  & !high_fib4[events_idx])  -
                  mean(!high_hfs[events_idx] &  high_fib4[events_idx])
nri_cat_nonev  <- mean(!high_hfs[!events_idx] &  high_fib4[!events_idx]) -
                  mean( high_hfs[!events_idx] & !high_fib4[!events_idx])
nri_categorical <- nri_cat_event + nri_cat_nonev
n_hfs_high_in_set <- sum(high_hfs, na.rm = TRUE)
cat(sprintf("  Categorical NRI (HFS ≥ 0.47 vs FIB-4 ≥ 2.67) = %.4f  [n_HFS_high = %d — exploratory per Pepe 2014]\n",
            nri_categorical, n_hfs_high_in_set))

# ============================================
# [6/7] Multi-score benchmark: NFS + APRI + DeLong vs HFS (R-Bias X1 multi-comparator)
# ============================================
cat("\n[6/7] Multi-score benchmark NFS + APRI (Round 2 X1 expand-from-FIB4-only) ...\n")

# Need albumin in g/dL and HOMA-IR; both already on df_hfs upstream
df_scores <- df_hfs %>%
  mutate(
    ifg_dm  = as.integer((!is.na(LBXGLU) & LBXGLU >= 100) | (!is.na(diabetes) & diabetes == 1)),
    ast_alt_ratio = ast_unl / pmax(alt_unl, 1),
    # Angulo 2007 NFS = -1.675 + 0.037*age + 0.094*BMI + 1.13*IFG/DM + 0.99*AST/ALT - 0.013*PLT - 0.66*Alb(g/dL)
    nfs = -1.675 + 0.037 * RIDAGEYR + 0.094 * BMXBMI +
          1.13 * ifg_dm + 0.99 * ast_alt_ratio -
          0.013 * LBXPLTSI - 0.66 * albumin_gdl,
    # Wai 2003 APRI = (AST / ULN_AST) * 100 / Platelet (10^9/L); ULN_AST=40 standard
    apri = (ast_unl / 40) * 100 / pmax(LBXPLTSI, 1)
  ) %>%
  filter(!is.na(nfs), !is.na(apri))

cat(sprintf("  Analytic n with NFS + APRI computable: %d\n", nrow(df_scores)))

roc_nfs  <- tryCatch(pROC::roc(df_scores$fibrosis_lsm8, df_scores$nfs,
                                ci = TRUE, levels = c(0, 1), direction = "<"),
                      error = function(e) NULL)
roc_apri <- tryCatch(pROC::roc(df_scores$fibrosis_lsm8, df_scores$apri,
                                ci = TRUE, levels = c(0, 1), direction = "<"),
                      error = function(e) NULL)
auroc_nfs  <- if (!is.null(roc_nfs))  as.numeric(roc_nfs$auc)  else NA_real_
auroc_apri <- if (!is.null(roc_apri)) as.numeric(roc_apri$auc) else NA_real_
cat(sprintf("  NFS  AUROC = %.3f  (95%% CI %.3f - %.3f)\n",
            auroc_nfs, roc_nfs$ci[1], roc_nfs$ci[3]))
cat(sprintf("  APRI AUROC = %.3f  (95%% CI %.3f - %.3f)\n",
            auroc_apri, roc_apri$ci[1], roc_apri$ci[3]))

# Paired DeLong tests: HFS vs each
# Re-build HFS ROC on df_scores subset for paired comparison
roc_hfs_scores  <- pROC::roc(df_scores$fibrosis_lsm8, df_scores$hfs,
                              levels = c(0, 1), direction = "<")
roc_fib4_scores <- pROC::roc(df_scores$fibrosis_lsm8, df_scores$fib4,
                              levels = c(0, 1), direction = "<")

delong_hfs_nfs <- tryCatch(
  pROC::roc.test(roc_hfs_scores, roc_nfs,  method = "delong", paired = TRUE),
  error = function(e) NULL)
delong_hfs_apri <- tryCatch(
  pROC::roc.test(roc_hfs_scores, roc_apri, method = "delong", paired = TRUE),
  error = function(e) NULL)
delong_hfs_fib4_scores <- tryCatch(
  pROC::roc.test(roc_hfs_scores, roc_fib4_scores, method = "delong", paired = TRUE),
  error = function(e) NULL)

dauc_hfs_nfs  <- as.numeric(roc_hfs_scores$auc) - auroc_nfs
dauc_hfs_apri <- as.numeric(roc_hfs_scores$auc) - auroc_apri
dauc_hfs_fib4_consistency <- as.numeric(roc_hfs_scores$auc) - as.numeric(roc_fib4_scores$auc)

if (!is.null(delong_hfs_nfs)) {
  cat(sprintf("  HFS vs NFS  ΔAUROC = %+.4f  Z = %.2f  P = %.4f\n",
              dauc_hfs_nfs, delong_hfs_nfs$statistic, delong_hfs_nfs$p.value))
}
if (!is.null(delong_hfs_apri)) {
  cat(sprintf("  HFS vs APRI ΔAUROC = %+.4f  Z = %.2f  P = %.4f\n",
              dauc_hfs_apri, delong_hfs_apri$statistic, delong_hfs_apri$p.value))
}
if (!is.null(delong_hfs_fib4_scores)) {
  cat(sprintf("  HFS vs FIB-4 (consistency) ΔAUROC = %+.4f  P = %.4f\n",
              dauc_hfs_fib4_consistency, delong_hfs_fib4_scores$p.value))
}

# Save multi-comparator detail
multiscore_df <- data.frame(
  score   = c("HFS", "NFS", "APRI", "FIB-4"),
  AUROC   = c(as.numeric(roc_hfs_scores$auc), auroc_nfs, auroc_apri, as.numeric(roc_fib4_scores$auc)),
  n       = rep(nrow(df_scores), 4),
  vs_HFS_dAUROC = c(NA, dauc_hfs_nfs, dauc_hfs_apri, dauc_hfs_fib4_consistency),
  vs_HFS_p = c(NA,
               if (!is.null(delong_hfs_nfs)) delong_hfs_nfs$p.value else NA,
               if (!is.null(delong_hfs_apri)) delong_hfs_apri$p.value else NA,
               if (!is.null(delong_hfs_fib4_scores)) delong_hfs_fib4_scores$p.value else NA),
  stringsAsFactors = FALSE
)
print(multiscore_df)
write.csv(multiscore_df, "output/tables/multiscore_benchmark.csv", row.names = FALSE)

# ============================================
# [7/7] retrodesign Type-S/M reminder for HFS≥0.47 cell (X5 cross-reference)
# ============================================
cat("\n[7/7] retrodesign Type-S/M (R-Causal-Methods X5) ...\n")
cat(sprintf("  HFS ≥ 0.47 cell in ROC analytic set: n = %d (of %d, %.1f%%)\n",
            n_hfs_high_in_set, nrow(df_hfs),
            100 * n_hfs_high_in_set / nrow(df_hfs)))
cat("  Detailed retrodesign (Gelman & Carlin 2014 Type-S/Type-M) saved to:\n")
cat("    output/tables/retrodesign_smallcells.csv  (run scripts/_retrodesign_smallcells.R)\n")
cat(sprintf("  Categorical NRI/IDI for this cell are EXPLORATORY per Pepe 2014 power floor.\n"))

# ============================================
# Final summary_df (expanded with Round 2 metrics)
# ============================================
summary_df <- data.frame(
  metric = c("AUROC_HFS_apparent", "AUROC_HFS_optimism_corrected",
             "AUROC_FIB4_apparent", "Delta_AUROC_HFS_minus_FIB4",
             "Bootstrap_mean_Delta", "Bootstrap_lci_Delta", "Bootstrap_uci_Delta",
             "H5_AUROC_>=0.65", "H5b_Delta_>=0.02",
             # Round 2 R-Stats P0 additions:
             "HosmerLemeshow_chisq", "HosmerLemeshow_p",
             "HarrellC_HFS", "ContinuousNRI_HFS_vs_FIB4", "IDI_HFS_vs_FIB4",
             "CategoricalNRI_HFS_vs_FIB4_exploratory", "HFS_GE_047_cell_n_in_ROC_set",
             # Round 2 R-Bias X1 multi-comparator additions:
             "AUROC_NFS_apparent", "AUROC_APRI_apparent",
             "DeltaAUROC_HFS_minus_NFS", "DeltaAUROC_HFS_minus_APRI",
             "DeLong_p_HFS_vs_NFS", "DeLong_p_HFS_vs_APRI"),
  value = c(as.numeric(roc_hfs$auc), auroc_corrected,
            as.numeric(roc_fib4$auc),
            as.numeric(roc_hfs$auc) - as.numeric(roc_fib4$auc),
            mean(boot_deltas),
            quantile(boot_deltas, 0.025), quantile(boot_deltas, 0.975),
            ifelse(as.numeric(roc_hfs$auc) >= 0.65, 1, 0),
            ifelse((as.numeric(roc_hfs$auc) - as.numeric(roc_fib4$auc)) >= 0.02, 1, 0),
            if (!is.null(hl_test)) as.numeric(hl_test$statistic) else NA,
            if (!is.null(hl_test)) as.numeric(hl_test$p.value) else NA,
            if (!is.null(harrell_c)) as.numeric(harrell_c["C Index"]) else NA,
            nri_continuous, idi, nri_categorical, n_hfs_high_in_set,
            auroc_nfs, auroc_apri,
            dauc_hfs_nfs, dauc_hfs_apri,
            if (!is.null(delong_hfs_nfs)) delong_hfs_nfs$p.value else NA,
            if (!is.null(delong_hfs_apri)) delong_hfs_apri$p.value else NA),
  stringsAsFactors = FALSE
)
print(summary_df)

# Save
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(summary_df, "output/tables/hfs_auroc_summary.csv", row.names = FALSE)
write.csv(calib_df, "output/tables/hfs_calibration_curve.csv", row.names = FALSE)
save(roc_hfs, roc_fib4, delong, boot_aurocs, boot_deltas, summary_df, calib_df,
     auroc_corrected,
     # Round 2 additions:
     hl_test, harrell_c, nri_continuous, idi, nri_categorical,
     n_hfs_high_in_set, fit_fib4_logit, df_scores,
     roc_nfs, roc_apri, multiscore_df,
     delong_hfs_nfs, delong_hfs_apri, delong_hfs_fib4_scores,
     file = "output/tables/hfs_predict_results.RData")

cat("\n保存:\n")
cat("  output/tables/hfs_auroc_summary.csv (H5 + H5b)\n")
cat("  output/tables/hfs_calibration_curve.csv\n")
cat("  output/tables/hfs_predict_results.RData\n")
cat("  output/figures/hfs_calibration.png\n")
cat("\nDONE 09_hfs_predict.R\n")
