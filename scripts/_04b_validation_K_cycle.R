# ============================================
# 006v2 / _04b_validation_K_cycle.R
# 在 NHANES K cycle (2021-2023) 复制 P_ 主分析 3 个核心结果:
#   1. Se→CAP RCS U-shape 非线性 LR-test (P_: 1.27e-11)
#   2. Se/Zn T3 vs T1 OR for CAP>=275 (P_: OR 1.22 (1.10-1.35) P=0.0025)
#   3. HFS AUROC vs LSM>=8 (P_: 0.731)
#
# Round 1 R-NHANES P0 FIX (2026-05-23, GAMMA):
#   K cycle 是 complex multistage probability sample, 必须用 WTMEC2YR +
#   SDMVPSU + SDMVSTRA. 上一版用 unweighted lm/glm/pROC 违反 project
#   CLAUDE.md "all statistical analyses must use survey weighted functions".
#   重写 — H2: svyglm + regTermTest cubic (rcs 5-knot/ns 在 design df=15 下
#   失败因 variance matrix degenerate); H3: svyglm quasibinomial;
#   H5: WeightedROC + PSU-cluster bootstrap; Hg/Se: svyquantile.
#
# 输出: output/tables/k_cycle_validation.csv (含 unweighted/weighted 双列)
# ============================================
library(dplyr); library(survey); library(pROC); library(WeightedROC)
set.seed(20260523)

options(survey.lonely.psu = "adjust")  # 与 P_ 主分析 04_survey_design.R 一致

cat("==== K cycle 复制验证分析 (weighted, R-NHANES P0 fix) ====\n")
load("data/processed/nhanes_k_final.RData")
d <- nhanes_k_final
cat(sprintf("K cohort n=%d (raw)\n", nrow(d)))

# Factor 化
d$SMQ020_f <- factor(d$SMQ020)
d$ALQ111_f <- factor(d$ALQ111)

# Hg/Se molar ratio
d$hg_se_molar <- ifelse(!is.na(d$hg_ugl) & !is.na(d$LBXBSE) & d$LBXBSE > 0,
                        (d$hg_ugl / 200.59) / (d$LBXBSE / 78.96), NA)

# Se centred + scaled, plus polynomial columns (避免 rcs/ns 在 svyglm regTermTest
# variance matrix degenerate 的问题)
se_med <- median(d$LBXBSE, na.rm = TRUE)
se_sd  <- sd(d$LBXBSE, na.rm = TRUE)
d$Se_z  <- (d$LBXBSE - se_med) / se_sd
d$Se_z2 <- d$Se_z^2
d$Se_z3 <- d$Se_z^3

# 共变量字符串
covs_str <- "age + sex_male + race + DMDEDUC2 + pir + SMQ020_f + ALQ111_f + kcal"

# === Survey design 构造 (K cycle 用 WTMEC2YR) ===
# 与 P_ 主分析对齐: ids = ~SDMVPSU, strata = ~SDMVSTRA, nest=TRUE
k_dsn <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~WTMEC2YR, data = d, nest = TRUE
)
cat(sprintf("k_dsn (full K cohort, WTMEC2YR): %d rows, design df=%d\n",
            nrow(d), degf(k_dsn)))

# Cohort subset (subset 保留 design 结构)
k_dsn_h2 <- subset(k_dsn, !is.na(LBXBSE) & !is.na(cap))            # H2: Se → CAP
k_dsn_h3 <- subset(k_dsn, !is.na(se_zn_ratio_diet) & !is.na(steatosis_cap275))  # H3
k_dsn_h5 <- subset(k_dsn, !is.na(hfs) & !is.na(fibrosis_lsm8))    # H5
cat(sprintf("k_dsn_h2 (Se+CAP):           n=%d\n", nrow(k_dsn_h2$variables)))
cat(sprintf("k_dsn_h3 (Se/Zn+CAP275):     n=%d\n", nrow(k_dsn_h3$variables)))
cat(sprintf("k_dsn_h5 (HFS+LSM8):         n=%d\n\n", nrow(k_dsn_h5$variables)))

# === 1. H2 Se → CAP cubic U-shape (svyglm + regTermTest, weighted) ===
cat("--- 1. H2 Se → CAP cubic non-linear (svyglm, weighted) ---\n")
# Note: rcs(5)/ns(4) 在 design df=15 下 regTermTest 变量协方差矩阵 degenerate
# (Lumley svyglm + 5+ knot multicollinearity 已知限制); 用 z-centred 立方多项式
# 替代是保守做法 — 仍可测严格的"非线性 (drop Se^2 + Se^3)" Wald F-test.
m_full_w <- svyglm(
  as.formula(paste("cap ~ Se_z + Se_z2 + Se_z3 +", covs_str)),
  design = k_dsn_h2
)

rt_nl <- regTermTest(m_full_w, ~Se_z2 + Se_z3,
                      method = "Wald", df = degf(k_dsn_h2))
p_nonlin_w <- as.numeric(rt_nl$p)
f_stat_w   <- as.numeric(rt_nl$Ftest)
cat(sprintf("  K cycle weighted Se cubic non-linear (drop Se^2 + Se^3) F=%.2f P_nl=%.3e\n",
            f_stat_w, p_nonlin_w))
cat(sprintf("  (P_ 主分析 RCS 5-knot: 1.27e-11; unweighted K cubic 等价: 4.6e-05)\n"))

# Nadir estimation — predict on Se 80-250 µg/L grid with weighted median/mode covars
new_grid <- data.frame(Se_z = seq((80 - se_med)/se_sd,
                                  (250 - se_med)/se_sd, length.out = 171))
new_grid$Se_z2 <- new_grid$Se_z^2
new_grid$Se_z3 <- new_grid$Se_z^3
for (cv in c("age","sex_male","pir","kcal","DMDEDUC2")) {
  qv <- tryCatch({
    q <- svyquantile(as.formula(paste0("~", cv)), k_dsn_h2, 0.5, na.rm = TRUE)
    if (is.list(q)) q <- q[[1]]
    # svyquantile 4.5 returns 1x4 matrix even for single quantile
    if (is.matrix(q)) as.numeric(q[1, "quantile"]) else as.numeric(q)
  }, error = function(e) median(d[[cv]], na.rm = TRUE))
  new_grid[[cv]] <- qv
}
for (cv in c("race","SMQ020_f","ALQ111_f")) {
  v <- d[[cv]]
  tab <- tryCatch(svytable(as.formula(paste0("~", cv)), k_dsn_h2),
                  error = function(e) table(v))
  mode_val <- names(sort(tab, decreasing = TRUE))[1]
  new_grid[[cv]] <- factor(mode_val, levels = levels(factor(v)))
}
pred_w <- predict(m_full_w, newdata = new_grid)
pred_vec <- as.numeric(as.matrix(pred_w))
new_grid$cap_pred <- pred_vec
nadir_idx <- which.min(new_grid$cap_pred)
nadir_se_w <- new_grid$Se_z[nadir_idx] * se_sd + se_med
cat(sprintf("  K cycle weighted CAP nadir at Se=%.1f µg/L (P_ nadir window: 130-170 µg/L)\n",
            nadir_se_w))

# === 2. H3 Se/Zn T3 vs T1 logistic for CAP>=275 (svyglm, weighted) ===
cat("\n--- 2. H3 Se/Zn T3 vs T1 logistic (svyglm, quasibinomial weighted) ---\n")

# Tertile 划分 (用未加权 ntile, 与主分析一致)
d$se_zn_t <- factor(ntile(d$se_zn_ratio_diet, 3), levels = c(1, 2, 3),
                    labels = c("T1","T2","T3"))
d$se_zn_t <- relevel(d$se_zn_t, ref = "T1")

# 重建 design (新加 se_zn_t 列)
k_dsn2 <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA,
  weights = ~WTMEC2YR, data = d, nest = TRUE
)
k_dsn2_h3 <- subset(k_dsn2, !is.na(se_zn_t) & !is.na(steatosis_cap275))
n_h3 <- nrow(k_dsn2_h3$variables)
cat(sprintf("  k_dsn2_h3 (Se/Zn tertile + CAP275 non-NA): n=%d\n", n_h3))

m_zn_w <- svyglm(
  as.formula(paste("steatosis_cap275 ~ se_zn_t +", covs_str)),
  design = k_dsn2_h3, family = quasibinomial()
)
cf_w <- summary(m_zn_w)$coefficients
or_t3_w <- ci_lo_w <- ci_hi_w <- p_t3_w <- NA_real_
if ("se_zn_tT3" %in% rownames(cf_w)) {
  cf_t3_w <- cf_w["se_zn_tT3", ]
  or_t3_w <- exp(cf_t3_w[1])
  ci_lo_w <- exp(cf_t3_w[1] - 1.96 * cf_t3_w[2])
  ci_hi_w <- exp(cf_t3_w[1] + 1.96 * cf_t3_w[2])
  p_t3_w  <- cf_t3_w[4]
  cat(sprintf("  K cycle weighted Se/Zn T3 vs T1: OR=%.2f (95%% CI %.2f-%.2f), P=%.4f\n",
              or_t3_w, ci_lo_w, ci_hi_w, p_t3_w))
  cat(sprintf("  P_ 主分析: 1.22 (1.10-1.35), P=0.0025\n"))
}

# === 3. H5 HFS AUROC vs LSM>=8 (WeightedROC, weighted) ===
cat("\n--- 3. H5 HFS AUROC vs LSM >= 8 (WeightedROC, survey-weighted) ---\n")
d_hfs_full <- d %>% filter(!is.na(hfs), !is.na(fibrosis_lsm8))
cat(sprintf("  HFS subset n=%d (with WTMEC2YR)\n", nrow(d_hfs_full)))

roc_w <- WeightedROC(
  guess  = d_hfs_full$hfs,
  label  = d_hfs_full$fibrosis_lsm8,
  weight = d_hfs_full$WTMEC2YR
)
auroc_w <- WeightedAUC(roc_w)
cat(sprintf("  K cycle weighted HFS AUROC (point estimate)=%.3f\n", auroc_w))

# PSU-cluster bootstrap 95% CI
set.seed(20260523)
B <- 500
boot_aucs <- numeric(B)
psu_keys  <- paste0(d_hfs_full$SDMVSTRA, "_", d_hfs_full$SDMVPSU)
psu_unique <- unique(psu_keys)
for (b in seq_len(B)) {
  picks <- sample(psu_unique, length(psu_unique), replace = TRUE)
  idx_b <- unlist(lapply(picks, function(p) which(psu_keys == p)))
  rb <- tryCatch({
    WeightedROC(d_hfs_full$hfs[idx_b],
                d_hfs_full$fibrosis_lsm8[idx_b],
                d_hfs_full$WTMEC2YR[idx_b])
  }, error = function(e) NULL)
  boot_aucs[b] <- if (!is.null(rb)) WeightedAUC(rb) else NA_real_
}
boot_aucs <- boot_aucs[!is.na(boot_aucs)]
ci_lo_auc <- as.numeric(quantile(boot_aucs, 0.025, na.rm = TRUE))
ci_hi_auc <- as.numeric(quantile(boot_aucs, 0.975, na.rm = TRUE))
cat(sprintf("  K cycle weighted AUROC=%.3f (95%% CI %.3f-%.3f, PSU-cluster boot B=%d)\n",
            auroc_w, ci_lo_auc, ci_hi_auc, length(boot_aucs)))
cat(sprintf("  P_ 主分析: 0.731 (0.703-0.759)\n"))

# === 4. ln(Se) → GGT_high (svyglm, weighted) ===
cat("\n--- 4. ln(Se) → GGT_high logistic (svyglm, weighted) ---\n")
k_dsn_ggt <- subset(k_dsn, !is.na(ggt_high) & !is.na(LBXBSE))
n_ggt <- nrow(k_dsn_ggt$variables)
m_ggt_w <- svyglm(
  as.formula(paste("ggt_high ~ ln_se_blood +", covs_str)),
  design = k_dsn_ggt, family = quasibinomial()
)
cf_ggt_w <- summary(m_ggt_w)$coefficients
or_ggt_w <- ci_lo_ggt <- ci_hi_ggt <- p_ggt_w <- NA_real_
if ("ln_se_blood" %in% rownames(cf_ggt_w)) {
  cf_se_w <- cf_ggt_w["ln_se_blood", ]
  or_ggt_w  <- exp(cf_se_w[1])
  ci_lo_ggt <- exp(cf_se_w[1] - 1.96 * cf_se_w[2])
  ci_hi_ggt <- exp(cf_se_w[1] + 1.96 * cf_se_w[2])
  p_ggt_w   <- cf_se_w[4]
  cat(sprintf("  K cycle weighted ln(Se) → GGT_high: OR=%.2f (95%% CI %.2f-%.2f) P=%.4f (n=%d)\n",
              or_ggt_w, ci_lo_ggt, ci_hi_ggt, p_ggt_w, n_ggt))
}

# === 5. Hg/Se molar ratio (svyquantile, weighted) ===
cat("\n--- 5. Hg/Se molar ratio (svyquantile, weighted) ---\n")
k_dsn_hgse <- subset(k_dsn, !is.na(hg_se_molar))
n_hgse <- nrow(k_dsn_hgse$variables)
hgse_q <- svyquantile(~hg_se_molar, k_dsn_hgse,
                      quantiles = c(0.25, 0.5, 0.75), na.rm = TRUE)
# survey 4.5 svyquantile returns list of 3x4 matrix (rows=quantile, cols=quantile/ci/se)
if (is.list(hgse_q)) hgse_q <- hgse_q[[1]]
hgse_med_w <- as.numeric(hgse_q["0.5",  "quantile"])
hgse_lo_w  <- as.numeric(hgse_q["0.25", "quantile"])
hgse_hi_w  <- as.numeric(hgse_q["0.75", "quantile"])
cat(sprintf("  K cycle weighted Hg/Se molar median=%.3f IQR=%.3f-%.3f (n=%d)\n",
            hgse_med_w, hgse_lo_w, hgse_hi_w, n_hgse))
cat(sprintf("  P_ 主分析: median 0.001, IQR 0.001-0.003\n"))

# === 输出 — schema (unweighted vs weighted 双列) ===
validation_summary <- data.frame(
  Test = c("1. Se→CAP non-linear P_nl",
           "2. Se/Zn T3 vs T1 OR for CAP≥275",
           "3. HFS AUROC vs LSM≥8",
           "4. ln(Se)→GGT_high OR",
           "5. Hg/Se molar median"),
  P_cycle_main = c("1.27e-11 (highly sig)",
                   "1.22 (1.10-1.35), P=0.0025",
                   "0.731 (0.703-0.759)",
                   "新分析 (v2 待 K 复制)",
                   "0.001 (IQR 0.001-0.003)"),
  K_cycle_unweighted_OLD = c("4.640e-05 (F=7.60)",
                             "0.95 (0.80-1.13), P=0.5706",
                             "0.642 (0.602-0.683)",
                             "OR=1.07, P=0.8268",
                             "0.001 (IQR 0.001-0.003)"),
  K_cycle_weighted_NEW = c(
    sprintf("%.3e (F=%.2f)", p_nonlin_w, f_stat_w),
    sprintf("%.2f (%.2f-%.2f), P=%.4f", or_t3_w, ci_lo_w, ci_hi_w, p_t3_w),
    sprintf("%.3f (%.3f-%.3f)", auroc_w, ci_lo_auc, ci_hi_auc),
    sprintf("OR=%.2f (%.2f-%.2f), P=%.4f", or_ggt_w, ci_lo_ggt, ci_hi_ggt, p_ggt_w),
    sprintf("%.3f (%.3f-%.3f)", hgse_med_w, hgse_lo_w, hgse_hi_w)),
  Note = c("svyglm cubic + Wald F (df=15)",
           "svyglm quasibinomial (WTMEC2YR)",
           "WeightedROC + PSU-cluster bootstrap (B=500)",
           "svyglm quasibinomial (WTMEC2YR)",
           "svyquantile (WTMEC2YR)")
)

write.csv(validation_summary, "output/tables/k_cycle_validation.csv", row.names = FALSE)
cat("\n=== K cycle weighted vs unweighted (R-NHANES P0 fix) ===\n")
print(validation_summary, row.names = FALSE)
cat("\n[OK] saved output/tables/k_cycle_validation.csv\n")
cat(sprintf("\n[OK] _04b_validation_K_cycle.R (weighted) done at %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# Save weighted nadir for downstream use
weighted_summary <- list(
  n_full = nrow(d),
  n_h2   = nrow(k_dsn_h2$variables),
  n_h3   = n_h3,
  n_h5   = nrow(d_hfs_full),
  n_ggt  = n_ggt,
  n_hgse = n_hgse,
  h2_p_nl = p_nonlin_w, h2_F = f_stat_w, h2_nadir = nadir_se_w,
  h3_or = or_t3_w, h3_ci_lo = ci_lo_w, h3_ci_hi = ci_hi_w, h3_p = p_t3_w,
  h5_auroc = auroc_w, h5_ci_lo = ci_lo_auc, h5_ci_hi = ci_hi_auc,
  ggt_or = or_ggt_w, ggt_p = p_ggt_w,
  hgse_med = hgse_med_w, hgse_lo = hgse_lo_w, hgse_hi = hgse_hi_w
)
saveRDS(weighted_summary, "output/tables/k_cycle_weighted_summary.rds")
