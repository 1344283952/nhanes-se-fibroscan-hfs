# ============================================================
# 006v2 / 20_hgse_antagonism.R — Hg-Se 拮抗专项 (经典 redox 故事)
# Ralston 2017 Selenium-Health-Benefit Value: Hg/Se 摩尔比 > 1 = 净毒,< 1 = Se 保护
# 2026-05-18 v2 升级
# ============================================================
library(dplyr); library(survey); library(rms); library(splines)
set.seed(20260518)

cat("==== 20_hgse_antagonism ====\n")
load("data/processed/nhanes_final.RData")

d <- nhanes_final %>%
  filter(!is.na(hg_ugl), !is.na(se_ugl), !is.na(hg_se_molar_ratio),
         !is.na(RIDAGEYR), !is.na(sex_male))
cat(sprintf("Analytic n: %d\n", nrow(d)))

# 描述
cat(sprintf("\nHg µg/L: median=%.2f IQR=%.2f-%.2f\n",
            median(d$hg_ugl), quantile(d$hg_ugl,0.25), quantile(d$hg_ugl,0.75)))
cat(sprintf("Se µg/L: median=%.1f IQR=%.1f-%.1f\n",
            median(d$se_ugl), quantile(d$se_ugl,0.25), quantile(d$se_ugl,0.75)))
cat(sprintf("Hg/Se molar ratio: median=%.3f IQR=%.3f-%.3f\n",
            median(d$hg_se_molar_ratio),
            quantile(d$hg_se_molar_ratio,0.25),
            quantile(d$hg_se_molar_ratio,0.75)))
cat(sprintf("  Hg/Se > 1 (净毒, Ralston 2017): %d (%.1f%%)\n",
            sum(d$hg_se_molar_ratio > 1, na.rm=TRUE),
            100*mean(d$hg_se_molar_ratio > 1, na.rm=TRUE)))

# 1. Hg × Se 二元交互项 (logistic for CAP>=275, LSM>=8, HFS>=0.47)
cat("\n=== 1. Hg × Se 二元交互项 (logistic) ===\n")
covs <- c("RIDAGEYR","sex_male","race","education","pir","SMQ020","ALQ111","DR1TKCAL")
results_int <- list()
for (oc in c("steatosis_cap275","fibrosis_lsm8","hfs_high","ggt_high")) {
  if (!(oc %in% colnames(d))) next
  # Main effects + interaction
  fml <- as.formula(paste(oc, "~ ln_hg * ln_se_blood +",
                          paste(covs, collapse = "+")))
  fit <- tryCatch(glm(fml, data = d, family = binomial()),
                  error = function(e) NULL)
  if (is.null(fit)) next
  cf <- summary(fit)$coefficients
  if ("ln_hg:ln_se_blood" %in% rownames(cf)) {
    row <- cf["ln_hg:ln_se_blood", ]
    results_int[[oc]] <- data.frame(
      outcome = oc,
      beta_interaction = row[1],
      se = row[2],
      z = row[3],
      p = row[4],
      OR_per_logHgSe = exp(row[1]),
      ci_lo = exp(row[1] - 1.96*row[2]),
      ci_hi = exp(row[1] + 1.96*row[2]))
  }
}
int_df <- if (length(results_int) > 0) do.call(rbind, results_int) else data.frame()
if (nrow(int_df) > 0) {
  int_df$p_BH <- p.adjust(int_df$p, "BH")
  print(int_df, row.names = FALSE)
}

# 2. Hg/Se 摩尔比 RCS (5-knot)
cat("\n=== 2. Hg/Se 摩尔比 RCS (5-knot) → 多 outcomes ===\n")
covs_str <- paste(c("RIDAGEYR","sex_male","DR1TKCAL"), collapse = "+")
results_rcs <- list()
for (oc in c("cap","lsm","ln_ggt","hfs")) {
  if (!(oc %in% colnames(d))) next
  fml <- as.formula(sprintf("%s ~ rcs(hg_se_molar_ratio, 5) + %s", oc, covs_str))
  m_full <- tryCatch(lm(fml, data = d), error = function(e) NULL)
  if (is.null(m_full)) next
  fml_lin <- as.formula(sprintf("%s ~ hg_se_molar_ratio + %s", oc, covs_str))
  m_lin <- lm(fml_lin, data = d)
  # ANOVA 非线性检验
  anv <- anova(m_lin, m_full)
  p_overall_F <- summary(m_full)$fstatistic
  p_overall <- pf(p_overall_F[1], p_overall_F[2], p_overall_F[3], lower.tail = FALSE)
  results_rcs[[oc]] <- data.frame(
    outcome = oc,
    F_stat_nonlinear = anv$F[2],
    p_nonlinear = anv$`Pr(>F)`[2],
    p_overall = unname(p_overall),
    n = nrow(model.frame(m_full)))
}
rcs_df <- if (length(results_rcs) > 0) do.call(rbind, results_rcs) else data.frame()
if (nrow(rcs_df) > 0) {
  rcs_df$p_BH <- p.adjust(rcs_df$p_nonlinear, "BH")
  print(rcs_df, row.names = FALSE)
}

# 3. Stratified analysis: Hg/Se median split
cat("\n=== 3. Median Hg/Se split: 高 Hg/Se 组 vs 低组 ===\n")
median_hgse <- median(d$hg_se_molar_ratio, na.rm = TRUE)
d$hgse_high <- as.integer(d$hg_se_molar_ratio > median_hgse)
cat(sprintf("Median Hg/Se = %.3f (split at 50th percentile)\n", median_hgse))
for (oc in c("steatosis_cap275","fibrosis_lsm8","ggt_high")) {
  if (!(oc %in% colnames(d))) next
  tab <- table(d[[oc]], d$hgse_high, useNA = "no")
  cat(sprintf("\n%s × Hg/Se high:\n", oc))
  print(tab)
  ct <- chisq.test(tab)
  cat(sprintf("  Chi-sq: %.2f (df=%d), P=%.4f\n", ct$statistic, ct$parameter, ct$p.value))
}

# 4. Save
out <- list(interaction = int_df, rcs = rcs_df,
            median_hgse_ratio = median_hgse,
            n_analytic = nrow(d))
save(out, file = "output/tables/hgse_antagonism.RData")
write.csv(int_df, "output/tables/hgse_interaction.csv", row.names = FALSE)
write.csv(rcs_df, "output/tables/hgse_rcs.csv",         row.names = FALSE)
cat("\n[OK] saved hgse_antagonism.RData / hgse_interaction.csv / hgse_rcs.csv\n")
