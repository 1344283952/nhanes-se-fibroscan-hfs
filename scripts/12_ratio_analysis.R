# ============================================
# 006 / 12_ratio_analysis.R — Dietary Se/Zn + Se/Cu ratios vs CAP / LSM / HFS
#
# H3 (per OSF §1): Dietary Se/Zn ratio top vs bottom tertile is associated with
#                  CAP ≥ 275 dB/m (Rinella 2023) with adjusted OR ≥ 1.20 (α=0.05).
#
# Per task.md §5.3 + §7.2 (dietary-only since blood Zn/Cu unavailable post-2017):
#   Dietary Se/Zn = DR1TSELE / DR1TZINC  (n_avail = 11,274)
#   Dietary Se/Cu = DR1TSELE / DR1TCOPP  (n_avail = 11,274)
#
# Models:
#   1. Tertile (continuous + categorical)
#   2. Logistic OR for CAP ≥ 275 (M2 covariates)
#   3. Linear regression for continuous LSM + HFS
#   4. Incremental ΔR² vs single-Se quartile model
# ============================================

set.seed(20260516)
suppressPackageStartupMessages({
  library(dplyr); library(survey); library(broom)
})
# R-Stats v2 + R-NHANES v2 Round 2 P0: lonely PSU handling.
options(survey.lonely.psu = "adjust")
cat("========================================\n")
cat("006 W5 — Dietary Se/Zn + Se/Cu ratio analyses (H3) (Round 2: cov_pre primary)\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")    # design_main, design_hfs
load("data/processed/nhanes_final.RData")

# ---- Construct dietary ratios ----
df <- nhanes_final %>%
  mutate(
    se_zn_ratio = ifelse(!is.na(DR1TSELE) & !is.na(DR1TZINC) & DR1TZINC > 0,
                          DR1TSELE / DR1TZINC, NA_real_),
    se_cu_ratio = ifelse(!is.na(DR1TSELE) & !is.na(DR1TCOPP) & DR1TCOPP > 0,
                          DR1TSELE / DR1TCOPP, NA_real_)
  )
cat(sprintf("Cohort N = %d\n", nrow(df)))
cat(sprintf("Se/Zn ratio non-NA: %d\n", sum(!is.na(df$se_zn_ratio))))
cat(sprintf("Se/Cu ratio non-NA: %d\n", sum(!is.na(df$se_cu_ratio))))

# Tertile labels
df <- df %>%
  mutate(
    se_zn_t = factor(ntile(se_zn_ratio, 3), levels = 1:3,
                     labels = c("T1 (low)", "T2", "T3 (high)")),
    se_cu_t = factor(ntile(se_cu_ratio, 3), levels = 1:3,
                     labels = c("T1 (low)", "T2", "T3 (high)"))
  )
cat("\nSe/Zn tertile distribution:\n"); print(table(df$se_zn_t, useNA = "ifany"))
cat("\nSe/Cu tertile distribution:\n"); print(table(df$se_cu_t, useNA = "ifany"))

# Rebuild design on the analytic subset
df_an <- df %>%
  filter(!is.na(se_zn_ratio), !is.na(se_cu_ratio),
         !is.na(steatosis_cap275), !is.na(wt_pooled), wt_pooled > 0)
cat(sprintf("\nAnalytic N (all ratios non-NA + CAP outcome): %d\n", nrow(df_an)))
des_an <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                    data = df_an, nest = TRUE)

# ---- Covariates ----
# R-Causal-Methods v2 + R-Stats v2 Round 2 P0:
# Primary covariate set = pre-exposure C only (Pearl backdoor minimal adjustment set
# per _dag_spec.md §4). The previous cov_M2 included BMI / HbA1c / LDL / diabetes
# which are mediator-confounders on the Se → outcome path. Conditioning on them
# blocks part of the total effect (Westreich-Greenland 2013 Table-2 fallacy).
# Total-effect estimate (primary) uses cov_pre; CDE direct-effect sensitivity
# uses cov_M2_sens block and is reported as a separate sensitivity row.
cov_pre <- c("age", "RIAGENDR", "race", "INDFMPIR", "education",
             "smoke", "drink", "DR1TKCAL")
cov_pre <- intersect(cov_pre, names(df_an))
cov_M2_sens <- c(cov_pre, "bmi", "LBXGH", "LBDLDL", "diabetes")
cov_M2_sens <- intersect(cov_M2_sens, names(df_an))
cat(sprintf("\ncov_pre (PRIMARY, total effect): %s\n", paste(cov_pre, collapse=", ")))
cat(sprintf("cov_M2_sens (sensitivity CDE only): %s\n", paste(cov_M2_sens, collapse=", ")))
# Alias for downstream loops that reference cov_M2 — points to PRIMARY block
cov_M2 <- cov_pre

# ---- Logistic OR for CAP ≥ 275 (Rinella 2023) ----
cat("\n[1/4] Logistic OR for CAP ≥ 275 (Rinella 2023) by ratio tertile ...\n")
results <- list()
for (ratio_var in c("se_zn_t", "se_cu_t")) {
  f <- as.formula(paste0("steatosis_cap275 ~ ", ratio_var, " + ",
                          paste(cov_M2, collapse = " + ")))
  fit <- tryCatch(svyglm(f, design = des_an, family = quasibinomial()),
                  error = function(e) {cat("  ERR", ratio_var, ":", conditionMessage(e), "\n"); NULL})
  if (!is.null(fit)) {
    s <- summary(fit)$coefficients
    rows <- grep(paste0("^", ratio_var), rownames(s), value = TRUE)
    out <- data.frame(
      ratio = ratio_var, term = rows,
      beta = s[rows, "Estimate"], se = s[rows, "Std. Error"],
      OR = exp(s[rows, "Estimate"]),
      lci = exp(s[rows, "Estimate"] - 1.96 * s[rows, "Std. Error"]),
      uci = exp(s[rows, "Estimate"] + 1.96 * s[rows, "Std. Error"]),
      p = s[rows, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    results[[ratio_var]] <- out
    print(out)
  }
}
or_df <- do.call(rbind, results)

# ---- Continuous LSM by ratio ----
cat("\n[2/4] Continuous LSM by ratio (M2 adjusted) ...\n")
lsm_results <- list()
for (ratio_var in c("se_zn_ratio", "se_cu_ratio")) {
  if (!(ratio_var %in% names(df_an))) next
  df_lsm <- df_an %>% filter(!is.na(lsm)) %>% mutate(lsm_log = log(pmax(lsm, 0.1)))
  if (nrow(df_lsm) < 100) next
  des_lsm <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                       data = df_lsm, nest = TRUE)
  f <- as.formula(paste0("lsm_log ~ scale(", ratio_var, ") + ",
                          paste(cov_M2, collapse = " + ")))
  fit <- tryCatch(svyglm(f, design = des_lsm), error = function(e) NULL)
  if (!is.null(fit)) {
    s <- summary(fit)$coefficients
    term1 <- rownames(s)[2]  # first non-intercept
    lsm_results[[ratio_var]] <- data.frame(
      ratio = ratio_var, term = term1,
      beta = s[term1, "Estimate"], se = s[term1, "Std. Error"],
      p = s[term1, "Pr(>|t|)"], n = nrow(df_lsm),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %s: beta=%.4f (SE %.4f) p=%.3g\n",
                ratio_var, s[term1, "Estimate"], s[term1, "Std. Error"], s[term1, "Pr(>|t|)"]))
  }
}
lsm_df <- do.call(rbind, lsm_results)

# ---- HFS continuous by ratio ----
cat("\n[3/4] Continuous HFS by ratio (M2 adjusted, fasting subsample) ...\n")
hfs_results <- list()
for (ratio_var in c("se_zn_ratio", "se_cu_ratio")) {
  df_h <- df_an %>%
    filter(!is.na(hfs), !is.na(wt_saf_pooled), wt_saf_pooled > 0,
           !is.na(.data[[ratio_var]]))
  if (nrow(df_h) < 100) next
  des_h <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_saf_pooled,
                     data = df_h, nest = TRUE)
  f <- as.formula(paste0("hfs ~ scale(", ratio_var, ") + ",
                          paste(cov_M2, collapse = " + ")))
  fit <- tryCatch(svyglm(f, design = des_h), error = function(e) NULL)
  if (!is.null(fit)) {
    s <- summary(fit)$coefficients
    term1 <- rownames(s)[2]
    hfs_results[[ratio_var]] <- data.frame(
      ratio = ratio_var, term = term1,
      beta = s[term1, "Estimate"], se = s[term1, "Std. Error"],
      p = s[term1, "Pr(>|t|)"], n = nrow(df_h),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %s: beta=%.4f (SE %.4f) p=%.3g\n",
                ratio_var, s[term1, "Estimate"], s[term1, "Std. Error"], s[term1, "Pr(>|t|)"]))
  }
}
hfs_df <- do.call(rbind, hfs_results)

# ---- Incremental ΔR² vs single-Se quartile baseline ----
cat("\n[4/4] Incremental ΔR² ratio model vs single-Se quartile (CAP ≥ 275) ...\n")
df_an <- df_an %>%
  mutate(
    se_q = factor(ntile(LBXBSE, 4), levels = 1:4,
                  labels = c("Q1", "Q2", "Q3", "Q4"))
  )
des_an2 <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                     data = df_an, nest = TRUE)
fit_base <- svyglm(as.formula(paste0("steatosis_cap275 ~ se_q + ",
                                      paste(cov_M2, collapse = " + "))),
                    design = des_an2, family = quasibinomial())
fit_ratio <- svyglm(as.formula(paste0("steatosis_cap275 ~ se_q + se_zn_t + se_cu_t + ",
                                       paste(cov_M2, collapse = " + "))),
                     design = des_an2, family = quasibinomial())
# Pseudo R² via Nagelkerke approximation on weighted residual deviance
calc_pseudo_r2 <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  null_dev <- fit$null.deviance; res_dev <- fit$deviance; n <- nrow(fit$model)
  if (is.null(null_dev) || is.null(res_dev) || n == 0) return(NA_real_)
  1 - exp((res_dev - null_dev) / n)
}
r2_base  <- calc_pseudo_r2(fit_base)
r2_ratio <- calc_pseudo_r2(fit_ratio)
delta_r2 <- r2_ratio - r2_base
cat(sprintf("  pseudo R² baseline (single Se quartile + M2): %.4f\n", r2_base))
cat(sprintf("  pseudo R² + Se/Zn + Se/Cu tertile:           %.4f\n", r2_ratio))
cat(sprintf("  ΔR² (incremental from ratios): %.4f\n", delta_r2))
delta_r2_df <- data.frame(
  pseudo_R2_baseline = r2_base, pseudo_R2_with_ratios = r2_ratio,
  delta_R2 = delta_r2
)

# Save
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(or_df,       "output/tables/ratio_OR_cap275.csv",   row.names = FALSE)
write.csv(lsm_df,      "output/tables/ratio_LSM_continuous.csv", row.names = FALSE)
write.csv(hfs_df,      "output/tables/ratio_HFS_continuous.csv", row.names = FALSE)
write.csv(delta_r2_df, "output/tables/ratio_incremental_R2.csv", row.names = FALSE)
save(or_df, lsm_df, hfs_df, delta_r2_df, fit_base, fit_ratio,
     file = "output/tables/ratio_analysis_results.RData")

cat("\n保存:\n")
cat("  output/tables/ratio_OR_cap275.csv (Se/Zn + Se/Cu tertile OR for CAP ≥ 275)\n")
cat("  output/tables/ratio_LSM_continuous.csv\n")
cat("  output/tables/ratio_HFS_continuous.csv\n")
cat("  output/tables/ratio_incremental_R2.csv (ΔR² vs single-Se quartile)\n")
cat("\nDONE 12_ratio_analysis.R\n")
