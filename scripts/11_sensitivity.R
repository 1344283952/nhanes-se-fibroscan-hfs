# ============================================
# 006 / 11_sensitivity.R — S1-S7 sensitivity analyses
#
# Per task.md §7.7 + OSF §6 (Round 5 R-Bias NK-8 added S7):
#   S1: retain baseline HBV/HCV (no exclusion)
#   S2: single-exposure only (only DR1TSELE OR only LBXBSE)
#   S3: MetALD subset standalone GAM
#   S4: MICE m=20 — deferred to manuscript
#   S5: single-cycle-J (2017-2018) single-Se within-cohort robustness
#   S6: LUX_IQR/LUXSMED ≤30% strict-validity (Boursier 2013)
#   S7 (NK-8): trimmed-exposure GAM (LBXBSE + DR1TSELE ∈ [P1, P99]) for non-rare CAP 44.3%
# ============================================

set.seed(20260516)
suppressPackageStartupMessages({
  library(dplyr); library(survey); library(mgcv)
})
# R-Stats v2 + R-NHANES v2 Round 2 P0: lonely PSU handling.
options(survey.lonely.psu = "adjust")
cat("========================================\n")
cat("006 W7 — Sensitivity S1-S7 (Round 2 — cov_pre primary; S4 deferred)\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")
load("data/processed/nhanes_final.RData")

# R-Causal-Methods v2 + R-Stats v2 Round 2 P0:
# Primary covariate set = pre-exposure C only (Pearl backdoor minimal adjustment set per
# _dag_spec.md §4). Sensitivity analyses are reported against the PRIMARY covariate
# set so estimates are comparable to the main analysis. The previous cov_M2 with
# BMI/HbA1c/LDL/diabetes conditions on mediator-confounders and blocks part of the
# Se → outcome total effect (Westreich-Greenland 2013 Table-2 fallacy).
cov_pre <- c("age", "RIAGENDR", "race", "INDFMPIR", "education",
             "smoke", "drink", "DR1TKCAL")
cov_pre_use <- intersect(cov_pre, names(nhanes_final))
# Alias preserved for code paths that reference cov_M2_use — semantically = cov_pre
cov_M2_use <- cov_pre_use
cat(sprintf("cov_pre PRIMARY (pre-exposure C only): %s\n",
            paste(cov_pre_use, collapse = ", ")))

# Helper: logistic OR for CAP ≥ 275 with M2 covariates, weighted
logistic_cap <- function(df, exp_var) {
  if (nrow(df) < 50 || !(exp_var %in% names(df))) return(NULL)
  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                   data = df, nest = TRUE)
  f <- as.formula(paste0("steatosis_cap275 ~ scale(", exp_var, ") + ",
                          paste(cov_M2_use, collapse = " + ")))
  fit <- tryCatch(svyglm(f, design = des, family = quasibinomial()),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coefficients
  data.frame(exposure = exp_var,
             beta = s[2, "Estimate"], se = s[2, "Std. Error"],
             OR = exp(s[2, "Estimate"]),
             lci = exp(s[2, "Estimate"] - 1.96 * s[2, "Std. Error"]),
             uci = exp(s[2, "Estimate"] + 1.96 * s[2, "Std. Error"]),
             p = s[2, "Pr(>|t|)"], n = nrow(df), stringsAsFactors = FALSE)
}

# ---- S1: Retain HBV/HCV ----
# Main cohort filtered them already; can't recover. Doc-only note.
cat("\n[S1] Retain HBV/HCV: DOC-ONLY (main cohort already filtered; reverse requires re-merge)\n")

# ---- S2a: Single-exposure DR1TSELE only ----
cat("\n[S2a] Single-exposure: DR1TSELE only (vs main dual)\n")
s2a <- logistic_cap(nhanes_final, "DR1TSELE")
if (!is.null(s2a)) print(s2a)

# ---- S2b: Single-exposure LBXBSE only ----
cat("\n[S2b] Single-exposure: LBXBSE only (vs main dual)\n")
s2b <- logistic_cap(nhanes_final, "LBXBSE")
if (!is.null(s2b)) print(s2b)

# ---- S3: MetALD subset only ----
cat("\n[S3] MetALD subset standalone (Se → CAP within MetALD class)\n")
df_metald <- nhanes_final %>% filter(metald_group == "MetALD")
cat(sprintf("  MetALD subset N = %d\n", nrow(df_metald)))
s3_diet <- logistic_cap(df_metald, "DR1TSELE")
s3_blood <- logistic_cap(df_metald, "LBXBSE")
if (!is.null(s3_diet)) print(s3_diet)
if (!is.null(s3_blood)) print(s3_blood)

# ---- S4: MICE m=20 (deferred) ----
cat("\n[S4] MICE m=20: DEFERRED to manuscript draft (Rtools required)\n")

# ---- S5: single-cycle-J 2017-2018 single-Se within-cohort robustness ----
# Wang 2021 ENV RES used NHANES 2011-2016 (predating FibroScan); we restrict
# to J 2017-2018 only as the closest single-cycle robustness comparison.
cat("\n[S5] Single-cycle J 2017-2018 single-Se robustness (LBXBSE → CAP)\n")
df_j <- nhanes_final %>% filter(!is_prepandemic)  # only J cycle (2017-2018), not P_
if (nrow(df_j) >= 50) {
  s5 <- logistic_cap(df_j, "LBXBSE")
  cat(sprintf("  Cycle-J only N = %d\n", nrow(df_j)))
  if (!is.null(s5)) print(s5)
} else {
  cat("  WARN: no cycle-J only rows available (006 cohort already P_ only)\n")
  s5 <- NULL
}

# ---- S6: LUX IQR ≤ 30% strict-validity ----
cat("\n[S6] LUX IQR/LSM ≤ 30% strict validity (Boursier 2013)\n")
df_iqr <- nhanes_final %>% filter(lsm_iqr_valid)
cat(sprintf("  IQR ≤ 30 cohort N = %d\n", nrow(df_iqr)))
s6_diet <- logistic_cap(df_iqr, "DR1TSELE")
s6_blood <- logistic_cap(df_iqr, "LBXBSE")
if (!is.null(s6_diet)) print(s6_diet)
if (!is.null(s6_blood)) print(s6_blood)

# ---- S7 (Round 5 R-Bias NK-8): Trimmed-exposure ----
cat("\n[S7] Trimmed-exposure GAM (LBXBSE + DR1TSELE in [P1, P99]) for non-rare CAP 44.3%\n")
p1_b <- quantile(nhanes_final$LBXBSE, 0.01, na.rm = TRUE)
p99_b <- quantile(nhanes_final$LBXBSE, 0.99, na.rm = TRUE)
p1_d <- quantile(nhanes_final$DR1TSELE, 0.01, na.rm = TRUE)
p99_d <- quantile(nhanes_final$DR1TSELE, 0.99, na.rm = TRUE)
df_trim <- nhanes_final %>%
  filter(LBXBSE >= p1_b, LBXBSE <= p99_b,
          DR1TSELE >= p1_d, DR1TSELE <= p99_d)
cat(sprintf("  Trim N = %d (removed %d extreme rows)\n",
            nrow(df_trim), nrow(nhanes_final) - nrow(df_trim)))
s7_diet <- logistic_cap(df_trim, "DR1TSELE")
s7_blood <- logistic_cap(df_trim, "LBXBSE")
if (!is.null(s7_diet)) print(s7_diet)
if (!is.null(s7_blood)) print(s7_blood)

# ---- Save ----
sens_all <- bind_rows(
  if (!is.null(s2a))   s2a %>% mutate(sensitivity = "S2a_diet_only"),
  if (!is.null(s2b))   s2b %>% mutate(sensitivity = "S2b_blood_only"),
  if (!is.null(s3_diet)) s3_diet %>% mutate(sensitivity = "S3_MetALD_diet"),
  if (!is.null(s3_blood)) s3_blood %>% mutate(sensitivity = "S3_MetALD_blood"),
  if (!is.null(s5))    s5 %>% mutate(sensitivity = "S5_single_cycle_J_robustness"),
  if (!is.null(s6_diet)) s6_diet %>% mutate(sensitivity = "S6_IQR30_diet"),
  if (!is.null(s6_blood)) s6_blood %>% mutate(sensitivity = "S6_IQR30_blood"),
  if (!is.null(s7_diet)) s7_diet %>% mutate(sensitivity = "S7_trim_diet"),
  if (!is.null(s7_blood)) s7_blood %>% mutate(sensitivity = "S7_trim_blood")
)
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(sens_all, "output/tables/sensitivity_S1-S7.csv", row.names = FALSE)
save(sens_all, file = "output/tables/sensitivity_results.RData")

cat("\n保存:\n")
cat("  output/tables/sensitivity_S1-S7.csv\n")
cat("\nDONE 11_sensitivity.R\n")
