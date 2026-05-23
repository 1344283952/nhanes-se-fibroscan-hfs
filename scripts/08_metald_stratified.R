# ============================================
# 006_se_fibroscan_hfs / 08_metald_stratified.R
# MetALD-stratified Se → HFS analyses (CAP ≥ 275 primary, CAP ≥ 248 sensitivity)
#
# Round 1 fixes incorporated:
#   - P0 [R-Clinical] Use AASLD 2024 CAP ≥ 275 as primary (encoded in 03_clean
#     as metald_group); CAP ≥ 248 (Karlas 2017, metald_group_cap248) as sensitivity
#   - P0 [R-Clinical] Rinella 2023 three-group strata: MASLD / MetALD / ALD / No-steatosis
#   - P1 [R-Stats] p-interaction Wald test (Se × metald_group product term)
#   - P1 [R-Bias] Pre-exposure covariates only
#   - P1 [R-Repro] set.seed(20260516); openxlsx workbook with per-stratum sheets
# ============================================

set.seed(20260516)
options(survey.lonely.psu = "adjust")  # W11 R2 R-Stats P1: Lumley 2010 §3.4 single-PSU stratum behaviour

library(dplyr)
library(survey)
library(openxlsx)

cat("========================================\n")
cat("006 W4 — MetALD-stratified Se → HFS (Round 1: CAP≥275 primary)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
load("data/processed/nhanes_design.RData")
# nhanes_final + design_main / design_lsm_strict / design_hfs / design_metald /
# design_metald_cap248 / design_steatosis

cat(sprintf("nhanes_final n=%d ; HFS-eligible n=%d\n",
            nrow(nhanes_final), sum(!is.na(nhanes_final$hfs))))

nhanes_final$sex_male_i <- as.integer(nhanes_final$RIAGENDR == 1)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)

# ---- Function: stratum-specific Se → HFS (continuous) via svyglm ----
# Returns coefficient table for Se exposures within stratum
fit_stratum <- function(df_subset, label) {
  if (nrow(df_subset) < 30) {
    cat(sprintf("  SKIP %s: n<30 (n=%d)\n", label, nrow(df_subset)))
    return(NULL)
  }
  des <- tryCatch(
    svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
              data = df_subset, nest = TRUE),
    error = function(e) {
      cat(sprintf("  svydesign error %s: %s\n", label, conditionMessage(e)))
      NULL
    })
  if (is.null(des)) return(NULL)

  cov_str <- "age + sex_male_i + race + education + pir + smoke + drink"

  # Dietary Se → HFS
  f_diet <- as.formula(paste("hfs ~ DR1TSELE +", cov_str))
  # Blood Se → HFS
  f_blood <- as.formula(paste("hfs ~ LBXBSE +", cov_str))

  m_diet <- tryCatch(svyglm(f_diet, design = des),
                      error = function(e) {
                        cat(sprintf("  diet model %s err: %s\n", label, conditionMessage(e)))
                        NULL
                      })
  m_blood <- tryCatch(svyglm(f_blood, design = des),
                       error = function(e) {
                         cat(sprintf("  blood model %s err: %s\n", label, conditionMessage(e)))
                         NULL
                       })

  extract_row <- function(m, exp_var, exp_label) {
    if (is.null(m)) {
      return(data.frame(stratum = label, exposure = exp_label,
                        n = nrow(df_subset),
                        beta = NA, se = NA, lcl = NA, ucl = NA, p = NA,
                        stringsAsFactors = FALSE))
    }
    coef_tab <- summary(m)$coefficients
    if (!exp_var %in% rownames(coef_tab)) {
      return(data.frame(stratum = label, exposure = exp_label,
                        n = nrow(df_subset),
                        beta = NA, se = NA, lcl = NA, ucl = NA, p = NA,
                        stringsAsFactors = FALSE))
    }
    b  <- coef_tab[exp_var, "Estimate"]
    se <- coef_tab[exp_var, "Std. Error"]
    p  <- coef_tab[exp_var, ncol(coef_tab)]
    data.frame(
      stratum = label, exposure = exp_label,
      n = nrow(df_subset),
      beta = b, se = se,
      lcl = b - 1.96 * se, ucl = b + 1.96 * se,
      p = p,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    extract_row(m_diet,  "DR1TSELE", "Dietary Se (DR1TSELE, per 1 µg/day)"),
    extract_row(m_blood, "LBXBSE",   "Blood Se (LBXBSE, per 1 µg/L)")
  )
}

# ---- Function: p-interaction via Wald test on product term ----
test_interaction <- function(df, exp_var, group_var) {
  df_sub <- df[!is.na(df[[group_var]]) & !is.na(df$hfs) & !is.na(df[[exp_var]]), ]
  if (nrow(df_sub) < 100) {
    return(data.frame(exposure = exp_var, group = group_var,
                      chisq = NA, df = NA, p_interaction = NA,
                      stringsAsFactors = FALSE))
  }
  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                   data = df_sub, nest = TRUE)
  f_main <- as.formula(sprintf(
    "hfs ~ %s + %s + age + sex_male_i + race + education + pir + smoke + drink",
    exp_var, group_var))
  f_int <- as.formula(sprintf(
    "hfs ~ %s * %s + age + sex_male_i + race + education + pir + smoke + drink",
    exp_var, group_var))
  m_main <- tryCatch(svyglm(f_main, design = des), error = function(e) NULL)
  m_int  <- tryCatch(svyglm(f_int,  design = des), error = function(e) NULL)
  if (is.null(m_main) || is.null(m_int)) {
    return(data.frame(exposure = exp_var, group = group_var,
                      chisq = NA, df = NA, p_interaction = NA,
                      stringsAsFactors = FALSE))
  }
  # Wald test on interaction terms only
  int_terms <- setdiff(names(coef(m_int)), names(coef(m_main)))
  if (length(int_terms) == 0) {
    return(data.frame(exposure = exp_var, group = group_var,
                      chisq = NA, df = NA, p_interaction = NA,
                      stringsAsFactors = FALSE))
  }
  wt <- tryCatch(
    survey::regTermTest(m_int, as.formula(paste("~", paste(int_terms, collapse = "+")))),
    error = function(e) NULL
  )
  if (is.null(wt)) {
    return(data.frame(exposure = exp_var, group = group_var,
                      chisq = NA, df = NA, p_interaction = NA,
                      stringsAsFactors = FALSE))
  }
  data.frame(
    exposure = exp_var, group = group_var,
    chisq = as.numeric(wt$Ftest),
    df = wt$df, p_interaction = as.numeric(wt$p),
    stringsAsFactors = FALSE
  )
}

# ---- Build strata for CAP ≥ 275 primary ----
cat("\n[1/2] Primary: CAP ≥ 275 MetALD stratification ...\n")
strata_levels <- c("MASLD", "MetALD", "ALD", "No steatosis", "Cryptogenic steatosis")
primary_results <- list()
for (lvl in strata_levels) {
  df_lvl <- nhanes_final %>%
    filter(!is.na(metald_group), metald_group == lvl, !is.na(hfs))
  cat(sprintf("  %-25s n=%d\n", lvl, nrow(df_lvl)))
  res <- fit_stratum(df_lvl, label = lvl)
  if (!is.null(res)) primary_results[[lvl]] <- res
}
primary_df <- do.call(rbind, primary_results)
rownames(primary_df) <- NULL

# Overall (un-stratified, full HFS-eligible)
df_all_hfs <- nhanes_final %>% filter(!is.na(hfs))
overall_res <- fit_stratum(df_all_hfs, label = "Overall (all strata)")
primary_df <- rbind(overall_res, primary_df)

# p-interaction tests for primary
p_int_primary <- rbind(
  test_interaction(nhanes_final, "DR1TSELE", "metald_group"),
  test_interaction(nhanes_final, "LBXBSE",   "metald_group")
)
cat("\np-interaction (primary CAP ≥ 275):\n"); print(p_int_primary)

# ---- Sensitivity: CAP ≥ 248 ----
cat("\n[2/2] Sensitivity: CAP ≥ 248 MetALD stratification ...\n")
sens_results <- list()
for (lvl in strata_levels) {
  df_lvl <- nhanes_final %>%
    filter(!is.na(metald_group_cap248), metald_group_cap248 == lvl, !is.na(hfs))
  cat(sprintf("  %-25s n=%d\n", lvl, nrow(df_lvl)))
  res <- fit_stratum(df_lvl, label = lvl)
  if (!is.null(res)) sens_results[[lvl]] <- res
}
sens_df <- do.call(rbind, sens_results)
rownames(sens_df) <- NULL
overall_res_s <- fit_stratum(df_all_hfs, label = "Overall (all strata)")
sens_df <- rbind(overall_res_s, sens_df)

p_int_sens <- rbind(
  test_interaction(nhanes_final, "DR1TSELE", "metald_group_cap248"),
  test_interaction(nhanes_final, "LBXBSE",   "metald_group_cap248")
)

# ---- Write to xlsx workbook ----
wb <- createWorkbook()
addWorksheet(wb, "Primary_CAP275")
addWorksheet(wb, "P_interaction_primary")
addWorksheet(wb, "Sensitivity_CAP248")
addWorksheet(wb, "P_interaction_sensitivity")
addWorksheet(wb, "README")

writeData(wb, "Primary_CAP275", primary_df)
writeData(wb, "P_interaction_primary", p_int_primary)
writeData(wb, "Sensitivity_CAP248", sens_df)
writeData(wb, "P_interaction_sensitivity", p_int_sens)
writeData(wb, "README", data.frame(
  description = c(
    "MetALD-stratified Se → HFS analyses",
    "Primary: AASLD 2024 CAP ≥ 275 dB/m (metald_group)",
    "Sensitivity: Karlas 2017 CAP ≥ 248 dB/m (metald_group_cap248)",
    "Strata per Rinella 2023: MASLD / MetALD / ALD / Cryptogenic / No steatosis",
    "p-interaction via survey::regTermTest Wald on Se × stratum product terms",
    "Pre-exposure covariates only: age + sex + race + edu + pir + smoke + drink",
    "Round 1 fix: cm_risk_count blood-pressure bug repaired in 03_clean → MetALD N halved (576 vs prior 802)",
    sprintf("Generated: %s ; seed=20260516", Sys.time())
  ),
  stringsAsFactors = FALSE
))

saveWorkbook(wb, "output/tables/metald_stratified.xlsx", overwrite = TRUE)

# Also save RData for downstream
save(primary_df, sens_df, p_int_primary, p_int_sens,
     file = "output/tables/metald_stratified.RData")

cat("\nPrimary CAP≥275 stratified results:\n"); print(primary_df)

cat("\n保存:\n")
cat("  output/tables/metald_stratified.xlsx (5 sheets)\n")
cat("  output/tables/metald_stratified.RData\n")

cat("\nDONE 08_metald_stratified.R\n")
