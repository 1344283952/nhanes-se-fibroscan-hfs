# ============================================
# 006 / 10_subgroup_forest.R — 9-subgroup × 4 Se exposures × 3 outcomes
#
# Per task.md §7.5: 9 subgroups (sex/age/race/edu/DM/HTN/PIR/smoke/drink)
# Outcomes: CAP ≥ 275 (steatosis), LSM ≥ 8 (F2+ fibrosis), HFS continuous
# Exposures: DR1TSELE, LBXBSE, se_zn_ratio, se_cu_ratio
# FDR family: 9 strata × 3 outcomes × 4 exposures = 108 tests at BH q=0.05
# (separately from primary 4 exposures × 3 outcomes = 12 tests)
# ============================================

set.seed(20260516)
suppressPackageStartupMessages({
  library(dplyr); library(survey)
})
# R-Stats v2 + R-NHANES v2 Round 2 P0: lonely PSU handling for stratified subset designs.
options(survey.lonely.psu = "adjust")
cat("========================================\n")
cat("006 W7 — 9-subgroup × 4 Se × 3 outcomes forest (Round 2: cov_pre primary)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")

df <- nhanes_final %>%
  mutate(
    sex_strata = factor(ifelse(RIAGENDR == 1, "Male", "Female")),
    age_strata = factor(ifelse(age < 60, "<60", ">=60")),
    race_strata = factor(case_when(
      race == "Non-Hispanic White" ~ "NH White",
      race == "Non-Hispanic Black" ~ "NH Black",
      race == "Mexican American" ~ "Mexican-American",
      TRUE ~ "Other"
    )),
    edu_strata = factor(case_when(
      education %in% c("Less than HS", "HS/GED") ~ "<=HS",
      TRUE ~ ">HS"
    )),
    dm_strata = factor(ifelse(diabetes == 1, "DM", "no DM")),
    htn_strata = factor(ifelse(hypertension == 1, "HTN", "no HTN")),
    pir_strata = factor(case_when(
      INDFMPIR < 1.30 ~ "Low (<1.30)",
      INDFMPIR < 3.50 ~ "Mid (1.30-3.50)",
      TRUE ~ "High (>=3.50)"
    )),
    smoke_strata = factor(ifelse(smoke == "Ever", "ever-smoker", "never-smoker")),
    drink_strata = factor(ifelse(drink == "Yes", "drinker", "non-drinker")),
    se_zn_ratio = ifelse(!is.na(DR1TSELE) & !is.na(DR1TZINC) & DR1TZINC > 0,
                          DR1TSELE / DR1TZINC, NA_real_),
    se_cu_ratio = ifelse(!is.na(DR1TSELE) & !is.na(DR1TCOPP) & DR1TCOPP > 0,
                          DR1TSELE / DR1TCOPP, NA_real_)
  )

subgroups <- c("sex_strata", "age_strata", "race_strata", "edu_strata",
               "dm_strata", "htn_strata", "pir_strata", "smoke_strata", "drink_strata")
exposures <- c("DR1TSELE", "LBXBSE", "se_zn_ratio", "se_cu_ratio")
outcomes <- list(
  cap275 = list(var = "steatosis_cap275", family = "binomial"),
  lsm8   = list(var = "fibrosis_lsm8",    family = "binomial"),
  hfs    = list(var = "hfs",              family = "gaussian", subset_var = "hfs",
                weight_var = "wt_saf_pooled")
)
# R-Causal-Methods v2 + R-Stats v2 Round 2 P0:
# Primary covariate set = pre-exposure C only (Pearl backdoor minimal adjustment set
# per _dag_spec.md §4). BMI / HbA1c / LDL / diabetes are mediator-confounders
# (descendant of Se per DAG), conditioning on them blocks part of total Se → outcome
# path — Westreich-Greenland 2013 Table-2 fallacy. cov_M2 (with mediator block) kept
# as commented sensitivity variant for direct-effect (CDE) interpretation only.
cov_pre <- c("age", "RIAGENDR", "race", "INDFMPIR", "education",
             "smoke", "drink", "DR1TKCAL")
# cov_M2_sens (post-exposure mediator-confounder block) — sensitivity CDE only,
# NOT used as primary; uncomment + re-run only for CDE sensitivity reporting.
# cov_M2_sens <- c(cov_pre, "bmi", "LBXGH", "LBDLDL", "diabetes")
cov_M2 <- cov_pre  # alias preserved for downstream code paths; semantically = cov_pre primary

results <- list()
for (sg in subgroups) {
  cat(sprintf("\n[%s] ...\n", sg))
  levs <- levels(df[[sg]])
  for (lev in levs) {
    df_sub <- df %>% filter(.data[[sg]] == lev)
    if (nrow(df_sub) < 100) next
    for (out_name in names(outcomes)) {
      o <- outcomes[[out_name]]
      df_o <- df_sub
      if (!is.null(o$subset_var)) df_o <- df_o %>% filter(!is.na(.data[[o$subset_var]]))
      w_var <- if (!is.null(o$weight_var) && o$weight_var %in% names(df_o)) o$weight_var else "wt_pooled"
      df_o <- df_o %>% filter(!is.na(.data[[w_var]]), .data[[w_var]] > 0)
      if (nrow(df_o) < 50) next
      des_o <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
                         weights = as.formula(paste0("~", w_var)),
                         data = df_o, nest = TRUE)
      cov_keep <- setdiff(cov_M2, sg) ; cov_keep <- intersect(cov_keep, names(df_o))
      for (exp_var in exposures) {
        f <- as.formula(paste0(o$var, " ~ scale(", exp_var, ") + ",
                                paste(cov_keep, collapse = " + ")))
        fit <- tryCatch(if (o$family == "binomial")
                          svyglm(f, design = des_o, family = quasibinomial())
                        else svyglm(f, design = des_o),
                        error = function(e) NULL)
        if (!is.null(fit)) {
          s <- summary(fit)$coefficients
          if (nrow(s) >= 2) {
            term1 <- rownames(s)[2]
            p_col <- if ("Pr(>|t|)" %in% colnames(s)) "Pr(>|t|)" else colnames(s)[ncol(s)]
            p_val <- as.numeric(s[term1, p_col])
            if (is.na(p_val) || is.nan(p_val)) {
              z <- s[term1, "Estimate"] / s[term1, "Std. Error"]
              p_val <- 2 * pnorm(abs(z), lower.tail = FALSE)
            }
            results[[paste(sg, lev, exp_var, out_name, sep="|")]] <- data.frame(
              subgroup = sg, stratum = lev, exposure = exp_var, outcome = out_name,
              beta = s[term1, "Estimate"], se = s[term1, "Std. Error"],
              effect = if (o$family == "binomial") exp(s[term1, "Estimate"]) else s[term1, "Estimate"],
              lci = if (o$family == "binomial") exp(s[term1, "Estimate"] - 1.96 * s[term1, "Std. Error"]) else s[term1, "Estimate"] - 1.96 * s[term1, "Std. Error"],
              uci = if (o$family == "binomial") exp(s[term1, "Estimate"] + 1.96 * s[term1, "Std. Error"]) else s[term1, "Estimate"] + 1.96 * s[term1, "Std. Error"],
              p = p_val, n = nrow(df_o), stringsAsFactors = FALSE)
          }
        }
      }
    }
  }
}
res_df <- do.call(rbind, results)

if (!is.null(res_df) && nrow(res_df) > 0) {
  res_df <- res_df %>% mutate(p_BH_subgroup = p.adjust(p, method = "BH"))
}
cat(sprintf("\nTotal subgroup rows: %d\n", nrow(res_df)))

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(res_df, "output/tables/subgroup_forest_data.csv", row.names = FALSE)
save(res_df, file = "output/tables/subgroup_forest_results.RData")

cat("\n保存:\n")
cat("  output/tables/subgroup_forest_data.csv\n")
cat("\nDONE 10_subgroup_forest.R\n")
