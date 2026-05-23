# ============================================================
# 006v2 / 21_cmaverse_redox.R — CMAverse 4-way decomposition
# Se → GGT-mediator → HFS / steatosis
# VanderWeele 2014 4-way: CDE + INTREF + PIE + INTMED
# 2026-05-18 v2 升级
# ============================================================
library(CMAverse); library(dplyr)
set.seed(20260518)

cat("==== 21_cmaverse_redox: Se → GGT → HFS / Steatosis ====\n")
load("data/processed/nhanes_final.RData")

# 分析子集: 5,885 主 cohort
d <- nhanes_final %>%
  filter(!is.na(LBXBSE), !is.na(ggt_iul),
         !is.na(RIDAGEYR), !is.na(sex_male))

# CMAverse 要求 covariates 数值 / factor 干净
d$age <- d$RIDAGEYR
d$sex <- factor(d$sex_male, levels = c(0, 1), labels = c("F","M"))
d$race_f <- as.factor(d$race)
d$smk <- as.factor(d$SMQ020)
d$drk <- as.factor(d$ALQ111)
d$bmi <- d$BMXBMI
d$kcal <- d$DR1TKCAL

# 主分析变量
d$x_se      <- d$LBXBSE                    # exposure (continuous)
d$m_ggt     <- d$ggt_iul                   # mediator (continuous)
d$y_steat   <- d$steatosis_cap275          # primary binary outcome
d$y_hfs     <- d$hfs_high                  # secondary
d$y_lsm     <- d$fibrosis_lsm8

# Reference / target levels for exposure contrast
# 用 IQR 对比: Q1 (低硒) vs Q3 (高硒) 模拟"病理性高硒"假设
q1_se <- quantile(d$x_se, 0.25, na.rm = TRUE)
q3_se <- quantile(d$x_se, 0.75, na.rm = TRUE)
cat(sprintf("\nSe Q1 = %.1f µg/L; Q3 = %.1f µg/L\n", q1_se, q3_se))

run_cmaverse <- function(outcome_var, outcome_type = "binary") {
  required_cols <- c(outcome_var, "x_se", "m_ggt",
                     "age", "sex", "bmi", "kcal", "smk", "drk", "race_f")
  d_use <- d %>%
    select(all_of(required_cols)) %>%
    filter(complete.cases(.))
  cat(sprintf("\nOutcome %s, n=%d (complete-case across all required cols)\n",
              outcome_var, nrow(d_use)))
  fit <- tryCatch(
    CMAverse::cmest(
      data = d_use,
      model = "rb",                   # regression-based (gold standard)
      outcome = outcome_var,
      exposure = "x_se",
      mediator = "m_ggt",
      basec = c("age","sex","bmi","kcal","smk","drk","race_f"),
      EMint = TRUE,
      mreg = list("linear"),
      yreg = if (outcome_type == "binary") "logistic" else "linear",
      astar = q1_se, a = q3_se,        # Se Q3 vs Q1 contrast
      mval = list(median(d_use$m_ggt, na.rm = TRUE)),
      estimation = "imputation",
      inference = "bootstrap",
      nboot = 500
    ),
    error = function(e) { cat("  cmaverse failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)
  fit_summary <- summary(fit)
  list(outcome = outcome_var, fit = fit, summary_text = capture.output(print(fit_summary)))
}

results <- list()
results$steat_v2 <- run_cmaverse("y_steat", "binary")
results$hfs_v2   <- run_cmaverse("y_hfs",   "binary")
results$lsm_v2   <- run_cmaverse("y_lsm",   "binary")

save(results, file = "output/tables/cmaverse_redox.RData")

# 抽 effect 表
extract_effects <- function(fit, name) {
  s <- summary(fit$fit)
  if (is.null(s)) return(NULL)
  ef <- s$effect.pe
  ci <- s$effect.ci.low
  ch <- s$effect.ci.high
  pv <- s$effect.pval
  data.frame(outcome = name,
             effect = names(ef), pe = ef,
             ci_lo = ci, ci_hi = ch, p = pv,
             row.names = NULL)
}
df_list <- list()
for (nm in names(results)) {
  if (!is.null(results[[nm]])) {
    df_list[[nm]] <- extract_effects(results[[nm]], nm)
  }
}
if (length(df_list) > 0) {
  effect_df <- do.call(rbind, df_list)
  write.csv(effect_df, "output/tables/cmaverse_redox_effects.csv", row.names = FALSE)
  cat("\n=== CMAverse 4-way 效应分解 ===\n")
  print(effect_df, row.names = FALSE)
  cat("\n[OK] saved cmaverse_redox.RData / cmaverse_redox_effects.csv\n")
} else {
  cat("[ERROR] 0 successful cmaverse fits\n")
}
