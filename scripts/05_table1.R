# ============================================
# 006 / 05_table1.R — Baseline characteristics (Round 2 rewrite)
#
# Round 2 P0 fixes (2026-05-16):
#   - 不再用 `dplyr::coalesce(WTMEC2YR, WTMECPRP) / 2` (违反 NCHS Series 2 No. 190)
#     直接读 03_clean 已算好的 `wt_pooled` 列 (per templates/_shared/pooled_weight.R)
#   - Sheet 3 MetALD 用 CAP ≥ 275 primary (`metald_group`, AASLD 2024)
#   - Sheet 4 改用 LSM ≥ 8 binary (`fibrosis_lsm8`)
#
# 4 sheets:
#   Sheet 1 "by_Se_blood_Q4"   — `se_blood_q_lbl`
#   Sheet 2 "by_Se_diet_Q4"    — `se_diet_q_lbl`
#   Sheet 3 "by_MetALD_primary"— `metald_group` (CAP ≥ 275 primary)
#   Sheet 4 "by_LSM_F2plus"    — `fibrosis_lsm8` (LSM ≥ 8 binary)
# ============================================

library(survey); library(tableone); library(openxlsx); library(dplyr)

set.seed(20260516)  # Round 1 R-Repro P0: 全 pipeline 一致 seed
options(survey.lonely.psu = "adjust")  # W11 R2 R-Stats P1: Lumley 2010 §3.4 single-PSU stratum behaviour

cat("========================================\n")
cat("006 / 05_table1 — Round 2 rewrite (Se Q4 × 2 + MetALD primary + LSM F2+)\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")
load("data/processed/nhanes_final.RData")

# --- 卫生检查 ---
required_cols <- c("wt_pooled", "SDMVPSU", "SDMVSTRA",
                   "se_blood_q", "se_diet_q",
                   "metald_group", "fibrosis_lsm8")
missing_cols <- setdiff(required_cols, names(nhanes_final))
if (length(missing_cols) > 0) {
  stop(sprintf("nhanes_final 缺少必备列 (Round 1/2 fix 后): %s",
               paste(missing_cols, collapse = ", ")))
}

# --- 构建 strata 标签 ---
nhanes_final <- nhanes_final %>%
  mutate(
    se_blood_q_lbl = factor(
      paste0("Q", se_blood_q),
      levels = c("Q1", "Q2", "Q3", "Q4")
    ),
    se_diet_q_lbl = factor(
      paste0("Q", se_diet_q),
      levels = c("Q1", "Q2", "Q3", "Q4")
    ),
    lsm_f2plus_lbl = factor(
      ifelse(fibrosis_lsm8 == 1, "LSM>=8 (F2+)", "LSM<8"),
      levels = c("LSM<8", "LSM>=8 (F2+)")
    )
  )

# --- 主 design (直接读 wt_pooled) ---
design_t1 <- svydesign(
  ids     = ~SDMVPSU,
  strata  = ~SDMVSTRA,
  weights = ~wt_pooled,
  data    = nhanes_final,
  nest    = TRUE
)

# --- Table 1 变量集 ---
t1_vars <- c("age", "RIAGENDR", "race", "education",
             "pir", "bmi", "smoke", "drink",
             "LBXSATSI", "LBXSASSI", "LBDSALSI", "LBXPLTSI",
             "fib4", "fib4_advanced", "fib4_advanced_aged",
             "hfs", "hfs_high",
             "LBXGH", "LBXGLU", "homa_ir", "diabetes", "hypertension",
             "LBXBSE", "DR1TSELE", "DR1TZINC", "DR1TCOPP",
             "se_zn_ratio_diet", "se_cu_ratio_diet",
             "lsm", "cap",
             "fibrosis_lsm8", "fibrosis_lsm12",
             "steatosis_cap275", "steatosis_cap248",
             "alcohol_gwk", "selfreport_liver")
factor_vars <- c("RIAGENDR", "race", "education", "smoke", "drink",
                 "fib4_advanced", "fib4_advanced_aged", "hfs_high",
                 "diabetes", "hypertension",
                 "fibrosis_lsm8", "fibrosis_lsm12",
                 "steatosis_cap275", "steatosis_cap248",
                 "selfreport_liver")
t1_vars     <- intersect(t1_vars, names(nhanes_final))
factor_vars <- intersect(factor_vars, t1_vars)

# ============================================
# Sheet 1: by Se blood Q4 (LBXBSE)
# ============================================
cat("--- Sheet 1: by_Se_blood_Q4 (LBXBSE quartile) ---\n")
t1a <- svyCreateTableOne(vars = t1_vars, factorVars = factor_vars,
                         strata = "se_blood_q_lbl", data = design_t1,
                         test = TRUE, smd = TRUE)
t1a_mat <- print(t1a, showAllLevels = TRUE, smd = TRUE, printToggle = FALSE,
                 contDigits = 2, catDigits = 1, missing = TRUE)

# ============================================
# Sheet 2: by Se dietary Q4 (DR1TSELE)
# ============================================
cat("\n--- Sheet 2: by_Se_diet_Q4 (DR1TSELE quartile) ---\n")
t1b <- svyCreateTableOne(vars = t1_vars, factorVars = factor_vars,
                         strata = "se_diet_q_lbl", data = design_t1,
                         test = TRUE, smd = TRUE)
t1b_mat <- print(t1b, showAllLevels = TRUE, smd = TRUE, printToggle = FALSE,
                 contDigits = 2, catDigits = 1, missing = TRUE)

# ============================================
# Sheet 3: by MetALD primary (CAP ≥ 275, AASLD 2024)
# ============================================
cat("\n--- Sheet 3: by_MetALD_primary (CAP ≥ 275, Rinella 2023 4-group) ---\n")
nhanes_metald <- nhanes_final %>% filter(!is.na(metald_group))
design_metald <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  data = nhanes_metald, nest = TRUE
)
t1c <- svyCreateTableOne(vars = t1_vars, factorVars = factor_vars,
                         strata = "metald_group", data = design_metald,
                         test = TRUE, smd = TRUE)
t1c_mat <- print(t1c, showAllLevels = TRUE, smd = TRUE, printToggle = FALSE,
                 contDigits = 2, catDigits = 1, missing = TRUE)

# ============================================
# Sheet 4: by LSM ≥ 8 binary (significant fibrosis F2+)
# ============================================
cat("\n--- Sheet 4: by_LSM_F2plus (LSM ≥ 8 binary) ---\n")
nhanes_fib <- nhanes_final %>% filter(!is.na(fibrosis_lsm8))
design_fib <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  data = nhanes_fib, nest = TRUE
)
t1d <- svyCreateTableOne(vars = t1_vars, factorVars = factor_vars,
                         strata = "lsm_f2plus_lbl", data = design_fib,
                         test = TRUE, smd = TRUE)
t1d_mat <- print(t1d, showAllLevels = TRUE, smd = TRUE, printToggle = FALSE,
                 contDigits = 2, catDigits = 1, missing = TRUE)

# ============================================
# 写 xlsx (4 sheets)
# ============================================
wb <- createWorkbook()
addWorksheet(wb, "by_Se_blood_Q4")
writeData(wb, "by_Se_blood_Q4", t1a_mat, rowNames = TRUE)
addWorksheet(wb, "by_Se_diet_Q4")
writeData(wb, "by_Se_diet_Q4", t1b_mat, rowNames = TRUE)
addWorksheet(wb, "by_MetALD_primary")
writeData(wb, "by_MetALD_primary", t1c_mat, rowNames = TRUE)
addWorksheet(wb, "by_LSM_F2plus")
writeData(wb, "by_LSM_F2plus", t1d_mat, rowNames = TRUE)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
saveWorkbook(wb, "output/tables/table1.xlsx", overwrite = TRUE)
cat("\n已保存 output/tables/table1.xlsx (4 sheets: Se_blood_Q4 / Se_diet_Q4 / MetALD_primary / LSM_F2plus)\n")

# ============================================
# Round 2 R-NHANES-Domain v2 + R-Stats v2 P0:
# Design effect (DEFF) + degrees of freedom report — required to demonstrate
# the survey weights are actually doing work (PLOS-Bio NHANES red flag).
# ============================================
cat("\n--- DEFF + df report (Round 2 R-NHANES-Domain + R-Stats P0) ---\n")
# Round 2 bug fix: naive `nPSU - nStrata` is wrong when subset reduces PSU count
# below stratum count (yields negative df). Use survey::degf() which correctly
# computes the design-based effective df via PSU clustering within strata.
df_survey <- as.integer(survey::degf(design_t1))
# Also report raw counts for transparency
n_psu_raw   <- length(unique(nhanes_final$SDMVPSU[!is.na(nhanes_final$SDMVPSU)]))
n_strat_raw <- length(unique(nhanes_final$SDMVSTRA[!is.na(nhanes_final$SDMVSTRA)]))
cat(sprintf("Survey design df (via survey::degf) = %d  [raw n_PSU = %d, n_strata = %d]\n",
            df_survey, n_psu_raw, n_strat_raw))

# DEFF for key primary outcomes + exposures
deff_vars <- c("steatosis_cap275", "fibrosis_lsm8", "fibrosis_lsm12",
               "LBXBSE", "DR1TSELE", "hfs", "fib4")
deff_tbl <- data.frame(variable = character(0), DEFF = numeric(0),
                       N_unweighted = integer(0), stringsAsFactors = FALSE)
for (v in deff_vars) {
  if (!(v %in% names(nhanes_final))) next
  sv <- tryCatch(svymean(as.formula(paste0("~", v)),
                          design = design_t1, na.rm = TRUE, deff = TRUE),
                  error = function(e) NULL)
  if (!is.null(sv)) {
    deff_val <- as.numeric(attr(sv, "deff")[1])
    deff_tbl <- rbind(deff_tbl, data.frame(
      variable = v, DEFF = deff_val,
      N_unweighted = sum(!is.na(nhanes_final[[v]])),
      stringsAsFactors = FALSE))
  }
}
deff_tbl$design_df <- df_survey
print(deff_tbl)
write.csv(deff_tbl, "output/tables/table1_deff.csv", row.names = FALSE)
cat("\n已保存 output/tables/table1_deff.csv (DEFF for 7 key primaries + design df)\n")
cat("========================================\n")
