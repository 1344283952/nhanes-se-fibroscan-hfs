# ============================================
# 006 / 04_survey_design.R (Round 1 fixes)
#
# Round 1 P0 fixes:
#   - 主 design 用 nhanes_final (P_ only, J ⊂ P_ dedupe done in 03_clean)
#   - Weight wt_pooled 已在 03_clean 算好 (= WTMECPRP * 3.2/3.2 = WTMECPRP)
#   - design_lsm_strict 加 sensitivity (LUXSIQRM ≤ 30%, S2)
#   - design_cap248 加 sensitivity (Karlas 2017 cut-off, S3)
# ============================================

library(survey); library(dplyr)

# R-NHANES-Domain v2 Round 2 P0: lonely PSU handling for subset designs.
# Without this, single-PSU strata in metald / hfs / cap248 subsets either error
# or silently inflate SE (Lumley 2010 §3.4). "adjust" applies a center-at-grand-mean
# correction conservative for inference.
options(survey.lonely.psu = "adjust")

cat("========================================\n")
cat("006 / 04_survey_design (Round 2 — lonely PSU adjust + P_ only main)\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")

# 1. design_main (P_ only, all LSM-valid, 包含 IQR > 30% 的)
design_main <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  data = nhanes_final, nest = TRUE
)
cat(sprintf("design_main (P_ only): %d rows\n", nrow(nhanes_final)))

# 2. design_lsm_strict (S2 sensitivity: 仅 LSM IQR ≤ 30%, Castera 2019)
# Round 1 R-Bias P0-6: 这个 filter 是 outcome-dependent selection, 不能做主筛
df_strict <- nhanes_final %>% filter(lsm_iqr_valid)
design_lsm_strict <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  data = df_strict, nest = TRUE
)
cat(sprintf("design_lsm_strict (S2 sensitivity, IQR ≤ 30%%): %d rows\n",
            nrow(df_strict)))

# 3. design_hfs (fasting subsample with HFS)
# Round 4 R-NHANES P0: HFS analysis uses the SAF fasting-subsample weight
# (WTSAFPRP on P_), NOT the MEC main weight. Add defensive filter
# `!is.na(wt_saf_pooled)` so any HFS-eligible row with missing SAF weight is
# dropped before svydesign() — Lumley 2010 §3.4 forbids zero/NA weights in
# the design column (otherwise svyglm silently demotes them to a self-rep PSU
# with bogus DF).
df_hfs <- nhanes_final %>%
  filter(!is.na(hfs), !is.na(wt_saf_pooled), wt_saf_pooled > 0)
design_hfs <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_saf_pooled,
  data = df_hfs, nest = TRUE
)
cat(sprintf("design_hfs (fasting subsample, SAF weight non-NA & > 0): %d rows\n",
            nrow(df_hfs)))

# 4. design_metald (with MASLD/MetALD/ALD classified by CAP ≥ 275 primary)
df_metald <- nhanes_final %>% filter(!is.na(metald_group),
                                      metald_group != "No steatosis")
design_metald <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  data = df_metald, nest = TRUE
)
cat(sprintf("design_metald (steatosis-having only, CAP ≥ 275): %d rows\n",
            nrow(df_metald)))

# 5. design_metald_cap248 (S3 sensitivity, Karlas 2017 cut-off)
df_metald248 <- nhanes_final %>% filter(!is.na(metald_group_cap248),
                                         metald_group_cap248 != "No steatosis")
design_metald_cap248 <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  data = df_metald248, nest = TRUE
)
cat(sprintf("design_metald_cap248 (S3 sensitivity, CAP ≥ 248): %d rows\n",
            nrow(df_metald248)))

# 6. design_steatosis_any (all steatosis-having, for fibrosis progression analyses)
df_steatosis <- nhanes_final %>%
  filter(metald_group %in% c("MASLD", "MetALD", "ALD", "Cryptogenic steatosis"))
design_steatosis <- svydesign(
  ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
  data = df_steatosis, nest = TRUE
)
cat(sprintf("design_steatosis: %d rows\n", nrow(df_steatosis)))

save(design_main, design_lsm_strict, design_hfs, design_metald,
     design_metald_cap248, design_steatosis,
     file = "data/processed/nhanes_design.RData")
cat("\n已保存 data/processed/nhanes_design.RData (6 个 design)\n")
cat("========================================\n")
