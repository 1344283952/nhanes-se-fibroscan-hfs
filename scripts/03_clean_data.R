# ============================================
# 006_se_fibroscan_hfs / 03_clean_data.R (Round 1 fixes)
# 输入: data/processed/nhanes_raw_merged.RData
# 输出: data/processed/nhanes_final.RData
#
# Round 1 fixes (2026-05-16, 10 reviewer):
#   - P0-1 [NHANES] 主分析仅用 P_ rows (J ⊂ P_, 删 _J 避免双重计数)
#   - P0-2 [Clinical] AASLD 2024 CAP cut-off 改 275 primary (248 sensitivity)
#   - P0-3 [Clinical] LSM MASLD-specific cut-off 8/12 primary (8/9.5/12.5 sensitivity)
#   - P0-4 [Clinical] Hepamet HFS albumin continuous fix (templates/_shared/)
#   - P0-5 [Clinical+Repro] MetALD 血压 bug fix
#   - P0-6 [Bias] LUX IQR ≤30% 改 sensitivity (不做主筛, 避免 BMI-dep selection)
#   - P0-7 [Stats] FIB-4 加 age-adjusted (McPherson 2017)
#   - P0-8 [Repro] set.seed(20260516) + scale attrs saved
#   - P0-9 [NHANES] Pooled weight 用 helper
# ============================================

library(dplyr); library(tidyr); library(purrr)

source("../../templates/_shared/fib4.R")
source("../../templates/_shared/hepamet_hfs.R")
source("../../templates/_shared/metald_classify.R")
source("../../templates/_shared/pooled_weight.R")
# Round 3 Stage A: Wave A shared modules wire-up
source("../../templates/_shared/var_aliases.R")     # apply_var_aliases() + coalesce_cols()
source("../../templates/_shared/cycle_config.R")    # make_cycle_spec() + cycle_years_vector() + validate_no_J_in_PP()
source("../../templates/_shared/diabetes_define.R") # diabetes_comprehensive() + detect_antidiabetic_rx()

set.seed(20260516)

cat("========================================\n")
cat("006 P12: 清洗 + Se 双暴露 + 救援切角 A (Round 3 Stage A: Wave A modules wired)\n")
cat("========================================\n\n")

load("data/processed/nhanes_raw_merged.RData")
cat(sprintf("原始合并: %d 行 × %d 列\n\n", nrow(nhanes_all), ncol(nhanes_all)))

flow <- list()
log_flow <- function(label, n) {
  cat(sprintf("  [流程] %-58s n = %d\n", label, n))
  flow[[length(flow) + 1]] <<- data.frame(step = label, n = n, stringsAsFactors = FALSE)
}
log_flow("Raw merged table (J + P_)", nrow(nhanes_all))

# Util (na_codes / zero_to_na 局部保留; coalesce_cols 由 var_aliases.R 提供)
na_codes <- function(x, codes) ifelse(x %in% codes, NA, x)
zero_to_na <- function(x) ifelse(!is.na(x) & x == 0, NA, x)

# Step 0: 跨周期变量统一 → 用 templates/_shared/var_aliases.R 集中维护
# (R-NHANES P1-2 / R-TechDebt P0-B)
nhanes_all <- apply_var_aliases(nhanes_all)
cat(sprintf("apply_var_aliases: canonical 列已生成 (hdl_mgdl/tg_mgdl/ast_unl/alt_unl/alb_gl/...)\n\n"))

# Step 1: 主分析 cohort 选择 — Plan C: 仅 P_ rows (J ⊂ P_)
# Round 1 R-NHANES P0-1: 不能同时纳入 J 和 P_, 因为 P_ 包含 J 的所有参与者
df <- nhanes_all %>% filter(is_prepandemic == TRUE)
log_flow("Pre-pandemic cohort only (J subset of P_, dedupe)", nrow(df))

# Step 2: age + 非妊娠 (Round 5 R-Clinical P0 fix: 仅排除 RIDEXPRG==1 确认妊娠;
# RIDEXPRG==3 "无法判定" 是 missing-status 不是 pregnancy, 排除会引入 selection
# bias. 标准 NHANES 公开发表 convention 仅排除 confirmed pregnancies)
df <- df %>% filter(RIDAGEYR >= 20)
log_flow("Age >= 20", nrow(df))
df <- df %>% filter(is.na(RIDEXPRG) | RIDEXPRG != 1)
log_flow("+ not pregnant (excl RIDEXPRG==1 only)", nrow(df))

# Step 3: FibroScan valid (LUXSMED + LUXCAPM 非 NA)
# NOTE Round 1 R-Bias P0-6: IQR/median ≤30% filter 改 sensitivity (BMI-dependent selection)
# 主分析保留所有 valid LUXSMED+LUXCAPM, IQR ≤ 30% 作 S2 sensitivity
df <- df %>% filter(!is.na(LUXSMED), !is.na(LUXCAPM))
log_flow("+ LUX valid (LUXSMED + LUXCAPM non-NA)", nrow(df))
df <- df %>% mutate(
  lsm_iqr_valid = ifelse(!is.na(LUXSIQRM), LUXSIQRM <= 30, TRUE)
)
cat(sprintf("\n  LSM IQR/median ≤ 30%% (Castera 2019 strict, for sensitivity S2): %d\n",
            sum(df$lsm_iqr_valid)))

# Step 4: Se 双暴露
df <- df %>% filter(!is.na(LBXBSE), !is.na(DR1TSELE))
log_flow("+ LBXBSE (blood Se) + DR1TSELE (dietary Se)", nrow(df))

# Step 5: 肝功能 + CBC PLT
df <- df %>% filter(
  !is.na(ast_unl), !is.na(alt_unl), !is.na(alb_gl), !is.na(LBXPLTSI)
)
log_flow("+ BIOPRO (AST/ALT/Alb) + CBC PLT", nrow(df))

# Step 6: Se 双暴露 log + standardize (save scale attrs)
df <- df %>% mutate(
  ln_se_blood = log(pmax(LBXBSE, 0.01)),
  ln_se_diet  = log(pmax(DR1TSELE, 0.01))
)
scale_attrs <- list(
  se_blood = list(center = mean(df$ln_se_blood), scale = sd(df$ln_se_blood)),
  se_diet  = list(center = mean(df$ln_se_diet),  scale = sd(df$ln_se_diet))
)
df <- df %>% mutate(
  z_se_blood = (ln_se_blood - scale_attrs$se_blood$center) / scale_attrs$se_blood$scale,
  z_se_diet  = (ln_se_diet  - scale_attrs$se_diet$center)  / scale_attrs$se_diet$scale,
  se_blood_q = ntile(LBXBSE, 4),
  se_diet_q  = ntile(DR1TSELE, 4)
)

# Step 7: Dietary 比值 (救援切角 A; blood Zn/Cu 2017+ NHANES 停测)
df <- df %>% mutate(
  se_zn_ratio_diet = ifelse(!is.na(DR1TSELE) & !is.na(DR1TZINC) & DR1TZINC > 0,
                            DR1TSELE / DR1TZINC, NA),
  se_cu_ratio_diet = ifelse(!is.na(DR1TSELE) & !is.na(DR1TCOPP) & DR1TCOPP > 0,
                            DR1TSELE / DR1TCOPP, NA)
) %>% mutate(
  se_zn_q = ntile(se_zn_ratio_diet, 4),
  se_cu_q = ntile(se_cu_ratio_diet, 4)
)

# ============================================================
# Step 7b: v2 UPGRADE — 重金属共暴露 (P_PBCD 全集)
# 2026-05-18 升级: 加 Pb/Cd/Hg/Mn 作 co-exposures
# 用途: qgcomp/WQS/BKMR mixture; Hg-Se 拮抗专项
# 不作主筛, NA 保留, 下游 mixture 脚本 complete-case 自取
# 来源: var_aliases.R lead_ugdl/cadmium_ugl/hg_total_ugl/manganese_ugl
# ============================================================
df <- df %>% mutate(
  pb_ugdl  = lead_ugdl,
  cd_ugl   = cadmium_ugl,
  hg_ugl   = hg_total_ugl,
  mn_ugl   = manganese_ugl,
  se_ugl   = LBXBSE,                       # alias for symmetry with metal name
  ln_pb    = ifelse(!is.na(pb_ugdl)  & pb_ugdl  > 0, log(pb_ugdl), NA),
  ln_cd    = ifelse(!is.na(cd_ugl)   & cd_ugl   > 0, log(cd_ugl),  NA),
  ln_hg    = ifelse(!is.na(hg_ugl)   & hg_ugl   > 0, log(hg_ugl),  NA),
  ln_mn    = ifelse(!is.na(mn_ugl)   & mn_ugl   > 0, log(mn_ugl),  NA),
  # Hg/Se 摩尔比 (经典 redox 拮抗指标; Ralston 2017 Selenium-Health-Benefit Value)
  # Mw: Hg 200.59 g/mol, Se 78.96 g/mol
  # LBXTHG µg/L (= ng/mL), LBXBSE µg/L: 摩尔数 = mass(µg/L) / Mw(g/mol)
  hg_se_molar_ratio = ifelse(!is.na(hg_ugl) & !is.na(se_ugl) & se_ugl > 0,
                             (hg_ugl / 200.59) / (se_ugl / 78.96), NA),
  # Quartile bins (对 mixture analyses)
  pb_q = ntile(pb_ugdl, 4),
  cd_q = ntile(cd_ugl,  4),
  hg_q = ntile(hg_ugl,  4),
  mn_q = ntile(mn_ugl,  4)
)
cat(sprintf(
  "  v2 金属共暴露 (派生变量, 不作主筛):\n    Pb (LBXBPB) valid: %d  Cd (LBXBCD): %d  Hg (LBXTHG): %d  Mn (LBXBMN): %d\n    Hg/Se 摩尔比 valid: %d\n",
  sum(!is.na(df$pb_ugdl)), sum(!is.na(df$cd_ugl)),
  sum(!is.na(df$hg_ugl)),  sum(!is.na(df$mn_ugl)),
  sum(!is.na(df$hg_se_molar_ratio))
))

# ============================================================
# Step 7c: v2 UPGRADE — Vitamin D 撤回 (NCHS P_ 重映射 SEQN, J 行 VID 无法 join 进 P_)
# 2026-05-18 实查: J SEQN 93703-102956, P_ SEQN 109263-124822, 交集=0
# NCHS Pre-pandemic 文件公开发布时主动 re-SEQN 防 J 重识别. Vit D 只能做 J-only 敏感性
# 不进 main P_ analysis. 字段保留为 sensitivity 但不算主升级.
# ============================================================

# Step 8: FibroScan outcomes — Round 1 P0-2/P0-3 AASLD 2024 cut-offs
df <- df %>% mutate(
  lsm = LUXSMED,
  cap = LUXCAPM,
  # CAP: Primary AASLD 2024 ≥ 275 (Karlas 2017 ≥ 248 sensitivity)
  steatosis_cap275 = as.integer(cap >= 275),       # PRIMARY (AASLD 2024)
  steatosis_cap248 = as.integer(cap >= 248),       # sensitivity (Karlas 2017)
  # LSM: MASLD-specific AASLD 2024 (8/12 kPa) for primary
  fibrosis_lsm8  = as.integer(lsm >= 8.0),         # significant fibrosis ≥ F2
  fibrosis_lsm12 = as.integer(lsm >= 12.0),        # MASLD-specific F3-F4 (AASLD 2024)
  # Generic Castera 2019 mixed-etiology cut-offs (sensitivity)
  fibrosis_lsm9_5  = as.integer(lsm >= 9.5),
  fibrosis_lsm12_5 = as.integer(lsm >= 12.5),
  lsm_cat = factor(case_when(
    lsm < 8.0  ~ "F0-F1",
    lsm < 12.0 ~ "F2 (significant)",
    lsm >= 12.0 ~ "F3-F4 (advanced)"
  ), levels = c("F0-F1", "F2 (significant)", "F3-F4 (advanced)"))
)
cat(sprintf("\n  CAP ≥ 275 (steatosis, AASLD 2024 primary): %d (%.1f%%)\n",
            sum(df$steatosis_cap275, na.rm=TRUE), 100*mean(df$steatosis_cap275, na.rm=TRUE)))
cat(sprintf("  CAP ≥ 248 (sensitivity): %d (%.1f%%)\n",
            sum(df$steatosis_cap248, na.rm=TRUE), 100*mean(df$steatosis_cap248, na.rm=TRUE)))
cat(sprintf("  LSM ≥ 8.0 (sig fibrosis): %d (%.1f%%)\n",
            sum(df$fibrosis_lsm8, na.rm=TRUE), 100*mean(df$fibrosis_lsm8, na.rm=TRUE)))
cat(sprintf("  LSM ≥ 12.0 (MASLD-specific advanced): %d (%.1f%%)\n",
            sum(df$fibrosis_lsm12, na.rm=TRUE), 100*mean(df$fibrosis_lsm12, na.rm=TRUE)))

# Step 9: FIB-4 + age-adjusted (Round 3 Stage A — R-Clinical P0-B: age-adj 为 PRIMARY)
df$fib4 <- calc_fib4(df$RIDAGEYR, df$ast_unl, df$alt_unl, df$LBXPLTSI)
df$fib4_advanced <- fib4_binary_advanced(df$fib4)
df$fib4_advanced_aged <- fib4_advanced_age_adj(df$fib4, df$RIDAGEYR)
# Round 3 Stage A: 给 downstream 暴露 fib4_primary alias (McPherson 2017 age-adj)
# Legacy 1.30 cut-off 保留为 sensitivity (transparency)
df$fib4_primary <- df$fib4_advanced_aged   # PRIMARY (age-adj: <65 用 1.30, ≥65 用 2.0)
df$fib4_legacy  <- df$fib4_advanced        # SENSITIVITY (legacy 1.30 cut-off, all-ages)
cat(sprintf("  FIB-4 legacy ≥ 1.30 (sensitivity): %d (%.1f%%)\n",
            sum(df$fib4_legacy, na.rm=TRUE),
            100*mean(df$fib4_legacy, na.rm=TRUE)))
cat(sprintf("  FIB-4 age-adj advanced (PRIMARY, McPherson 2017): %d (%.1f%%)\n",
            sum(df$fib4_primary, na.rm=TRUE),
            100*mean(df$fib4_primary, na.rm=TRUE)))

# Step 10: 血压均值 (MOVED EARLIER — Round 1 P0-5)
df <- df %>% mutate(
  BPXSY1 = zero_to_na(BPXSY1), BPXSY2 = zero_to_na(BPXSY2),
  BPXSY3 = zero_to_na(BPXSY3), BPXSY4 = zero_to_na(BPXSY4),
  BPXDI1 = zero_to_na(BPXDI1), BPXDI2 = zero_to_na(BPXDI2),
  BPXDI3 = zero_to_na(BPXDI3), BPXDI4 = zero_to_na(BPXDI4)
) %>%
  rowwise() %>%
  mutate(
    sbp = mean(c(BPXSY2, BPXSY3, BPXSY4), na.rm = TRUE),
    dbp = mean(c(BPXDI2, BPXDI3, BPXDI4), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    sbp = ifelse(is.nan(sbp), NA, sbp),
    dbp = ifelse(is.nan(dbp), NA, dbp)
  )

# Step 11: Hepamet HFS (fasting subsample, albumin continuous fix)
df <- df %>% mutate(
  homa_ir = ifelse(!is.na(LBXIN) & !is.na(LBXGLU),
                   (LBXIN * LBXGLU) / 405, NA),
  sex_male = as.integer(RIAGENDR == 1),
  albumin_gdl = alb_gl / 10
)
# Diabetes: Round 3 Stage A — Wire diabetes_comprehensive() (ADA 2025, 7 criteria)
# (R-Clinical P0-A: 旧 5-criteria inline 缺 OGTT + RXQ_RX antidiabetic drug detection)
# Step A: 从 RXQ_RX 长表抽 antidiabetic drug hits 并 left_join 到 df (CLAUDE.md pitfall #2)
if (exists("rx_all") && nrow(rx_all) > 0 && "RXDDRUG" %in% names(rx_all)) {
  rx_dm <- detect_antidiabetic_rx(rx_all)
  df <- df %>% left_join(rx_dm, by = "SEQN") %>%
    mutate(rx_diabetes_yes = ifelse(is.na(rx_diabetes_yes), FALSE, rx_diabetes_yes))
  cat(sprintf("  RXQ_RX antidiabetic drug detected (rx_diabetes_yes=TRUE): %d\n",
              sum(df$rx_diabetes_yes, na.rm = TRUE)))
} else {
  cat("  [警告] rx_all 缺失或无 RXDDRUG 列, 跳过 RXQ_RX drug detection\n")
}
# Step B: 7-criteria composite (DIQ010/050/070 + HbA1c + FPG + OGTT 2h + RXQ_RX)
df$diabetes_raw <- diabetes_comprehensive(df)     # 含 NA, sensitivity 保留
n_dm_na <- sum(is.na(df$diabetes_raw))
cat(sprintf("  Diabetes raw composite: yes=%d / no=%d / NA=%d\n",
            sum(df$diabetes_raw == 1, na.rm = TRUE),
            sum(df$diabetes_raw == 0, na.rm = TRUE),
            n_dm_na))
# Round 4 R-NHANES P0: NA→0 默认行为 logged 显式 (下游 cm_risk_count / HFS 需 0/1
# 二值, 不接受 NA; ADA 2025 实操中"all-NA inputs" 等价 no diabetes signal → 视为
# not-diabetic, but 必须 log 该假设以便 S-sensitivity 复核)
df$diabetes <- df$diabetes_raw
df$diabetes[is.na(df$diabetes)] <- 0L
df$diabetes <- as.integer(df$diabetes)
cat(sprintf("  Diabetes (ADA 2025 7-criteria composite, NA→0 recoded %d rows): %d (%.1f%%)\n",
            n_dm_na, sum(df$diabetes, na.rm = TRUE), 100 * mean(df$diabetes, na.rm = TRUE)))
cat(sprintf("    [sensitivity reserve] df$diabetes_raw 保留 NA (供 S-sensitivity 复核)\n"))
dm_brk <- diabetes_criteria_breakdown(df)
print(dm_brk)
df$hfs <- with(df, ifelse(
  !is.na(homa_ir) & !is.na(albumin_gdl),
  calc_hepamet_hfs(sex_male, RIDAGEYR, ast_unl, albumin_gdl,
                   LBXPLTSI, homa_ir, diabetes),
  NA
))
df$hfs_cat <- hepamet_category(df$hfs)
df$hfs_high <- as.integer(!is.na(df$hfs) & df$hfs >= 0.47)

# ============================================================
# Step 11b: v2 UPGRADE — Redox biomarker outcomes
# 2026-05-18 升级: GGT/Iron/Bilirubin/Uric acid 作 secondary outcomes + mediators
# 用途: CMAverse 4-way (Se→GGT→HFS); 直接 outcome 测试; effect modifier
# 不作主筛, NA 保留
# 来源: var_aliases.R ggt_unl/total_bili_mgdl/uric_acid + 原始 LBXSIR
# ============================================================
df <- df %>% mutate(
  ggt_iul    = ggt_unl,                              # GGT U/L (氧化应激敏感)
  total_bili = total_bili_mgdl,                      # 总胆红素 mg/dL (内源抗氧化)
  uric_mgdl  = uric_acid,                            # 尿酸 mg/dL (双面抗氧化)
  iron_ugdl  = LBXSIR,                               # 血清铁 µg/dL (氧化驱动 Fenton)
  # GGT 异常阈 (sex-specific, ACG 2017 guidelines):
  # men > 40 U/L, women > 32 U/L
  ggt_high   = ifelse(!is.na(ggt_iul) & !is.na(sex_male),
                      as.integer((sex_male == 1 & ggt_iul > 40) |
                                 (sex_male == 0 & ggt_iul > 32)), NA),
  # log-transform for skewed distributions
  ln_ggt        = ifelse(!is.na(ggt_iul)    & ggt_iul    > 0, log(ggt_iul),    NA),
  ln_total_bili = ifelse(!is.na(total_bili) & total_bili > 0, log(total_bili), NA),
  ln_iron       = ifelse(!is.na(iron_ugdl)  & iron_ugdl  > 0, log(iron_ugdl),  NA),
  # Iron overload: ferritin not in NHANES Pre-pandemic; 用 iron_ugdl > sex-specific ULN
  # men > 160 µg/dL, women > 145 µg/dL (Mayo Clinic reference)
  iron_high  = ifelse(!is.na(iron_ugdl) & !is.na(sex_male),
                      as.integer((sex_male == 1 & iron_ugdl > 160) |
                                 (sex_male == 0 & iron_ugdl > 145)), NA),
  # 高尿酸血症 (Hyperuricemia, ACR 2020): >7 mg/dL 男, >6 mg/dL 女
  hyperuricemia = ifelse(!is.na(uric_mgdl) & !is.na(sex_male),
                         as.integer((sex_male == 1 & uric_mgdl > 7) |
                                    (sex_male == 0 & uric_mgdl > 6)), NA)
)
cat(sprintf(
  "  v2 Redox outcomes:\n    GGT (LBXSGTSI) valid: %d (high: %d)  Iron (LBXSIR): %d (high: %d)\n    Total bili (LBXSTB): %d  Uric acid (LBXSUA): %d (hyperuricemia: %d)\n",
  sum(!is.na(df$ggt_iul)),  sum(df$ggt_high == 1, na.rm = TRUE),
  sum(!is.na(df$iron_ugdl)), sum(df$iron_high == 1, na.rm = TRUE),
  sum(!is.na(df$total_bili)),
  sum(!is.na(df$uric_mgdl)), sum(df$hyperuricemia == 1, na.rm = TRUE)
))
# Note: intermediate HFS computation point — not the canonical analytic-cohort
# count. Canonical HFS chain (within 5,885 analytic cohort) appended at end.
# log_flow call deliberately omitted here to avoid misleading intermediate row;
# canonical 5,885 / 2,908 / 2,810 / 5 cascade added after final-cohort filter.

# Step 12: Alcohol + MetALD (用 CAP ≥ 275 主, 血压 fix)
df$ALQ130 <- na_codes(df$ALQ130, c(77, 99, 777, 999))
df$alcohol_gwk <- alcohol_gwk_from_alq130(df$ALQ130)
df$sex_chr <- ifelse(df$sex_male == 1, "Male", "Female")
df$BPQ020_c <- na_codes(df$BPQ020, c(7, 9))
df$htn_med <- as.integer(!is.na(df$BPQ020_c) & df$BPQ020_c == 1)
df$cm_risk_count <- count_cardiometabolic_risk(
  bmi = df$BMXBMI, waist = df$BMXWAIST, sex = df$sex_chr,
  fbg = df$LBXGLU, hba1c = df$LBXGH,
  t2d_yes = df$diabetes, t2d_med_yes = NA,
  sbp = df$sbp, dbp = df$dbp, htn_med_yes = df$htn_med,
  tg = df$tg_mgdl, lipid_med_yes = NA,
  hdl = df$hdl_mgdl
)
# 用 AASLD 2024 primary CAP ≥ 275
df$metald_group <- metald_classify(
  steatosis_yes = df$steatosis_cap275,
  alcohol_gwk   = df$alcohol_gwk,
  sex           = df$sex_chr,
  cm_risk_count = df$cm_risk_count
)
# Sensitivity: 用 CAP ≥ 248 (Karlas)
df$metald_group_cap248 <- metald_classify(
  steatosis_yes = df$steatosis_cap248,
  alcohol_gwk   = df$alcohol_gwk,
  sex           = df$sex_chr,
  cm_risk_count = df$cm_risk_count
)

# Step 13: Self-reported liver
df$MCQ160L <- na_codes(df$MCQ160L, c(7, 9))
df$selfreport_liver <- ifelse(df$MCQ160L %in% c(1, 2),
                               as.integer(df$MCQ160L == 1), NA)

# Step 14: Mortality (006 不做主分析, 字段保留为 sensitivity only)
# Round 1 R-Bias note: P_ 2019-March 2020 部分死亡未跟踪到 NCHS 2019 release, 仅 J 部分有 mortality
# 此处仅 J 部分 (2017-2018) 有 ELIGSTAT==1 with follow-up; 仍保留字段供后续 sensitivity (S4)
df <- df %>% mutate(
  mort_allcause_sens = ifelse(!is.na(MORTSTAT), as.integer(MORTSTAT == 1), NA),
  permth = ifelse(!is.na(PERMTH_EXM), PERMTH_EXM, PERMTH_INT)
)

# Step 15: 协变量 (drink 用 ALQ111, ALQ101 2017+ 停用)
df <- df %>% mutate(
  DMDEDUC2 = na_codes(DMDEDUC2, c(7, 9)),
  DMDMARTL = na_codes(DMDMARTL, c(77, 99)),
  SMQ020   = na_codes(SMQ020, c(7, 9)),
  ALQ111   = na_codes(ALQ111, c(7, 9))
)
df <- df %>% mutate(
  age = RIDAGEYR,
  age_group = factor(case_when(
    age < 40 ~ "20-39", age < 60 ~ "40-59", TRUE ~ ">=60"
  ), levels = c("20-39", "40-59", ">=60")),
  race = factor(recode(as.character(RIDRETH1),
                       "1" = "Mexican American", "2" = "Other Hispanic",
                       "3" = "Non-Hispanic White", "4" = "Non-Hispanic Black",
                       "5" = "Other Race"),
                levels = c("Non-Hispanic White", "Non-Hispanic Black",
                           "Mexican American", "Other Hispanic", "Other Race")),
  education = factor(case_when(
    DMDEDUC2 %in% c(1, 2) ~ "Less than HS",
    DMDEDUC2 == 3 ~ "High school",
    DMDEDUC2 %in% c(4, 5) ~ "College or above"
  ), levels = c("Less than HS", "High school", "College or above")),
  pir = INDFMPIR,
  pir_group = factor(case_when(
    pir <= 1.3 ~ "<=1.3", pir <= 3.5 ~ "1.3-3.5", pir > 3.5 ~ ">3.5"
  ), levels = c("<=1.3", "1.3-3.5", ">3.5")),
  bmi = BMXBMI,
  bmi_cat = factor(case_when(
    bmi < 25 ~ "<25", bmi < 30 ~ "25-29.9", bmi >= 30 ~ ">=30"
  ), levels = c("<25", "25-29.9", ">=30")),
  smoke = factor(case_when(SMQ020 == 2 ~ "Never", SMQ020 == 1 ~ "Ever"),
                 levels = c("Never", "Ever")),
  drink = factor(case_when(ALQ111 == 1 ~ "Yes", ALQ111 == 2 ~ "No"),
                 levels = c("No", "Yes")),
  hypertension = df$htn_med
)

# Step 16: Pooled weight — Round 3 Stage A wire cycle_config (R-NHANES P1-1 / R-TechDebt P0-3)
# 006 主 cohort = P_ only (J ⊂ P_ dedupe), weight = WTMECPRP, 但保留 helper 调用保持架构一致
cycles_006_main <- make_cycle_spec(year_range = integer(0), include_prepandemic = TRUE)
validate_no_J_in_PP(cycles_006_main)  # 防御 — 仅 P_, validator 应静默通过
cycle_years_006_main <- cycle_years_vector(cycles_006_main)
stopifnot("PrePandemic_2017_March2020" %in% names(cycle_years_006_main))
cat(sprintf("  cycle_years_006_main (from cycle_config.R): %s\n",
            paste(names(cycle_years_006_main), cycle_years_006_main, sep = "=", collapse = ", ")))
df$wt_pooled     <- pooled_mec_weight(df, cycle_years_006_main)
df$wt_diet_pooled <- pooled_diet_weight(df, cycle_years_006_main)
df$wt_saf_pooled <- pooled_saf_weight(df, cycle_years_006_main)

# Step 17: 核心协变量完整
core_cov_strict <- c("age", "race", "education", "pir", "bmi",
                     "RIAGENDR", "SDMVPSU", "SDMVSTRA")
n_before <- nrow(df)
df <- df %>% filter(if_all(all_of(core_cov_strict), ~ !is.na(.)))
log_flow(sprintf("Core covariate-complete (main analytic cohort; %d excluded)", n_before - nrow(df)), nrow(df))

# Save + 汇总
nhanes_final <- df
cat("\n========================================\n")
cat(sprintf("最终样本 (主分析 = P_ only): %d 行 × %d 列\n",
            nrow(nhanes_final), ncol(nhanes_final)))
cat(sprintf("  Se blood (LBXBSE) 中位: %.1f µg/L\n",
            median(nhanes_final$LBXBSE, na.rm=TRUE)))
cat(sprintf("  Se dietary (DR1TSELE) 中位: %.1f µg/day\n",
            median(nhanes_final$DR1TSELE, na.rm=TRUE)))
cat(sprintf("\n  CAP ≥ 275 (steatosis, PRIMARY AASLD 2024): %d (%.1f%%)\n",
            sum(nhanes_final$steatosis_cap275, na.rm=TRUE),
            100*mean(nhanes_final$steatosis_cap275, na.rm=TRUE)))
cat(sprintf("  LSM ≥ 8.0 (sig fibrosis): %d (%.1f%%)\n",
            sum(nhanes_final$fibrosis_lsm8, na.rm=TRUE),
            100*mean(nhanes_final$fibrosis_lsm8, na.rm=TRUE)))
cat(sprintf("  LSM ≥ 12.0 (MASLD-specific advanced): %d (%.1f%%)\n",
            sum(nhanes_final$fibrosis_lsm12, na.rm=TRUE),
            100*mean(nhanes_final$fibrosis_lsm12, na.rm=TRUE)))
cat(sprintf("  HFS-eligible (fasting): %d / HFS ≥ 0.47: %d\n",
            sum(!is.na(nhanes_final$hfs)),
            sum(nhanes_final$hfs_high, na.rm=TRUE)))
cat(sprintf("  FIB-4 legacy ≥ 1.30 (sensitivity): %d (%.1f%%)\n",
            sum(nhanes_final$fib4_legacy, na.rm=TRUE),
            100*mean(nhanes_final$fib4_legacy, na.rm=TRUE)))
cat(sprintf("  FIB-4 age-adj advanced (PRIMARY, McPherson 2017): %d (%.1f%%)\n",
            sum(nhanes_final$fib4_primary, na.rm=TRUE),
            100*mean(nhanes_final$fib4_primary, na.rm=TRUE)))
cat(sprintf("  Diabetes (ADA 2025 composite): %d (%.1f%%)\n",
            sum(nhanes_final$diabetes, na.rm=TRUE),
            100*mean(nhanes_final$diabetes, na.rm=TRUE)))
cat(sprintf("\n  MetALD groups (CAP ≥ 275 primary, 血压 fix):\n"))
print(table(nhanes_final$metald_group, useNA = "ifany"))
cat(sprintf("\n  MetALD groups (CAP ≥ 248 sensitivity):\n"))
print(table(nhanes_final$metald_group_cap248, useNA = "ifany"))
cat(sprintf("\n  Self-reported liver: %d (%.1f%%)\n",
            sum(nhanes_final$selfreport_liver, na.rm=TRUE),
            100*mean(nhanes_final$selfreport_liver, na.rm=TRUE)))
cat(sprintf("========================================\n"))

if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
save(nhanes_final, mort_all, scale_attrs,
     file = "data/processed/nhanes_final.RData")
cat("已保存 data/processed/nhanes_final.RData\n")

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)

# W17 canonical HFS chain in the 5,885 analytic cohort
# (avoids the pre-W11 misleading intermediate "HFS-eligible: 3349 / HFS ≥ 0.47: 7"
# which was computed on the 6,747 pre-final-filter cohort).
n_hfs_eligible_final  <- sum(!is.na(nhanes_final$hfs))
n_hfs_auroc_subset    <- sum(!is.na(nhanes_final$hfs) &
                              !is.na(nhanes_final$fib4) &
                              !is.na(nhanes_final$wt_saf_pooled) &
                              nhanes_final$wt_saf_pooled > 0)
n_hfs_high_in_auroc   <- sum(!is.na(nhanes_final$hfs) &
                              !is.na(nhanes_final$fib4) &
                              !is.na(nhanes_final$wt_saf_pooled) &
                              nhanes_final$wt_saf_pooled > 0 &
                              nhanes_final$hfs_high == 1)
log_flow("HFS-eligible (fasting insulin+glucose, within 5,885)", n_hfs_eligible_final)
log_flow("HFS AUROC analytic subset (HFS + FIB-4 + WTSAFPRP>0)", n_hfs_auroc_subset)
log_flow("HFS >= 0.47 high-risk cell (within AUROC subset)", n_hfs_high_in_auroc)

flow_df <- do.call(rbind, flow)
write.csv(flow_df, "output/tables/flow_counts.csv", row.names = FALSE)
cat("Saved output/tables/flow_counts.csv\n")
