# ============================================
# 006v2 / _03b_clean_K_cycle.R
# 复制主分析清洗 pipeline 到 NHANES K cycle (2021-2023)
# 输入: data/k_cycle_raw/*.xpt
# 输出: data/processed/nhanes_k_final.RData
# ============================================
library(dplyr); library(haven); library(tidyr); library(purrr)
source("../../templates/_shared/fib4.R")
source("../../templates/_shared/hepamet_hfs.R")
set.seed(20260518)

cat("==== K cycle 清洗 (复制主分析 pipeline) ====\n")
raw <- "data/k_cycle_raw"

# 读 + merge 所有需要的模块
read_x <- function(f) tryCatch(haven::read_xpt(file.path(raw, paste0(f, ".xpt"))),
                                error = function(e) NULL)
demo  <- read_x("DEMO_L")
bmx   <- read_x("BMX_L")
bpxo  <- read_x("BPXO_L")   # 注意 K cycle 用 oscillometric BPXO 不是 BPX
smq   <- read_x("SMQ_L")
alq   <- read_x("ALQ_L")
diq   <- read_x("DIQ_L")
bpq   <- read_x("BPQ_L")
mcq   <- read_x("MCQ_L")
paq   <- read_x("PAQ_L")
pbcd  <- read_x("PBCD_L")
lux   <- read_x("LUX_L")
bio   <- read_x("BIOPRO_L")
cbc   <- read_x("CBC_L")
ins   <- read_x("INS_L")
glu   <- read_x("GLU_L")
ghb   <- read_x("GHB_L")
hdl   <- read_x("HDL_L")
tch   <- read_x("TCHOL_L")
trg   <- read_x("TRIGLY_L")
dr1   <- read_x("DR1TOT_L")
dr2   <- read_x("DR2TOT_L")

# Left join 全部到 DEMO 上
df <- demo
join_safe <- function(d, x) {
  if (is.null(x)) return(d)
  overlap <- intersect(colnames(d), colnames(x))
  overlap <- setdiff(overlap, "SEQN")
  if (length(overlap) > 0) x <- x %>% select(-all_of(overlap))
  d %>% left_join(x, by = "SEQN")
}
for (x in list(bmx, bpxo, smq, alq, diq, bpq, mcq, paq, pbcd, lux,
               bio, cbc, ins, glu, ghb, hdl, tch, trg, dr1, dr2)) {
  df <- join_safe(df, x)
}
cat(sprintf("Merged K cycle: %d 行 × %d 列\n", nrow(df), ncol(df)))

flow <- list()
log_flow <- function(label, n) {
  cat(sprintf("  [流程] %-58s n = %d\n", label, n))
  flow[[length(flow) + 1]] <<- data.frame(step = label, n = n, stringsAsFactors = FALSE)
}
log_flow("Raw merged K cycle", nrow(df))

# Step 1: Age >= 20 + 非妊娠
df <- df %>% filter(RIDAGEYR >= 20)
log_flow("Age >= 20", nrow(df))
if ("RIDEXPRG" %in% colnames(df)) {
  df <- df %>% filter(is.na(RIDEXPRG) | RIDEXPRG != 1)
}
log_flow("+ not pregnant", nrow(df))

# Step 2: FibroScan valid
df <- df %>% filter(!is.na(LUXSMED), !is.na(LUXCAPM))
log_flow("+ LUX valid (FibroScan)", nrow(df))

# Step 3: Se 双暴露
df <- df %>% filter(!is.na(LBXBSE), !is.na(DR1TSELE))
log_flow("+ Se dual (LBXBSE + DR1TSELE)", nrow(df))

# Step 4: BIOPRO basic + CBC PLT
df <- df %>% filter(!is.na(LBXSASSI), !is.na(LBXSATSI), !is.na(LBXSAL), !is.na(LBXPLTSI))
log_flow("+ BIOPRO + CBC PLT", nrow(df))

# Step 5: 衍生 — Se 双 + 比值 + 5 metals + GGT outcomes (与 P_ 主分析一致)
df <- df %>% mutate(
  ln_se_blood = log(pmax(LBXBSE, 0.01)),
  ln_se_diet  = log(pmax(DR1TSELE, 0.01)),
  se_blood_q  = ntile(LBXBSE, 4),
  se_diet_q   = ntile(DR1TSELE, 4),
  se_zn_ratio_diet = ifelse(!is.na(DR1TZINC) & DR1TZINC > 0, DR1TSELE / DR1TZINC, NA),
  se_cu_ratio_diet = ifelse(!is.na(DR1TCOPP) & DR1TCOPP > 0, DR1TSELE / DR1TCOPP, NA),
  cap = LUXCAPM, lsm = LUXSMED,
  steatosis_cap275 = as.integer(cap >= 275),
  fibrosis_lsm8    = as.integer(lsm >= 8.0),
  fibrosis_lsm12   = as.integer(lsm >= 12.0),
  sex_male = as.integer(RIAGENDR == 1),
  ast_unl  = LBXSASSI,
  alt_unl  = LBXSATSI,
  alb_gl   = LBXSAL,
  albumin_gdl = alb_gl / 10,
  homa_ir  = ifelse(!is.na(LBXIN) & !is.na(LBXGLU), (LBXIN * LBXGLU) / 405, NA),
  # 4 redox biomarkers
  ggt_iul  = LBXSGTSI,
  uric_mgdl = LBXSUA,
  iron_ugdl = LBXSIR,
  total_bili = LBXSTB,
  ggt_high = ifelse(!is.na(ggt_iul) & !is.na(sex_male),
                    as.integer((sex_male == 1 & ggt_iul > 40) |
                               (sex_male == 0 & ggt_iul > 32)), NA),
  # 4 co-exposure metals
  pb_ugdl  = LBXBPB,
  cd_ugl   = LBXBCD,
  hg_ugl   = LBXTHG,
  mn_ugl   = LBXBMN
)

# Step 6: 排除 (与 P_ 主分析一致 — age + race + edu + smoke + alcohol + bmi + kcal complete)
df <- df %>% mutate(
  age = RIDAGEYR,
  race = factor(RIDRETH1),
  DMDEDUC2 = ifelse(DMDEDUC2 %in% c(7, 9), NA, DMDEDUC2),
  SMQ020 = ifelse(SMQ020 %in% c(7, 9), NA, SMQ020),
  ALQ111 = ifelse(ALQ111 %in% c(7, 9), NA, ALQ111),
  bmi = BMXBMI,
  pir = INDFMPIR,
  kcal = DR1TKCAL
)
core_cov <- c("age", "race", "DMDEDUC2", "bmi", "pir", "kcal", "SMQ020", "ALQ111")
df_main <- df %>% filter(complete.cases(.[, core_cov]))
log_flow("+ Core covariate complete (main analytic K cohort)", nrow(df_main))

# Step 7: HFS subset + AUROC subset
df_main$hfs <- with(df_main, ifelse(
  !is.na(homa_ir) & !is.na(albumin_gdl),
  calc_hepamet_hfs(sex_male, age, ast_unl, albumin_gdl, LBXPLTSI, homa_ir, 0),  # diabetes=0 简化
  NA))
hfs_sub <- df_main %>% filter(!is.na(hfs))
log_flow("HFS-eligible subset", nrow(hfs_sub))

# Step 8: 权重 (K cycle 用 WTMEC2YR)
df_main$wt_mec <- df_main$WTMEC2YR
df_main$wt_diet <- if("WTDRD1" %in% colnames(df_main)) df_main$WTDRD1 else df_main$WTMEC2YR

# Save
nhanes_k_final <- df_main
save(nhanes_k_final, file = "data/processed/nhanes_k_final.RData")
cat(sprintf("\n[OK] saved nhanes_k_final.RData (n=%d)\n", nrow(nhanes_k_final)))
write.csv(do.call(rbind, flow), "output/tables/k_cycle_flow.csv", row.names = FALSE)
