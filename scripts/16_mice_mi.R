# 006 — Multiple Imputation Sensitivity (mice m=20, OSF v1.5 S9)
# Predictive mean matching m=20 + Rubin pool
# Re-fit primary CAP/LSM/HFS logistic + compare to complete-case

suppressPackageStartupMessages({
  library(mice)
  library(dplyr)
  library(broom)
})

set.seed(20260516)
cat("==========================================================\n")
cat("006 Multiple Imputation Sensitivity (mice m=20)\n")
cat("OSF v1.5 S9 (2026-05-18)\n")
cat("==========================================================\n\n")

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE)

log_file <- file.path("logs", paste0("mi_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

# Load main cohort
data_path <- "data/processed/nhanes_final.RData"
if (!file.exists(data_path)) {
  cat("!! data/processed/nhanes_final.RData missing\n")
  cat("!! Run 03_clean_data.R first\n")
  quit(save = "no", status = 1)
}
load(data_path)

# Try common naming conventions
df <- if (exists("nhanes_final")) nhanes_final else
      if (exists("nhanes_main")) nhanes_main else
      stop("Could not find cohort data frame in nhanes_final.RData")

cat("Cohort N =", nrow(df), "\n")

# ============================================================
# 1. Variables for imputation
# ============================================================
mi_vars <- c("BMXBMI", "DR1TKCAL", "INDFMPIR", "DMDEDUC2",
             "SMQ020", "ALQ111", "LBXBSE", "DR1TSELE",
             "LUXCAPM", "LUXSMED", "LBXSATSI", "LBXSASSI", "LBDSALSI", "LBXPLTSI",
             "RIDAGEYR", "RIAGENDR", "RIDRETH3")
mi_vars <- intersect(mi_vars, names(df))
sub <- df[, mi_vars]

cat("\n[1] Imputation variables:", length(mi_vars), "vars\n")
miss_pct <- round(colMeans(is.na(sub)) * 100, 2)
print(data.frame(var = names(miss_pct), miss_pct = miss_pct))

# ============================================================
# 2. mice imputation m=20
# ============================================================
cat("\n[2] Running mice m=20 (PMM, maxit=10)\n")
ptm <- proc.time()
imp <- mice(sub, m = 20, method = "pmm", maxit = 10, printFlag = FALSE, seed = 20260516)
elapsed <- (proc.time() - ptm)[3]
cat("  mice finished in", round(elapsed/60, 1), "min\n")
saveRDS(imp, "output/tables/mi_imp_m20.rds")

# ============================================================
# 3. Re-fit primary: CAP >= 275 ~ LBXBSE + LBXBSE^2 + DR1TSELE + adjust
# ============================================================
cat("\n[3] Re-fit primary CAP>=275 logistic on m=20 imputed sets\n")
fits_cap <- with(imp, glm(I(LUXCAPM >= 275) ~ LBXBSE + I(LBXBSE^2) + DR1TSELE +
                          RIDAGEYR + RIAGENDR + factor(RIDRETH3) + factor(DMDEDUC2) + INDFMPIR +
                          factor(SMQ020) + factor(ALQ111) + BMXBMI + DR1TKCAL,
                        family = quasibinomial))
pool_cap <- pool(fits_cap)
sum_cap <- summary(pool_cap, conf.int = TRUE, exponentiate = TRUE)
sum_cap$endpoint <- "CAP>=275"
print(sum_cap)
write.csv(sum_cap, "output/tables/mi_cap_logistic.csv", row.names = FALSE)

# ============================================================
# 4. Re-fit primary: LSM >= 8 ~ ...
# ============================================================
cat("\n[4] Re-fit primary LSM>=8 logistic on m=20 imputed sets\n")
fits_lsm <- with(imp, glm(I(LUXSMED >= 8) ~ LBXBSE + I(LBXBSE^2) + DR1TSELE +
                          RIDAGEYR + RIAGENDR + factor(RIDRETH3) + factor(DMDEDUC2) + INDFMPIR +
                          factor(SMQ020) + factor(ALQ111) + BMXBMI + DR1TKCAL,
                        family = quasibinomial))
pool_lsm <- pool(fits_lsm)
sum_lsm <- summary(pool_lsm, conf.int = TRUE, exponentiate = TRUE)
sum_lsm$endpoint <- "LSM>=8"
write.csv(sum_lsm, "output/tables/mi_lsm_logistic.csv", row.names = FALSE)

# ============================================================
# 5. Compare to complete-case
# ============================================================
cat("\n[5] Complete-case fits for comparison\n")
cc_cap <- glm(I(LUXCAPM >= 275) ~ LBXBSE + I(LBXBSE^2) + DR1TSELE +
              RIDAGEYR + RIAGENDR + factor(RIDRETH3) + factor(DMDEDUC2) + INDFMPIR +
              factor(SMQ020) + factor(ALQ111) + BMXBMI + DR1TKCAL,
            data = df, family = quasibinomial)
cc_cap_tidy <- broom::tidy(cc_cap, conf.int = TRUE, exponentiate = TRUE)
cc_cap_tidy$endpoint <- "CAP>=275"
write.csv(cc_cap_tidy, "output/tables/mi_cap_completecase.csv", row.names = FALSE)

# Side-by-side comparison table
key_term <- "LBXBSE"
mi_se_or <- sum_cap[sum_cap$term == key_term, c("estimate", "conf.low", "conf.high", "p.value")]
cc_se_or <- cc_cap_tidy[cc_cap_tidy$term == key_term, c("estimate", "conf.low", "conf.high", "p.value")]

side_by_side <- rbind(
  data.frame(method = "Complete-case", as.list(cc_se_or)),
  data.frame(method = "mice m=20", as.list(mi_se_or))
)
write.csv(side_by_side, "output/tables/mi_vs_cc_comparison.csv", row.names = FALSE)
cat("\n=== Side-by-side comparison: LBXBSE effect on CAP>=275 ===\n")
print(side_by_side)

# Convergence check
plot_path <- "output/figures/mi_convergence.png"
tryCatch({
  png(plot_path, width = 1200, height = 800, res = 100)
  plot(imp, layout = c(2, 4))
  dev.off()
  cat("\n  Convergence plot saved:", plot_path, "\n")
}, error = function(e) cat("  plot fail:", e$message, "\n"))

cat("\n==========================================================\n")
cat("MI sensitivity complete (elapsed", round(elapsed/60, 1), "min)\n")
cat("Output: output/tables/mi_cap_logistic.csv + mi_lsm_logistic.csv +\n")
cat("        mi_vs_cc_comparison.csv + mi_imp_m20.rds\n")
cat("==========================================================\n")
