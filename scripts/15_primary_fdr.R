# ============================================
# 006 / 15_primary_fdr.R — Primary-family BH-FDR
# Round 2 R-Stats v2 P0 fix: Methods §2.5 claims BH-FDR on primary family
# (4 Se exposures × 3 outcomes = 12 tests) but original 10_subgroup_forest.R
# only adjusted the 108-test subgroup family. This script collects the 12
# primary p-values from RCS (07) + ratio (12) + GAM main smooth (06)
# and applies BH q=0.05.
#
# Primary family (12 tests, per OSF §5):
#   - 6 RCS non-linearity tests (3 outcomes × 2 exposures): from rcs_pvalues.csv
#   - 4 ratio tertile tests (Se/Zn T2,T3 + Se/Cu T2,T3 against CAP≥275): from
#     ratio_OR_cap275.csv (T2,T3 each ratio = 4 tests)
#   - 2 GAM tensor-interaction tests (CAP, LSM): from gam_dual_se.RData
# Total = 12 tests
# ============================================

set.seed(20260516)
suppressPackageStartupMessages({
  library(dplyr); library(mgcv)
})
cat("========================================\n")
cat("006 / 15_primary_fdr — primary 12-test BH-FDR (Round 2 R-Stats P0)\n")
cat("========================================\n\n")

p_primary <- data.frame(
  test_id  = character(0),
  family   = character(0),
  exposure = character(0),
  outcome  = character(0),
  p_raw    = numeric(0),
  stringsAsFactors = FALSE
)

# ---- 1. RCS non-linearity (6 tests: 3 outcomes × 2 exposures) ----
if (file.exists("output/tables/rcs_pvalues.csv")) {
  rcs <- read.csv("output/tables/rcs_pvalues.csv", stringsAsFactors = FALSE)
  for (i in seq_len(nrow(rcs))) {
    p_primary <- rbind(p_primary, data.frame(
      test_id  = sprintf("RCS_nonlin_%s_%s", rcs$exposure[i], rcs$outcome[i]),
      family   = "RCS_nonlinearity",
      exposure = rcs$exposure[i],
      outcome  = rcs$outcome[i],
      p_raw    = as.numeric(rcs$p_nonlinear_anova[i]),
      stringsAsFactors = FALSE))
  }
  cat(sprintf("  RCS family: %d tests loaded\n", nrow(rcs)))
} else cat("  WARN: rcs_pvalues.csv missing — skip RCS family\n")

# ---- 2. Ratio tertile (4 tests: Se/Zn T2,T3 + Se/Cu T2,T3) ----
if (file.exists("output/tables/ratio_OR_cap275.csv")) {
  rat <- read.csv("output/tables/ratio_OR_cap275.csv", stringsAsFactors = FALSE)
  # Each ratio_var has 2 rows (T2, T3); 2 ratios × 2 tertiles = 4 tests
  for (i in seq_len(nrow(rat))) {
    p_primary <- rbind(p_primary, data.frame(
      test_id  = sprintf("ratio_%s", rat$term[i]),
      family   = "ratio_tertile",
      exposure = rat$ratio[i],
      outcome  = "steatosis_cap275",
      p_raw    = as.numeric(rat$p[i]),
      stringsAsFactors = FALSE))
  }
  cat(sprintf("  Ratio family: %d tests loaded\n", nrow(rat)))
} else cat("  WARN: ratio_OR_cap275.csv missing — skip ratio family\n")

# ---- 3. GAM tensor-interaction (2 tests: CAP, LSM) ----
if (file.exists("output/tables/gam_dual_se.RData")) {
  load("output/tables/gam_dual_se.RData")
  for (oname in names(results)) {
    if (oname == "hfs") next  # only CAP + LSM primary
    s <- results[[oname]]$summary
    if (!is.null(s) && "s.table" %in% names(s)) {
      st <- as.data.frame(s$s.table)
      ti_row <- grep("ti\\(", rownames(st), value = TRUE)
      if (length(ti_row) >= 1) {
        pti <- as.numeric(st[ti_row[1], "p-value"])
        p_primary <- rbind(p_primary, data.frame(
          test_id  = sprintf("GAM_tensor_%s", oname),
          family   = "GAM_tensor",
          exposure = "DR1TSELE_x_LBXBSE",
          outcome  = oname,
          p_raw    = pti,
          stringsAsFactors = FALSE))
      }
    }
  }
  cat(sprintf("  GAM tensor family: %d tests loaded\n",
              sum(p_primary$family == "GAM_tensor")))
} else cat("  WARN: gam_dual_se.RData missing — skip GAM tensor family\n")

# ---- BH-FDR ----
cat(sprintf("\nTotal primary tests: %d\n", nrow(p_primary)))
if (nrow(p_primary) == 0) stop("No primary p-values collected — re-run upstream scripts")
p_primary$p_BH_primary <- p.adjust(p_primary$p_raw, method = "BH")
p_primary <- p_primary[order(p_primary$p_raw), ]
p_primary$signif_BH_005 <- as.integer(p_primary$p_BH_primary < 0.05)

cat("\nPrimary-family BH-FDR result (q=0.05):\n")
print(p_primary)
cat(sprintf("\nN signif at BH q<0.05: %d / %d\n",
            sum(p_primary$signif_BH_005), nrow(p_primary)))

# Save
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(p_primary, "output/tables/primary_fdr.csv", row.names = FALSE)
save(p_primary, file = "output/tables/primary_fdr.RData")
cat("\n保存:\n")
cat("  output/tables/primary_fdr.csv\n")
cat("  output/tables/primary_fdr.RData\n")
cat("\nDONE 15_primary_fdr.R\n")
