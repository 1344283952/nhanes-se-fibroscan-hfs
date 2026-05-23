# ============================================
# 006 / 13_cross_classification.R — 16-cell Q1-Q4 DR1TSELE × Q1-Q4 LBXBSE
#
# Per task.md §7.1.3: cross-classification 16-cell heatmap for CAP ≥ 275
#   (and LSM ≥ 8 as supplementary).
# Each cell = adjusted prevalence (M2) of outcome.
# ============================================

set.seed(20260516)
options(survey.lonely.psu = "adjust")  # W11 R2 R-Stats P1: Lumley 2010 §3.4 single-PSU stratum behaviour
suppressPackageStartupMessages({
  library(dplyr); library(survey)
})
cat("========================================\n")
cat("006 W5 — Q4 × Q4 cross-classification heatmap\n")
cat("========================================\n\n")

load("data/processed/nhanes_design.RData")
load("data/processed/nhanes_final.RData")

df <- nhanes_final %>%
  filter(!is.na(DR1TSELE), !is.na(LBXBSE),
         !is.na(steatosis_cap275), !is.na(fibrosis_lsm8),
         !is.na(wt_pooled), wt_pooled > 0) %>%
  mutate(
    diet_q  = factor(ntile(DR1TSELE, 4), levels = 1:4,
                     labels = c("D1", "D2", "D3", "D4")),
    blood_q = factor(ntile(LBXBSE, 4), levels = 1:4,
                     labels = c("B1", "B2", "B3", "B4"))
  )
cat(sprintf("Analytic N = %d\n", nrow(df)))

des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pooled,
                 data = df, nest = TRUE)

# Adjusted prevalence per cell via svyby
cell_cap275 <- svyby(~steatosis_cap275, ~ diet_q + blood_q,
                     design = des, FUN = svymean, na.rm = TRUE)
cell_lsm8 <- svyby(~fibrosis_lsm8, ~ diet_q + blood_q,
                   design = des, FUN = svymean, na.rm = TRUE)

cat("\nCAP >= 275 prevalence per Se 4x4 cell:\n")
print(cell_cap275)
cat("\nLSM >= 8 prevalence per Se 4x4 cell:\n")
print(cell_lsm8)

# Heat-map data prep
cap_heat <- cell_cap275 %>%
  rename(prev = steatosis_cap275, se = se) %>%
  mutate(diet_q = factor(diet_q, levels = c("D1","D2","D3","D4")),
         blood_q = factor(blood_q, levels = c("B4","B3","B2","B1")))  # reverse for plot
lsm_heat <- cell_lsm8 %>%
  rename(prev = fibrosis_lsm8, se = se) %>%
  mutate(diet_q = factor(diet_q, levels = c("D1","D2","D3","D4")),
         blood_q = factor(blood_q, levels = c("B4","B3","B2","B1")))

# Quick PNG heatmap
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
png("output/figures/cross_cap275_heatmap.png", width = 700, height = 600, res = 110)
par(mar = c(5, 5, 3, 1))
mat <- matrix(NA, 4, 4, dimnames = list(c("B4","B3","B2","B1"), c("D1","D2","D3","D4")))
for (i in seq_len(nrow(cap_heat))) {
  mat[as.character(cap_heat$blood_q[i]), as.character(cap_heat$diet_q[i])] <- cap_heat$prev[i]
}
image(1:4, 1:4, t(mat[nrow(mat):1, ]),
      xaxt = "n", yaxt = "n",
      xlab = "Dietary Se quartile (DR1TSELE)",
      ylab = "Blood Se quartile (LBXBSE)",
      col = hcl.colors(20, "viridis"),
      main = "CAP >= 275 prevalence (Rinella 2023): 16-cell Se cross-classification")
axis(1, at = 1:4, labels = c("D1 (low)", "D2", "D3", "D4 (high)"))
axis(2, at = 1:4, labels = c("B1 (low)", "B2", "B3", "B4 (high)"))
for (i in 1:4) for (j in 1:4) {
  v <- mat[5 - i, j]
  if (!is.na(v)) text(j, i, sprintf("%.1f%%", 100 * v), cex = 0.9,
                       col = ifelse(v > 0.5, "white", "black"))
}
dev.off()

# Save tables
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(cell_cap275, "output/tables/cross_cap275_cells.csv", row.names = FALSE)
write.csv(cell_lsm8,   "output/tables/cross_lsm8_cells.csv", row.names = FALSE)
save(cell_cap275, cell_lsm8, cap_heat, lsm_heat,
     file = "output/tables/cross_classification_results.RData")

cat("\n保存:\n")
cat("  output/tables/cross_cap275_cells.csv (16 cells with SE)\n")
cat("  output/tables/cross_lsm8_cells.csv\n")
cat("  output/figures/cross_cap275_heatmap.png\n")
cat("\nDONE 13_cross_classification.R\n")
