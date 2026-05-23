# ============================================
# 006 / _retrodesign_smallcells.R — Type S/M post-hoc detectable effect (X5)
# Round 2 R-Causal-Methods + R-Stats v2 P0:
# Per Gelman & Carlin 2014, when an analysis is underpowered, the conditional
# probabilities of (a) wrong sign (Type S) and (b) effect-size inflation
# (Type M, also called magnitude) are non-trivial. Required for:
#   - MetALD substratum n=576 (planned stratified GAM smooths)
#   - HFS ≥ 0.47 cell n=5 (in AUROC analytic subset n=2,810; categorical NRI/IDI
#     underpowered per Pepe 2014). W11 R2 R-DataChain P0 fix: was n=6 stale.
# ============================================

set.seed(20260516)
# retrodesign() (Gelman & Carlin 2014) helper — no CRAN package universally
# available, so reproduce the 8-line implementation here
retrodesign <- function(A, s, alpha = 0.05, df = Inf, n.sims = 10000) {
  # A = hypothesized true effect size; s = standard error; alpha = 2-sided
  z <- qt(1 - alpha / 2, df)
  p.hi  <- 1 - pt(z - A / s, df)
  p.lo  <-     pt(-z - A / s, df)
  power <- p.hi + p.lo
  type.s <- p.lo / power
  estimate <- A + s * rt(n.sims, df)
  significant <- abs(estimate) > s * z
  exaggeration <- mean(abs(estimate)[significant]) / abs(A)
  list(power = power, type.s = type.s, exaggeration = exaggeration)
}

cat("========================================\n")
cat("006 / retrodesign Type S/M (X5) — MetALD + HFS≥0.47\n")
cat("========================================\n\n")

# ---- Case 1: MetALD substratum n=576, stratified GAM detectable OR ----
# Assume binary outcome at base prevalence 12% (LSM≥8 in MetALD cell)
# A hypothesized minimal clinically relevant OR = 1.30 → log(OR)= 0.262
# Approximate SE for logistic OR on n=576 with p≈0.12:
#   var(log_OR) ≈ 1/(n * p * (1-p)) ≈ 1/(576 * 0.12 * 0.88) = 1/60.83 → SE≈ 0.128
metald_res <- retrodesign(A = log(1.30), s = 0.128, alpha = 0.05, df = Inf)
cat("MetALD n=576 (LSM≥8 outcome, hypothesized OR=1.30):\n")
cat(sprintf("  Statistical power (alpha=0.05): %.3f\n", metald_res$power))
cat(sprintf("  Type S (P[wrong sign | significant]): %.3f\n", metald_res$type.s))
cat(sprintf("  Type M (exaggeration ratio): %.2fx\n\n", metald_res$exaggeration))

# Worse-power scenario: hypothesized smaller OR=1.15
metald_res_small <- retrodesign(A = log(1.15), s = 0.128, alpha = 0.05, df = Inf)
cat("MetALD n=576 with smaller hypothesized OR=1.15:\n")
cat(sprintf("  Power: %.3f / Type S: %.3f / Type M: %.2fx\n\n",
            metald_res_small$power, metald_res_small$type.s,
            metald_res_small$exaggeration))

# ---- Case 2: HFS≥0.47 cell n=5, categorical NRI ----
# Categorical NRI for binary high-risk vs low-risk vs FIB-4. Underpowered
# with n_high = 5 (W11 R2 P0: was n=6 stale; the AUROC analytic subset
# n=2,810 contains exactly 5 HFS≥0.47 cases). Hypothesized NRI ~ 0.20.
# Pepe 2014 power floor: log_NRI ≈ 0.20, SE scales as 1/sqrt(n_high).
# SE for n=5 ≈ 0.40 × sqrt(6/5) = 0.438 (slightly larger than the n=6 case)
nri_res <- retrodesign(A = 0.20, s = 0.438, alpha = 0.05, df = Inf)
cat("HFS≥0.47 cell n=5 (categorical NRI vs FIB-4, hypothesized NRI=0.20):\n")
cat(sprintf("  Power: %.3f\n", nri_res$power))
cat(sprintf("  Type S: %.3f\n", nri_res$type.s))
cat(sprintf("  Type M: %.2fx\n\n", nri_res$exaggeration))

# Continuous AUROC by comparison (powered on n=2,810 with ΔAUROC=0.057)
# SE(ΔAUROC) ≈ 0.012 per DeLong/Hanley
auroc_res <- retrodesign(A = 0.057, s = 0.012, alpha = 0.05, df = Inf)
cat("Continuous ΔAUROC (n=2,810, observed Δ=0.057, SE≈0.012):\n")
cat(sprintf("  Power: %.3f / Type S: %.4f / Type M: %.3fx\n\n",
            auroc_res$power, auroc_res$type.s, auroc_res$exaggeration))

# Save
res_tbl <- data.frame(
  analysis = c("MetALD_n576_OR1.30", "MetALD_n576_OR1.15",
               "HFS_high_n5_NRI0.20", "Continuous_dAUROC_n2810"),
  hypothesized_effect = c(log(1.30), log(1.15), 0.20, 0.057),
  SE = c(0.128, 0.128, 0.438, 0.012),
  power = c(metald_res$power, metald_res_small$power,
            nri_res$power, auroc_res$power),
  type_S = c(metald_res$type.s, metald_res_small$type.s,
             nri_res$type.s, auroc_res$type.s),
  type_M = c(metald_res$exaggeration, metald_res_small$exaggeration,
             nri_res$exaggeration, auroc_res$exaggeration)
)
print(res_tbl)
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(res_tbl, "output/tables/retrodesign_smallcells.csv", row.names = FALSE)

cat("\nInterpretation (Gelman & Carlin 2014):\n")
cat("  Power > 0.80 = adequately powered\n")
cat("  Type S > 5% = high risk of getting the sign wrong if significant\n")
cat("  Type M > 2x = significant findings inflate the true effect by >2x\n")
cat("\n保存:\n")
cat("  output/tables/retrodesign_smallcells.csv\n")
cat("\nDONE _retrodesign_smallcells.R\n")
