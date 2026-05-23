# 006 — Two-sample Mendelian Randomization: Blood Selenium → MASLD/NAFLD
# OSF v1.5 amendment H6 (2026-05-18, v2)
# Backend: ieugwasr (OpenGWAS API) + MendelianRandomization (CRAN) [+ MRPRESSO if available]
# Bypasses TwoSampleMR (GitHub install unreliable on Windows R 4.6 ucrt)
#
# IVW + MR-Egger + Weighted Median + Simple Median
# + Egger intercept (pleiotropy) + Cochran Q (heterogeneity)
# + leave-one-out + MR-PRESSO (optional)

suppressPackageStartupMessages({
  library(ieugwasr)
  library(MendelianRandomization)
})
has_presso <- requireNamespace("MRPRESSO", quietly = TRUE)

set.seed(20260516)

cat("==========================================================\n")
cat("006 Two-sample MR (v2): Blood Selenium -> MASLD/NAFLD\n")
cat("OSF H6 (v1.5 2026-05-18)\n")
cat("Backend: ieugwasr + MendelianRandomization", if (has_presso) "+ MRPRESSO" else "(no MRPRESSO)", "\n")
cat("==========================================================\n\n")

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE)

log_file <- file.path("logs", paste0("mr_v2_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

# ============================================================
# 1. Exposure: Cornelis 2015 blood-selenium GWAS (ieu-a-1077)
# ============================================================
cat("[1] Fetching blood-Se instruments (Cornelis 2015 ieu-a-1077)\n")
exp_dat <- tryCatch({
  th <- tophits(id = "ieu-a-1077", clump = TRUE, r2 = 0.001, kb = 10000)
  if (is.null(th) || nrow(th) == 0) stop("empty")
  th
}, error = function(e) {
  cat("  ! tophits API fail:", e$message, "\n")
  cat("  ! Fallback: Cornelis 2015 hand-coded 2 lead SNPs (Table 2)\n")
  data.frame(rsid = c("rs921943", "rs6859667"),
             beta = c(0.094, 0.110), se = c(0.011, 0.014),
             ea = c("T", "A"), nea = c("C", "G"),
             eaf = c(0.27, 0.10), p = c(2e-19, 1e-15),
             trait = "Selenium",
             stringsAsFactors = FALSE)
})

cat("  Got", nrow(exp_dat), "Se instruments\n")
cols_to_print <- intersect(c("rsid","beta","se","ea","nea","eaf","p"), names(exp_dat))
print(exp_dat[, cols_to_print, drop = FALSE])
write.csv(exp_dat, "output/tables/mr_instruments_se.csv", row.names = FALSE)

F_stat <- mean((exp_dat$beta / exp_dat$se)^2, na.rm = TRUE)
cat("\n  Mean F-statistic:", round(F_stat, 1), "(>10 = strong instruments)\n")

# ============================================================
# 2. Outcome: NAFLD/MASLD GWAS
# ============================================================
cat("\n[2] Fetching NAFLD/MASLD outcome via ieugwasr::associations\n")
out_ids <- c("ebi-a-GCST90091033", "ieu-b-7", "finn-b-K11_NAFLD")
out_dat <- NULL; chosen <- NA
for (oid in out_ids) {
  cat("  Trying:", oid, "\n")
  od <- tryCatch(associations(variants = exp_dat$rsid, id = oid),
                 error = function(e) { cat("    ! err:", e$message, "\n"); NULL })
  if (!is.null(od) && nrow(od) >= 1) {
    cat("    ✓ matched", nrow(od), "SNPs\n")
    out_dat <- od; chosen <- oid; break
  }
}

if (is.null(out_dat)) {
  cat("\n!! No outcome GWAS data could be fetched.\n")
  writeLines(c(
    "MR Step 2 failed: ieugwasr::associations returned no matches.",
    paste("Tried:", paste(out_ids, collapse = ", ")),
    paste("Timestamp:", Sys.time()),
    "Likely cause: OpenGWAS server unreachable from this workstation or SNPs not present.",
    "Fallback: manually download Anstee 2020 GWAS summary stats or use hand-coded effect sizes for Liu 2024 Sci Rep reference values."
  ), "output/tables/mr_FAILED_note.txt")
  quit(save = "no", status = 0)
}
cat("\n  Selected outcome GWAS:", chosen, "\n")
write.csv(out_dat, "output/tables/mr_outcome_masld.csv", row.names = FALSE)

# ============================================================
# 3. Harmonise
# ============================================================
cat("\n[3] Harmonising exposure-outcome SNPs by effect allele\n")
har <- merge(exp_dat, out_dat, by = "rsid", suffixes = c("_exp", "_out"))
cat("  Merged on rsid:", nrow(har), "SNPs\n")

ea_exp <- toupper(if ("ea_exp" %in% names(har)) har$ea_exp else har$ea)
ea_out_col <- intersect(c("ea_out", "effect_allele", "ea.outcome"), names(har))[1]
beta_out_col <- intersect(c("beta_out", "beta.outcome"), names(har))[1]
se_out_col <- intersect(c("se_out", "se.outcome"), names(har))[1]

if (!is.na(ea_out_col)) {
  ea_out <- toupper(har[[ea_out_col]])
  flip <- !is.na(ea_out) & ea_exp != ea_out
  cat("  Effect-allele flips:", sum(flip), "of", nrow(har), "\n")
  har$beta_y_aligned <- ifelse(flip, -1 * har[[beta_out_col]], har[[beta_out_col]])
} else {
  cat("  ! Outcome effect-allele column not found — assuming aligned\n")
  har$beta_y_aligned <- har[[beta_out_col]]
}

beta_x <- if ("beta_exp" %in% names(har)) har$beta_exp else har$beta
se_x <- if ("se_exp" %in% names(har)) har$se_exp else har$se
beta_y <- har$beta_y_aligned
se_y <- har[[se_out_col]]

ok <- complete.cases(beta_x, se_x, beta_y, se_y)
beta_x <- beta_x[ok]; se_x <- se_x[ok]
beta_y <- beta_y[ok]; se_y <- se_y[ok]
snp_keep <- har$rsid[ok]
cat("  Final SNPs after QC:", length(beta_x), "\n")

write.csv(data.frame(rsid = snp_keep, beta_x, se_x, beta_y, se_y),
          "output/tables/mr_harmonised.csv", row.names = FALSE)

if (length(beta_x) < 2) {
  cat("!! <2 SNPs — cannot run MR\n")
  quit(save = "no", status = 0)
}

# ============================================================
# 4. MR estimates
# ============================================================
cat("\n[4] MR estimates (IVW + MR-Egger + Weighted Median + Simple Median)\n")
mr_obj <- mr_input(bx = beta_x, bxse = se_x, by = beta_y, byse = se_y,
                   exposure = "Blood Se", outcome = chosen, snps = snp_keep)

ivw <- mr_ivw(mr_obj)
egger <- mr_egger(mr_obj)
wmed <- mr_median(mr_obj, weighting = "weighted", iterations = 1000)
smed <- tryCatch(mr_median(mr_obj, weighting = "simple"), error = function(e) NULL)

build_row <- function(method, est, se, p, lo, hi, n) {
  data.frame(method = method, nsnp = n, beta = est, se = se, pval = p,
             CIlower = lo, CIupper = hi,
             or = exp(est), or_lci = exp(lo), or_uci = exp(hi))
}

res <- rbind(
  build_row("IVW", ivw@Estimate, ivw@StdError, ivw@Pvalue, ivw@CILower, ivw@CIUpper, ivw@SNPs),
  build_row("MR-Egger", egger@Estimate, egger@StdError.Est, egger@Pvalue.Est,
            egger@CILower.Est, egger@CIUpper.Est, egger@SNPs),
  build_row("Weighted Median", wmed@Estimate, wmed@StdError, wmed@Pvalue,
            wmed@CILower, wmed@CIUpper, wmed@SNPs)
)
if (!is.null(smed)) {
  res <- rbind(res, build_row("Simple Median", smed@Estimate, smed@StdError, smed@Pvalue,
                               smed@CILower, smed@CIUpper, smed@SNPs))
}
print(res)
write.csv(res, "output/tables/mr_estimates.csv", row.names = FALSE)

# ============================================================
# 5. Sensitivity
# ============================================================
cat("\n[5] Egger intercept (directional pleiotropy)\n")
ei <- data.frame(intercept = egger@Intercept, se = egger@StdError.Int,
                 pval = egger@Pvalue.Int,
                 CIlower = egger@CILower.Int, CIupper = egger@CIUpper.Int)
print(ei)
write.csv(ei, "output/tables/mr_egger_intercept.csv", row.names = FALSE)

cat("\n[6] Cochran Q heterogeneity (IVW)\n")
het <- data.frame(method = "IVW",
                  Q = ivw@Heter.Stat[1], dof = ivw@Heter.Stat[2], Qpval = ivw@Heter.Stat[3])
print(het)
write.csv(het, "output/tables/mr_heterogeneity.csv", row.names = FALSE)

if (has_presso && length(beta_x) >= 4) {
  cat("\n[7] MR-PRESSO (outlier + distortion)\n")
  presso_input <- data.frame(bx = beta_x, by = beta_y, sebx = se_x, seby = se_y)
  presso <- tryCatch({
    MRPRESSO::mr_presso(BetaOutcome = "by", BetaExposure = "bx",
                        SdOutcome = "seby", SdExposure = "sebx",
                        data = presso_input, OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                        NbDistribution = 1000, SignifThreshold = 0.05)
  }, error = function(e) { cat("  ! err:", e$message, "\n"); NULL })
  if (!is.null(presso)) {
    saveRDS(presso, "output/tables/mr_presso.rds")
    cat("  PRESSO Global Test P =", presso$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
  }
} else cat("\n[7] Skip MR-PRESSO (need >=4 SNPs AND MRPRESSO installed)\n")

if (length(beta_x) >= 3) {
  cat("\n[8] Leave-one-out\n")
  loo <- data.frame()
  for (i in seq_along(beta_x)) {
    obj_i <- mr_input(bx = beta_x[-i], bxse = se_x[-i],
                      by = beta_y[-i], byse = se_y[-i])
    ivw_i <- mr_ivw(obj_i)
    loo <- rbind(loo, data.frame(dropped_idx = i, dropped_snp = snp_keep[i],
                                  beta = ivw_i@Estimate, se = ivw_i@StdError,
                                  pval = ivw_i@Pvalue))
  }
  print(loo)
  write.csv(loo, "output/tables/mr_leaveoneout.csv", row.names = FALSE)
}

# ============================================================
# 6. Summary
# ============================================================
cat("\n==========================================================\n")
cat("MR analysis complete\n")
cat("Primary IVW: OR =", round(res$or[1], 3),
    "(", round(res$or_lci[1], 3), "-", round(res$or_uci[1], 3), ")",
    "P =", signif(res$pval[1], 3), "\n")
cat("Egger intercept P =", signif(ei$pval, 3), "(P>0.05 = no directional pleiotropy)\n")
cat("Cochran Q P =", signif(het$Qpval, 3), "\n")
cat("Mean F-stat =", round(F_stat, 1), "\n")
cat("Outcome GWAS used:", chosen, "\n")
cat("==========================================================\n")
