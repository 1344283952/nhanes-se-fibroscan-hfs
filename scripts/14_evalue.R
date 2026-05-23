# ============================================
# 006 / 14_evalue.R — E-value (VanderWeele 2017) for ratio OR + MetALD effects
#
# Per task.md §7.8 + OSF §5:
#   Logistic OR (CAP ≥ 275 prevalence 44.3%, NON-rare): convert OR → RR via
#   Zhang & Yu 1998 DOI 10.1001/jama.280.19.1690 before E-value (rare=FALSE).
#   PR_zhang_yu = OR / (1 - p0 + p0 * OR) where p0 = baseline outcome prob
# ============================================

set.seed(20260516)
suppressPackageStartupMessages({
  library(dplyr); library(EValue)
})
cat("========================================\n")
cat("006 — E-value for ratio OR + HFS AUROC results\n")
cat("========================================\n\n")

load("output/tables/ratio_analysis_results.RData")    # or_df from ratios
load("data/processed/nhanes_final.RData")

# Baseline CAP ≥ 275 prevalence
p0 <- mean(nhanes_final$steatosis_cap275 == 1, na.rm = TRUE)
cat(sprintf("Baseline CAP ≥ 275 prevalence p0 = %.3f (non-rare; OR → RR conversion needed)\n", p0))

# Zhang & Yu 1998 OR → RR conversion
or_to_rr <- function(or, p0) {
  if (is.na(or) || or <= 0) return(NA_real_)
  or / (1 - p0 + p0 * or)
}

# ---- E-value for ratio ORs ----
cat("\n[1/2] Ratio tertile OR E-values (after Zhang & Yu 1998 OR→RR conversion) ...\n")
ratio_ev <- or_df %>%
  rowwise() %>%
  mutate(RR = or_to_rr(OR, p0),
         RR_lci = or_to_rr(lci, p0),
         RR_uci = or_to_rr(uci, p0)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(ev = list(
    if (is.na(RR) || RR <= 0) c(EV_point = NA, EV_lcl = NA) else
      tryCatch({
        e <- evalues.RR(est = RR, lo = RR_lci, hi = RR_uci)
        c(EV_point = e["E-values", "point"], EV_lcl = e["E-values", "lower"])
      }, error = function(e) c(EV_point = NA, EV_lcl = NA))
  )) %>%
  ungroup() %>%
  mutate(EV_point = sapply(ev, function(x) as.numeric(x["EV_point"])),
         EV_lcl   = sapply(ev, function(x) as.numeric(x["EV_lcl"]))) %>%
  select(-ev)
print(ratio_ev %>% select(ratio, term, OR, lci, uci, RR, EV_point, EV_lcl))

# ---- E-value for HFS prediction (AUROC is not HR/OR; this section reports ----
#       only that the AUROC of 0.731 is well above the AUROC=0.50 null,
#       so unmeasured confounding doesn't apply directly to discrimination.
#       Document this limitation rather than computing a pseudo E-value.)
cat("\n[2/2] HFS AUROC interpretation:\n")
cat("  AUROC is a discrimination metric, not an exposure-outcome effect estimate.\n")
cat("  VanderWeele 2017 E-value is for confounding-adjusted HR/OR/RR.\n")
cat("  HFS H5 confirmation (AUROC=0.731 ≥ 0.65) is not subject to E-value adjustment.\n")
cat("  Documented as Limitation in manuscript Methods §2.6.\n")

# Save
if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
write.csv(ratio_ev %>% select(ratio, term, OR, lci, uci, RR, RR_lci, RR_uci, EV_point, EV_lcl),
          "output/tables/evalue_ratio_OR.csv", row.names = FALSE)
save(ratio_ev, p0, file = "output/tables/evalue_results.RData")

cat("\n保存:\n")
cat("  output/tables/evalue_ratio_OR.csv (Se/Zn + Se/Cu tertile OR + E-value)\n")
cat("\nDONE 14_evalue.R\n")
