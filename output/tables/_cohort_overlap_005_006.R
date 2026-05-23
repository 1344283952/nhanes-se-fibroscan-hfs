# ============================================
# Cohort overlap analysis: 005 (metal_fib4_liver) vs 006 (se_fibroscan_hfs)
# Computes SEQN overlap between the two analytic cohorts, with explicit
# handling of CDC's P_ SEQN re-randomization policy (J ↔ P_ mapping is
# RDC-only, not publicly disclosed).
#
# Output: _cohort_overlap_005_006.csv
# ============================================

cat("========================================\n")
cat("Cohort overlap analysis: 005 vs 006\n")
cat("========================================\n\n")

# --- Load 005 ---
load("D:/临时学习/nhanes-pipeline/nhanes-pipeline/projects/005_metal_fib4_liver/data/processed/nhanes_final.RData")
df005_main <- nhanes_final           # 2011-2018 G/H/I/J main analytic cohort
df005_mort <- df_mort                # mortality-linked subset of main
df005_pp   <- nhanes_pp_sens         # P_ pre-pandemic sensitivity cohort (S5)
rm(nhanes_final, df_mort, nhanes_pp_sens)
cat(sprintf("005 main (2011-2018 G/H/I/J): N = %d\n", nrow(df005_main)))
cat(sprintf("005 mort cohort:              N = %d\n", nrow(df005_mort)))
cat(sprintf("005 P_ sensitivity (S5):      N = %d\n", nrow(df005_pp)))

# --- Load 006 ---
load("D:/临时学习/nhanes-pipeline/nhanes-pipeline/projects/006_se_fibroscan_hfs/data/processed/nhanes_final.RData")
df006 <- nhanes_final                # P_ only (J biologically ⊂ P_, but SEQN re-randomized)
rm(nhanes_final)
cat(sprintf("006 main (P_ only):           N = %d\n\n", nrow(df006)))

# SEQN ranges (sanity)
cat("===== SEQN ranges =====\n")
cat(sprintf("005 main G (2011-12) SEQN: %d - %d\n",
            min(df005_main$SEQN[df005_main$cycle_tag == "NHANES_2011_2012"]),
            max(df005_main$SEQN[df005_main$cycle_tag == "NHANES_2011_2012"])))
cat(sprintf("005 main H (2013-14) SEQN: %d - %d\n",
            min(df005_main$SEQN[df005_main$cycle_tag == "NHANES_2013_2014"]),
            max(df005_main$SEQN[df005_main$cycle_tag == "NHANES_2013_2014"])))
cat(sprintf("005 main I (2015-16) SEQN: %d - %d\n",
            min(df005_main$SEQN[df005_main$cycle_tag == "NHANES_2015_2016"]),
            max(df005_main$SEQN[df005_main$cycle_tag == "NHANES_2015_2016"])))
cat(sprintf("005 main J (2017-18) SEQN: %d - %d\n",
            min(df005_main$SEQN[df005_main$cycle_tag == "NHANES_2017_2018"]),
            max(df005_main$SEQN[df005_main$cycle_tag == "NHANES_2017_2018"])))
cat(sprintf("005 P_ sens SEQN:          %d - %d\n",
            min(df005_pp$SEQN), max(df005_pp$SEQN)))
cat(sprintf("006 main (P_) SEQN:        %d - %d\n",
            min(df006$SEQN), max(df006$SEQN)))

# --- Direct SEQN intersection ---
seqn_005_main <- unique(df005_main$SEQN)
seqn_006_main <- unique(df006$SEQN)
shared_main_direct <- intersect(seqn_005_main, seqn_006_main)
cat(sprintf("\n===== Direct SEQN intersection =====\n"))
cat(sprintf("005 main vs 006 main shared SEQN: %d\n", length(shared_main_direct)))
cat("(Expected 0 — CDC re-randomized SEQN in P_ files; J↔P_ mapping is RDC-only)\n")

# --- 005 PP sensitivity vs 006 (both pull from P_ → real overlap) ---
seqn_005_pp <- unique(df005_pp$SEQN)
shared_pp_vs_006 <- intersect(seqn_005_pp, seqn_006_main)
cat(sprintf("\n===== 005 P_ sensitivity vs 006 (both from P_ files) =====\n"))
cat(sprintf("005 PP sens N:                %d\n", length(seqn_005_pp)))
cat(sprintf("006 main N:                   %d\n", length(seqn_006_main)))
cat(sprintf("Shared SEQN (both pull P_):   %d\n", length(shared_pp_vs_006)))

# --- Estimate biological participant overlap ---
# The J cycle (Aug 2017-2018) is BIOLOGICALLY a subset of P_ (Aug 2017-March 2020).
# Under CDC's P_ design, J participants account for ~9,254 of 15,560 P_ rows
# (raw merged level). After analytic exclusions:
#   - 005 keeps 4,154 J participants (PBCD metals + BIOPRO + PLT complete + age>=20, etc.)
#   - 006 keeps 5,885 P_ participants (FibroScan valid + Se + BIOPRO + PLT, etc.)
# A J participant is BIOLOGICALLY in 006 iff:
#   - they survived 006's exclusion criteria (FibroScan valid, age>=20, Se, etc.)
# We cannot identify which specific 006 SEQN corresponds to which 005 J SEQN
# (CDC RDC-only), but we can compute the EXPECTED biological overlap as:
#   006 P_ × Pr(J | P_) = 5,885 × (J_raw / P_raw) = 5,885 × 9,254 / 15,560 ≈ 3,500
# This is an UPPER BOUND on the biological cohort overlap.

n_J_raw <- 9254     # NHANES_2017_2018 raw rows (from 005 raw_merged earlier)
n_P_raw <- 15560    # PrePandemic_2017_March2020 raw rows
n_006 <- nrow(df006)
n_005_main <- nrow(df005_main)
n_005_J <- sum(df005_main$cycle_tag == "NHANES_2017_2018")

frac_J_in_P <- n_J_raw / n_P_raw
biol_overlap_upper_006side <- n_006 * frac_J_in_P
biol_overlap_upper_005side <- n_005_J  # all 4,154 J participants are biologically in P_

cat(sprintf("\n===== Estimated BIOLOGICAL participant overlap =====\n"))
cat(sprintf("(CDC RDC-only mapping prevents exact match; this is upper-bound estimate)\n"))
cat(sprintf("J raw rows / P_ raw rows = %d / %d = %.4f (fraction of P_ that is J)\n",
            n_J_raw, n_P_raw, frac_J_in_P))
cat(sprintf("006 main × Pr(J|P_): %.0f participants are J-cycle biologically\n",
            biol_overlap_upper_006side))
cat(sprintf("005 main cycle J N (all biologically in P_): %d\n",
            biol_overlap_upper_005side))
cat(sprintf("Biological overlap = min of (5,885 × fraction-passing-006-exclusions among J)\n"))
cat(sprintf("                          (4,154 × fraction-passing-006-exclusions)\n"))

# --- Build CSV output ---
overlap_rows <- data.frame(
  comparison = c(
    "005_main_total",
    "006_main_total",
    "005_cycle_2011_2012_only",
    "005_cycle_2013_2014_only",
    "005_cycle_2015_2016_only",
    "005_cycle_2017_2018_only",
    "005_PP_sensitivity_total",
    "005_mortality_cohort_total",
    "shared_SEQN_005_main_vs_006_main",
    "shared_SEQN_005_PP_sens_vs_006_main",
    "estimated_biological_overlap_upper_bound",
    "raw_NHANES_J_rows",
    "raw_NHANES_P_rows",
    "Pr_J_given_Prepandemic"
  ),
  n_seqn = c(
    n_005_main,
    n_006,
    sum(df005_main$cycle_tag == "NHANES_2011_2012"),
    sum(df005_main$cycle_tag == "NHANES_2013_2014"),
    sum(df005_main$cycle_tag == "NHANES_2015_2016"),
    n_005_J,
    nrow(df005_pp),
    nrow(df005_mort),
    length(shared_main_direct),
    length(shared_pp_vs_006),
    round(biol_overlap_upper_006side),
    n_J_raw,
    n_P_raw,
    round(frac_J_in_P, 4)
  ),
  notes = c(
    "005 main analytic cohort, NHANES 2011-2018 (G/H/I/J) post-exclusion",
    "006 main analytic cohort, NHANES P_ pre-pandemic (Aug 2017-March 2020) post-exclusion",
    "005 cycle G subset (no FibroScan available; no overlap with 006)",
    "005 cycle H subset",
    "005 cycle I subset",
    "005 cycle J subset — biologically all in P_, but CDC SEQN re-randomized",
    "005 P_ sensitivity (S5, separate cohort drawn from same P_ files as 006)",
    "005 mortality-linked subset of main (ELIGSTAT==1); P_ has no mortality linkage",
    "Direct SEQN match = 0. CDC re-randomized SEQN in P_; J↔P_ mapping is NCHS RDC-only. SEQN ranges: J = 93,703-102,956; P_ = 109,263-124,822.",
    "005 P_ sensitivity cohort and 006 main BOTH pull from P_ files, so SEQN match works here. Difference (6,599 vs 5,885) reflects 006's additional FibroScan-valid + Se filter.",
    "Estimated biological participant overlap (upper bound): 006 main × (J raw / P_ raw) = 5,885 × (9,254/15,560). True overlap is smaller because 006 applies FibroScan-valid filter that excludes a fraction of J participants.",
    "Raw NHANES_J rows in 005 merged data (before exclusions)",
    "Raw P_ rows in 005/006 merged data (before exclusions)",
    "Fraction of P_ raw rows that were originally measured in NHANES cycle J"
  ),
  stringsAsFactors = FALSE
)

out_path <- "D:/临时学习/nhanes-pipeline/nhanes-pipeline/projects/006_se_fibroscan_hfs/output/tables/_cohort_overlap_005_006.csv"
write.csv(overlap_rows, out_path, row.names = FALSE)
cat(sprintf("\n========== CSV WRITTEN ==========\n%s\n", out_path))

# --- Cover letter friendly disclosure ---
cat("\n========== COVER LETTER DISCLOSURE ==========\n\n")

cat("OPTION A (technically precise, SEQN-level):\n\n")
cat(sprintf(paste0(
  "Of the analytic cohorts (006: n=%d from NHANES pre-pandemic 2017-March 2020; ",
  "005: n=%d from NHANES 2011-2018), no participants share SEQN identifiers ",
  "because CDC re-randomized SEQN in the pre-pandemic (P_) release. However, ",
  "biologically the two cohorts may share up to ~%d participants who took part ",
  "in NHANES cycle J (2017-2018), which forms a subset of the P_ pre-pandemic ",
  "dataset. Cycle-level mapping between J and P_ SEQN is restricted to NCHS RDC ",
  "and not publicly available.\n"),
  n_006, n_005_main, round(biol_overlap_upper_006side)))

cat("\nOPTION B (concise, for cover letter):\n\n")
cat(sprintf(paste0(
  "Of the n=%d (006) and n=%d (005) analytic cohorts, the two studies draw ",
  "from non-overlapping NHANES public-use SEQN ranges (005: 2011-2018 cycles ",
  "G/H/I/J; 006: P_ pre-pandemic, Aug 2017-March 2020), so SEQN-level overlap ",
  "is zero. Biologically, however, an estimated ~%d participants who took part ",
  "in NHANES cycle J (2017-2018) appear in both datasets under different ",
  "CDC-assigned SEQN (J↔P_ mapping is NCHS Research Data Center restricted). ",
  "We declare this biological overlap for transparency; outcomes (FIB-4 vs ",
  "FibroScan), exposures (5-metal mixture vs Se dual-exposure), and analytic ",
  "models differ between the two studies.\n"),
  n_006, n_005_main, round(biol_overlap_upper_006side)))

cat("\n========== DONE ==========\n")
