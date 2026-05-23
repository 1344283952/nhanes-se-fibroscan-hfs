# ============================================
# 006 / run_all.R — End-to-end pipeline runner
# Round 2 R-DataChain P0: file was absent (process-standard violation vs 001-005);
# this is the canonical entry point per CLAUDE.md.
# Execute from project root: Rscript scripts/run_all.R
# Total wall time: ~30-60 min depending on bootstrap reps.
# ============================================

# Ensure run from project dir
if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("run_all.R must be executed from project root (projects/006_se_fibroscan_hfs/).")
}

set.seed(20260516)
options(survey.lonely.psu = "adjust")  # global lonely-PSU policy

START_T <- Sys.time()
cat(sprintf("\n[run_all] start %s\n\n", format(START_T, "%Y-%m-%d %H:%M:%S")))

scripts <- c(
  "scripts/01_download_data.R",          # skipped if data already present
  "scripts/02_merge_data.R",
  "scripts/03_clean_data.R",
  "scripts/04_survey_design.R",
  "scripts/05_table1.R",
  "scripts/06_gam_dual_exposure.R",
  "scripts/07_rcs_ushape.R",
  "scripts/08_metald_stratified.R",
  "scripts/09_dag.R",                    # DAG render (independent)
  "scripts/09_hfs_predict.R",            # HFS prediction + HL + Harrell C + NRI/IDI
  "scripts/10_subgroup_forest.R",
  "scripts/11_sensitivity.R",
  "scripts/12_ratio_analysis.R",
  "scripts/13_cross_classification.R",
  "scripts/14_consort.R",
  "scripts/14_evalue.R",
  "scripts/15_primary_fdr.R",            # Round 2 R-Stats P0: 12-test BH-FDR
  "scripts/_ipw_selection_sens.R",       # X3 IPW selection-bias sensitivity (S8)
  "scripts/_retrodesign_smallcells.R",   # X5 Type S/M for small cells
  "scripts/_render_figures.R",
  "scripts/_render_fig2_dag.R"
)

# Optionally skip download if data exists
if (file.exists("data/processed/nhanes_raw_merged.RData")) {
  scripts <- scripts[!grepl("01_download_data", scripts)]
  cat("[run_all] data/processed/nhanes_raw_merged.RData present — skipping 01_download_data.R\n\n")
}

# W17 R-DataChain fresh-test fix: loop variable renamed from `s` to `.script_file`
# because sourced scripts (e.g., 12_ratio_analysis.R) rebind `s` to a large
# summary() coefficient matrix, then sprintf("%-45s", s) tries to format the
# matrix and overflows R's 8192-char sprintf limit.
for (.script_file in scripts) {
  if (!file.exists(.script_file)) {
    cat(sprintf("[run_all] SKIP %s (file not found)\n", .script_file)); next
  }
  .t0 <- Sys.time()
  cat(sprintf("[run_all] --- source %s ---\n", .script_file))
  .res <- tryCatch({
    source(.script_file, echo = FALSE, local = new.env(parent = globalenv()))
    list(ok = TRUE, err = NULL)
  }, error = function(e) list(ok = FALSE, err = conditionMessage(e)))
  .dt <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  if (.res$ok) {
    cat(sprintf("[run_all] OK  %-45s (%.1f s)\n\n", .script_file, .dt))
  } else {
    cat(sprintf("[run_all] FAIL %-45s (%.1f s)\n  ERROR: %s\n\n",
                .script_file, .dt, .res$err))
  }
}

END_T <- Sys.time()
cat(sprintf("\n[run_all] done %s — total %.1f min\n",
            format(END_T, "%Y-%m-%d %H:%M:%S"),
            as.numeric(difftime(END_T, START_T, units = "mins"))))
