# ============================================
# 006 / _rerun_w2.R — sequential re-run of modified scripts after Round 2 fixes
# Run from project root via:
#   Rscript scripts/_rerun_w2.R
# Generates fresh output/tables/*.csv + .RData with cov_pre primary,
# lonely-PSU adjust, DEFF, HL + Harrell C + NRI/IDI + NFS + APRI + multiscore + retrodesign.
# ============================================

if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("_rerun_w2.R must be executed from project root (projects/006_se_fibroscan_hfs/).")
}

set.seed(20260516)
options(survey.lonely.psu = "adjust")

t0 <- Sys.time()
cat(sprintf("\n[rerun_w2] START %s\n\n", format(t0, "%Y-%m-%d %H:%M:%S")))

scripts_in_order <- c(
  "scripts/04_survey_design.R",          # lonely PSU + design regenerate
  "scripts/05_table1.R",                  # DEFF + df + new design
  "scripts/12_ratio_analysis.R",          # cov_pre primary
  "scripts/14_evalue.R",                  # new OR → new E-values
  "scripts/10_subgroup_forest.R",         # cov_pre
  "scripts/15_primary_fdr.R",             # new: 12-test primary BH-FDR
  "scripts/09_hfs_predict.R",             # HL + Harrell C + NRI/IDI + NFS + APRI + multiscore
  "scripts/_ipw_selection_sens.R",        # X3 IPW selection sensitivity (S8)
  "scripts/_retrodesign_smallcells.R"     # X5 Type S/M
)

for (.rerun_script_path in scripts_in_order) {
  if (!file.exists(.rerun_script_path)) {
    cat(sprintf("[rerun_w2] SKIP missing %s\n", .rerun_script_path)); next
  }
  .rerun_ts <- Sys.time()
  cat(sprintf("[rerun_w2] >>> %s\n", .rerun_script_path))
  flush.console()
  # source() in a fresh local() environment to prevent the sourced script from
  # rebinding any of the loop variables (e.g., 12_ratio_analysis.R reassigns `s`
  # to summary(fit)$coefficients, which corrupted the original runner). The
  # sourced scripts save outputs to disk, so we don't need to preserve their
  # in-memory bindings between iterations.
  .rerun_res <- tryCatch({
    local(source(.rerun_script_path, echo = FALSE), envir = new.env(parent = globalenv()))
    list(ok = TRUE, err = NULL)
  }, error = function(e) list(ok = FALSE, err = conditionMessage(e)))
  .rerun_dt <- as.numeric(difftime(Sys.time(), .rerun_ts, units = "secs"))
  if (isTRUE(.rerun_res$ok)) {
    cat(sprintf("[rerun_w2] OK   %-45s (%.1f s) @ %s\n\n",
                .rerun_script_path, .rerun_dt,
                format(Sys.time(), "%H:%M:%S")))
  } else {
    cat(sprintf("[rerun_w2] FAIL %-45s (%.1f s)\n  ERROR: %s\n\n",
                .rerun_script_path, .rerun_dt, .rerun_res$err))
  }
  flush.console()
}

cat(sprintf("\n[rerun_w2] DONE %s — total %.1f min\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
