# ============================================
# 006 / _rerun_w2_part2.R — continuation of _rerun_w2.R after halt at 15
# Scripts not yet run: 09_hfs_predict, 11_sensitivity, _ipw_selection, _retrodesign
# Plus re-run 05 (DEFF df bug fix)
# ============================================

if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("_rerun_w2_part2.R must be executed from project root (projects/006_se_fibroscan_hfs/).")
}

set.seed(20260516)
options(survey.lonely.psu = "adjust")

t0 <- Sys.time()
cat(sprintf("\n[rerun_w2_part2] START %s\n\n", format(t0, "%Y-%m-%d %H:%M:%S")))

scripts_in_order <- c(
  "scripts/05_table1.R",                 # DEFF df bug fix
  "scripts/09_hfs_predict.R",            # HL + Harrell C + NRI/IDI + NFS + APRI
  "scripts/11_sensitivity.R",            # cov_pre primary
  "scripts/_ipw_selection_sens.R",       # X3
  "scripts/_retrodesign_smallcells.R"    # X5
)

for (.script_path in scripts_in_order) {
  if (!file.exists(.script_path)) {
    cat(sprintf("[rerun_w2_part2] SKIP missing %s\n", .script_path)); next
  }
  .ts <- Sys.time()
  cat(sprintf("[rerun_w2_part2] >>> %s\n", .script_path))
  flush.console()
  .res <- tryCatch({
    local(source(.script_path, echo = FALSE), envir = new.env(parent = globalenv()))
    list(ok = TRUE, err = NULL)
  }, error = function(e) list(ok = FALSE, err = conditionMessage(e)))
  .dt <- as.numeric(difftime(Sys.time(), .ts, units = "secs"))
  if (isTRUE(.res$ok)) {
    cat(sprintf("[rerun_w2_part2] OK   %-45s (%.1f s) @ %s\n\n",
                .script_path, .dt, format(Sys.time(), "%H:%M:%S")))
  } else {
    cat(sprintf("[rerun_w2_part2] FAIL %-45s (%.1f s)\n  ERROR: %s\n\n",
                .script_path, .dt, .res$err))
  }
  flush.console()
}

cat(sprintf("\n[rerun_w2_part2] DONE %s — total %.1f min\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
