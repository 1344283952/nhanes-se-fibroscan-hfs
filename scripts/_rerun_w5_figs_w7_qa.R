# ============================================
# 006 / _rerun_w5_figs_w7_qa.R
# W5 figure re-render (after part2 data updates) + W7 QA pipeline reruns
# ============================================

if (!dir.exists("scripts") || !dir.exists("data")) {
  stop("_rerun_w5_figs_w7_qa.R must be executed from project root.")
}

set.seed(20260516)
options(survey.lonely.psu = "adjust")

t0 <- Sys.time()
cat(sprintf("\n[w5_w7] START %s\n\n", format(t0, "%Y-%m-%d %H:%M:%S")))

scripts_seq <- c(
  "scripts/09_dag.R",                   # Fig 2 DAG (new 3-color palette + spread positions)
  "scripts/14_consort.R",                # Fig 1 CONSORT (numbers from updated flow_counts)
  "scripts/_render_fig2_dag.R",          # alternative Fig 2 DAG render via figures_helpers
  "scripts/_render_figures.R"            # Fig 1/3/4/5/6/7/8 batch
)

for (.sp in scripts_seq) {
  if (!file.exists(.sp)) { cat(sprintf("[w5_w7] SKIP %s\n", .sp)); next }
  .ts <- Sys.time()
  cat(sprintf("[w5_w7] >>> %s\n", .sp)); flush.console()
  .res <- tryCatch({
    local(source(.sp, echo = FALSE), envir = new.env(parent = globalenv()))
    list(ok = TRUE, err = NULL)
  }, error = function(e) list(ok = FALSE, err = conditionMessage(e)))
  .dt <- as.numeric(difftime(Sys.time(), .ts, units = "secs"))
  if (isTRUE(.res$ok)) {
    cat(sprintf("[w5_w7] OK   %-45s (%.1f s) @ %s\n\n",
                .sp, .dt, format(Sys.time(), "%H:%M:%S")))
  } else {
    cat(sprintf("[w5_w7] FAIL %-45s (%.1f s)\n  ERROR: %s\n\n",
                .sp, .dt, .res$err))
  }
  flush.console()
}

# W7 QA: run consistency check (extended; ≥40 numbers)
cat("\n[w5_w7] >>> templates/_consistency_check.R\n"); flush.console()
.ts <- Sys.time()
.res <- tryCatch({
  local(source("../../templates/_consistency_check.R", echo = FALSE),
        envir = new.env(parent = globalenv()))
  list(ok = TRUE, err = NULL)
}, error = function(e) list(ok = FALSE, err = conditionMessage(e)))
.dt <- as.numeric(difftime(Sys.time(), .ts, units = "secs"))
if (isTRUE(.res$ok)) {
  cat(sprintf("[w5_w7] OK   _consistency_check.R (%.1f s)\n\n", .dt))
} else {
  cat(sprintf("[w5_w7] FAIL _consistency_check.R (%.1f s)\n  ERROR: %s\n\n", .dt, .res$err))
}

cat(sprintf("\n[w5_w7] DONE %s — total %.1f min\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
