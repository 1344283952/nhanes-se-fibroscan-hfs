# Re-diagnose 006 BKMR 4-metal (drop Mn) with bkmrhat::kmbayes_diagnose
# Vehtari 2021 rank-normalized split-R̂ at multiple burn-in proportions
# Plus PIP chain-half consistency check (Research 8 substantive-stability argument)
#
# Author: Claude Code automation
# Date: 2026-05-23
# Trigger: Research 8 finding — classical Gelman-Rubin (coda::gelman.diag) gave r_Se = 1.13
#          Hypothesis: Vehtari rank-normalized split R̂ + aggressive burn-in may drop r_Se ≤ 1.10
#          without re-running BKMR.
suppressPackageStartupMessages({
  library(bkmr)
  library(bkmrhat)
  library(coda)
})

# Optional rstan for classic Rhat (only if available)
HAS_RSTAN <- requireNamespace("rstan", quietly = TRUE)

CHECKPOINT <- "data/processed/bkmr_4metal_drop_Mn_006_checkpoint.rds"
OUT_DIR    <- "output/tables"
N_ITER     <- 10000
METAL_LABS <- c("Pb", "Cd", "Hg", "Se")

stopifnot(file.exists(CHECKPOINT))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Loading checkpoint ===\n")
cat("File:", CHECKPOINT, "\n")
cp <- readRDS(CHECKPOINT)
stopifnot(length(cp$accumulated_fits) == 2)

cat("done_iter_per_chain:", cp$done_iter_per_chain, "\n")
cat("target_iter:", cp$target_iter, "\n")
cat("n chains:", length(cp$accumulated_fits), "\n")

fitlist <- cp$accumulated_fits
# Apply the bkmrfit_list class that bkmrhat expects
class(fitlist) <- c("bkmrfit.list", class(fitlist))

# ---------------------------------------------------------------
# Helper: manual rank-normalized split-R̂ if kmbayes_diagnose fails
# ---------------------------------------------------------------
manual_rhat_table <- function(fitlist, burnin_iter, n_iter = N_ITER) {
  sel_post <- seq(burnin_iter + 1, n_iter)
  rows <- list()
  # r_m (kernel widths) — one per metal
  for (j in seq_along(METAL_LABS)) {
    chains <- lapply(fitlist, function(f) f$r[sel_post, j])
    M <- do.call(cbind, chains)
    rhat_classic_split <- if (HAS_RSTAN) {
      tryCatch(rstan::Rhat(M), error = function(e) NA_real_)
    } else NA_real_
    # ESS
    ess_bulk <- if (HAS_RSTAN) {
      tryCatch(rstan::ess_bulk(M), error = function(e) NA_real_)
    } else NA_real_
    ess_tail <- if (HAS_RSTAN) {
      tryCatch(rstan::ess_tail(M), error = function(e) NA_real_)
    } else NA_real_
    rows[[paste0("r_", METAL_LABS[j])]] <- c(
      rhat_classic_split = rhat_classic_split,
      ess_bulk = ess_bulk,
      ess_tail = ess_tail
    )
  }
  # lambda + sigsq.eps (scalar per iter)
  for (par in c("lambda", "sigsq.eps")) {
    chains <- lapply(fitlist, function(f) f[[par]][sel_post])
    M <- do.call(cbind, chains)
    rhat_classic_split <- if (HAS_RSTAN) {
      tryCatch(rstan::Rhat(M), error = function(e) NA_real_)
    } else NA_real_
    ess_bulk <- if (HAS_RSTAN) {
      tryCatch(rstan::ess_bulk(M), error = function(e) NA_real_)
    } else NA_real_
    ess_tail <- if (HAS_RSTAN) {
      tryCatch(rstan::ess_tail(M), error = function(e) NA_real_)
    } else NA_real_
    rows[[par]] <- c(
      rhat_classic_split = rhat_classic_split,
      ess_bulk = ess_bulk,
      ess_tail = ess_tail
    )
  }
  out <- do.call(rbind, rows)
  out <- as.data.frame(out)
  out$param <- rownames(out)
  out[, c("param", "rhat_classic_split", "ess_bulk", "ess_tail")]
}

# ---------------------------------------------------------------
# Sweep burn-ins: 50%, 60%, 70%
# ---------------------------------------------------------------
all_results <- list()
for (burn_pct in c(0.50, 0.60, 0.70)) {
  burnin_iter <- floor(N_ITER * burn_pct)
  cat(sprintf("\n=== burnin = %d  (%.0f%% of %d) ===\n",
              burnin_iter, burn_pct * 100, N_ITER))

  # Try bkmrhat::kmbayes_diagnose first
  diag_obj <- tryCatch({
    bkmrhat::kmbayes_diagnose(fitlist, warmup = burnin_iter)
  }, error = function(e) {
    cat("  kmbayes_diagnose error:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(diag_obj)) {
    cat("  kmbayes_diagnose succeeded — capturing output:\n")
    # Capture printed simsummary as text (don't try as.data.frame — simsummary [<- breaks)
    diag_text <- capture.output(print(diag_obj))
    # Save raw text for inspection
    writeLines(diag_text, file.path(OUT_DIR,
      sprintf("bkmr_4metal_bkmrhat_diag_burnin%d.txt", burnin_iter)))
    cat("  (saved full diag text)\n")
    # Extract the Q5/Q50/Q95/Mean/SD/Rhat/Bulk_ESS/Tail_ESS table for r1-r4 + lambda + sigsq.eps
    keep_pat <- "^(r[1-4]|lambda|sigsq\\.eps)\\s"
    keep_lines <- diag_text[grepl(keep_pat, diag_text)]
    if (length(keep_lines) > 0) {
      cat("  Key parameter R̂ rows (rank-normalized split, Vehtari 2021):\n")
      # Header
      hdr_idx <- grep("^\\s+Q5\\s+Q50", diag_text)
      if (length(hdr_idx) > 0) cat("   ", diag_text[hdr_idx[1]], "\n")
      for (ln in keep_lines) cat("   ", ln, "\n")
      # Parse to data.frame manually — tokenize whitespace
      df_rows <- list()
      for (ln in keep_lines) {
        toks <- strsplit(trimws(ln), "\\s+")[[1]]
        # toks: param Q5 Q50 Q95 Mean SD Rhat Bulk_ESS Tail_ESS
        if (length(toks) >= 9) {
          df_rows[[length(df_rows) + 1]] <- data.frame(
            param = toks[1],
            Q5 = as.numeric(toks[2]),
            Q50 = as.numeric(toks[3]),
            Q95 = as.numeric(toks[4]),
            Mean = as.numeric(toks[5]),
            SD = as.numeric(toks[6]),
            Rhat = as.numeric(toks[7]),
            Bulk_ESS = as.numeric(toks[8]),
            Tail_ESS = as.numeric(toks[9]),
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(df_rows) > 0) {
        diag_df <- do.call(rbind, df_rows)
        # Rename r1-r4 → r_Pb, r_Cd, r_Hg, r_Se
        for (k in seq_along(METAL_LABS)) {
          diag_df$param[diag_df$param == paste0("r", k)] <- paste0("r_", METAL_LABS[k])
        }
        diag_df$burnin_pct  <- burn_pct
        diag_df$burnin_iter <- burnin_iter
        out_csv <- file.path(OUT_DIR,
          sprintf("bkmr_4metal_bkmrhat_diag_burnin%d.csv", burnin_iter))
        write.csv(diag_df, out_csv, row.names = FALSE)
        cat("  wrote:", out_csv, "\n")
        all_results[[as.character(burn_pct)]] <- diag_df
      }
    }
  }

  # Always also do manual rank-normalized split-R̂ via rstan::Rhat for cross-validation
  cat("\n  Manual split-Rhat table (rstan::Rhat, rank-normalized split):\n")
  mtab <- manual_rhat_table(fitlist, burnin_iter)
  print(mtab, row.names = FALSE)
  mtab$burnin_pct  <- burn_pct
  mtab$burnin_iter <- burnin_iter
  out_csv2 <- file.path(OUT_DIR,
    sprintf("bkmr_4metal_bkmrhat_manual_rhat_burnin%d.csv", burnin_iter))
  write.csv(mtab, out_csv2, row.names = FALSE)
  cat("  wrote:", out_csv2, "\n")
  if (is.null(all_results[[as.character(burn_pct)]])) {
    all_results[[as.character(burn_pct)]] <- mtab
  }
}

# ---------------------------------------------------------------
# PIP chain-half consistency check (Research 8 substantive stability)
# ---------------------------------------------------------------
cat("\n=== PIP chain-half consistency ===\n")
pip_tab <- list()
for (ci in seq_along(fitlist)) {
  fit <- fitlist[[ci]]
  burnin <- 5000  # 50% burn-in for PIP (standard)
  post_iter <- (burnin + 1):N_ITER
  half_pt <- floor(length(post_iter) / 2)
  idx1 <- post_iter[1:half_pt]
  idx2 <- post_iter[(half_pt + 1):length(post_iter)]
  if (!is.null(fit$delta)) {
    d1 <- fit$delta[idx1, , drop = FALSE]
    d2 <- fit$delta[idx2, , drop = FALSE]
    pip1 <- colMeans(d1)
    pip2 <- colMeans(d2)
    diff <- abs(pip1 - pip2)
    cat(sprintf("\nChain %d  (burnin=%d, n_post=%d, half=%d each)\n",
                ci, burnin, length(post_iter), half_pt))
    cat(sprintf("  PIP_half1: %s\n", paste0(METAL_LABS, "=", sprintf("%.3f", pip1), collapse = "  ")))
    cat(sprintf("  PIP_half2: %s\n", paste0(METAL_LABS, "=", sprintf("%.3f", pip2), collapse = "  ")))
    cat(sprintf("  |diff|:    %s\n", paste0(METAL_LABS, "=", sprintf("%.3f", diff), collapse = "  ")))
    pip_tab[[paste0("chain", ci)]] <- data.frame(
      chain = ci, metal = METAL_LABS,
      pip_half1 = pip1, pip_half2 = pip2, abs_diff = diff
    )
  } else {
    cat(sprintf("Chain %d  (no delta indicator stored; skip)\n", ci))
  }
}
if (length(pip_tab) > 0) {
  pip_df <- do.call(rbind, pip_tab)
  rownames(pip_df) <- NULL
  out_csv <- file.path(OUT_DIR, "bkmr_4metal_pip_chainhalf.csv")
  write.csv(pip_df, out_csv, row.names = FALSE)
  cat("\nwrote:", out_csv, "\n")
}

# ---------------------------------------------------------------
# Summary: best r_Se R̂ across burn-in sweeps
# ---------------------------------------------------------------
cat("\n=== SUMMARY ===\n")
cat("Comparing classic split-Rhat (rank-normalized via rstan::Rhat) at each burn-in:\n\n")
for (burn_pct in names(all_results)) {
  df <- all_results[[burn_pct]]
  cat(sprintf("--- burnin = %s (iter = %d) ---\n",
              burn_pct, floor(N_ITER * as.numeric(burn_pct))))
  # Detect column for Rhat value (bkmrhat: "Rhat"; manual: "rhat_classic_split")
  rhat_col <- intersect(c("Rhat", "rhat_classic_split", "rhat", "PSRF", "psrf",
                          "point_est", "Point.est."), names(df))
  if (length(rhat_col) > 0) {
    rcol <- rhat_col[1]
    show <- df[, c("param", rcol), drop = FALSE]
    show[[rcol]] <- round(as.numeric(show[[rcol]]), 4)
    print(show, row.names = FALSE)
  } else {
    print(df)
  }
  cat("\n")
}

cat("\n=== DONE ===\n")
cat("CSVs written to:", OUT_DIR, "\n")
cat("  - bkmr_4metal_bkmrhat_diag_burnin*.csv (if kmbayes_diagnose succeeded)\n")
cat("  - bkmr_4metal_bkmrhat_manual_rhat_burnin*.csv (rstan::Rhat manual)\n")
cat("  - bkmr_4metal_pip_chainhalf.csv\n")
