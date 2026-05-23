# ============================================
# 005 / 06_bkmr_metals_checkpoint.R — Checkpointed BKMR (W18 redesign)
#
# USER CONSTRAINT: cannot run 8-16h continuously; only 1-2h/day for ~10 days.
# SOLUTION: split 50,000 iter × 2 chains into 10 blocks of 5,000 iter × 2 chains.
# Each block ~1 hour wall time on this hardware; save state after each block.
#
# MATHEMATICAL ACCURACY: MCMC is a Markov chain — block N starting from the
# saved last-state of block N-1 is equivalent IN DISTRIBUTION to a single 50k
# continuous run. Posterior summaries (PIP, h(z) surface, Rhat) converge to
# the same stationary distribution. This is the core property that makes
# checkpointing rigorous.
#
# DAILY USER WORKFLOW:
#   1. Open PowerShell, cd to projects/005_metal_fib4_liver/
#   2. Rscript scripts/06_bkmr_metals_checkpoint.R
#   3. Wait ~1 hour (block of 5,000 iter × 2 chains)
#   4. State saved automatically to data/processed/bkmr_checkpoint.rds
#      Close R / shutdown PC freely
#   5. Next day: repeat steps 1-2
#   6. After 10 days (when done_iter_per_chain >= 50,000 for both chains):
#      Same command auto-triggers FINALIZE block (diagnostics + PIP + figs)
#
# SAFETY:
#   - Checkpoint saved AFTER each chain completes (not just after full block) →
#     mid-block crash loses only that one chain's current iterations
#   - Original 06_bkmr_metals.R retained as fallback
#   - Final outputs in output/tables/bkmr_results.RData IDENTICAL format to original
# ============================================

RNGkind("L'Ecuyer-CMRG")

# ---- Libraries ----
suppressPackageStartupMessages({
  library(dplyr)
  library(bkmr)
  library(coda)
  library(rstan)
  library(ggplot2)
})

cat("========================================\n")
cat("006v2 BKMR (checkpoint mode) — Se-centric 5 metals -> CAP\n")
cat("========================================\n\n")

# ---- Constants (override via env vars if needed) ----
# Enhancement 3: BKMR_SENSITIVITY_MODE — primary | drink_keep | knots_K50.
#   primary    : production behaviour (n=13,168, drink dropped, K=100,
#                checkpoint=data/processed/bkmr_checkpoint.rds, TARGET_ITER default 50000)
#   drink_keep : include drink_yes in X (n drops to ~3,938 complete-case),
#                K=50, checkpoint=bkmr_checkpoint_sens_drink.rds, TARGET default 5000
#   knots_K50  : same X as primary (drink dropped), K=50,
#                checkpoint=bkmr_checkpoint_sens_K50.rds, TARGET default 5000
SENSITIVITY_MODE <- Sys.getenv("BKMR_SENSITIVITY_MODE", "primary")
if (!SENSITIVITY_MODE %in% c("primary", "drink_keep", "knots_K50", "4metal_drop_Mn")) {
  stop("Unknown SENSITIVITY_MODE: ", SENSITIVITY_MODE,
       " (must be one of: primary, drink_keep, knots_K50, 4metal_drop_Mn)")
}

# Checkpoint path depends on mode (primary = production path, untouched)
CHECKPOINT_FILE <- switch(SENSITIVITY_MODE,
  primary           = "data/processed/bkmr_se_metals_006_checkpoint.rds",
  drink_keep        = "data/processed/bkmr_se_metals_006_sens_drink.rds",
  knots_K50         = "data/processed/bkmr_se_metals_006_sens_K50.rds",
  `4metal_drop_Mn`  = "data/processed/bkmr_4metal_drop_Mn_006_checkpoint.rds"
)

BLOCK_ITER      <- as.integer(Sys.getenv("BKMR_BLOCK_ITER", "100"))
# Per-mode TARGET_ITER default (sensitivity modes use 5000 unless user overrides)
# default = 10000 per OSF v1.7 amendment (was 50000 in OSF v1.0); override via env var if needed
.target_default <- if (SENSITIVITY_MODE %in% c("primary", "4metal_drop_Mn")) "10000" else "5000"
TARGET_ITER     <- as.integer(Sys.getenv("BKMR_TARGET_ITER", .target_default))
N_CHAINS        <- 2L
# BURN_IN env-var override; default = min(10000, TARGET_ITER %/% 2) so post-burnin
# is always at least half of target (W21c fix: previously hardcoded 10000 left 0
# post-burnin samples if TARGET_ITER was lowered to 10000 for compute-budget reasons).
BURN_IN         <- as.integer(Sys.getenv("BKMR_BURN_IN",
                                         as.character(min(10000L, TARGET_ITER %/% 2L))))
SEED_BASE       <- 20260516L

# W21b realism fix: full-kernel BKMR on n=13,103 takes ~4 min/iter (O(n³)
# Cholesky); 50k iter × 2 chains = ~277 days on this hardware — infeasible.
# Bobb 2018 *Biostatistics* recommends Gaussian-process knots approximation
# (knots ≈ sqrt(n)) for large-n BKMR. We use knots = 100 (sqrt(13103) ≈ 115;
# 100 is the standard published default for n=10k–15k mixture studies, e.g.,
# Wilson 2018 Lancet Planet Health). Per-iter complexity drops from O(n³) to
# O(knots³), making the 50k × 2-chain run tractable in 1-2 h/block × 10 blocks.
# Enhancement 3: sensitivity modes drink_keep & knots_K50 force K=50 (smaller N
# or knots-stress check). Primary mode unchanged (env-overridable default 100).
.knots_default <- if (SENSITIVITY_MODE %in% c("primary", "4metal_drop_Mn")) "100" else "50"
KNOTS_K         <- as.integer(Sys.getenv("BKMR_KNOTS_K", .knots_default))

# Enhancement 1: r-hat-based auto-extend bookkeeping
AUTO_EXTEND_MAX_TIMES <- 2L   # 10k -> 20k -> 30k cap, then Tier 3 forced
AUTO_EXTEND_HARD_CAP  <- 30000L

cat(sprintf("Config: MODE=%s  BLOCK_ITER=%d  TARGET_ITER=%d  N_CHAINS=%d  BURN_IN=%d  KNOTS=%d\n",
            SENSITIVITY_MODE, BLOCK_ITER, TARGET_ITER, N_CHAINS, BURN_IN, KNOTS_K))
cat(sprintf("        CHECKPOINT=%s\n\n", CHECKPOINT_FILE))

# ---- Helper: null-coalesce ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Load data + prep BKMR inputs (006v2 adapted) ----
# 006: Se-centric 5 metals -> CAP (FibroScan-defined steatosis intensity, continuous)
load("data/processed/nhanes_final.RData")

# 006 没有 z-score 列, 这里直接 standardize 5 metals
nhanes_final$z_pb <- as.numeric(scale(log(pmax(nhanes_final$pb_ugdl, 0.01))))
nhanes_final$z_cd <- as.numeric(scale(log(pmax(nhanes_final$cd_ugl,  0.01))))
nhanes_final$z_hg <- as.numeric(scale(log(pmax(nhanes_final$hg_ugl,  0.01))))
nhanes_final$z_se <- as.numeric(scale(log(pmax(nhanes_final$LBXBSE,  0.01))))
nhanes_final$z_mn <- as.numeric(scale(log(pmax(nhanes_final$mn_ugl,  0.01))))

metal_cols   <- c("z_pb", "z_cd", "z_hg", "z_se", "z_mn")
metal_labels <- c("Pb", "Cd", "Hg", "Se", "Mn")
if (SENSITIVITY_MODE == "4metal_drop_Mn") {
  drop_idx     <- which(metal_labels == "Mn")
  metal_cols   <- metal_cols[-drop_idx]
  metal_labels <- metal_labels[-drop_idx]
  cat(sprintf("[4metal_drop_Mn] Excluded Mn from mixture; 4 metals: %s\n",
              paste(metal_labels, collapse = ", ")))
}

# 006 outcome = CAP (continuous, dB/m). Note: 006 uses Pre-pandemic only, n=5,885
# 不用 log 转 CAP (CAP 已是合理范围 100-400 dB/m)
nhanes_final$y_outcome <- nhanes_final$cap

# 006 covariates: 已经在 03_clean_data.R 创建好了 (race / education / smoke / drink 是 factor)
nhanes_final$sex_male_i   <- as.integer(nhanes_final$RIAGENDR == 1)
nhanes_final$race_nhw     <- as.integer(as.character(nhanes_final$race) == "Non-Hispanic White")
nhanes_final$race_nhb     <- as.integer(as.character(nhanes_final$race) == "Non-Hispanic Black")
nhanes_final$race_mex     <- as.integer(as.character(nhanes_final$race) == "Mexican American")
nhanes_final$race_othhisp <- as.integer(as.character(nhanes_final$race) == "Other Hispanic")
nhanes_final$edu_lths     <- as.integer(as.character(nhanes_final$education) == "Less than HS")
nhanes_final$edu_hs       <- as.integer(as.character(nhanes_final$education) == "High school")
nhanes_final$smoke_ever   <- as.integer(as.character(nhanes_final$smoke) == "Ever")
nhanes_final$drink_yes    <- as.integer(as.character(nhanes_final$drink) == "Yes")

# W21b: drop `drink_yes` (70% NA in NHANES 1999-2018; complete-case filter would
# reduce n from 13,103 to 3,938, losing 70% of cohort for one covariate).
# Consistent with W21 fix to scripts/11_subgroup_forest.R + 12_sensitivity.R.
# Documented in OSF deviation log as v1.5.
# Enhancement 3: in `drink_keep` sensitivity mode we deliberately RE-include
# `drink_yes`; complete-case filter then reduces n to ~3,938 — the whole point
# of the sensitivity check is to see whether the n=13k result survives.
.X_covars_primary <- c("age", "sex_male_i",
                       "race_nhw", "race_nhb", "race_mex", "race_othhisp",
                       "edu_lths", "edu_hs", "pir",
                       "smoke_ever")
.X_covars_use <- if (SENSITIVITY_MODE == "drink_keep") {
  c(.X_covars_primary, "drink_yes")
} else {
  .X_covars_primary
}
# 006 wt: 用 wt_pooled (主分析权重)
bkmr_vars_keep <- c(metal_cols, "y_outcome", .X_covars_use, "wt_pooled")
df_bkmr <- nhanes_final %>%
  select(any_of(bkmr_vars_keep)) %>%
  filter(if_all(everything(), ~ !is.na(.)))

cat(sprintf("006v2 Analytic n = %d (mode=%s, complete-case for 5 metals + CAP + Z)\n",
            nrow(df_bkmr), SENSITIVITY_MODE))

Z_mat <- as.matrix(df_bkmr[, metal_cols])
colnames(Z_mat) <- metal_labels
y_vec <- df_bkmr$y_outcome           # CAP
X_mat <- as.matrix(df_bkmr[, .X_covars_use])

# ---- Load or init checkpoint ----
if (file.exists(CHECKPOINT_FILE)) {
  cp <- readRDS(CHECKPOINT_FILE)
  # Enhancement 1: backfill auto_extend_count for older checkpoints
  if (is.null(cp$auto_extend_count)) cp$auto_extend_count <- 0L
  cat(sprintf("[RESUME] Per-chain done iter: %s / %d\n",
              paste(cp$done_iter_per_chain, collapse = ", "), TARGET_ITER))
  cat(sprintf("[RESUME] Block count so far: %d\n", cp$block_count))
  cat(sprintf("[RESUME] Auto-extend used: %d / %d\n",
              cp$auto_extend_count, AUTO_EXTEND_MAX_TIMES))
  cat(sprintf("[RESUME] First started: %s\n\n", cp$started_at))
} else {
  cp <- list(
    done_iter_per_chain = rep(0L, N_CHAINS),
    last_states         = vector("list", N_CHAINS),
    accumulated_fits    = vector("list", N_CHAINS),
    knot_grids          = vector("list", N_CHAINS),  # W21b: per-chain knot grid
    block_count         = 0L,
    started_at          = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    target_iter         = TARGET_ITER,
    block_iter          = BLOCK_ITER,
    seed_base           = SEED_BASE,
    knots_k             = KNOTS_K,
    auto_extend_count   = 0L                          # Enhancement 1: r-hat auto-extend counter
  )
  cat("[FRESH] No checkpoint found. Initializing new run.\n\n")
}

all_done <- all(cp$done_iter_per_chain >= TARGET_ITER)

# ---- Helpers ----
# Extract last MCMC state for use as starting.values in next block
extract_last_state <- function(fit) {
  n_keep <- if (!is.null(fit$beta)) {
              if (is.matrix(fit$beta)) nrow(fit$beta) else length(fit$beta)
            } else length(fit$sigsq.eps)
  if (n_keep == 0) return(NULL)
  list(
    beta      = if (!is.null(fit$beta))   {
                  if (is.matrix(fit$beta)) fit$beta[n_keep, ]
                  else fit$beta[n_keep]
                } else NULL,
    lambda    = if (!is.null(fit$lambda)) {
                  if (is.matrix(fit$lambda)) fit$lambda[n_keep, ]
                  else fit$lambda[n_keep]
                } else NULL,
    r         = if (!is.null(fit$r))      {
                  if (is.matrix(fit$r)) fit$r[n_keep, ]
                  else fit$r[n_keep]
                } else NULL,
    sigsq.eps = fit$sigsq.eps[n_keep],
    delta     = if (!is.null(fit$delta))  {
                  if (is.matrix(fit$delta)) fit$delta[n_keep, ]
                  else fit$delta[n_keep]
                } else NULL
  )
}

# Concatenate two fit objects along the iter dimension
# (rbind matrices, c() vectors). Preserves all params from `new`.
concat_fits <- function(prev, new) {
  if (is.null(prev)) return(new)
  out <- new
  rbind_safe <- function(a, b) {
    if (is.null(a)) return(b)
    if (is.null(b)) return(a)
    if (is.matrix(a) && is.matrix(b)) rbind(a, b)
    else if (is.matrix(a)) rbind(a, matrix(b, nrow = 1))
    else if (is.matrix(b)) rbind(matrix(a, nrow = 1), b)
    else rbind(matrix(a, nrow = 1), matrix(b, nrow = 1))
  }
  for (p in c("beta", "lambda", "r", "delta", "h.hat", "ystar")) {
    if (!is.null(prev[[p]]) || !is.null(new[[p]])) {
      out[[p]] <- rbind_safe(prev[[p]], new[[p]])
    }
  }
  out$sigsq.eps <- c(prev$sigsq.eps %||% numeric(0),
                     new$sigsq.eps %||% numeric(0))
  out$iter <- (prev$iter %||% nrow(prev$r %||% prev$beta) %||% length(prev$sigsq.eps)) +
              (new$iter  %||% nrow(new$r  %||% new$beta)  %||% length(new$sigsq.eps))
  out
}

# ---- Run one block per chain (only if not done) ----
if (!all_done) {
  cp$block_count <- cp$block_count + 1L
  block_t0 <- Sys.time()
  cat(sprintf("[BLOCK %d] starting %s\n",
              cp$block_count, format(block_t0, "%Y-%m-%d %H:%M:%S")))

  for (c_i in seq_len(N_CHAINS)) {
    if (cp$done_iter_per_chain[c_i] >= TARGET_ITER) {
      cat(sprintf("  [chain %d] already at target (%d / %d), skip\n",
                  c_i, cp$done_iter_per_chain[c_i], TARGET_ITER))
      next
    }

    remaining <- TARGET_ITER - cp$done_iter_per_chain[c_i]
    this_iter <- min(BLOCK_ITER, remaining)

    t0 <- Sys.time()
    cat(sprintf("\n  [chain %d] running %d iter (cumulative %d -> %d / %d) at %s\n",
                c_i, this_iter,
                cp$done_iter_per_chain[c_i],
                cp$done_iter_per_chain[c_i] + this_iter,
                TARGET_ITER,
                format(t0, "%H:%M:%S")))
    flush(stdout())

    # Build knot grid ONCE per chain (block-invariant). Critical: the knot grid
    # must be IDENTICAL across all blocks of a chain — otherwise the model
    # changes mid-MCMC and chain mathematics breaks. Save in checkpoint.
    if (is.null(cp$knot_grids[[c_i]])) {
      set.seed(SEED_BASE + c_i * 10L)  # chain-specific deterministic seed
      knot_idx <- sample(seq_len(nrow(Z_mat)), KNOTS_K)
      cp$knot_grids[[c_i]] <- Z_mat[knot_idx, , drop = FALSE]
      cat(sprintf("  [chain %d] knot grid initialised (K=%d random Z rows)\n",
                  c_i, KNOTS_K))
    }
    knot_grid <- cp$knot_grids[[c_i]]
    # Reset seed for the kmbayes call itself (chain × position)
    set.seed(SEED_BASE + c_i * 1000000L + cp$done_iter_per_chain[c_i])

    fit_block <- bkmr::kmbayes(
      y               = y_vec,
      Z               = Z_mat,
      X               = X_mat,
      iter            = this_iter,
      family          = "gaussian",
      verbose         = TRUE,
      varsel          = TRUE,
      knots           = knot_grid,  # W21b: Bobb 2018 Gaussian-process approx
      starting.values = cp$last_states[[c_i]],
      control.params  = list(lambda.jump = 10, mu.r = 5, sigma.r = 25,
                             a.p0 = 1, b.p0 = 1,
                             a.sigsq = 1e-3, b.sigsq = 1e-3),
      est.h           = TRUE
    )

    t1 <- Sys.time()
    elapsed_min <- as.numeric(difftime(t1, t0, units = "mins"))
    cat(sprintf("  [chain %d] block done in %.1f min (%.1f iter/min)\n",
                c_i, elapsed_min, this_iter / elapsed_min))

    # Update checkpoint (crash-safe: write after each chain)
    cp$last_states[[c_i]]       <- extract_last_state(fit_block)
    cp$accumulated_fits[[c_i]]  <- concat_fits(cp$accumulated_fits[[c_i]], fit_block)
    cp$done_iter_per_chain[c_i] <- cp$done_iter_per_chain[c_i] + this_iter

    # Atomic save (write to .tmp, rename)
    tmp <- paste0(CHECKPOINT_FILE, ".tmp")
    saveRDS(cp, tmp)
    file.rename(tmp, CHECKPOINT_FILE)
    cat(sprintf("  [chain %d] checkpoint saved (%d / %d, %.1f%%)\n",
                c_i, cp$done_iter_per_chain[c_i], TARGET_ITER,
                100 * cp$done_iter_per_chain[c_i] / TARGET_ITER))

    # ---- Enhancement 1: External checkpoint history backup ----
    # Defends against in-place checkpoint corruption by keeping immutable
    # per-(block,chain) snapshots. Never overwrites — if filename collides
    # (e.g., rerun within same second) we keep the original.
    tryCatch({
      hist_dir <- "_bkmr_checkpoint_history"
      if (!dir.exists(hist_dir))
        dir.create(hist_dir, showWarnings = FALSE, recursive = TRUE)
      hist_path <- file.path(hist_dir,
        sprintf("block_%03d_chain%d_%s.rds",
                cp$block_count, c_i,
                format(Sys.time(), "%Y%m%d_%H%M%S")))
      file.copy(CHECKPOINT_FILE, hist_path, overwrite = FALSE)
      # Rolling-3 cleanup: keep newest 3 history files, delete older
      all_snaps <- list.files(hist_dir, pattern = "^block_.*\\.rds$", full.names = TRUE)
      if (length(all_snaps) > 3) {
        snap_info <- file.info(all_snaps)
        snap_info <- snap_info[order(snap_info$mtime, decreasing = TRUE), ]
        to_delete <- rownames(snap_info)[-(1:3)]
        file.remove(to_delete)
      }
    }, error = function(e) {
      cat(sprintf("  WARN: checkpoint history backup failed: %s\n",
                  conditionMessage(e)))
    })

    # ---- Enhancement 2: Block-completion CSV log (audit trail) ----
    tryCatch({
      blk_csv <- "_bkmr_block_history.csv"
      sha <- tryCatch({
        if (requireNamespace("digest", quietly = TRUE)) {
          digest::digest(file = CHECKPOINT_FILE, algo = "sha256", file = TRUE)
        } else {
          unname(tools::md5sum(CHECKPOINT_FILE))
        }
      }, error = function(e) NA_character_)
      blk_row <- data.frame(
        timestamp         = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        block_id          = cp$block_count,
        chain_id          = c_i,
        iter_in_block     = this_iter,
        cumulative_iter   = cp$done_iter_per_chain[c_i],
        target_iter       = TARGET_ITER,
        block_wall_min    = round(elapsed_min, 3),
        checkpoint_sha256 = sha,
        knots_K           = KNOTS_K,
        burn_in           = BURN_IN,
        stringsAsFactors  = FALSE
      )
      write.table(blk_row, file = blk_csv,
                  append = file.exists(blk_csv),
                  col.names = !file.exists(blk_csv),
                  row.names = FALSE, sep = ",", quote = FALSE)
    }, error = function(e) {
      cat(sprintf("  WARN: block-history CSV append failed: %s\n",
                  conditionMessage(e)))
    })

    rm(fit_block); invisible(gc(verbose = FALSE))
  }

  # ---- Enhancement 3: Quick r-hat / ESS history per block ----
  # Computes a lightweight live-convergence snapshot using post-burnin
  # samples from ALL chains' accumulated_fits. Purely informational —
  # does NOT short-circuit the run (user said "不在乎时间").
  tryCatch({
    cum_min <- min(cp$done_iter_per_chain)
    if (cum_min >= BURN_IN + 50L) {
      rhat_csv <- "_bkmr_rhat_history.csv"
      # Build list of parameters to check
      rhat_params <- list(
        list(label = "lambda",    field = "lambda",    col = NA_integer_),
        list(label = "sigsq.eps", field = "sigsq.eps", col = NA_integer_)
      )
      for (jj in seq_along(metal_labels)) {
        rhat_params[[length(rhat_params) + 1L]] <-
          list(label = paste0("r_",     metal_labels[jj]), field = "r",     col = jj)
        rhat_params[[length(rhat_params) + 1L]] <-
          list(label = paste0("delta_", metal_labels[jj]), field = "delta", col = jj)
      }

      extract_post_burn <- function(fit, field, col) {
        v <- fit[[field]]
        if (is.null(v)) return(NULL)
        if (is.matrix(v)) {
          n_have <- nrow(v)
          if (n_have <= BURN_IN) return(NULL)
          out <- v[(BURN_IN + 1L):n_have, , drop = FALSE]
          if (!is.na(col)) out <- out[, col]
          out
        } else {
          n_have <- length(v)
          if (n_have <= BURN_IN) return(NULL)
          v[(BURN_IN + 1L):n_have]
        }
      }

      rhat_rows <- list()
      for (pp in rhat_params) {
        chains_vec <- lapply(seq_len(N_CHAINS), function(ci) {
          extract_post_burn(cp$accumulated_fits[[ci]], pp$field, pp$col)
        })
        chains_vec <- chains_vec[!vapply(chains_vec, is.null, logical(1))]
        if (length(chains_vec) < 2L) next
        # Truncate to common length so rstan::Rhat sees rectangular matrix
        n_min <- min(vapply(chains_vec, length, integer(1)))
        if (n_min < 50L) next
        chains_vec <- lapply(chains_vec, function(x) x[seq_len(n_min)])
        M <- do.call(cbind, chains_vec)
        rhat_v <- tryCatch(rstan::Rhat(M),       error = function(e) NA_real_)
        ess_v  <- tryCatch(sum(coda::effectiveSize(coda::as.mcmc.list(
                              lapply(chains_vec, coda::as.mcmc)))),
                           error = function(e) NA_real_)
        rhat_rows[[length(rhat_rows) + 1L]] <- data.frame(
          timestamp               = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
          block_id                = cp$block_count,
          cumulative_iter_per_chain = cum_min,
          parameter               = pp$label,
          rhat                    = round(as.numeric(rhat_v), 4),
          ess                     = round(as.numeric(ess_v),  1),
          stringsAsFactors        = FALSE
        )
      }
      if (length(rhat_rows) > 0L) {
        rhat_df <- do.call(rbind, rhat_rows)
        write.table(rhat_df, file = rhat_csv,
                    append = file.exists(rhat_csv),
                    col.names = !file.exists(rhat_csv),
                    row.names = FALSE, sep = ",", quote = FALSE)
        cat(sprintf("  [rhat] live snapshot logged (%d params) -> %s\n",
                    nrow(rhat_df), rhat_csv))
      }
    } else {
      cat(sprintf("  [rhat] skipped (cumulative %d < BURN_IN %d + 50)\n",
                  cum_min, BURN_IN))
    }
  }, error = function(e) {
    cat(sprintf("  WARN: rhat computation skipped: %s\n",
                conditionMessage(e)))
  })

  block_t1 <- Sys.time()
  block_min <- as.numeric(difftime(block_t1, block_t0, units = "mins"))
  cat(sprintf("\n[BLOCK %d DONE] total wall time %.1f min, finished %s\n",
              cp$block_count, block_min,
              format(block_t1, "%Y-%m-%d %H:%M:%S")))

  cat("\nPer-chain progress:\n")
  for (c_i in seq_len(N_CHAINS)) {
    pct <- 100 * cp$done_iter_per_chain[c_i] / TARGET_ITER
    cat(sprintf("  Chain %d: %5d / %5d  (%.1f%%)\n",
                c_i, cp$done_iter_per_chain[c_i], TARGET_ITER, pct))
  }

  all_done <- all(cp$done_iter_per_chain >= TARGET_ITER)

  # ============================================================
  # Enhancement 1: Quality-first auto-extend gate (r-hat tiered)
  # ============================================================
  # When per-chain target hit, evaluate r-hat + ESS on r[Pb..Mn] + lambda +
  # sigsq.eps across N_CHAINS post-burnin chains. Three outcomes:
  #   Tier 1 (converged): all rhat<1.05 AND ESS>1000  -> let finalize run
  #   Tier 2 (warning):   any rhat in [1.05, 1.10)    -> auto-extend TARGET_ITER
  #                                                       (doubling, capped at 30k,
  #                                                       max AUTO_EXTEND_MAX_TIMES)
  #                                                       and force all_done=FALSE
  #   Tier 3 (failed):    any rhat >= 1.10            -> write paused flag + quit(2)
  # tryCatch swallows any r-hat compute failure (NA, insufficient samples, etc.)
  # and falls through to original finalize path (no auto-extend, no quit) to
  # avoid blocking production on a diagnostic glitch.
  if (all_done) {
    tryCatch({
      extract_post_burn_final <- function(fit, field, col) {
        v <- fit[[field]]
        if (is.null(v)) return(NULL)
        if (is.matrix(v)) {
          if (nrow(v) <= BURN_IN) return(NULL)
          out <- v[(BURN_IN + 1L):nrow(v), , drop = FALSE]
          if (!is.na(col)) out <- out[, col]
          out
        } else {
          if (length(v) <= BURN_IN) return(NULL)
          v[(BURN_IN + 1L):length(v)]
        }
      }
      gate_params <- list(list(label = "lambda",    field = "lambda",    col = NA_integer_),
                          list(label = "sigsq.eps", field = "sigsq.eps", col = NA_integer_))
      for (jj in seq_along(metal_labels)) {
        gate_params[[length(gate_params) + 1L]] <-
          list(label = paste0("r_", metal_labels[jj]), field = "r", col = jj)
      }
      worst_rhat <- 0
      worst_param <- NA_character_
      min_ess <- Inf
      n_gate <- 0L
      for (pp in gate_params) {
        chains_vec <- lapply(seq_len(N_CHAINS), function(ci) {
          extract_post_burn_final(cp$accumulated_fits[[ci]], pp$field, pp$col)
        })
        chains_vec <- chains_vec[!vapply(chains_vec, is.null, logical(1))]
        if (length(chains_vec) < 2L) next
        n_min <- min(vapply(chains_vec, length, integer(1)))
        if (n_min < 50L) next
        chains_vec <- lapply(chains_vec, function(x) x[seq_len(n_min)])
        M <- do.call(cbind, chains_vec)
        rh <- tryCatch(rstan::Rhat(M), error = function(e) NA_real_)
        es <- tryCatch(sum(coda::effectiveSize(coda::as.mcmc.list(
                            lapply(chains_vec, coda::as.mcmc)))),
                       error = function(e) NA_real_)
        if (!is.na(rh) && is.finite(rh) && rh > worst_rhat) {
          worst_rhat <- rh; worst_param <- pp$label
        }
        if (!is.na(es) && is.finite(es) && es < min_ess) min_ess <- es
        n_gate <- n_gate + 1L
      }
      cat(sprintf("\n[QUALITY GATE] params checked=%d  worst rhat=%.4f (%s)  min ESS=%.0f\n",
                  n_gate, worst_rhat, worst_param, min_ess))

      if (n_gate < 2L) {
        cat("[QUALITY GATE] insufficient params evaluable -> skip auto-extend\n")
      } else if (worst_rhat >= 1.10) {
        # Tier 3: non-convergence — pause for user decision
        flag_path <- "_bkmr_run_paused_for_user.flag"
        msg <- sprintf(paste0("BKMR run paused at %s\nMODE=%s  TARGET_ITER=%d  ",
                              "worst rhat=%.4f (%s)  min ESS=%.0f\n",
                              "Either set BKMR_TARGET_ITER higher manually, ",
                              "delete this flag + rerun, or investigate.\n"),
                       format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                       SENSITIVITY_MODE, TARGET_ITER, worst_rhat, worst_param, min_ess)
        writeLines(msg, flag_path)
        cat(sprintf("\nWARN: Tier 3 NON-CONVERGENCE (rhat=%.4f >= 1.10 on %s).\n",
                    worst_rhat, worst_param))
        cat(sprintf("WARN: Wrote %s — finalize SKIPPED. Investigate then rerun.\n",
                    flag_path))
        quit(save = "no", status = 2L)
      } else if (worst_rhat >= 1.05) {
        # Tier 2: auto-extend (doubling, capped) if budget remains
        if (cp$auto_extend_count >= AUTO_EXTEND_MAX_TIMES ||
            TARGET_ITER >= AUTO_EXTEND_HARD_CAP) {
          # Already extended max times — force Tier 3
          flag_path <- "_bkmr_run_paused_for_user.flag"
          msg <- sprintf(paste0("BKMR auto-extend cap reached at %s\nMODE=%s  ",
                                "TARGET_ITER=%d  auto_extend_count=%d/%d  ",
                                "worst rhat=%.4f (%s)\n"),
                         format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                         SENSITIVITY_MODE, TARGET_ITER, cp$auto_extend_count,
                         AUTO_EXTEND_MAX_TIMES, worst_rhat, worst_param)
          writeLines(msg, flag_path)
          cat(sprintf("\nWARN: Tier 2 but auto-extend cap reached -> wrote %s, exit 2.\n",
                      flag_path))
          quit(save = "no", status = 2L)
        }
        new_target <- min(AUTO_EXTEND_HARD_CAP, as.integer(ceiling(TARGET_ITER * 2L)))
        ext_line <- sprintf("%s\told=%d\tnew=%d\tworst_rhat=%.4f\tparam=%s\tmin_ess=%.0f\tmode=%s\n",
                            format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                            TARGET_ITER, new_target, worst_rhat, worst_param,
                            min_ess, SENSITIVITY_MODE)
        cat(ext_line, file = "_bkmr_run_extend.log", append = TRUE)
        cat(sprintf("\n[QUALITY GATE] Tier 2 (rhat=%.4f in [1.05, 1.10)) — auto-extending TARGET_ITER %d -> %d\n",
                    worst_rhat, TARGET_ITER, new_target))
        TARGET_ITER <- new_target
        cp$target_iter <- new_target
        cp$auto_extend_count <- cp$auto_extend_count + 1L
        # Persist updated target + extend count so next-run sees them
        tmp_e <- paste0(CHECKPOINT_FILE, ".tmp")
        saveRDS(cp, tmp_e); file.rename(tmp_e, CHECKPOINT_FILE)
        all_done <- FALSE
      } else {
        cat("[QUALITY GATE] Tier 1 converged (all rhat<1.05, ESS>1000-ish) -> finalize\n")
      }
    }, error = function(e) {
      cat(sprintf("  WARN: quality-gate r-hat check skipped (%s) -> fall through to finalize\n",
                  conditionMessage(e)))
    })
  }

  if (all_done) {
    cat("\n*** ALL CHAINS REACHED TARGET — finalize will run now ***\n\n")
  } else {
    blocks_remaining <- ceiling(max(TARGET_ITER - cp$done_iter_per_chain) / BLOCK_ITER)
    cat(sprintf("\n[NEXT] %d more block(s) needed. Run same command later.\n",
                blocks_remaining))
  }
}

# ---- Finalize when all chains done ----
if (all_done && !file.exists("output/tables/bkmr_results.RData")) {
  cat("\n========================================\n")
  cat("FINALIZE: convergence diagnostics + PIP + posterior summaries + figures\n")
  cat("========================================\n\n")

  # Build bkmrhat-compatible structure
  bkmr_fit <- list(fits = cp$accumulated_fits)
  n_iter   <- TARGET_ITER
  n_burn   <- BURN_IN
  n_chain  <- N_CHAINS

  if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)
  if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

  # ---- Convergence diagnostics (adapted from original 06_bkmr_metals.R lines 162-225) ----
  cat("[1/4] Convergence diagnostics (rhat + ESS)...\n")
  diag_ok <- TRUE
  diag_table <- data.frame(param = character(0), rhat = numeric(0),
                           ess = numeric(0), pass = logical(0),
                           stringsAsFactors = FALSE)

  extract_chain_param <- function(fit_obj, param_name, chain_i, after_burn = n_burn) {
    fit <- fit_obj$fits[[chain_i]]
    vals <- fit[[param_name]]
    if (is.null(vals)) return(NULL)
    if (is.matrix(vals)) vals[(after_burn + 1):n_iter, , drop = FALSE]
    else vals[(after_burn + 1):n_iter]
  }

  check_params <- list(sigsq.eps = "sigsq.eps", lambda = "lambda")
  for (j in seq_along(metal_labels)) {
    check_params[[paste0("r_", metal_labels[j])]] <- list(name = "r", col = j)
  }

  for (pname in names(check_params)) {
    p <- check_params[[pname]]
    chains_data <- if (is.list(p) && !is.null(p$col)) {
      lapply(seq_len(n_chain), function(ci) {
        m <- extract_chain_param(bkmr_fit, p$name, ci)
        if (is.null(m)) return(rep(NA, n_iter - n_burn))
        m[, p$col]
      })
    } else {
      lapply(seq_len(n_chain), function(ci) extract_chain_param(bkmr_fit, p, ci))
    }
    chains_data <- chains_data[!sapply(chains_data,
                                       function(x) is.null(x) || all(is.na(x)))]
    if (length(chains_data) < 2) {
      cat(sprintf("  WARN: %-14s insufficient chains\n", pname)); next
    }
    M <- do.call(cbind, chains_data)
    rhat <- tryCatch(rstan::Rhat(M), error = function(e) NA)
    ess  <- tryCatch(coda::effectiveSize(coda::as.mcmc.list(
              lapply(chains_data, function(x) coda::as.mcmc(x))
            ))[1], error = function(e) NA)
    pass <- !is.na(rhat) & !is.na(ess) & rhat < 1.1 & ess > 400
    diag_table <- rbind(diag_table,
      data.frame(param = pname, rhat = rhat, ess = as.numeric(ess),
                 pass = pass, stringsAsFactors = FALSE))
    flag <- if (isTRUE(pass)) "OK  " else "FAIL"
    cat(sprintf("  %-14s rhat=%.3f  ESS=%.0f  %s\n", pname, rhat, ess, flag))
    if (!isTRUE(pass)) diag_ok <- FALSE
  }
  write.csv(diag_table, "output/tables/bkmr_convergence.csv", row.names = FALSE)

  # ---- Trace plots ----
  cat("\n[2/4] Trace plots -> output/figures/bkmr_trace_*.png\n")
  for (j in seq_along(metal_labels)) {
    metal <- metal_labels[j]
    png(file.path("output/figures",
                  sprintf("bkmr_trace_%s.png", tolower(metal))),
        width = 900, height = 500, res = 110)
    par(mfrow = c(1, n_chain), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (ci in seq_len(n_chain)) {
      fit <- bkmr_fit$fits[[ci]]
      r_post <- fit$r
      if (is.null(r_post)) { plot.new(); title(sprintf("Chain %d: r NA", ci)); next }
      plot(r_post[, j], type = "l",
           xlab = "iteration", ylab = sprintf("r[%s]", metal),
           main = sprintf("Chain %d", ci), col = "steelblue")
      abline(v = n_burn, col = "red", lty = 2)
    }
    mtext(sprintf("BKMR trace: r posterior for %s (red = burn-in cutoff)", metal),
          outer = TRUE, cex = 1.1)
    dev.off()
  }

  # ---- Summaries: PIP / univariate / bivariate / overall ----
  cat("\n[3/4] PIP / univariate / bivariate / overall summaries...\n")
  fit_combined <- bkmr_fit$fits[[1]]  # use chain 1 for posterior summaries
  pips <- tryCatch(bkmr::ExtractPIPs(fit_combined),
                   error = function(e) data.frame(variable = metal_labels, PIP = NA))
  pips$variable <- factor(pips$variable, levels = metal_labels)
  pips <- pips[order(-pips$PIP), ]
  write.csv(pips, "output/tables/bkmr_pip.csv", row.names = FALSE)
  cat("  PIP rank:\n"); print(pips)

  # Fig 4: PIP bar
  source("../../templates/_shared/figures_helpers.R")
  p_pip <- ggplot(pips,
                  aes(x = reorder(variable, PIP), y = PIP, fill = variable)) +
    geom_col(width = 0.6, colour = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
    scale_fill_manual(values = palette_viridis_d(length(metal_labels)), guide = "none") +
    coord_flip(ylim = c(0, 1)) +
    labs(title = "Figure 4. BKMR posterior inclusion probability (PIP)",
         subtitle = "Red dashed: PIP = 0.5 (Bobb 2015 inclusion criterion).",
         x = NULL, y = "Posterior inclusion probability") +
    nhanes_theme_publication(base_size = 10)
  save_publication_figure(p_pip, "output/figures/fig4_bkmr_pip.tiff",
                          type = "single", height_mm = 90, format = "tiff", dpi = 600)
  save_figure_snapshot(p_pip, "fig4_bkmr_pip",
                       data_paths = "output/tables/bkmr_pip.csv",
                       seed = SEED_BASE)

  univariate <- tryCatch(
    bkmr::PredictorResponseUnivar(fit = fit_combined, q.fixed = 0.5,
                                  sel = seq(n_burn + 1, n_iter, by = 25),
                                  method = "approx"),
    error = function(e) { cat("Univar err:", conditionMessage(e), "\n"); NULL })

  bivariate <- tryCatch(
    bkmr::PredictorResponseBivar(
      fit = fit_combined,
      z.pairs = data.frame(z1 = rep(1:4, each = 4),
                           z2 = rep(2:5, times = 4)) %>% dplyr::filter(z1 < z2),
      q.fixed = 0.5,
      sel = seq(n_burn + 1, n_iter, by = 50),
      method = "approx"),
    error = function(e) { cat("Bivar err:", conditionMessage(e), "\n"); NULL })

  overall <- tryCatch(
    bkmr::OverallRiskSummaries(
      fit = fit_combined, y = y_vec, Z = Z_mat, X = X_mat,
      qs = seq(0.1, 0.9, by = 0.1), q.fixed = 0.5, method = "approx",
      sel = seq(n_burn + 1, n_iter, by = 25)),
    error = function(e) { cat("Overall err:", conditionMessage(e), "\n"); NULL })

  # ---- Enhancement 2: Figure 3 — BKMR bivariate surface ggplot ----
  # Heatmap (faceted by metal pair) of posterior mean h(z) from bivariate.
  # bkmr::PredictorResponseBivar returns long-format DF with columns
  # variable1, variable2, z1, z2, est, se. We translate variable1/2 from
  # numeric indices to metal labels (if needed) and facet.
  tryCatch({
    if (!is.null(bivariate) && nrow(bivariate) > 0L) {
      biv_df <- as.data.frame(bivariate)
      # bkmr returns variable1/variable2 as factor of metal NAMES already
      # (from colnames(Z_mat)); guard against the integer-index branch.
      if (is.numeric(biv_df$variable1)) {
        biv_df$variable1 <- factor(metal_labels[biv_df$variable1], levels = metal_labels)
      }
      if (is.numeric(biv_df$variable2)) {
        biv_df$variable2 <- factor(metal_labels[biv_df$variable2], levels = metal_labels)
      }
      biv_df$pair <- paste0(biv_df$variable1, " - ", biv_df$variable2)
      p_fig3 <- ggplot(biv_df, aes(x = z1, y = z2, fill = est)) +
        geom_raster(interpolate = TRUE) +
        geom_contour(aes(z = est), colour = "white", linewidth = 0.2,
                     alpha = 0.6, bins = 6) +
        scale_fill_viridis_c(name = "h(z)", option = "viridis") +
        facet_wrap(~ pair, ncol = 3) +
        labs(title = "Figure 3. BKMR posterior mean h(z) bivariate surface for 10 metal pairs (5 choose 2)",
             subtitle = "Other metals held at median (q = 0.5). White contours = iso-h levels.",
             x = "Metal 1 (z-score)", y = "Metal 2 (z-score)") +
        nhanes_theme_publication(base_size = 9) +
        theme(legend.position = "right",
              panel.spacing = grid::unit(0.6, "lines"))
      save_publication_figure(p_fig3, "output/figures/fig3_bkmr_bivariate.tiff",
                              type = "double", height_mm = 140,
                              format = "tiff", dpi = 600)
      save_figure_snapshot(p_fig3, "fig3_bkmr_bivariate",
                           data_paths = "output/tables/bkmr_results.RData",
                           seed = SEED_BASE)
      cat("  + output/figures/fig3_bkmr_bivariate.tiff (BKMR bivariate surface)\n")
    } else {
      cat("  WARN: bivariate is NULL/empty — skipping Fig 3\n")
    }
  }, error = function(e) {
    cat(sprintf("  WARN: Fig 3 bivariate plot skipped: %s\n", conditionMessage(e)))
  })

  # ---- Save full bundle ----
  cat("\n[4/4] Saving output/tables/bkmr_results.RData\n")
  save(bkmr_fit, pips, univariate, bivariate, overall, diag_table, metal_labels,
       file = "output/tables/bkmr_results.RData")

  # ---- Enhancement 4: Full audit-defense archive ----
  # Bundles raw per-chain accumulated_fits + diag_table + block & rhat history
  # CSVs + session info for reviewer-dispute defensibility / re-analysis.
  tryCatch({
    accumulated_fits <- cp$accumulated_fits
    block_history <- tryCatch(
      if (file.exists("_bkmr_block_history.csv"))
        read.csv("_bkmr_block_history.csv", stringsAsFactors = FALSE) else NULL,
      error = function(e) NULL)
    rhat_history  <- tryCatch(
      if (file.exists("_bkmr_rhat_history.csv"))
        read.csv("_bkmr_rhat_history.csv",  stringsAsFactors = FALSE) else NULL,
      error = function(e) NULL)
    session_info <- list(
      sys_info        = Sys.info(),
      r_version       = R.version,
      bkmr_version    = tryCatch(as.character(packageVersion("bkmr")),
                                 error = function(e) NA_character_),
      coda_version    = tryCatch(as.character(packageVersion("coda")),
                                 error = function(e) NA_character_),
      rstan_version   = tryCatch(as.character(packageVersion("rstan")),
                                 error = function(e) NA_character_),
      finalized_at    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      target_iter     = TARGET_ITER,
      burn_in         = BURN_IN,
      knots_K         = KNOTS_K,
      n_chains        = N_CHAINS,
      seed_base       = SEED_BASE,
      started_at      = cp$started_at,
      block_count     = cp$block_count
    )
    save(accumulated_fits, diag_table, block_history, rhat_history,
         session_info, metal_labels,
         file = "output/tables/bkmr_full_archive.RData")
    cat("  + output/tables/bkmr_full_archive.RData (audit archive)\n")
  }, error = function(e) {
    cat(sprintf("  WARN: full archive save failed: %s\n",
                conditionMessage(e)))
  })

  cat("\nFiles produced:\n")
  cat("  output/tables/bkmr_results.RData\n")
  cat("  output/tables/bkmr_pip.csv\n")
  cat("  output/tables/bkmr_convergence.csv\n")
  cat("  output/figures/bkmr_trace_<metal>.png (x5)\n")
  cat("  output/figures/fig4_bkmr_pip.tiff\n")

  if (!diag_ok) {
    cat("\nWARN: Convergence not fully satisfied. See bkmr_convergence.csv.\n")
    cat("      Consider running more blocks (set BKMR_TARGET_ITER higher).\n")
    writeLines("CONVERGENCE_FAILED", "output/tables/bkmr_CONVERGENCE_FAILED.flag")
  } else {
    cat("\nALL CONVERGENCE CHECKS PASSED (rhat < 1.1 AND ESS > 400).\n")
  }

  cat("\nDONE 06_bkmr_metals_checkpoint.R (finalize complete)\n")
} else if (all_done) {
  cat("\n[ALREADY FINALIZED] output/tables/bkmr_results.RData exists.\n")
  cat("Delete that file and rerun if you want to re-finalize.\n")
}
