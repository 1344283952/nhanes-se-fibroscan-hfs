# ============================================
# 006 / _render_bkmr_4metal_figures.R
#
# Post-process only: 4-metal BKMR (Pb / Cd / Hg / Se, Mn dropped) finalize 后
# 独立渲染 PIP CSV + trace + convergence + fig13 PIP bar + univariate +
# Hg×Se bivariate surface。**不重启 chainer**, 仅从已存 checkpoint 读取
# accumulated_fits 做 post-burnin posterior summaries / figures。
#
# 输入:
#   data/processed/bkmr_4metal_drop_Mn_006_checkpoint.rds
#   (由 scripts/23_bkmr_se_metals_checkpoint.R 在 SENSITIVITY_MODE=4metal_drop_Mn
#    下生成; 包含 cp$accumulated_fits[[1..2]] + cp$done_iter_per_chain)
#
# 输出:
#   output/tables/bkmr_4metal_pip.csv
#   output/tables/bkmr_4metal_convergence.csv
#   output/figures/bkmr_4metal_trace_<metal>.png   (×4 个 metal)
#   output/figures/fig13_bkmr_4metal_pip.tiff
#   output/figures/fig13_bkmr_4metal_pip.pdf
#   output/figures/bkmr_4metal_univariate.png
#   output/figures/bkmr_4metal_hg_se_surface.png
#
# 调用:
#   cd projects/006_se_fibroscan_hfs/
#   Rscript scripts/_render_bkmr_4metal_figures.R
#
# 设计:
#   * 每个 figure section 用 tryCatch 包裹, 失败时落 stub (zero-byte placeholder
#     + .stub 标记), 不阻塞下一个 section
#   * 中文注释 + 英文 identifier
#   * 末尾 cat() 报告输出文件清单
#
# 已 hardcode 的 assumption (开新场景请改这里):
#   - metal_labels = c("Pb", "Cd", "Hg", "Se")  顺序固定
#   - BURN_IN  = 5000
#   - TARGET_ITER = 10000
#   - N_CHAINS = 2
#   - CHECKPOINT_FILE = "data/processed/bkmr_4metal_drop_Mn_006_checkpoint.rds"
#   - Bivariate surface 固定看 Hg × Se (索引 3 vs 4 in metal_labels)
# ============================================

RNGkind("L'Ecuyer-CMRG")

suppressPackageStartupMessages({
  library(dplyr)
  library(bkmr)
  library(coda)
  library(rstan)
  library(ggplot2)
})

cat("========================================\n")
cat("006 BKMR 4-metal (Pb/Cd/Hg/Se) post-process renderer\n")
cat("========================================\n\n")

# ---- Hardcoded config ----
CHECKPOINT_FILE <- "data/processed/bkmr_4metal_drop_Mn_006_checkpoint.rds"
BURN_IN         <- 5000L
TARGET_ITER     <- 10000L
N_CHAINS        <- 2L
metal_labels    <- c("Pb", "Cd", "Hg", "Se")
SEED_BASE       <- 20260516L

# ---- Output dirs ----
if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# ---- 工具: 输出 stub 文件 (在 section 失败时落地, 便于上游 inventory check) ----
write_stub <- function(path, reason) {
  con <- file(path, open = "wb"); close(con)            # zero-byte placeholder
  writeLines(reason, paste0(path, ".stub"))
  cat(sprintf("  STUB: %s  (%s)\n", path, reason))
}

# ---- Step 1: Read checkpoint ----
cat("[1] Reading checkpoint...\n")
if (!file.exists(CHECKPOINT_FILE)) {
  stop("Checkpoint not found: ", CHECKPOINT_FILE,
       "\n  请先用 BKMR_SENSITIVITY_MODE=4metal_drop_Mn 跑 23_bkmr_se_metals_checkpoint.R")
}
cp <- readRDS(CHECKPOINT_FILE)

# ---- Step 2: 验证 BKMR fully done ----
cat("[2] Verifying done_iter_per_chain == c(", TARGET_ITER, ",", TARGET_ITER, ")\n")
done_iter <- cp$done_iter_per_chain
if (is.null(done_iter)) {
  stop("cp$done_iter_per_chain missing — checkpoint structure unrecognized.")
}
cat(sprintf("  done_iter_per_chain = c(%s)\n", paste(done_iter, collapse = ", ")))
if (!all(done_iter == TARGET_ITER)) {
  stop(sprintf("BKMR not fully done: have c(%s), expected c(%d, %d). 先跑完再调用本脚本。",
               paste(done_iter, collapse = ", "), TARGET_ITER, TARGET_ITER))
}
if (length(cp$accumulated_fits) < N_CHAINS) {
  stop("accumulated_fits 数量不足: 期望 ", N_CHAINS, " 个 chain, 实得 ",
       length(cp$accumulated_fits))
}
cat("  PASS: BKMR fully done across both chains.\n\n")

# ---- Step 3: 提取 post-burnin samples ----
cat("[3] Extracting post-burnin samples (drop first ", BURN_IN, " iter)\n")
bkmr_fit <- list(fits = cp$accumulated_fits)
n_iter   <- TARGET_ITER
n_burn   <- BURN_IN
n_chain  <- N_CHAINS
sel_post <- seq(n_burn + 1L, n_iter)             # post-burnin indices

extract_chain_param <- function(fit_obj, param_name, chain_i,
                                after_burn = n_burn) {
  fit  <- fit_obj$fits[[chain_i]]
  vals <- fit[[param_name]]
  if (is.null(vals)) return(NULL)
  if (is.matrix(vals)) vals[(after_burn + 1):n_iter, , drop = FALSE]
  else vals[(after_burn + 1):n_iter]
}
cat("  PASS: post-burnin window =", length(sel_post), "samples / chain\n\n")

# ---- Step 4: PIP table ----
cat("[4] PIP table -> output/tables/bkmr_4metal_pip.csv\n")
fit_combined <- bkmr_fit$fits[[1]]   # bkmr::ExtractPIPs 单 chain 用即可
pips <- tryCatch(
  {
    p <- bkmr::ExtractPIPs(fit_combined)
    p$variable <- factor(p$variable, levels = metal_labels)
    p <- p[order(-p$PIP), ]
    write.csv(p, "output/tables/bkmr_4metal_pip.csv", row.names = FALSE)
    cat("  PIP rank:\n"); print(p)
    p
  },
  error = function(e) {
    cat(sprintf("  WARN: ExtractPIPs failed: %s\n", conditionMessage(e)))
    stub <- data.frame(variable = metal_labels, PIP = NA_real_)
    write.csv(stub, "output/tables/bkmr_4metal_pip.csv", row.names = FALSE)
    writeLines(conditionMessage(e), "output/tables/bkmr_4metal_pip.csv.stub")
    stub
  }
)
cat("\n")

# ---- Step 5: Trace plots ----
cat("[5] Trace plots -> output/figures/bkmr_4metal_trace_<metal>.png\n")
for (j in seq_along(metal_labels)) {
  metal     <- metal_labels[j]
  out_path  <- file.path("output/figures",
                         sprintf("bkmr_4metal_trace_%s.png", tolower(metal)))
  tryCatch({
    png(out_path, width = 900, height = 500, res = 110)
    par(mfrow = c(1, n_chain), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
    for (ci in seq_len(n_chain)) {
      fit    <- bkmr_fit$fits[[ci]]
      r_post <- fit$r
      if (is.null(r_post)) {
        plot.new(); title(sprintf("Chain %d: r NA", ci)); next
      }
      plot(r_post[, j], type = "l",
           xlab = "iteration", ylab = sprintf("r[%s]", metal),
           main = sprintf("Chain %d", ci), col = "steelblue")
      abline(v = n_burn, col = "red", lty = 2)
    }
    mtext(sprintf("BKMR trace: r posterior for %s (red = burn-in cutoff)",
                  metal), outer = TRUE, cex = 1.1)
    dev.off()
    cat(sprintf("  + %s\n", out_path))
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    write_stub(out_path, sprintf("trace %s err: %s", metal, conditionMessage(e)))
  })
}
cat("\n")

# ---- Step 6: Convergence diagnostics (rhat + ESS) ----
cat("[6] Convergence diagnostics -> output/tables/bkmr_4metal_convergence.csv\n")
diag_table <- data.frame(param = character(0), rhat = numeric(0),
                         ess = numeric(0), pass = logical(0),
                         stringsAsFactors = FALSE)
diag_ok <- TRUE
tryCatch({
  # 监控参数: lambda + sigsq.eps + r_<metal> 4 个
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
      lapply(seq_len(n_chain),
             function(ci) extract_chain_param(bkmr_fit, p, ci))
    }
    chains_data <- chains_data[!sapply(chains_data,
                                       function(x) is.null(x) || all(is.na(x)))]
    if (length(chains_data) < 2) {
      cat(sprintf("  WARN: %-14s insufficient chains\n", pname)); next
    }
    M    <- do.call(cbind, chains_data)
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
  write.csv(diag_table, "output/tables/bkmr_4metal_convergence.csv",
            row.names = FALSE)
}, error = function(e) {
  cat(sprintf("  WARN: convergence diag failed: %s\n", conditionMessage(e)))
  write_stub("output/tables/bkmr_4metal_convergence.csv",
             conditionMessage(e))
})
cat("\n")

# ---- Step 7: Fig13 PIP bar chart (TIFF + PDF) ----
cat("[7] Fig13 PIP bar chart -> output/figures/fig13_bkmr_4metal_pip.{tiff,pdf}\n")
fig13_tiff <- "output/figures/fig13_bkmr_4metal_pip.tiff"
fig13_pdf  <- "output/figures/fig13_bkmr_4metal_pip.pdf"
tryCatch({
  pips_plot <- pips
  if (is.null(pips_plot) || nrow(pips_plot) == 0L)
    stop("pips empty; cannot draw fig13.")
  pips_plot$variable <- factor(pips_plot$variable, levels = metal_labels)

  p_fig13 <- ggplot(pips_plot,
                    aes(x = reorder(variable, PIP), y = PIP,
                        fill = variable)) +
    geom_col(width = 0.6, colour = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    coord_flip(ylim = c(0, 1)) +
    labs(title = "Figure 13. BKMR posterior inclusion probability (PIP), 4-metal model (Mn dropped)",
         subtitle = "Red dashed: PIP = 0.5 (Bobb 2015 inclusion criterion). Metals: Pb / Cd / Hg / Se.",
         x = NULL, y = "Posterior inclusion probability") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey25"),
          axis.text     = element_text(colour = "black"))

  # 优先用 publication helper, 不存在则 fall back 到原生 device
  helper_path <- "../../templates/_shared/figures_helpers.R"
  used_helper <- FALSE
  if (file.exists(helper_path)) {
    tryCatch({
      source(helper_path)
      if (exists("save_publication_figure")) {
        save_publication_figure(p_fig13, fig13_tiff,
                                type = "single", height_mm = 90,
                                format = "tiff", dpi = 600)
        save_publication_figure(p_fig13, fig13_pdf,
                                type = "single", height_mm = 90,
                                format = "pdf")
        used_helper <- TRUE
      }
    }, error = function(e) {
      cat(sprintf("  NOTE: helper failed (%s); falling back\n",
                  conditionMessage(e)))
    })
  }
  if (!used_helper) {
    ggsave(fig13_tiff, p_fig13, width = 6, height = 4, dpi = 600,
           compression = "lzw")
    ggsave(fig13_pdf,  p_fig13, width = 6, height = 4)
  }
  cat(sprintf("  + %s\n  + %s\n", fig13_tiff, fig13_pdf))
}, error = function(e) {
  write_stub(fig13_tiff, sprintf("fig13 err: %s", conditionMessage(e)))
  write_stub(fig13_pdf,  sprintf("fig13 err: %s", conditionMessage(e)))
})
cat("\n")

# ---- Step 8: Univariate response (h(z_metal) at median other metals) ----
cat("[8] Univariate response -> output/figures/bkmr_4metal_univariate.png\n")
univar_out <- "output/figures/bkmr_4metal_univariate.png"
tryCatch({
  univariate <- bkmr::PredictorResponseUnivar(
    fit     = fit_combined,
    q.fixed = 0.5,
    sel     = seq(n_burn + 1, n_iter, by = 25),
    method  = "approx"
  )
  if (is.null(univariate) || nrow(univariate) == 0L)
    stop("PredictorResponseUnivar returned empty.")

  uni_df <- as.data.frame(univariate)
  # bkmr returns "variable" column as factor of metal names by default
  if (is.numeric(uni_df$variable)) {
    uni_df$variable <- factor(metal_labels[uni_df$variable],
                              levels = metal_labels)
  } else {
    uni_df$variable <- factor(as.character(uni_df$variable),
                              levels = metal_labels)
  }

  p_uni <- ggplot(uni_df, aes(x = z, y = est)) +
    geom_ribbon(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se),
                fill = "steelblue", alpha = 0.25) +
    geom_line(colour = "steelblue", linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ variable, ncol = 2, scales = "free_x") +
    labs(title = "BKMR univariate posterior mean h(z_metal), 4-metal model",
         subtitle = "Other 3 metals held at median (q = 0.5). Ribbon = 95% credible interval.",
         x = "Metal exposure (z-score, log-standardized)",
         y = "h(z) (relative log-scale on outcome)") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey25"),
          strip.background = element_rect(fill = "grey90"))

  ggsave(univar_out, p_uni, width = 7, height = 5, dpi = 300)
  cat(sprintf("  + %s\n", univar_out))
}, error = function(e) {
  write_stub(univar_out, sprintf("univar err: %s", conditionMessage(e)))
})
cat("\n")

# ---- Step 9: Bivariate interaction surface (Hg × Se) ----
cat("[9] Hg × Se bivariate surface -> output/figures/bkmr_4metal_hg_se_surface.png\n")
biv_out <- "output/figures/bkmr_4metal_hg_se_surface.png"
tryCatch({
  hg_idx <- which(metal_labels == "Hg")
  se_idx <- which(metal_labels == "Se")
  if (length(hg_idx) != 1L || length(se_idx) != 1L)
    stop("metal_labels 中找不到 Hg 或 Se 索引")

  bivariate <- bkmr::PredictorResponseBivar(
    fit     = fit_combined,
    z.pairs = data.frame(z1 = hg_idx, z2 = se_idx),
    q.fixed = 0.5,
    sel     = seq(n_burn + 1, n_iter, by = 50),
    method  = "approx"
  )
  if (is.null(bivariate) || nrow(bivariate) == 0L)
    stop("PredictorResponseBivar returned empty.")

  biv_df <- as.data.frame(bivariate)
  if (is.numeric(biv_df$variable1))
    biv_df$variable1 <- metal_labels[biv_df$variable1]
  if (is.numeric(biv_df$variable2))
    biv_df$variable2 <- metal_labels[biv_df$variable2]

  p_biv <- ggplot(biv_df, aes(x = z1, y = z2, fill = est)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = est), colour = "white",
                 linewidth = 0.25, alpha = 0.7, bins = 8) +
    scale_fill_viridis_c(name = "h(z)", option = "viridis") +
    labs(title = "BKMR posterior mean h(z) surface: Hg × Se",
         subtitle = "Other 2 metals (Pb, Cd) at median (q = 0.5). White contours = iso-h levels.",
         x = "Hg exposure (z-score, log-standardized)",
         y = "Se exposure (z-score, log-standardized)") +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, colour = "grey25"),
          legend.position = "right")

  ggsave(biv_out, p_biv, width = 6, height = 5, dpi = 300)
  cat(sprintf("  + %s\n", biv_out))
}, error = function(e) {
  write_stub(biv_out, sprintf("bivariate err: %s", conditionMessage(e)))
})
cat("\n")

# ---- Final inventory ----
cat("========================================\n")
cat("DONE _render_bkmr_4metal_figures.R\n")
cat("========================================\n")
cat("Files produced (CSV + PNG + TIFF + PDF):\n")
cat("  output/tables/bkmr_4metal_pip.csv\n")
cat("  output/tables/bkmr_4metal_convergence.csv\n")
for (metal in metal_labels) {
  cat(sprintf("  output/figures/bkmr_4metal_trace_%s.png\n", tolower(metal)))
}
cat("  output/figures/fig13_bkmr_4metal_pip.tiff\n")
cat("  output/figures/fig13_bkmr_4metal_pip.pdf\n")
cat("  output/figures/bkmr_4metal_univariate.png\n")
cat("  output/figures/bkmr_4metal_hg_se_surface.png\n")
if (!diag_ok) {
  cat("\nWARN: convergence not fully satisfied — 检查 bkmr_4metal_convergence.csv\n")
}
cat("\nNote: 任何带 .stub 后缀的文件意味着该 section 失败, 文件本体为 0 字节占位\n")
