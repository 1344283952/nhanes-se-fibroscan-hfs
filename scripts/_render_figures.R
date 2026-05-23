# ============================================================
# scripts/_render_figures.R
# 渲染 publication-grade figures for 006_se_fibroscan_hfs
# 目标期刊: FRBM / Hepatology Communications (TIFF LZW)
#
# 依赖: ggplot2, viridis, viridisLite, patchwork, rsvg (可选), magick (可选)
# 输出: output/figures/fig{1,3,4,5,6,7,8}_*.tiff + .pdf preview
#       output/figures/_snapshot/fig{1,3,4,5,6,7,8}.RData
# 用法: 必须在 projects/006_se_fibroscan_hfs/ 工作目录下执行
#       Rscript scripts/_render_figures.R
# ============================================================

set.seed(20260516)

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(viridisLite)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# 加载共享 helpers
source("../../templates/_shared/figures_helpers.R")

# 确保输出目录存在
fig_dir   <- "output/figures"
snap_dir  <- "output/figures/_snapshot"
if (!dir.exists(fig_dir))  dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(snap_dir)) dir.create(snap_dir, recursive = TRUE)

# 记录运行结果
produced <- character()
errors   <- character()

safe_run <- function(fig_id, expr) {
  out <- tryCatch(expr, error = function(e) {
    msg <- sprintf("%s: %s", fig_id, conditionMessage(e))
    errors <<- c(errors, msg)
    message("ERROR ", msg)
    NULL
  })
  invisible(out)
}

# ============================================================
# Fig 1 — CONSORT (从已有 svg/png 升级 TIFF 600 dpi)
# ============================================================
safe_run("fig1", {
  svg_path  <- file.path(fig_dir, "fig1_consort.svg")
  png_path  <- file.path(fig_dir, "fig1_consort.png")
  tiff_path <- file.path(fig_dir, "fig1_consort.tiff")

  if (file.exists(svg_path) && requireNamespace("rsvg", quietly = TRUE)) {
    # rsvg 输出 TIFF 600 dpi (178 mm 宽 × 200 mm 高)
    # 178 mm = 7.01 in; @ 600 dpi -> 4206 px 宽
    rsvg::rsvg_pdf(svg_path, paste0(tools::file_path_sans_ext(tiff_path), ".pdf"))
    # rsvg 没有原生 TIFF，但可以输出 PNG @ 600 dpi 再转 TIFF
    rsvg_png_path <- tempfile(fileext = ".png")
    rsvg::rsvg_png(svg_path, rsvg_png_path, width = 4206)
    # 用 grDevices::tiff 转写 (LZW)
    img <- png::readPNG(rsvg_png_path)
    h_in <- 200 / 25.4
    w_in <- 178 / 25.4
    tiff(tiff_path, width = w_in, height = h_in, units = "in", res = 600,
         compression = "lzw", bg = "white")
    grid::grid.raster(img)
    dev.off()
    produced <<- c(produced, tiff_path)
    message("Saved ", tiff_path, " via rsvg")
  } else {
    # fallback: 把 PNG 重打包为 TIFF (用 grDevices::tiff)
    if (file.exists(png_path) && requireNamespace("png", quietly = TRUE)) {
      img <- png::readPNG(png_path)
      h_in <- 200 / 25.4
      w_in <- 178 / 25.4
      tiff(tiff_path, width = w_in, height = h_in, units = "in", res = 600,
           compression = "lzw", bg = "white")
      grid::grid.raster(img)
      dev.off()
      produced <<- c(produced, tiff_path)
      message("Saved ", tiff_path, " from PNG fallback")
    } else {
      stop("Neither rsvg+svg nor png fallback available for Fig 1")
    }
  }

  # 快照只存路径
  snap <- list(source_svg = svg_path, source_png = png_path,
               flow_csv = "output/tables/flow_counts.csv",
               seed = 20260516, timestamp = Sys.time())
  save(snap, file = file.path(snap_dir, "fig1.RData"))
})

# ============================================================
# Fig 3 — GAM bivariate heatmap (DR1TSELE × LBXBSE → CAP/HFS/LSM)
# 3 panel composite, viridis-D fill, TIFF 300 dpi half-tone
# ============================================================
safe_run("fig3", {
  gam_cap <- read.csv("output/tables/gam_grid_cap.csv", stringsAsFactors = FALSE)
  gam_hfs <- read.csv("output/tables/gam_grid_hfs.csv", stringsAsFactors = FALSE)
  gam_lsm <- read.csv("output/tables/gam_grid_lsm.csv", stringsAsFactors = FALSE)

  make_heatmap <- function(df, fill_label, panel_title) {
    # W13 R3 R-Figure P0 fix: the GAM prediction grid is irregular (19 unique
    # x and 19 unique y quantile-spaced values, gaps up to 38.6 vs median 7.75),
    # so geom_raster(interpolate=TRUE) produced vertical stripes (W11 attempt)
    # and geom_tile() without explicit width/height left white gaps (W12 fail).
    # stat_contour_filled() integrates the surface and produces seamless bands
    # that are robust to irregular sampling. Median cross-hairs and 8-band
    # binning preserve the visual gradient.
    df_sum <- df %>%
      group_by(DR1TSELE, LBXBSE) %>%
      summarise(pred = mean(pred, na.rm = TRUE), .groups = "drop")

    med_x <- stats::median(df_sum$DR1TSELE, na.rm = TRUE)
    med_y <- stats::median(df_sum$LBXBSE,   na.rm = TRUE)

    ggplot(df_sum, aes(x = DR1TSELE, y = LBXBSE, z = pred)) +
      stat_contour_filled(bins = 9, alpha = 0.95) +
      geom_hline(yintercept = med_y, linetype = "dotted",
                 color = "white", linewidth = 0.35) +
      geom_vline(xintercept = med_x, linetype = "dotted",
                 color = "white", linewidth = 0.35) +
      scale_fill_viridis_d(option = "D", name = fill_label) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(title = panel_title,
           x = expression("Dietary Se (DR1TSELE, "*mu*"g/day)"),
           y = expression("Blood Se (LBXBSE, "*mu*"g/L)")) +
      nhanes_theme_publication(base_size = 9) +
      theme(legend.position = "right",
            legend.key.width  = unit(3, "mm"),
            legend.key.height = unit(8, "mm"),
            legend.text       = element_text(size = 6),
            plot.title        = element_text(size = 8, face = "bold"),
            plot.margin       = margin(2, 8, 2, 2))
  }

  p1 <- make_heatmap(gam_cap, "CAP\n(dB/m)",  "A. CAP (steatosis)")
  p2 <- make_heatmap(gam_lsm, "LSM\n(kPa)",   "B. LSM (fibrosis)")
  p3 <- make_heatmap(gam_hfs, "HFS",          "C. Hepamet HFS")

  fig3 <- p1 + p2 + p3 + plot_layout(ncol = 3)

  save_publication_figure(
    fig3,
    file.path(fig_dir, "fig4_gam_heatmap.tiff"),
    type = "double",
    height_mm = 90,
    format = "tiff",
    dpi = 300
  )
  produced <<- c(produced, file.path(fig_dir, "fig4_gam_heatmap.tiff"))

  save_figure_snapshot(
    fig3, "fig3",
    data_paths = c("output/tables/gam_grid_cap.csv",
                   "output/tables/gam_grid_hfs.csv",
                   "output/tables/gam_grid_lsm.csv"),
    out_dir = snap_dir
  )
})

# ============================================================
# Fig 4 — RCS curves (从 rcs_ushape.RData$results_list[[*]]$pred 重渲染)
# 6 panel: cap/lsm/hfs × DR1TSELE/LBXBSE  —— 用 patchwork 拼
# TIFF 300 dpi
# ============================================================
safe_run("fig4", {
  load("output/tables/rcs_ushape.RData")  # results_list, pvals_df

  panel_labels <- c(
    cap_DR1TSELE = "A. CAP × Dietary Se",
    lsm_DR1TSELE = "B. LSM × Dietary Se",
    hfs_DR1TSELE = "C. HFS × Dietary Se",
    cap_LBXBSE   = "D. CAP × Blood Se",
    lsm_LBXBSE   = "E. LSM × Blood Se",
    hfs_LBXBSE   = "F. HFS × Blood Se"
  )
  y_labels <- c(
    cap_DR1TSELE = "Predicted CAP (dB/m)",
    lsm_DR1TSELE = "Predicted LSM (log-OR)",
    hfs_DR1TSELE = "Predicted HFS",
    cap_LBXBSE   = "Predicted CAP (dB/m)",
    lsm_LBXBSE   = "Predicted LSM (log-OR)",
    hfs_LBXBSE   = "Predicted HFS"
  )
  x_labels <- c(
    cap_DR1TSELE = expression("Dietary Se ("*mu*"g/day)"),
    lsm_DR1TSELE = expression("Dietary Se ("*mu*"g/day)"),
    hfs_DR1TSELE = expression("Dietary Se ("*mu*"g/day)"),
    cap_LBXBSE   = expression("Blood Se ("*mu*"g/L)"),
    lsm_LBXBSE   = expression("Blood Se ("*mu*"g/L)"),
    hfs_LBXBSE   = expression("Blood Se ("*mu*"g/L)")
  )

  panel_colors <- palette_viridis_d(2)  # dietary = 第 1, blood = 第 2

  build_panel <- function(name) {
    res  <- results_list[[name]]
    pred <- as.data.frame(res$pred)
    is_blood <- grepl("LBXBSE", name)
    color <- if (is_blood) panel_colors[2] else panel_colors[1]
    x_var <- if (is_blood) "LBXBSE" else "DR1TSELE"

    # 找对应 p-value
    parts <- strsplit(name, "_")[[1]]
    outc  <- parts[1]
    expo  <- parts[2]
    p_row <- pvals_df[pvals_df$outcome == outc & pvals_df$exposure == expo, ]
    p_nl  <- ifelse(nrow(p_row) == 1, p_row$p_nonlinear_anova, NA)
    p_lbl <- if (!is.na(p_nl)) {
      if (p_nl < 0.001) "italic(P)[nonlinear] < 0.001"
      else sprintf("italic(P)[nonlinear] == %.3f", p_nl)
    } else ""

    p <- ggplot(pred, aes_string(x = x_var, y = "yhat")) +
      geom_ribbon(aes_string(ymin = "lower", ymax = "upper"),
                  fill = color, alpha = 0.25) +
      geom_line(color = color, linewidth = 0.6) +
      labs(title = panel_labels[[name]],
           x = x_labels[[name]],
           y = y_labels[[name]]) +
      nhanes_theme_publication(base_size = 9)

    if (nzchar(p_lbl)) {
      p <- p + annotate("text",
                        x = -Inf, y = Inf,
                        label = p_lbl,
                        parse = TRUE,
                        hjust = -0.1, vjust = 1.4,
                        size = 2.6, color = "grey20")
    }
    p
  }

  panels <- lapply(names(panel_labels), build_panel)
  fig4   <- patchwork::wrap_plots(panels, ncol = 3)

  save_publication_figure(
    fig4,
    file.path(fig_dir, "fig3_rcs_se.tiff"),
    type = "double",
    height_mm = 140,
    format = "tiff",
    dpi = 300
  )
  produced <<- c(produced, file.path(fig_dir, "fig3_rcs_se.tiff"))

  save_figure_snapshot(fig4, "fig4",
                       data_paths = c("output/tables/rcs_ushape.RData",
                                      "output/tables/rcs_pvalues.csv"),
                       out_dir = snap_dir)
})

# ============================================================
# Fig 5 — HFS calibration (从 hfs_calibration_curve.csv)
# TIFF 600 dpi
# ============================================================
safe_run("fig5", {
  cal <- read.csv("output/tables/hfs_calibration_curve.csv", stringsAsFactors = FALSE)
  auroc <- read.csv("output/tables/hfs_auroc_summary.csv", stringsAsFactors = FALSE)
  v <- setNames(auroc$value, auroc$metric)

  auroc_lbl <- sprintf("AUROC = %.3f", v["AUROC_HFS_apparent"])
  delta_lbl <- sprintf("Delta AUROC (HFS - FIB-4) = %.3f",
                       v["Delta_AUROC_HFS_minus_FIB4"])

  fig5 <- ggplot(cal, aes(x = predicted, y = observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.4) +
    geom_line(color = palette_viridis_d(5)[1], linewidth = 0.7) +
    geom_point(color = palette_viridis_d(5)[1], size = 0.5, alpha = 0.6) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    annotate("text", x = 0.02, y = 0.98, label = auroc_lbl,
             hjust = 0, vjust = 1, size = 3) +
    annotate("text", x = 0.02, y = 0.92, label = delta_lbl,
             hjust = 0, vjust = 1, size = 3) +
    labs(x = "Predicted probability (HFS)",
         y = "Observed event rate") +
    nhanes_theme_publication(base_size = 10)

  save_publication_figure(
    fig5,
    file.path(fig_dir, "fig5_hfs_calibration.tiff"),
    type = "single",
    height_mm = 90,
    format = "tiff",
    dpi = 600
  )
  produced <<- c(produced, file.path(fig_dir, "fig5_hfs_calibration.tiff"))

  save_figure_snapshot(fig5, "fig5",
                       data_paths = c("output/tables/hfs_calibration_curve.csv",
                                      "output/tables/hfs_auroc_summary.csv"),
                       out_dir = snap_dir)
})

# ============================================================
# Fig 6 — Q1-Q4 cross-classification heatmap (CAP ≥ 275)
# 从 cross_cap275_cells.csv —— 4×4 heatmap with text annotation
# TIFF 600 dpi
# ============================================================
safe_run("fig6", {
  cross <- read.csv("output/tables/cross_cap275_cells.csv", stringsAsFactors = FALSE)
  cross$diet_q  <- factor(cross$diet_q,  levels = c("D1", "D2", "D3", "D4"))
  cross$blood_q <- factor(cross$blood_q, levels = c("B1", "B2", "B3", "B4"))
  cross$pct_lab <- sprintf("%.1f%%", cross$steatosis_cap275 * 100)

  fig6 <- ggplot(cross, aes(x = diet_q, y = blood_q, fill = steatosis_cap275)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = pct_lab), size = 3, color = "white", fontface = "bold") +
    scale_fill_viridis_c(option = "D",
                         name = "CAP >= 275\nprevalence",
                         labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = c("D1 (low)", "D2", "D3", "D4 (high)")) +
    scale_y_discrete(labels = c("B1 (low)", "B2", "B3", "B4 (high)")) +
    labs(x = "Dietary Se quartile (DR1TSELE)",
         y = "Blood Se quartile (LBXBSE)") +
    nhanes_theme_publication(base_size = 10) +
    coord_fixed()

  save_publication_figure(
    fig6,
    file.path(fig_dir, "fig6_cross_heatmap.tiff"),
    type = "single",
    height_mm = 100,
    format = "tiff",
    dpi = 600
  )
  produced <<- c(produced, file.path(fig_dir, "fig6_cross_heatmap.tiff"))

  save_figure_snapshot(fig6, "fig6",
                       data_paths = "output/tables/cross_cap275_cells.csv",
                       out_dir = snap_dir)
})

# ============================================================
# Fig 7 — Ratio tertile OR forest (Se/Zn + Se/Cu)
# 从 ratio_OR_cap275.csv —— T2, T3 OR (95% CI) ratio family colored
# TIFF 600 dpi line art
# ============================================================
safe_run("fig7", {
  ratio <- read.csv("output/tables/ratio_OR_cap275.csv", stringsAsFactors = FALSE)

  ratio$ratio_label <- ifelse(ratio$ratio == "se_zn_t", "Se/Zn", "Se/Cu")
  ratio$tertile     <- ifelse(grepl("T2", ratio$term), "T2 vs T1",
                              ifelse(grepl("T3", ratio$term), "T3 (high) vs T1", "T1"))
  ratio$label       <- paste(ratio$ratio_label, ratio$tertile, sep = ": ")
  ratio$label       <- factor(ratio$label, levels = rev(ratio$label))

  rdbu_cols <- palette_rdbu_diverging(11)
  ratio_pal <- c("Se/Zn" = rdbu_cols[10], "Se/Cu" = rdbu_cols[2])  # blue / red

  fig7 <- ggplot(ratio, aes(x = OR, y = label, color = ratio_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.2, linewidth = 0.5) +
    geom_point(size = 2.5, shape = 18) +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
                  limits = c(0.5, 2.5)) +
    scale_color_manual(values = ratio_pal, name = "Ratio family") +
    labs(x = "Odds ratio (95% CI), CAP >= 275 vs T1",
         y = NULL) +
    nhanes_theme_publication(base_size = 10) +
    theme(legend.position = "bottom")

  save_publication_figure(
    fig7,
    file.path(fig_dir, "fig7_ratio_forest.tiff"),
    type = "double",
    height_mm = 90,
    format = "tiff",
    dpi = 600
  )
  produced <<- c(produced, file.path(fig_dir, "fig7_ratio_forest.tiff"))

  save_figure_snapshot(fig7, "fig7",
                       data_paths = c("output/tables/ratio_OR_cap275.csv",
                                      "output/tables/evalue_ratio_OR.csv"),
                       out_dir = snap_dir)
})

# ============================================================
# Fig 8 — Subgroup forest (主图: CAP ≥ 275 × DR1TSELE; 9 subgroups)
# 从 subgroup_forest_data.csv 筛 outcome=cap275, exposure=DR1TSELE
# TIFF 600 dpi
# ============================================================
safe_run("fig8", {
  sub <- read.csv("output/tables/subgroup_forest_data.csv", stringsAsFactors = FALSE)
  sub_main <- sub %>%
    dplyr::filter(outcome == "cap275", exposure == "DR1TSELE") %>%
    dplyr::mutate(
      subgroup_clean = sub("_strata$", "", subgroup),
      label = paste0(subgroup_clean, ": ", stratum),
      OR    = effect,
      OR    = ifelse(is.na(OR) | OR <= 0, exp(beta), OR),
      lci_p = ifelse(is.na(lci) | lci <= 0, exp(beta - 1.96 * se), lci),
      uci_p = ifelse(is.na(uci) | uci <= 0, exp(beta + 1.96 * se), uci),
      p_lab = sprintf("P[BH]=%.2f", p_BH_subgroup)
    )
  sub_main$label <- factor(sub_main$label, levels = rev(sub_main$label))

  # 颜色按 subgroup family (9 个 subgroup 颜色块)
  unique_groups <- unique(sub_main$subgroup_clean)
  group_cols    <- setNames(palette_viridis_d(length(unique_groups)), unique_groups)

  fig8 <- ggplot(sub_main, aes(x = OR, y = label, color = subgroup_clean)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lci_p, xmax = uci_p), height = 0.25, linewidth = 0.4) +
    geom_point(size = 1.8, shape = 16) +
    scale_x_log10(breaks = c(0.3, 0.5, 0.75, 1, 1.5, 2, 3),
                  limits = c(0.3, 3.0)) +
    scale_color_manual(values = group_cols, name = "Subgroup", guide = "none") +
    labs(x = "OR per IQR Dietary Se (CAP >= 275)",
         y = NULL,
         title = "Subgroup analysis: dietary Se on CAP >= 275") +
    nhanes_theme_publication(base_size = 9) +
    theme(axis.text.y = element_text(size = 7))

  save_publication_figure(
    fig8,
    file.path(fig_dir, "fig8_subgroup_forest.tiff"),
    type = "double",
    height_mm = 180,
    format = "tiff",
    dpi = 600
  )
  produced <<- c(produced, file.path(fig_dir, "fig8_subgroup_forest.tiff"))

  save_figure_snapshot(fig8, "fig8",
                       data_paths = "output/tables/subgroup_forest_data.csv",
                       out_dir = snap_dir)
})

# ============================================================
# Summary
# ============================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("RENDERED FIGURES:\n")
if (length(produced) == 0) {
  cat("  (none)\n")
} else {
  for (f in produced) {
    info <- file.info(f)
    cat(sprintf("  %s  (%.1f KB)\n", f, info$size / 1024))
  }
}
cat("\nERRORS:\n")
if (length(errors) == 0) {
  cat("  (none)\n")
} else {
  for (e in errors) cat("  ", e, "\n", sep = "")
}
cat(strrep("=", 60), "\n", sep = "")
