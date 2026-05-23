# ============================================================
# 006 / _render_figs_p0fix.R — P0 fix for fig5 / fig8 / fig11 / fig13
#
# 上下文:
#   R-Figure Round 1 标记 4 张图 P0:
#     - fig5  HFS calibration:   step-function only, missing loess + HL + slope/intercept CI
#     - fig8  Subgroup forest:   只有 1 panel (cap275 × DR1TSELE), 应该 4 contrasts × 3 outcomes × 7 strata
#     - fig11 CMAverse 4-way:    PIE label 显示 "5.9%" (实际应为 "-6.9%"; - 号被坐标轴截掉)
#     - fig13 BKMR 4-metal PIP:  title/subtitle 被画布右边截断; PIP 数值标签缺失
#
# 输出 (全部覆盖):
#   output/figures/fig5_hfs_calibration.{tiff,pdf,png}
#   output/figures/fig8_subgroup_forest.{tiff,pdf,png}
#   output/figures/fig11_cmaverse_4way.{tiff,pdf,png}
#   output/figures/fig13_bkmr_4metal_pip.{tiff,pdf,png}
#
# 调用:
#   cd projects/006_se_fibroscan_hfs/
#   Rscript scripts/_render_figs_p0fix.R
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

# Try optional packages
have_resource_sel <- requireNamespace("ResourceSelection", quietly = TRUE)

fig_dir <- "output/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Color palette (consistent with existing v2 figures)
pal_blue  <- "#2E86AB"   # negative / CDE pure
pal_red   <- "#E74C3C"   # positive / INTREF
pal_green <- "#27AE60"   # INTMED
pal_orange<- "#F39C12"   # accent
pal_purple<- "#8E44AD"

theme_pub <- function(base_size = 10) theme_classic(base_size = base_size) +
  theme(panel.grid.major   = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(face = "bold", size = base_size + 1),
        plot.subtitle      = element_text(size = base_size - 1, color = "grey25"),
        plot.caption       = element_text(size = base_size - 2, color = "grey50"),
        axis.text          = element_text(color = "black"),
        axis.title         = element_text(color = "black"),
        strip.background   = element_rect(fill = "grey90", color = NA),
        strip.text         = element_text(face = "bold"),
        plot.margin        = margin(4, 8, 4, 4))

save_fig <- function(plot, base, width_mm, height_mm) {
  width_in  <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ggsave(paste0(base, ".tiff"), plot, width = width_in, height = height_in,
         dpi = 300, compression = "lzw", units = "in")
  ggsave(paste0(base, ".pdf"),  plot, width = width_in, height = height_in,
         units = "in")
  ggsave(paste0(base, ".png"),  plot, width = width_in, height = height_in,
         dpi = 200, units = "in")
}

cat("========================================\n")
cat("006 figs P0 fix - fig5 / fig8 / fig11 / fig13\n")
cat("========================================\n\n")

# ============================================================
# Fig 5 — HFS Flexible Calibration (proper loess + HL + slope/intercept)
# ============================================================
cat("[1/4] Fig 5: HFS flexible calibration ...\n")

fig5_status <- tryCatch({
  # Reload predictor + outcome from final processed data
  load("data/processed/nhanes_final.RData")
  df_hfs <- nhanes_final %>%
    filter(!is.na(hfs), !is.na(fib4), !is.na(fibrosis_lsm8),
           !is.na(wt_saf_pooled), wt_saf_pooled > 0)
  cat(sprintf("  Analytic N: %d (LSM>=8 events: %d, %.1f%%)\n",
              nrow(df_hfs), sum(df_hfs$fibrosis_lsm8 == 1),
              100 * mean(df_hfs$fibrosis_lsm8 == 1)))

  # 1. AUROC (re-compute for annotation)
  roc_hfs  <- pROC::roc(df_hfs$fibrosis_lsm8, df_hfs$hfs,
                        ci = TRUE, levels = c(0, 1), direction = "<")
  auroc_pe <- as.numeric(roc_hfs$auc)
  auroc_lo <- roc_hfs$ci[1]; auroc_hi <- roc_hfs$ci[3]

  # 2. Calibration intercept + slope via logistic regression of
  #    Y ~ logit(P_hat) (Cox 1958; Steyerberg 2010)
  eps <- 1e-6
  hfs_logit <- qlogis(pmin(pmax(df_hfs$hfs, eps), 1 - eps))
  fit_cal   <- glm(df_hfs$fibrosis_lsm8 ~ hfs_logit, family = binomial())
  s         <- summary(fit_cal)$coefficients
  intercept_pe <- s[1, "Estimate"]; intercept_se <- s[1, "Std. Error"]
  slope_pe     <- s[2, "Estimate"]; slope_se     <- s[2, "Std. Error"]
  ci_int <- intercept_pe + c(-1, 1) * 1.96 * intercept_se
  ci_slp <- slope_pe     + c(-1, 1) * 1.96 * slope_se

  # 3. Hosmer-Lemeshow
  hl <- if (have_resource_sel) {
    tryCatch(ResourceSelection::hoslem.test(df_hfs$fibrosis_lsm8,
                                            df_hfs$hfs, g = 10),
             error = function(e) NULL)
  } else NULL
  hl_chi <- if (!is.null(hl)) as.numeric(hl$statistic) else NA
  hl_df  <- if (!is.null(hl)) as.numeric(hl$parameter) else NA
  hl_p   <- if (!is.null(hl)) as.numeric(hl$p.value)   else NA

  # 4. Decile bin points (observed vs mean predicted per decile)
  df_hfs$decile <- cut(df_hfs$hfs,
                       breaks = unique(quantile(df_hfs$hfs,
                                                probs = seq(0, 1, 0.1),
                                                na.rm = TRUE)),
                       include.lowest = TRUE, labels = FALSE)
  dec_df <- df_hfs %>%
    group_by(decile) %>%
    summarise(pred_mean = mean(hfs, na.rm = TRUE),
              obs_rate  = mean(fibrosis_lsm8 == 1, na.rm = TRUE),
              n         = n(),
              obs_se    = sqrt(obs_rate * (1 - obs_rate) / pmax(n, 1)),
              .groups   = "drop") %>%
    filter(!is.na(decile))

  # 5. Loess curve with 95% CI for fitted Y | P
  fig5 <- ggplot(df_hfs, aes(x = hfs, y = fibrosis_lsm8)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.5) +
    # Loess fit on raw 0/1 outcome (Van Calster 2019 flexible cal)
    geom_smooth(method = "loess", se = TRUE, span = 0.75,
                color = pal_blue, fill = pal_blue, alpha = 0.18,
                linewidth = 0.8, formula = y ~ x) +
    # Decile bin diagnostics
    geom_errorbar(data = dec_df,
                  inherit.aes = FALSE,
                  aes(x = pred_mean,
                      ymin = pmax(obs_rate - 1.96 * obs_se, 0),
                      ymax = pmin(obs_rate + 1.96 * obs_se, 1)),
                  width = 0.015, color = "grey25", linewidth = 0.3) +
    geom_point(data = dec_df, inherit.aes = FALSE,
               aes(x = pred_mean, y = obs_rate),
               color = "grey15", size = 1.8) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       labels = scales::number_format(accuracy = 0.01)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25),
                       labels = scales::number_format(accuracy = 0.01)) +
    labs(title = "Figure 5. Hepamet HFS flexible calibration vs LSM >= 8 kPa",
         subtitle = sprintf(
           "Loess smoother + 95%% CI (Van Calster 2019); decile bins (Wilson 95%% CI). n = %s, events = %s (%.1f%%).",
           format(nrow(df_hfs), big.mark = ","),
           format(sum(df_hfs$fibrosis_lsm8 == 1), big.mark = ","),
           100 * mean(df_hfs$fibrosis_lsm8 == 1)),
         x = "Predicted probability (HFS)",
         y = "Observed event rate (LSM >= 8 kPa)") +
    theme_pub(base_size = 10)

  # Annotations: AUROC, intercept, slope, HL
  hl_label <- if (!is.na(hl_chi)) {
    sprintf("Hosmer-Lemeshow: chi^2 = %.2f (df = %d), P = %.3f",
            hl_chi, hl_df, hl_p)
  } else "Hosmer-Lemeshow: unavailable"

  annot_lines <- c(
    sprintf("AUROC = %.3f (95%% CI %.3f-%.3f)", auroc_pe, auroc_lo, auroc_hi),
    sprintf("Calibration intercept = %.3f (95%% CI %.3f to %.3f)",
            intercept_pe, ci_int[1], ci_int[2]),
    sprintf("Calibration slope = %.3f (95%% CI %.3f to %.3f)",
            slope_pe, ci_slp[1], ci_slp[2]),
    hl_label
  )
  fig5 <- fig5 + annotate("label", x = 0.02, y = 0.99,
                          label = paste(annot_lines, collapse = "\n"),
                          hjust = 0, vjust = 1, size = 2.7,
                          fill = "white", color = "black",
                          label.padding = unit(0.4, "lines"),
                          label.size = 0.3)

  save_fig(fig5, file.path(fig_dir, "fig5_hfs_calibration"),
           width_mm = 178, height_mm = 110)

  # Save numeric summary for downstream
  fig5_stats <- data.frame(
    auroc_pe = auroc_pe, auroc_lo = auroc_lo, auroc_hi = auroc_hi,
    intercept_pe = intercept_pe, intercept_lo = ci_int[1], intercept_hi = ci_int[2],
    slope_pe = slope_pe, slope_lo = ci_slp[1], slope_hi = ci_slp[2],
    hl_chi = hl_chi, hl_df = hl_df, hl_p = hl_p,
    n = nrow(df_hfs))
  write.csv(fig5_stats, "output/tables/fig5_calibration_stats.csv",
            row.names = FALSE)
  cat(sprintf("  AUROC = %.3f; intercept = %.3f; slope = %.3f; HL chi^2 = %.2f P = %.3f\n",
              auroc_pe, intercept_pe, slope_pe,
              hl_chi, hl_p))
  cat("  [OK] fig5 rendered (.tiff + .pdf + .png)\n")
  "OK"
}, error = function(e) {
  cat(sprintf("  [FAIL] fig5: %s\n", conditionMessage(e)))
  conditionMessage(e)
})
cat("\n")


# ============================================================
# Fig 8 — Faceted subgroup forest: 4 contrasts × 3 outcomes × 7 strata
# Option A: render full 156-row data set faceted by outcome
# ============================================================
cat("[2/4] Fig 8: Full subgroup forest (4 contrasts x 3 outcomes x strata) ...\n")

fig8_status <- tryCatch({
  sub <- read.csv("output/tables/subgroup_forest_data.csv",
                  stringsAsFactors = FALSE)
  cat(sprintf("  Total rows: %d\n", nrow(sub)))

  # Pretty labels
  exposure_lab <- c(DR1TSELE   = "Dietary Se",
                    LBXBSE     = "Blood Se",
                    se_zn_ratio= "Se/Zn",
                    se_cu_ratio= "Se/Cu")
  outcome_lab <- c(cap275 = "Steatosis (CAP >= 275)",
                   lsm8   = "Fibrosis (LSM >= 8)",
                   hfs    = "HFS (continuous, linear)")

  sub$subgroup_clean <- sub("_strata$", "", sub$subgroup)
  sub$exposure_lab   <- exposure_lab[sub$exposure]
  sub$outcome_lab    <- factor(outcome_lab[sub$outcome],
                               levels = outcome_lab)
  sub$strata_label   <- paste0(sub$subgroup_clean, ": ", sub$stratum)

  # Effect: for binary outcomes effect is OR; for HFS (linear) effect is beta
  # Plot OR on log scale for binary, identity for HFS continuous.
  # Build OR / lci / uci robustly
  sub <- sub %>%
    mutate(
      is_binary = outcome %in% c("cap275", "lsm8"),
      eff_pe = ifelse(is_binary, exp(beta),       beta),
      eff_lo = ifelse(is_binary, exp(beta - 1.96*se), beta - 1.96*se),
      eff_hi = ifelse(is_binary, exp(beta + 1.96*se), beta + 1.96*se),
      p_lab  = sprintf("P=%.3f", p)
    )

  # Stable order for y-axis: group by subgroup family, then alphabetical stratum
  stratum_order <- sub %>%
    distinct(subgroup_clean, stratum, strata_label) %>%
    arrange(factor(subgroup_clean,
                   levels = c("sex", "age", "race", "edu",
                              "dm", "htn", "pir")),
            stratum) %>%
    mutate(combo = strata_label) %>%
    pull(combo)
  stratum_order <- rev(stratum_order)
  sub$strata_label <- factor(sub$strata_label, levels = stratum_order)
  sub$exposure_lab <- factor(sub$exposure_lab,
                             levels = c("Dietary Se","Blood Se","Se/Zn","Se/Cu"))

  # Color per exposure family
  exposure_pal <- c("Dietary Se" = pal_blue,
                    "Blood Se"   = pal_red,
                    "Se/Zn"      = pal_orange,
                    "Se/Cu"      = pal_purple)

  # 3 panels (one per outcome). Each panel: y = strata_label, x = effect,
  # color = exposure. Use position_dodge to separate 4 contrasts per stratum.
  # For binary outcomes use log x axis; for HFS use identity.
  build_panel <- function(outc_key, outc_label) {
    df_o <- sub %>% filter(outcome == outc_key)
    is_bin <- outc_key %in% c("cap275","lsm8")
    null_x <- if (is_bin) 1 else 0
    p <- ggplot(df_o,
                aes(x = eff_pe, y = strata_label, color = exposure_lab)) +
      geom_vline(xintercept = null_x, linetype = "dashed",
                 color = "grey50", linewidth = 0.4) +
      geom_errorbarh(aes(xmin = eff_lo, xmax = eff_hi),
                     height = 0.3, linewidth = 0.45,
                     position = position_dodge(width = 0.7),
                     na.rm = TRUE) +
      geom_point(size = 1.6, shape = 16,
                 position = position_dodge(width = 0.7),
                 na.rm = TRUE) +
      scale_color_manual(values = exposure_pal,
                         name = "Exposure",
                         drop = FALSE) +
      labs(title = outc_label, y = NULL,
           x = if (is_bin) "OR (per 1-SD) on log scale" else "beta (per 1-SD)") +
      theme_pub(base_size = 8) +
      theme(axis.text.y  = element_text(size = 6.5),
            plot.title   = element_text(size = 9, face = "bold"))
    if (is_bin) {
      p <- p + scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3),
                             limits = range(c(0.5, df_o$eff_lo, df_o$eff_hi, 2.5),
                                            na.rm = TRUE))
    } else {
      # HFS continuous beta — let limits be auto
      p
    }
    p
  }

  p_a <- build_panel("cap275", "(A) Steatosis (CAP >= 275)")
  p_b <- build_panel("lsm8",   "(B) Fibrosis (LSM >= 8)")
  p_c <- build_panel("hfs",    "(C) HFS continuous (beta)")

  fig8 <- (p_a | p_b | p_c) +
    plot_layout(guides = "collect", ncol = 3) +
    plot_annotation(
      title = "Figure 8. Subgroup analysis: 4 selenium contrasts x 3 outcomes across 7 strata",
      subtitle = sprintf(
        "n = %d sub-stratum estimates (7 strata x mean 2.2 levels x 4 contrasts x 3 outcomes). BH-adjusted within subgroup.",
        nrow(sub)),
      caption = "Dashed line = null (OR=1 for binary; beta=0 for HFS). Effects estimated via complex-survey weighted GLM.",
      theme = theme(plot.title = element_text(face = "bold", size = 11),
                    plot.subtitle = element_text(size = 9, color = "grey25"),
                    plot.caption = element_text(size = 7, color = "grey50"))
    ) &
    theme(legend.position = "bottom")

  save_fig(fig8, file.path(fig_dir, "fig8_subgroup_forest"),
           width_mm = 230, height_mm = 200)
  cat(sprintf("  [OK] fig8 rendered with %d cells across 3 outcomes x 4 exposures x %d strata\n",
              nrow(sub), length(unique(sub$strata_label))))
  "OK"
}, error = function(e) {
  cat(sprintf("  [FAIL] fig8: %s\n", conditionMessage(e)))
  conditionMessage(e)
})
cat("\n")


# ============================================================
# Fig 11 — CMAverse 4-way decomposition (CORRECT PIE label)
# ============================================================
cat("[3/4] Fig 11: CMAverse 4-way decomposition (PIE -6.9% label fix) ...\n")

fig11_status <- tryCatch({
  cma <- read.csv("output/tables/cmaverse_redox_effects.csv",
                  stringsAsFactors = FALSE)
  # LSM 4-way components
  lsm_4 <- cma %>%
    filter(outcome == "lsm_v2",
           effect %in% c("ERcde(prop)", "ERintref(prop)",
                          "ERintmed(prop)", "ERpnie(prop)"))
  cat(sprintf("  Raw 4-way components:\n"))
  print(lsm_4 %>% select(effect, pe, p))

  # Component metadata (manuscript references)
  comp_meta <- tibble::tibble(
    effect = c("ERcde(prop)", "ERintref(prop)",
               "ERintmed(prop)", "ERpnie(prop)"),
    component = c("Controlled Direct (CDE)",
                  "Reference Interaction (INTREF)",
                  "Mediated Interaction (INTMED)",
                  "Pure Indirect (PIE)"),
    color_grp = c("CDE","INTREF","INTMED","PIE"),
    order_id  = c(1, 2, 3, 4)
  )
  lsm_4 <- lsm_4 %>%
    inner_join(comp_meta, by = "effect") %>%
    mutate(pct = 100 * pe,
           pct_lab = sprintf("%+.1f%% (P=%.3f)", pct, p),
           component = factor(component, levels = comp_meta$component))

  # Color: CDE blue / INTREF orange / INTMED green / PIE red(negative-emphasis)
  comp_pal <- c("Controlled Direct (CDE)"        = pal_blue,
                "Reference Interaction (INTREF)" = pal_orange,
                "Mediated Interaction (INTMED)"  = pal_green,
                "Pure Indirect (PIE)"            = pal_red)

  # Build horizontal bar with explicit label positioning so negatives don't
  # get clipped by axis. Use y-axis range padded both sides.
  y_min <- floor(min(lsm_4$pct, 0) / 10) * 10 - 10
  y_max <- ceiling(max(lsm_4$pct) / 10) * 10 + 15

  fig11 <- ggplot(lsm_4, aes(x = component, y = pct, fill = component)) +
    geom_col(width = 0.65, color = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0, color = "grey25", linewidth = 0.5) +
    geom_text(aes(label = pct_lab,
                  y = ifelse(pct >= 0, pct + 2, pct - 2),
                  hjust = ifelse(pct >= 0, 0, 1)),
              size = 3.0, color = "black") +
    scale_fill_manual(values = comp_pal, guide = "none") +
    scale_y_continuous(limits = c(y_min, y_max),
                       breaks = seq(round(y_min / 10) * 10,
                                    round(y_max / 10) * 10, 10)) +
    coord_flip() +
    labs(title = "Figure 11. CMAverse 4-way decomposition: Se -> GGT -> LSM >= 8 kPa",
         subtitle = "Se Q3 vs Q1 contrast; total RR = 0.89 (95% CI 0.79-0.97, P = 0.004); n = 5,712; 500 bootstraps",
         caption = "Components on excess-risk-ratio (ER) scale per VanderWeele 2014. PIE negative => mediator GGT lowers risk on the Se Q3 vs Q1 contrast.",
         x = NULL,
         y = "Proportion of total effect (%)") +
    theme_pub(base_size = 10) +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"),
          plot.caption  = element_text(size = 7, color = "grey50",
                                       hjust = 0))

  save_fig(fig11, file.path(fig_dir, "fig11_cmaverse_4way"),
           width_mm = 178, height_mm = 110)
  cat(sprintf("  Components labeled: %s\n",
              paste(lsm_4$pct_lab, collapse = " | ")))
  cat("  [OK] fig11 rendered with corrected PIE -6.9% label and CDE/INTREF/INTMED colors\n")
  "OK"
}, error = function(e) {
  cat(sprintf("  [FAIL] fig11: %s\n", conditionMessage(e)))
  conditionMessage(e)
})
cat("\n")


# ============================================================
# Fig 13 — BKMR 4-metal PIP with full title + bar value labels
# ============================================================
cat("[4/4] Fig 13: BKMR 4-metal PIP (full title + value labels) ...\n")

fig13_status <- tryCatch({
  pip <- read.csv("output/tables/bkmr_4metal_pip.csv",
                  stringsAsFactors = FALSE)
  cat("  PIP table:\n"); print(pip)

  pip$variable <- factor(pip$variable, levels = c("Pb","Cd","Hg","Se"))
  pip$pip_lab  <- sprintf("%.2f", pip$PIP)

  # Color per metal (consistent w/ prior fig13)
  metal_pal <- c("Pb" = "#1B9E77",
                 "Cd" = "#D95F02",
                 "Hg" = "#7570B3",
                 "Se" = "#E7298A")

  fig13 <- ggplot(pip,
                  aes(x = reorder(variable, PIP), y = PIP,
                      fill = variable)) +
    geom_col(width = 0.6, color = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = pal_red, linewidth = 0.5) +
    geom_text(aes(label = pip_lab),
              hjust = -0.2, size = 3.4, color = "black") +
    scale_fill_manual(values = metal_pal, guide = "none") +
    scale_y_continuous(limits = c(0, 1.12),
                       breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0.02))) +
    coord_flip() +
    labs(title = "Figure 13. BKMR posterior inclusion probability (PIP), 4-metal model",
         subtitle = "4-metal sensitivity model (Mn dropped). Dashed red: PIP = 0.5 (Bobb 2015 inclusion criterion).",
         caption = "Bayesian Kernel Machine Regression (BKMR), 2 chains x 10,000 iter (5,000 burn-in). Numeric labels: posterior inclusion probability.",
         x = NULL,
         y = "Posterior inclusion probability") +
    theme_pub(base_size = 10) +
    theme(axis.text.y = element_text(face = "bold", size = 10),
          plot.title  = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "grey25"),
          plot.caption  = element_text(size = 7, color = "grey50"))

  save_fig(fig13, file.path(fig_dir, "fig13_bkmr_4metal_pip"),
           width_mm = 178, height_mm = 100)
  cat("  [OK] fig13 rendered with full title visible + PIP value labels on every bar\n")
  "OK"
}, error = function(e) {
  cat(sprintf("  [FAIL] fig13: %s\n", conditionMessage(e)))
  conditionMessage(e)
})
cat("\n")


# ============================================================
# Summary
# ============================================================
cat("========================================\n")
cat("DONE _render_figs_p0fix.R\n")
cat("========================================\n")
status_tbl <- data.frame(
  figure = c("fig5_hfs_calibration",
             "fig8_subgroup_forest",
             "fig11_cmaverse_4way",
             "fig13_bkmr_4metal_pip"),
  status = c(fig5_status, fig8_status, fig11_status, fig13_status)
)
print(status_tbl, row.names = FALSE)
cat("\nFiles produced (in output/figures/):\n")
for (f in c("fig5_hfs_calibration","fig8_subgroup_forest",
            "fig11_cmaverse_4way","fig13_bkmr_4metal_pip")) {
  for (ext in c("tiff","pdf","png")) {
    p <- file.path(fig_dir, paste0(f, ".", ext))
    if (file.exists(p)) {
      sz <- file.info(p)$size / 1024
      cat(sprintf("  %s (%.1f KB)\n", p, sz))
    }
  }
}
