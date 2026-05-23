# ============================================
# 006v2 / _render_v2_figures.R — render new figures for v2 sections
# Output: output/figures/fig9..fig13 (qgcomp / WQS / Hg-Se / CMAverse / K cycle)
# ============================================
library(ggplot2); library(dplyr); library(patchwork); library(tidyr)
fig_dir <- "output/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(20260518)

theme_pub <- function() theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 9))

# === Fig 9: qgcomp + WQS combined ===
cat("==== Fig 9: qgcomp + WQS combined ====\n")
qg <- read.csv("output/tables/qgcomp_mixture.csv", stringsAsFactors = FALSE)
qg$label <- factor(qg$outcome, levels = qg$outcome)
qg$or_pe <- ifelse(qg$type == "binary", exp(qg$psi), qg$psi)
qg$or_lo <- ifelse(qg$type == "binary", exp(qg$ci_lo), qg$ci_lo)
qg$or_hi <- ifelse(qg$type == "binary", exp(qg$ci_hi), qg$ci_hi)

p9a <- ggplot(qg, aes(x = or_pe, y = reorder(label, or_pe), color = factor(sign(psi)))) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0.25) +
  geom_vline(xintercept = ifelse(qg$type[1] == "binary", 1, 0),
             linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("-1" = "#2E86AB", "1" = "#E74C3C"),
                     labels = c("-1" = "Negative", "1" = "Positive"),
                     name = "Direction") +
  labs(title = "(A) qgcomp 5-metal mixture (ψ effect / OR per quartile-shift)",
       x = "Effect estimate", y = NULL) +
  theme_pub()

wqs <- read.csv("output/tables/wqs_mixture.csv", stringsAsFactors = FALSE)
wqs$or <- exp(wqs$beta)
wqs$or_lo <- exp(wqs$beta - 1.96 * wqs$se)
wqs$or_hi <- exp(wqs$beta + 1.96 * wqs$se)
wqs$label <- paste0(wqs$outcome, " (", wqs$direction, ")")

p9b <- ggplot(wqs, aes(x = or, y = reorder(label, or), color = direction)) +
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("positive" = "#E74C3C", "negative" = "#2E86AB"),
                     name = "Direction") +
  scale_x_log10() +
  labs(title = "(B) WQS regression (OR per WQS-index increment)",
       x = "OR per index unit (log scale)", y = NULL) +
  theme_pub()

fig9 <- p9a / p9b + plot_layout(heights = c(1, 1))
ggsave(file.path(fig_dir, "fig9_mixture_qgcomp_wqs.tiff"), fig9,
       width = 8, height = 7, dpi = 300, compression = "lzw")
ggsave(file.path(fig_dir, "fig9_mixture_qgcomp_wqs.pdf"), fig9,
       width = 8, height = 7)
cat("[OK] fig9 saved\n")

# === Fig 10: Hg-Se 拮抗 (Hg/Se ratio histogram + interaction forest) ===
cat("==== Fig 10: Hg-Se 拮抗 ====\n")
load("data/processed/nhanes_final.RData")
d <- nhanes_final[!is.na(nhanes_final$hg_se_molar_ratio), ]

p10a <- ggplot(d, aes(x = hg_se_molar_ratio)) +
  geom_histogram(bins = 60, fill = "#2E86AB", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  annotate("text", x = 1.05, y = 50, label = "Ralston 2017 threshold (Hg/Se = 1)",
           hjust = 0, color = "#E74C3C", size = 2.7) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10),
                labels = c("10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "1", "10")) +
  labs(title = "(A) Hg/Se molar ratio distribution (NHANES Pre-pandemic, n=5,883)",
       x = "Hg/Se molar ratio (log)", y = "n") +
  theme_pub()

inter <- read.csv("output/tables/hgse_interaction.csv", stringsAsFactors = FALSE)
inter$or_pe <- inter$OR_per_logHgSe
inter <- inter[order(inter$OR_per_logHgSe), ]
p10b <- ggplot(inter, aes(x = or_pe, y = reorder(outcome, or_pe))) +
  geom_point(size = 2.5, color = "#E74C3C") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, color = "#E74C3C") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(x = ci_hi * 1.1, label = sprintf("P=%.3f", p)),
            hjust = 0, size = 2.7) +
  scale_x_log10() +
  labs(title = "(B) Hg × Se interaction (OR per log Hg/Se unit, BH-adj P)",
       x = "OR (log scale)", y = NULL) +
  theme_pub()

fig10 <- p10a / p10b + plot_layout(heights = c(1, 1))
ggsave(file.path(fig_dir, "fig10_hgse_antagonism.tiff"), fig10,
       width = 8, height = 7, dpi = 300, compression = "lzw")
ggsave(file.path(fig_dir, "fig10_hgse_antagonism.pdf"), fig10,
       width = 8, height = 7)
cat("[OK] fig10 saved\n")

# === Fig 11: CMAverse 4-way decomposition (LSM) ===
cat("==== Fig 11: CMAverse 4-way ====\n")
cma <- read.csv("output/tables/cmaverse_redox_effects.csv", stringsAsFactors = FALSE)
# Extract LSM 4-way components
lsm_4 <- cma[cma$outcome == "lsm_v2" & cma$effect %in% c("ERcde(prop)","ERintref(prop)","ERintmed(prop)","ERpnie(prop)"), ]
if (nrow(lsm_4) > 0) {
  lsm_4$component <- recode(lsm_4$effect,
    "ERcde(prop)" = "Controlled Direct\n(CDE)",
    "ERintref(prop)" = "Reference\nInteraction (INTREF)",
    "ERintmed(prop)" = "Mediated\nInteraction (INTMED)",
    "ERpnie(prop)" = "Pure Indirect\n(PIE)")
  lsm_4$pct <- 100 * lsm_4$pe
  lsm_4$pe_label <- sprintf("%.1f%% (P=%.3f)", lsm_4$pct, lsm_4$p)
  fig11 <- ggplot(lsm_4, aes(x = reorder(component, pct), y = pct, fill = pct > 0)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_text(aes(label = pe_label, y = ifelse(pct > 0, pct + 5, pct - 5)),
              size = 2.7) +
    scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#2E86AB"),
                      guide = "none") +
    coord_flip() +
    labs(title = "Figure 11. CMAverse 4-way decomposition: Se → GGT → LSM ≥ 8 kPa",
         subtitle = sprintf("Se Q3 vs Q1 contrast; Total RR=0.89 (0.79-0.97, P=0.048); n=5,712; 500 boot"),
         x = NULL,
         y = "Proportion of total effect (%)") +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 8, color = "grey30"))
  ggsave(file.path(fig_dir, "fig11_cmaverse_4way.tiff"), fig11,
         width = 8, height = 5, dpi = 300, compression = "lzw")
  ggsave(file.path(fig_dir, "fig11_cmaverse_4way.pdf"), fig11,
         width = 8, height = 5)
  cat("[OK] fig11 saved\n")
}

# === Fig 12: NHANES K cycle external validation comparison ===
# R4 P0 fix (2026-05-23): use design-weighted K-cycle numbers from
# output/tables/k_cycle_validation.csv column K_cycle_weighted_NEW
# (svyglm/svyquantile under WTMEC2YR, df=15)
cat("==== Fig 12: K cycle vs P_ comparison (weighted) ====\n")
kval <- data.frame(
  Test = c("U-shape\n(P-non-lin)", "Se/Zn T3 vs T1\n(OR for CAP≥275)",
           "HFS AUROC\n(vs LSM≥8)", "Hg/Se molar\nratio median"),
  P_cycle = c(1.27e-11, 1.22, 0.731, 0.001),
  P_cycle_lo = c(NA, 1.10, 0.699, NA),
  P_cycle_hi = c(NA, 1.35, 0.762, NA),
  K_cycle = c(2.77e-5, 1.08, 0.656, 0.001),
  K_cycle_lo = c(NA, 0.78, 0.606, NA),
  K_cycle_hi = c(NA, 1.49, 0.703, NA),
  stringsAsFactors = FALSE
)

# Use shape difference (P-value vs OR vs AUROC vs ratio) — make 4 small panels
make_compare <- function(label, p_val, k_val, p_lo = NA, p_hi = NA, k_lo = NA, k_hi = NA,
                         ylabel) {
  df <- data.frame(cycle = c("P_ (2017-2020)", "K (2021-2023)"),
                   est = c(p_val, k_val),
                   lo  = c(p_lo, k_lo),
                   hi  = c(p_hi, k_hi))
  ggplot(df, aes(x = cycle, y = est, color = cycle)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, na.rm = TRUE) +
    scale_color_manual(values = c("P_ (2017-2020)" = "#E74C3C",
                                  "K (2021-2023)" = "#2E86AB"),
                      guide = "none") +
    labs(title = label, x = NULL, y = ylabel) +
    theme_pub() +
    theme(axis.text.x = element_text(size = 8))
}

p12a <- make_compare("(A) U-shape RCS P-non-linear", 1.27e-11, 2.77e-5,
                    ylabel = "P-value") +
  scale_y_log10()
p12b <- make_compare("(B) Se/Zn T3 vs T1 OR (CAP≥275)",
                    1.22, 1.08, 1.10, 1.35, 0.78, 1.49,
                    ylabel = "OR") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50")
p12c <- make_compare("(C) HFS AUROC (vs LSM≥8)",
                    0.731, 0.656, 0.699, 0.762, 0.606, 0.703,
                    ylabel = "AUROC") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50")
p12d <- make_compare("(D) Hg/Se molar ratio (median)",
                    0.001, 0.001,
                    ylabel = "Hg/Se")

fig12 <- (p12a | p12b) / (p12c | p12d) +
  plot_annotation(title = "Figure 12. External validation: NHANES P_ (2017-2020) vs K cycle (2021-2023, design-weighted)",
                  subtitle = "P_ n=5,885; K n=3,310; HFS subset P_ n=2,810 / K n=1,835; K cycle WTMEC2YR svyglm/svyquantile, df=15",
                  theme = theme(plot.title = element_text(face = "bold", size = 11),
                                plot.subtitle = element_text(size = 8, color = "grey30")))
ggsave(file.path(fig_dir, "fig12_k_cycle_validation.tiff"), fig12,
       width = 9, height = 7, dpi = 300, compression = "lzw")
ggsave(file.path(fig_dir, "fig12_k_cycle_validation.pdf"), fig12,
       width = 9, height = 7)
ggsave(file.path(fig_dir, "fig12_k_cycle_validation.png"), fig12,
       width = 9, height = 7, dpi = 150)
cat("[OK] fig12 saved\n")

cat("\n==== v2 figures all rendered to output/figures/fig9-12 ====\n")
