# ============================================
# 006 / 09_dag.R — DAG render (Fig 2)
#
# Round 4 R-Figure-Planning P0:
#   - 之前 _dag_spec.md 描述了 DAG 但没有 R 渲染脚本; Fig 2 不可重现
#   - 本脚本对照 _dag_spec.md §3 ggdag 代码, 输出 publication-quality TIFF/PDF
#   - 使用 templates/_shared/figures_helpers.R 的 nhanes_theme_publication
#     + palette_dag_nodes + save_publication_figure(type="double", dpi=600)
#
# 渲染依赖:
#   - ggdag (DAG layout)
#   - dagitty (DAG declaration)
#   - figures_helpers.R (publication theme + TIFF LZW)
#
# 输出:
#   - output/figures/fig2_dag.tiff (双栏 178mm × 130mm, 600 dpi)
#   - output/figures/fig2_dag.pdf  (矢量预览)
# ============================================

set.seed(20260516)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggdag)
  library(dagitty)
  library(dplyr)
})
source("../../templates/_shared/figures_helpers.R")

cat("========================================\n")
cat("006 / 09_dag — Fig 2 DAG (dual-Se → CAP/LSM/HFS)\n")
cat("========================================\n\n")

# ---- DAG declaration (mirrors _dag_spec.md §3) ----
# Round 2 R-Figure P0: spread positions to eliminate text overlap;
# exposure (A1, A2) and outcome (Y) MUST use distinct colors (R-Figure said
# both were red in v1.0 → color-blind hazard). New palette: red=exposure,
# blue=outcome, yellow=mediator, white=confounder, grey=latent.
g <- dagitty('dag {
  C  [pos="-2.5,  -0.5"]
  A1 [pos="-0.8, -1.3", exposure]
  A2 [pos="-0.8,  0.3", exposure]
  L1 [pos=" 1.3, -1.6"]
  M  [pos=" 1.3,  0.3"]
  Y  [pos=" 3.6, -0.5", outcome]
  U  [pos=" 1.3,  1.8", latent]

  C  -> A1
  C  -> A2
  C  -> L1
  C  -> M
  C  -> Y
  A1 -> A2
  A1 -> M
  A1 -> Y
  A1 -> L1
  A2 -> M
  A2 -> Y
  A2 -> L1
  L1 -> M
  L1 -> Y
  M  -> Y
  U  -> A1
  U  -> A2
  U  -> M
  U  -> Y
}')

dag_df <- tidy_dagitty(g, seed = 20260516)

# ---- Node-type palette (Round 2 P0: exposure ≠ outcome color; color-blind safe) ----
# Override figures_helpers::palette_dag_nodes() locally — keep red for exposure
# but switch outcome to blue to satisfy deuteranopia / R-Figure v2 8/10 floor.
node_pal <- c(
  exposure   = "#E74C3C",  # red (A1, A2 — dual selenium source)
  outcome    = "#2E86AB",  # blue (Y — distinct from exposure)
  mediator   = "#F4D03F",  # yellow (M)
  confounder = "#FFFFFF",  # white with grey border (C, L1)
  latent     = "#BDC3C7"   # light grey (U)
)
node_role <- c(
  A1 = "exposure",
  A2 = "exposure",
  Y  = "outcome",
  M  = "mediator",
  L1 = "confounder",
  C  = "confounder",
  U  = "latent"
)
# Round 2 R-Figure P0: shorter labels to avoid overlap on the v1.0 layout.
node_lbl <- c(
  A1 = "A1\nDietary Se",
  A2 = "A2\nBlood Se",
  Y  = "Y\nCAP/LSM/HFS",
  M  = "M\nHOMA-IR\nhsCRP/NLR",
  L1 = "L1\nBMI/T2D\nHbA1c/LDL\nSBP/DBP",
  C  = "C\nAge/Sex/Race\nEdu/PIR\nSmoke/Drink\nkcal",
  U  = "U\nSELENOP\nGPx3"
)

dag_df$data$role <- node_role[dag_df$data$name]
dag_df$data$lbl  <- node_lbl[dag_df$data$name]

# ---- Plot ----
p_dag <- ggplot(dag_df$data,
                aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_colour = "grey25", edge_width = 0.4,
                 arrow_directed = grid::arrow(length = grid::unit(2.5, "mm"),
                                              type = "closed")) +
  geom_dag_point(aes(fill = role),
                 shape = 21, size = 20, stroke = 0.5, colour = "grey20") +
  geom_dag_text(aes(label = lbl), size = 2.4, lineheight = 0.85,
                family = "sans", colour = "grey15") +
  scale_fill_manual(values = node_pal,
                    breaks = c("exposure", "outcome", "mediator",
                               "confounder", "latent"),
                    labels = c("Exposure (A1, A2; red)",
                               "Outcome (Y; blue)",
                               "Mediator (M; yellow)",
                               "Confounder (C, L1; white)",
                               "Latent (U; grey)"),
                    name = NULL) +
  theme_dag() +
  nhanes_theme_publication(base_size = 9) +
  theme(axis.line  = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "bottom",
        legend.key.size = grid::unit(4, "mm")) +
  labs(
    title = "Figure 2. Directed acyclic graph: dual-source selenium → CAP / LSM / Hepamet HFS",
    subtitle = "Pre-registered (OSF v1.0, 2026-05-16). Dual-source A1 (dietary, DR1TSELE) + A2 (whole blood, LBXBSE) with tensor interaction ti(A1, A2); A_ratio (dietary Se/Zn, Se/Cu) entered separately. NHANES does not measure SELENOP or GPx3 (U) — limitation flagged in Discussion.",
    caption = "References: Pearl 2009 §3.2; VanderWeele 2014 DOI 10.1097/EDE.0000000000000121; Rayman 2012 DOI 10.1016/S0140-6736(11)61452-9; Rinella 2023 DOI 10.1097/HEP.0000000000000520."
  )

# ---- Save ----
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
save_publication_figure(
  plot      = p_dag,
  path      = "output/figures/fig2_dag.tiff",
  type      = "double",
  height_mm = 130,
  format    = "tiff",
  dpi       = 600
)

save_figure_snapshot(
  plot       = p_dag,
  fig_id     = "fig2_dag",
  data_paths = character(0),
  seed       = 20260516
)

cat("\nDONE 09_dag.R\n")
