# =====================================================================
# 006 / _render_fig2_dag.R — DAG (Fig 2)
#
# Per _dag_spec.md:
#   A1 = dietary Se (DR1TSELE); A2 = serum Se (LBXBSE);
#   M = HOMA-IR / hsCRP / NLR (insulin/inflammation mediators);
#   Y = LSM/CAP/HFS;
#   L1 = post-exposure mediator-confounders (BMI, T2D, HbA1c, LDL);
#   C = pre-exposure confounders (Age, Sex, Race, Edu, PIR, Smoke, Drink, kcal);
#   U = unmeasured (SELENOP, GPx3, occupational, microbiome).
# =====================================================================

set.seed(20260516)

suppressPackageStartupMessages({
  library(ggdag); library(ggplot2); library(dplyr); library(dagitty)
})
source("../../templates/_shared/figures_helpers.R")

cat("======================================== \n")
cat("006 / Fig 2 DAG (ggdag)\n")
cat("======================================== \n\n")

dag <- dagitty::dagitty('
dag {
  A1 [exposure, pos="2,2"]
  A2 [exposure, pos="2,2.7"]
  M  [pos="3.5,2"]
  Y  [outcome, pos="5,2"]
  C  [pos="0.5,2.7"]
  L1 [pos="3.5,3"]
  U  [pos="3.5,1"]

  C -> A1; C -> A2; C -> L1; C -> M; C -> Y
  A1 -> A2; A1 -> M; A1 -> Y; A1 -> L1
  A2 -> M; A2 -> Y; A2 -> L1
  L1 -> M; L1 -> Y
  M -> Y
  U -> A1; U -> A2; U -> M; U -> Y
}
')

tidy_dag <- ggdag::tidy_dagitty(dag)
tidy_dag$data$node_type <- dplyr::case_when(
  tidy_dag$data$name %in% c("A1", "A2") ~ "Exposure",
  tidy_dag$data$name == "Y"             ~ "Outcome",
  tidy_dag$data$name == "M"             ~ "Mediator",
  tidy_dag$data$name == "L1"            ~ "L1 (post-exposure)",
  tidy_dag$data$name == "U"             ~ "Unmeasured",
  TRUE                                  ~ "Confounder C"
)

palette_nodes <- c(
  "Exposure"            = "#E74C3C",  # red
  "Outcome"             = "#2E86AB",  # blue (W11 R2 P0: must NOT collide with Exposure)
  "Mediator"            = "#F4D03F",  # yellow
  "L1 (post-exposure)"  = "#85C1E9",  # light blue
  "Confounder C"        = "#FFFFFF",  # white
  "Unmeasured"          = "#999999"   # grey
)

tidy_dag$data$label_display <- dplyr::case_when(
  tidy_dag$data$name == "A1" ~ "Dietary Se",
  tidy_dag$data$name == "A2" ~ "Serum Se",
  tidy_dag$data$name == "M"  ~ "HOMA-IR/CRP/NLR",
  tidy_dag$data$name == "Y"  ~ "LSM/CAP/HFS",
  tidy_dag$data$name == "L1" ~ "L1 (BMI/T2D/HbA1c/LDL)",
  tidy_dag$data$name == "C"  ~ "C (Age/Sex/Race/Edu/PIR/Smoke/Drink)",
  tidy_dag$data$name == "U"  ~ "U (SELENOP/GPx3)",
  TRUE ~ tidy_dag$data$name
)

p_dag <- ggplot(tidy_dag$data,
                aes(x = x, y = y, xend = xend, yend = yend)) +
  ggdag::geom_dag_edges(edge_colour = "grey30", edge_width = 0.4,
                        arrow_directed = grid::arrow(length = grid::unit(2, "mm"),
                                                     type = "closed")) +
  ggdag::geom_dag_point(aes(fill = node_type), shape = 21,
                        size = 14, stroke = 0.4, colour = "black") +
  ggdag::geom_dag_text(aes(label = label_display),
                       colour = "black", size = 2.2, fontface = "bold") +
  scale_fill_manual(values = palette_nodes, name = NULL) +
  coord_cartesian(xlim = c(-0.3, 5.5), ylim = c(0.3, 3.5)) +
  labs(title = "Fig 2. DAG for dual-source Se (dietary + serum) on hepatic steatosis/stiffness/HFS",
       subtitle = "L1 are post-exposure mediator-confounders (VanderWeele 2014); U bounded by E-value sensitivity.") +
  nhanes_theme_publication(base_size = 9, base_family = "sans") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
if (!dir.exists("output/figures/_snapshot")) dir.create("output/figures/_snapshot", recursive = TRUE)

ggsave("output/figures/fig2_dag.tiff", p_dag,
       width = 178, height = 110, units = "mm", dpi = 600,
       device = grDevices::tiff, compression = "lzw")
ggsave("output/figures/fig2_dag.pdf", p_dag,
       width = 178, height = 110, units = "mm")

save(tidy_dag, p_dag, file = "output/figures/_snapshot/fig2.RData")

cat(sprintf("Fig 2 TIFF size: %d bytes\n",
            file.size("output/figures/fig2_dag.tiff")))
cat("\nDONE _render_fig2_dag.R\n")
