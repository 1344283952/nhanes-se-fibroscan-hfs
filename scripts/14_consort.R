# ============================================
# 006 / 14_consort.R — CONSORT-style flowchart (Fig 1)
#
# Uses templates/_shared/figures_helpers.R::consort_from_flow()
# ============================================

set.seed(20260516)
suppressPackageStartupMessages({
  library(dplyr); library(DiagrammeR); library(DiagrammeRsvg); library(rsvg)
})
source("../../templates/_shared/figures_helpers.R")

cat("========================================\n")
cat("006 / Fig 1 CONSORT flowchart\n")
cat("========================================\n\n")

flow_csv <- "output/tables/flow_counts.csv"
if (!file.exists(flow_csv)) stop("flow_counts.csv not found")
flow <- utils::read.csv(flow_csv, stringsAsFactors = FALSE)
if (!"label" %in% names(flow)) flow$label <- flow$step
if (!"step" %in% names(flow)) flow$step <- paste0("s", seq_len(nrow(flow)))
flow$step <- gsub("[^A-Za-z0-9_]", "_", flow$step)
flow$step <- ifelse(grepl("^[0-9]", flow$step), paste0("s_", flow$step), flow$step)

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
write.csv(flow, "output/tables/flow_counts_for_consort.csv", row.names = FALSE)
dot_str <- consort_from_flow("output/tables/flow_counts_for_consort.csv",
                              title = "Figure 1. Sample selection flow (NHANES 2017-March 2020 Pre-pandemic)")

g <- DiagrammeR::grViz(dot_str)
svg_str <- DiagrammeRsvg::export_svg(g)
writeLines(svg_str, "output/figures/fig1_consort.svg")
rsvg::rsvg_png("output/figures/fig1_consort.svg",
               "output/figures/fig1_consort.png",
               width = 1200, height = 900)

cat("\n保存:\n")
cat("  output/figures/fig1_consort.svg + .png\n")
cat("\nDONE 14_consort.R\n")
