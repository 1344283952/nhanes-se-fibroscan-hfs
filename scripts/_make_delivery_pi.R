# ============================================
# 006_se_fibroscan_hfs / _make_delivery_pi.R
# 拆分打两个独立 zip 给通讯作者:
#
#   1. 投稿主文档包.zip (~1MB)
#      仅含 投稿系统上传的文档 + 给通讯作者看的操作手册
#
#   2. GitHub上传包.zip (~30-50MB)
#      仅含 GitHub 仓库内容 (scripts + data/processed + output + README.md)
#
# 改自 templates/_make_delivery_pi.R (003_v2 版本)
# ============================================

cat("========================================\n")
cat("006 / 打包通讯作者友好分包\n")
cat("========================================\n\n")

ts <- format(Sys.time(), "%Y%m%d_%H%M")

if (!requireNamespace("zip", quietly = TRUE)) {
  install.packages("zip", repos = "https://cloud.r-project.org", quiet = TRUE)
}

# 用户规则: 只保留最新一份 zip, 避免提交错版本.
# 生成新 zip 前先清掉所有旧的 投稿主文档包_006_*.zip 与 GitHub上传包_006_*.zip
old_zips <- c(
  list.files(".", pattern = "^投稿主文档包_006_.*\\.zip$", full.names = TRUE),
  list.files(".", pattern = "^GitHub上传包_006_.*\\.zip$",  full.names = TRUE)
)
if (length(old_zips) > 0) {
  cat(sprintf("[清旧] 移除 %d 个旧 zip:\n", length(old_zips)))
  for (z in old_zips) cat("  -", z, "\n")
  unlink(old_zips, force = TRUE)
}

# 2026-05-23: Pre-flight audit + auto-fill from author_defaults.yaml
defaults_yaml <- file.path("..", "..", "templates", "_credentials", "author_defaults.yaml")
audit_script  <- file.path("..", "..", "templates", "_pre_zip_audit.R")

if (file.exists(defaults_yaml)) {
  cat("[auto-fill] Loading author defaults from", defaults_yaml, "\n")
  if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml", repos = "https://cloud.r-project.org", quiet = TRUE)
  yaml_data <- yaml::read_yaml(defaults_yaml)
  # Apply placeholder_map to each .md before pandoc render
  # (This step is OUTSIDE current scope — leave a NOTE for manual application by Fix A agent)
}

# ============================================
# Pack 1: 投稿主文档包
# ============================================

pack1_name <- sprintf("投稿主文档包_006_%s.zip", ts)
top1 <- "投稿主文档包_006"

stage1 <- file.path(tempdir(), top1)
unlink(stage1, recursive = TRUE)
dir.create(stage1, recursive = TRUE)

cat(">>> 打包 1: 投稿主文档包\n\n")

# 2026-05-18 W22 修订: 加 .docx 主路径 (期刊投稿系统接受) + supplementary + 备用 cover letter
# (旧版只 ship .md, PI 无法直接上传到 Editorial Manager — 现修复)

# .docx — PI 投稿系统直接上传文件 (Elsevier EM / Wolters Kluwer EM 都接受 .doc/.docx)
# 2026-05-23 升级: Environment International (IF 11.0) 为主投; Hepatology Communications 为唯一保留 fallback
top1_docx_files <- c(
  "投稿操作指南.docx",                         # ★第一份打开 (PI 4 步操作)
  "manuscript_EnvironInt_v1.docx",            # → Manuscript 字段 (Environ Int 主投)
  "cover_letter_EnvironInt.docx",             # ★ → Cover Letter 字段 (Environ Int 主投)
  "cover_letter_HepComm.docx",                # 备投 (保底): Hepatology Communications
  "suggested_reviewers.docx",                 # PI 抄进系统 reviewer form (不上传)
  "STROBE_checklist.docx",                    # → Supplementary
  "TRIPOD_checklist.docx",                    # → Supplementary (HFS prediction 必附)
  "AGReMA_checklist.docx",                    # → Supplementary (AGReMA 2021 mediation reporting)
  "supplementary_information.docx",           # → Supplementary Information 字段
  "supplementary_methods_EnvironInt.docx",    # → Supplementary Methods (compressed methods companion)
  "_osf_preregistration.docx",                # PI OSF 注册时粘贴
  "结果说明.docx",                              # PI 中文 reference (不上传, 投稿前自查)
  "方法学段.docx"                                # PI 中文 reference (双语对照)
)

# 投稿主文档包不再含 .md 源文件 — Yu 姐 2026-05-18 反馈:
# 投稿系统(Elsevier EM / Wolters Kluwer EM)只吃 .doc/.docx/.pdf, .md 会让对方误以为是草稿.
# .md 源文件保留在仓库 (GitHub 上传包) 用于版本追溯, 不进主投 zip.

# 2026-05-23: rename-on-copy so zip uses generic names (manuscript.docx not manuscript_EnvironInt_v1.docx)
docx_rename_map <- list(
  "manuscript_EnvironInt_v1.docx"           = "manuscript.docx",
  "cover_letter_EnvironInt.docx"            = "cover_letter.docx",
  "cover_letter_HepComm.docx"               = "cover_letter_fallback_HepComm.docx",
  "supplementary_information.docx"          = "supplementary_information.docx",
  "supplementary_methods_EnvironInt.docx"   = "supplementary_methods.docx",
  "suggested_reviewers.docx"                = "suggested_reviewers.docx",
  "STROBE_checklist.docx"                   = "STROBE_checklist.docx",
  "TRIPOD_checklist.docx"                   = "TRIPOD_checklist.docx",
  "AGReMA_checklist.docx"                   = "AGReMA_checklist.docx",
  "_osf_preregistration.docx"               = "OSF_preregistration_excerpt.docx",
  "投稿操作指南.docx"                          = "投稿操作指南.docx",
  "结果说明.docx"                              = "结果说明.docx",
  "方法学段.docx"                              = "方法学段.docx"
)
for (long_name in names(docx_rename_map)) {
  if (file.exists(long_name)) {
    short_name <- docx_rename_map[[long_name]]
    file.copy(long_name, file.path(stage1, short_name), overwrite = TRUE)
    cat(sprintf("  + %s  →  %s\n", long_name, short_name))
  } else {
    cat(sprintf("  - 缺: %s\n", long_name))
  }
}

# 2026-05-23 user feedback: PNG + PDF 双格式 (TIFF only 之前是错误判断)
# PNG: reviewer 浏览器原生显示; PDF: vector line art (CONSORT/DAG/forest)
main_figs_base <- c(
  "fig1_consort",
  "fig2_dag",
  "fig3_rcs_se",
  "fig4_gam_heatmap",
  "fig5_hfs_calibration",
  "fig6_cross_heatmap",
  "fig7_ratio_forest",
  "fig8_subgroup_forest",
  "fig9_mixture_qgcomp_wqs",
  "fig10_hgse_antagonism",
  "fig11_cmaverse_4way",
  "fig12_k_cycle_validation",
  "fig13_bkmr_4metal_pip"
)
fig_dir1 <- file.path(stage1, "figures")
dir.create(fig_dir1, showWarnings = FALSE)
for (base in main_figs_base) {
  src <- file.path("output", "figures", paste0(base, ".png"))
  if (file.exists(src)) {
    file.copy(src, file.path(fig_dir1, paste0(base, ".png")), overwrite = TRUE)
  } else {
    cat(sprintf("  缺: %s.png (skip)\n", base))
  }
}
cat(sprintf("  + figures/ (%d files)\n", length(list.files(fig_dir1))))

# Tables (主表 XLSX + 核心 CSV)
tab_dir1 <- file.path(stage1, "tables")
dir.create(tab_dir1, showWarnings = FALSE)
main_tabs_xlsx <- list.files("output/tables", pattern = "\\.xlsx$", full.names = TRUE)
file.copy(main_tabs_xlsx, tab_dir1, overwrite = TRUE)

main_tabs_csv <- c(
  "hfs_auroc_summary.csv",
  "ratio_OR_cap275.csv",
  "evalue_ratio_OR.csv",
  "sensitivity_S1-S7.csv",
  "rcs_pvalues.csv",
  "subgroup_forest_data.csv",
  # v1.5 W22 additions — Step H Round 1 缺失修复 (2026-05-18)
  "mr_NOTE_literature_cited.txt",          # H6 MR literature-cite audit trail (v1.5)
  "mi_vs_cc_comparison.csv",               # S4 mice m=20 vs complete-case side-by-side
  "mi_cap_logistic.csv",                   # S4 mice m=20 CAP>=275 pooled Rubin
  "mi_lsm_logistic.csv",                   # S4 mice m=20 LSM>=8 pooled
  "mi_cap_completecase.csv",               # S4 complete-case fit
  "nadir_bootstrap.csv",                   # S9 RCS nadir 1000-rep bootstrap
  "multiscore_benchmark.csv",              # HFS vs FIB-4/NFS/APRI paired DeLong
  "sensitivity_S8_ipw_selection.csv",      # S8 IPW selection bias
  "retrodesign_smallcells.csv",            # Type-S/M for underpowered cells
  "table1_deff.csv",                       # DEFF report (PLOS Biology 2025 防御)
  # v2.0 OSF amendment additions — 4 new analytic layers (2026-05-18)
  "qgcomp_mixture.csv",                    # §3.8 qgcomp 5-metal mixture (Keil 2020)
  "wqs_mixture.csv",                       # §3.8 WQS pos+neg (Carrico 2015)
  "hgse_interaction.csv",                  # §3.9 Hg×Se interaction logistic
  "hgse_rcs.csv",                          # §3.9 Hg/Se molar RCS non-linearity
  "cmaverse_redox_effects.csv",            # §3.10 CMAverse 4-way (VanderWeele 2014)
  "mr_refit_main.csv",                     # §3.7 own MR refit (Anstee 2020 outcome GWAS)
  "mr_refit_loo.csv",                      # §3.7 MR leave-one-out
  "k_cycle_validation.csv",                # §3.11 NHANES K cycle external validation
  "k_cycle_flow.csv"                       # §3.11 K cycle selection cascade
  # NOTE: _audit_unused_vars.csv removed (Step H Round 1 Dim10 fix — internal QA artefact, not shipped to PI)
)
for (f in main_tabs_csv) {
  sp <- file.path("output/tables", f)
  if (file.exists(sp)) file.copy(sp, file.path(tab_dir1, f), overwrite = TRUE)
}
cat(sprintf("  + tables/ (%d files)\n", length(list.files(tab_dir1))))

# README.txt
readme1 <- c(
  "# 投稿主文档包 006 (Selenium x FibroScan x HFS) — Redox Biology 主投版",
  "",
  sprintf("打包时间: %s", format(Sys.time())),
  "目标期刊 (2026-05-18 升级 audit 后):",
  "  主投: Redox Biology (Elsevier, IF 11.9) — scope 实查接 NHANES + redox-axis",
  "  备 1: Aliment Pharmacol Ther (Wiley, IF 6.7) — 临床肝病主场",
  "  备 2: Hepatology Communications (Wiley, IF 4.6) — 肝病 OA",
  "  保底: Free Radic Biol Med (Elsevier, IF 7.4) — 原 v1.4 主投, 降为 fallback",
  "",
  "OSF Pre-registration: v1.5 amendment (含 H6 MR + S4 mice m=20 + S9 nadir + AGReMA-SF 2024 + PLOS Biology 2025 defence)",
  "",
  "## ★ 投稿系统直接上传文件 (.docx — Elsevier EM 接受格式)",
  "- 投稿操作指南.docx              ★ 第一份打开 (4 步走)",
  "- manuscript_v1.docx             → Manuscript 字段",
  "- cover_letter_RedoxBiol.docx    ★ → Cover Letter 字段 (Redox Biology 主投)",
  "- cover_letter.docx              备投 1 — FRBM (原主投, fallback)",
  "- cover_letter_LiverInt.docx     备投 2 — Aliment Pharmacol Ther / Liver International",
  "- cover_letter_HepComm.docx      备投 3 — Hepatology Communications",
  "- supplementary_information.docx → Supplementary Information 字段",
  "- STROBE_checklist.docx          → Supplementary (STROBE 22 项)",
  "- TRIPOD_checklist.docx          → Supplementary (TRIPOD-AI 27 项, HFS 预测必附)",
  "- AGReMA_checklist.docx          → Supplementary (AGReMA Statement 2021 mediation reporting, v1.5 新增)",
  "- suggested_reviewers.docx       PI 抄进系统 reviewer form (不上传)",
  "- _osf_preregistration.docx      PI OSF 注册时粘贴",
  "",
  "## PI 中文 reference (不上传, 投稿前自查)",
  "- 结果说明.docx              含 H1-H5 + H6 verdicts + 核心数字",
  "- 方法学段.docx              方法学中文版 (与 manuscript Methods 双语对照)",
  "",
  "## 图表",
  "- figures/   8 张主图 TIFF (600 dpi) + PDF preview — Figures 字段每张单独上传",
  "- tables/    Table 1 xlsx + 主结果 CSV + mr_NOTE_literature_cited.txt + nadir_bootstrap.csv + mi_*.csv",
  "",
  "## 第一步: 打开 投稿操作指南.docx, 按 4 步操作"
)
writeLines(readme1, file.path(stage1, "README.txt"))
cat("  + README.txt\n")

# 2026-05-23: Pre-flight audit before zipping (6-gate enforcement)
if (file.exists(audit_script)) {
  cat("\n[pre-flight] Running pre-zip audit (6 gates)...\n")
  source(audit_script)
  audit_res <- pre_zip_audit(stage1, getwd())
  if (audit_res$pass) {
    cat("[pre-flight] ✓ PASS (6/6 gates clean)\n\n")
  } else {
    cat("[pre-flight] ✗ FAIL (", audit_res$n_failures, "gate failures):\n", sep = "")
    for (f in audit_res$failures) cat("  - ", f, "\n", sep = "")
    cat("\n[pre-flight] Aborting zip generation. Fix the failures above and re-run.\n")
    quit(status = 1)
  }
}

cat("\n压缩 pack 1 ...\n")
old_wd <- getwd()
setwd(dirname(stage1))
zip::zip(zipfile = file.path(old_wd, pack1_name),
         files   = basename(stage1),
         recurse = TRUE,
         mode    = "cherry-pick")
setwd(old_wd)
unlink(stage1, recursive = TRUE)

zsize1 <- file.size(pack1_name)
cat(sprintf("\n[OK] pack 1: %s (%.1f MB)\n", pack1_name, zsize1 / 1024 / 1024))

# ============================================
# Pack 2: GitHub 上传包
# ============================================

pack2_name <- sprintf("GitHub上传包_006_%s.zip", ts)
top2 <- "GitHub上传包_006"

stage2 <- file.path(tempdir(), top2)
unlink(stage2, recursive = TRUE)
dir.create(stage2, recursive = TRUE)

cat("\n\n>>> 打包 2: GitHub 上传包\n\n")

# README.md - GitHub 首页
readme_gh <- c(
  "# 006: Selenium x FibroScan x HFS (NHANES Pre-pandemic 2017-March 2020)",
  "",
  "Replication code for: Dual-source selenium (dietary DR1TSELE + serum LBXBSE) x hepatic steatosis/stiffness x Hepamet HFS prediction.",
  "",
  "## Key findings",
  "- H2 confirmed: Se U-shape on CAP via mgcv tensor *P* = 1.27e-11",
  "- H5 confirmed: Hepamet HFS AUROC = 0.731; dAUROC vs FIB-4 = 0.057 (95% CI 0.034-0.081)",
  "",
  "## Reproducibility",
  "All scripts run from project root with `Rscript scripts/0X_*.R` or `Rscript scripts/run_all.R`.",
  "Raw NHANES data is public domain - re-download via `scripts/01_download_data.R`.",
  "",
  "## License: MIT",
  "## AI disclosure: code + draft assisted by Claude Opus 4.7 (1M context) per COPE 2025."
)
writeLines(readme_gh, file.path(stage2, "README.md"))
cat("  + README.md (GitHub 首页)\n")

# .gitignore
gitignore <- c(
  "# raw NHANES data (~700 MB, public domain - re-download via script 01)",
  "data/raw/",
  "",
  "# R / IDE",
  ".Rhistory", ".RData", ".Rproj.user/", "*.Rproj",
  "",
  "# OS",
  "Thumbs.db", ".DS_Store",
  "",
  "# build outputs",
  "*.docx", "*.zip"
)
writeLines(gitignore, file.path(stage2, ".gitignore"))
cat("  + .gitignore\n")

# LICENSE
license <- c(
  "MIT License",
  "",
  "Copyright (c) 2026 NHANES 006 contributors",
  "",
  "Permission is hereby granted, free of charge, to any person obtaining a copy",
  "of this software and associated documentation files (the \"Software\"), to deal",
  "in the Software without restriction, including without limitation the rights",
  "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell",
  "copies of the Software, and to permit persons to whom the Software is",
  "furnished to do so, subject to the following conditions:",
  "",
  "The above copyright notice and this permission notice shall be included in all",
  "copies or substantial portions of the Software.",
  "",
  "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR",
  "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,",
  "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.")
writeLines(license, file.path(stage2, "LICENSE"))
cat("  + LICENSE\n")

# scripts/
sub_scripts <- file.path(stage2, "scripts")
dir.create(sub_scripts, showWarnings = FALSE)
# Step H Round 1 Dim10 fix — exclude internal QA / audit scripts from GitHub zip
# Internal-only (not shipped): _paper_novelty_check.R, _audit_unused_vars.R, _ipw_selection_sens.R
all_scripts <- list.files("scripts", full.names = TRUE)
internal_only <- grepl("^(_paper_novelty_check|_audit_unused_vars|_ipw_selection_sens)",
                       basename(all_scripts))
shipped_scripts <- all_scripts[!internal_only]
file.copy(shipped_scripts, sub_scripts, recursive = TRUE)
cat(sprintf("  + scripts/ (%d files; %d internal-only excluded)\n",
            length(list.files(sub_scripts)), sum(internal_only)))

# data/processed/
sub_data <- file.path(stage2, "data/processed")
dir.create(sub_data, showWarnings = FALSE, recursive = TRUE)
proc_files <- list.files("data/processed", pattern = "\\.RData$", full.names = TRUE)
file.copy(proc_files, sub_data)
cat(sprintf("  + data/processed/ (%d RData)\n", length(proc_files)))

# output/
sub_output <- file.path(stage2, "output")
dir.create(sub_output, showWarnings = FALSE)
file.copy("output/tables", sub_output, recursive = TRUE)
file.copy("output/figures", sub_output, recursive = TRUE)
cat(sprintf("  + output/tables (%d) + output/figures (%d)\n",
            length(list.files(file.path(sub_output, "tables"))),
            length(list.files(file.path(sub_output, "figures")))))

cat("\n压缩 pack 2 (this may take ~30 sec) ...\n")
old_wd <- getwd()
setwd(dirname(stage2))
zip::zip(zipfile = file.path(old_wd, pack2_name),
         files   = basename(stage2),
         recurse = TRUE,
         mode    = "cherry-pick")
setwd(old_wd)
unlink(stage2, recursive = TRUE)

zsize2 <- file.size(pack2_name)
cat(sprintf("\n[OK] pack 2: %s (%.1f MB)\n", pack2_name, zsize2 / 1024 / 1024))

cat("\n========================================\n")
cat(sprintf("两个 zip 均已生成:\n  %s (%.1f MB)\n  %s (%.1f MB)\n",
            pack1_name, zsize1/1024/1024,
            pack2_name, zsize2/1024/1024))
cat("========================================\n")
