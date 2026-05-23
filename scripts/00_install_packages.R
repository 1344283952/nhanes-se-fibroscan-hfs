# ============================================
# 00_install_packages.R
# 首次运行：安装所有依赖的 R 包
# 只需运行一次
# ============================================

cat("========================================\n")
cat("正在安装 R 包，首次运行可能需要 5-10 分钟...\n")
cat("========================================\n\n")

# 定义需要安装的包 (Round 2 R-DataChain P0: 补齐 06 scripts 实际使用的 14 个核心包)
packages <- c(
  # 基础
  "tidyverse",    # 数据处理全家桶（dplyr, tidyr, ggplot2, purrr, stringr 等）
  "survey",       # 复杂抽样设计与加权分析
  "haven",        # 读取 .xpt (SAS transport) 文件
  "broom",        # 提取回归模型结果为 tidy 格式
  "tableone",     # 快速生成基线特征表 (Table 1)
  "openxlsx",     # 输出 Excel 文件
  "DiagrammeR",   # 绘制流程图 (Figure 1)
  "DiagrammeRsvg",# CONSORT SVG → PNG/PDF/TIFF 渲染
  "rsvg",         # SVG → 高 dpi 输出
  "corrplot",     # 相关矩阵可视化（备用）
  # 统计模型 (R-Stats v2 必备)
  "mgcv",         # GAM tensor-product smooths (06_gam_dual_exposure.R)
  "rms",          # 限制立方样条 RCS + ols/lrm/cph (07_rcs_ushape.R)
  "pROC",         # AUROC + DeLong test (09_hfs_predict.R)
  "boot",         # bootstrap 包装 (09)
  "Hmisc",        # rcorr.cens for Harrell C-statistic (09)
  "ResourceSelection", # hoslem.test for Hosmer-Lemeshow calibration (09)
  "PredictABEL",  # reclassification NRI/IDI (09; nribin alternative)
  "EValue",       # E-value 不可测混杂敏感性 (14_evalue.R)
  "mgcv",         # repeat for clarity
  # 边际加强 (X3/X5/X10 等)
  # weightit + cobalt: declared but not invoked in current pipeline (_ipw_selection_sens.R uses manual 1/p_complete); kept for future Round-X bal.tab() addition
  "weightit",     # IPTW (planned)
  "cobalt",       # balance diagnostics (planned)
  "mice",         # multiple imputation (X10)
  "miceadds",     # mice + svydesign 集成
  "retrodesign",  # Type S/M error (Gelman-Carlin; X5)
  # IO / 工具
  "httr2",        # Crossref + OpenAlex API (_paper_novelty_check.R)
  "jsonlite",     # JSON parsing
  "magick",       # 图像处理 (CONSORT fallback)
  # 网络分析 / 其他备用
  "sandwich",     # robust SE
  "car",          # VIF / collinearity (R-Stats m2)
  # Round 1 v2 envint P0 fix (2026-05-23 R-DataChain): 补齐 scripts 18-23 实际使用包
  "qgcomp",       # quantile g-computation mixture (18_mixture_qgcomp.R)
  "gWQS",         # weighted quantile sum mixture (19_mixture_wqs.R)
  "CMAverse",     # 4-way decomposition mediation (21_cmaverse_redox.R)
  "ieugwasr",     # OpenGWAS API (22_mr_refit.R)
  "MendelianRandomization", # IVW/Egger/Median MR (22)
  "MRPRESSO",     # MR-PRESSO outlier (22 backup; CRAN may need archive install)
  "bkmr",         # Bayesian kernel machine regression (23_bkmr_se_metals_checkpoint.R)
  "rstan",        # Rhat diagnostic for BKMR convergence (_render_bkmr_4metal_figures.R)
  "coda",         # ESS diagnostic for BKMR (_render_bkmr_4metal_figures.R)
  "digest",       # md5 hash for processed data integrity (_verify_data_integrity.R)
  "WeightedROC",  # K-cycle weighted AUROC (_04b_validation_K_cycle.R after Round 1 P0 fix)
  "TwoSampleMR"   # alternative MR pipeline (19_mr_two_sample.R if used in 005-style fallback)
)
# Round 2 R-DataChain note: 重复条目自动去重，install.packages 容忍。

# 检查并安装缺失的包
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("正在安装: ", pkg, "\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(paste0("已安装: ", pkg, "\n"))
  }
}

cat("\n========================================\n")
cat("所有 R 包安装完成！\n")
cat("========================================\n")
