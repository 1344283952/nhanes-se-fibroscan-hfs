# ============================================
# 006_se_fibroscan_hfs / 02_merge_data.R
# 合并 NHANES 2017-2018 (_J) + Pre-pandemic 2017-March 2020 (P_)
# 截面 cohort, FibroScan + Se 双暴露 + Hepamet HFS
#
# 输入: data/raw/*.xpt + data/raw/mortality/*.dat (mortality 仅 2017-2018, secondary)
# 输出: data/processed/nhanes_raw_merged.RData
# ============================================

library(haven); library(dplyr); library(purrr)

cat("========================================\n")
cat("006 P12: NHANES 2017-March 2020 (J + P_) + LUX\n")
cat("========================================\n\n")

raw_dir  <- "data/raw"
mort_dir <- "data/raw/mortality"

xpt_files <- list.files(raw_dir, pattern = "\\.xpt$", full.names = TRUE, ignore.case = TRUE)
cat(sprintf("找到 %d 个 .xpt 文件\n\n", length(xpt_files)))

read_safe <- function(p) {
  tryCatch(read_xpt(p), error = function(e) {
    cat(sprintf("  [读失败] %s\n", basename(p))); NULL
  })
}
data_list <- map(xpt_files, read_safe)
names(data_list) <- gsub("\\.xpt$", "", basename(xpt_files), ignore.case = TRUE)
data_list <- data_list[!sapply(data_list, is.null)]
cat(sprintf("成功读 %d 个数据帧\n\n", length(data_list)))

# 仅 J + P_
cycles <- list(
  list(year = 2017, suffix = "_J", prefix = "", tag = "NHANES_2017_2018"),
  list(year = 2017, suffix = "",   prefix = "P_", tag = "PrePandemic_2017_March2020")
)

modules_uniform <- c(
  "DEMO", "BMX", "BPX", "SMQ", "ALQ", "DIQ", "BPQ", "MCQ", "PAQ",
  "PBCD",       # Se whole blood (LBXBSE)
  "LUX",        # FibroScan (LUXSMED + LUXCAPM)
  "BIOPRO",     # AST/ALT/Albumin/GGT
  "CBC",        # Platelet
  "INS", "GLU", "GHB",
  "HDL", "TCHOL", "TRIGLY",
  "DR1TOT", "DR2TOT",
  "VID", "ALB_CR"
)

safe_left_join <- function(x, y, key = "SEQN") {
  if (is.null(y) || nrow(y) == 0) return(x)
  if (!key %in% names(y)) return(x)
  dup <- intersect(setdiff(names(x), key), names(y))
  if (length(dup) > 0) y <- y[, setdiff(names(y), dup), drop = FALSE]
  left_join(x, y, by = key)
}

merge_cycle <- function(cyc) {
  cat(sprintf("--- %s ---\n", cyc$tag))
  base_name <- paste0(cyc$prefix, "DEMO", cyc$suffix)
  if (!base_name %in% names(data_list)) {
    cat(sprintf("  [警告] %s 不存在\n\n", base_name)); return(NULL)
  }
  result <- data_list[[base_name]]
  result$cycle_year     <- cyc$year
  result$cycle_tag      <- cyc$tag
  result$is_prepandemic <- (cyc$prefix == "P_")
  for (mod in modules_uniform[-1]) {
    key <- paste0(cyc$prefix, mod, cyc$suffix)
    if (key %in% names(data_list)) {
      result <- safe_left_join(result, data_list[[key]])
    }
  }
  cat(sprintf("  -> %d 行 × %d 列\n\n", nrow(result), ncol(result)))
  result
}

merged_list <- lapply(cycles, merge_cycle)
merged_list <- merged_list[!sapply(merged_list, is.null)]
nhanes_all <- bind_rows(merged_list)

cat(sprintf("========================================\n"))
cat(sprintf("J + P_ 合并: %d 行 × %d 列\n", nrow(nhanes_all), ncol(nhanes_all)))
cat(sprintf("========================================\n\n"))

# 长表 (006 不用 RXQ_RX 主分析，但下载了就保留)
rx_list <- list()
for (cyc in cycles) {
  key <- paste0(cyc$prefix, "RXQ_RX", cyc$suffix)
  if (key %in% names(data_list)) rx_list[[key]] <- data_list[[key]]
}
rx_all <- if (length(rx_list) > 0) bind_rows(rx_list) else data.frame(SEQN = integer(0))
cat(sprintf("rx_all: %d 行 × %d 列\n\n", nrow(rx_all), ncol(rx_all)))

# Mortality (secondary, 仅 2017-2018 NCHS 2019 release)
cat("--- Mortality (secondary, 2017-2018 only) ---\n")
mort_files <- list.files(mort_dir, pattern = "\\.dat$", full.names = TRUE)
read_mort <- function(path) {
  widths <- c(SEQN = 6, PADDING1 = 8, ELIGSTAT = 1, MORTSTAT = 1,
              UCOD_LEADING = 3, DIABETES = 1, HYPERTEN = 1, DODQTR = 1,
              DODYEAR = 4, WGT_NEW = 8, SA_WGT_NEW = 8,
              PERMTH_INT = 3, PERMTH_EXM = 3)
  df <- tryCatch(
    read.fwf(path, widths = widths, header = FALSE,
             na.strings = c("", "."), stringsAsFactors = FALSE,
             col.names = names(widths)),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)
  df$PADDING1 <- NULL
  df$SEQN <- as.integer(df$SEQN)
  for (c in c("ELIGSTAT","MORTSTAT","UCOD_LEADING","DIABETES",
              "HYPERTEN","PERMTH_INT","PERMTH_EXM")) {
    df[[c]] <- suppressWarnings(as.integer(df[[c]]))
  }
  df
}
mort_all <- bind_rows(lapply(mort_files, read_mort))
cat(sprintf("Mortality 链接: %d 行\n", nrow(mort_all)))
cat(sprintf("  ELIGSTAT == 1: %d\n", sum(mort_all$ELIGSTAT == 1, na.rm = TRUE)))
cat(sprintf("  Deaths: %d\n", sum(mort_all$MORTSTAT == 1, na.rm = TRUE)))

nhanes_all <- safe_left_join(nhanes_all, mort_all)
cat(sprintf("\n合并 mortality 后: %d 行 × %d 列\n", nrow(nhanes_all), ncol(nhanes_all)))

# Save
if (!dir.exists("data/processed")) dir.create("data/processed", recursive = TRUE)
save(nhanes_all, rx_all, mort_all, file = "data/processed/nhanes_raw_merged.RData")
cat("\n已保存 data/processed/nhanes_raw_merged.RData\n")
cat("========================================\n")
