# ============================================
# 006_se_fibroscan_hfs / 01_download_data.R
# 下载 NHANES 2017-2018 (_J) + Pre-pandemic 2017-March 2020 (P_)
# 主要 outcome: FibroScan (LUX) — LSM + CAP
# Hepamet HFS 需: AST + Albumin + Platelet + HOMA-IR + DM + age + sex
#
# 论文: Se 双暴露 (膳食 + 全血) + Se/Zn-Se/Cu 比值 × FibroScan × HFS
#       + MetALD 分层 + U-shape RCS
# ============================================

cat("========================================\n")
cat("006_se_fibroscan_hfs / NHANES 2017-March 2020\n")
cat("========================================\n\n")

raw_dir  <- "data/raw"
mort_dir <- "data/raw/mortality"
log_path <- "data/raw/_download.log"
if (!dir.exists(raw_dir))  dir.create(raw_dir,  recursive = TRUE)
if (!dir.exists(mort_dir)) dir.create(mort_dir, recursive = TRUE)

writeLines(c("# 006 NHANES download log", paste("# started:", Sys.time())),
           log_path)
log_line <- function(msg) {
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_path, append = TRUE)
}

cycles <- data.frame(
  suffix = c("_J", ""),
  prefix = c("",   "P_"),
  pathYear = c(2017, 2017),
  stringsAsFactors = FALSE
)

modules_uniform <- c(
  "DEMO", "BMX", "BPX", "SMQ", "ALQ", "DIQ", "BPQ", "MCQ", "PAQ",
  "RXQ_RX",
  # Se 双暴露 — PBCD 内含 LBXBSE 全血 Se
  "PBCD",
  # FibroScan
  "LUX",        # LSM + CAP
  # HFS 必需
  "BIOPRO",     # AST/ALT/Albumin/GGT
  "CBC",        # Platelet LBXPLTSI
  "INS",        # 胰岛素 (HOMA-IR fasting subsample)
  "GLU",        # 空腹血糖
  "GHB",        # HbA1c
  # 血脂 (cardiometabolic risk for MetALD)
  "HDL", "TCHOL", "TRIGLY",
  # 膳食 (DR1TSELE 膳食 Se + DR1TZINC + DR1TCOPP for ratios)
  "DR1TOT", "DR2TOT",
  "VID",
  "ALB_CR"
)

tasks <- list()
for (i in seq_len(nrow(cycles))) {
  for (mod in modules_uniform) {
    if (cycles$prefix[i] == "P_") {
      tasks[[length(tasks)+1]] <- list(
        filename = paste0("P_", mod),
        year = cycles$pathYear[i], kind = "xpt"
      )
    } else {
      tasks[[length(tasks)+1]] <- list(
        filename = paste0(mod, cycles$suffix[i]),
        year = cycles$pathYear[i], kind = "xpt"
      )
    }
  }
}

# Mortality 可选 (2017-2018 NHANES NCHS 2019 release)
mort_base <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality"
yr <- 2017
fn <- paste0("NHANES_", yr, "_", yr+1, "_MORT_2019_PUBLIC.dat")
tasks[[length(tasks)+1]] <- list(
  filename = fn, year = yr, kind = "mort",
  url = paste0(mort_base, "/", fn)
)

log_line(sprintf("总任务数: %d (xpt=%d, mort=%d)",
                 length(tasks),
                 sum(sapply(tasks, function(t) t$kind == "xpt")),
                 sum(sapply(tasks, function(t) t$kind == "mort"))))

download_task <- function(task, raw_dir, mort_dir) {
  is_real_xpt <- function(path) {
    if (!file.exists(path)) return(FALSE)
    if (file.size(path) < 1024) return(FALSE)
    con <- file(path, "rb"); on.exit(close(con))
    bytes <- readChar(con, 6, useBytes = TRUE)
    identical(bytes, "HEADER")
  }

  if (task$kind == "xpt") {
    destfile <- file.path(raw_dir, paste0(task$filename, ".xpt"))
    if (is_real_xpt(destfile)) {
      return(list(file = task$filename, status = "已存在"))
    }
    urls <- c(
      paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",
             task$year, "/DataFiles/", task$filename, ".xpt"),
      paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/",
             task$year, "/DataFiles/", task$filename, ".XPT")
    )
    for (url in urls) {
      ok <- tryCatch({
        download.file(url, destfile = destfile, mode = "wb",
                      quiet = TRUE, method = "libcurl")
        is_real_xpt(destfile)
      }, error = function(e) FALSE, warning = function(w) FALSE)
      if (isTRUE(ok)) return(list(file = task$filename, status = "ok"))
      if (file.exists(destfile)) file.remove(destfile)
    }
    return(list(file = task$filename, status = "FAIL"))
  } else {
    destfile <- file.path(mort_dir, task$filename)
    if (file.exists(destfile) && file.size(destfile) > 0) {
      return(list(file = task$filename, status = "已存在"))
    }
    ok <- tryCatch({
      download.file(task$url, destfile = destfile, mode = "wb",
                    quiet = TRUE, method = "libcurl")
      file.exists(destfile) && file.size(destfile) > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(list(file = task$filename, status = "ok"))
    return(list(file = task$filename, status = "FAIL"))
  }
}

library(parallel)
n_workers <- min(6, length(tasks))
log_line(sprintf("启动 %d worker 并行下载...", n_workers))

cl <- makeCluster(n_workers)
on.exit(stopCluster(cl), add = TRUE)

t0 <- Sys.time()
results <- parLapplyLB(cl, tasks, download_task,
                       raw_dir = raw_dir, mort_dir = mort_dir)
elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

stopCluster(cl)

status_tbl <- table(sapply(results, function(r) r$status))
log_line(sprintf("\n并行下载结束 (%.1f 秒):", elapsed))
for (s in names(status_tbl)) {
  log_line(sprintf("  %-10s : %d", s, status_tbl[[s]]))
}

fails <- sapply(results[sapply(results, function(r) r$status == "FAIL")],
                function(r) r$file)
if (length(fails) > 0) {
  log_line("\n失败文件:")
  for (f in fails) log_line(sprintf("  - %s", f))
}

n_xpt  <- length(list.files(raw_dir,  pattern = "\\.xpt$"))
n_mort <- length(list.files(mort_dir, pattern = "\\.dat$"))
log_line(sprintf("\ndata/raw/         %d 个 .xpt", n_xpt))
log_line(sprintf("data/raw/mortality/ %d 个 .dat", n_mort))
log_line("========================================")
