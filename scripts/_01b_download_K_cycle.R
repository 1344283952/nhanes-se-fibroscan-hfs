# ============================================
# 006v2 / _01b_download_K_cycle.R
# 下载 NHANES K cycle (Aug 2021 - Aug 2023) — 后缀 _L
# 用途: External validation in independent NHANES cycle
# 数据发布时间: Sep 2024 (LUX_L 含 FibroScan)
# ============================================

cat("========================================\n")
cat("006v2 / NHANES K cycle (2021-2023, _L) 外部验证数据下载\n")
cat("========================================\n\n")

raw_dir  <- "data/k_cycle_raw"
log_path <- "data/k_cycle_raw/_download.log"
if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)

writeLines(c("# 006v2 NHANES K cycle download log",
             paste("# started:", Sys.time())), log_path)
log_line <- function(msg) {
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_path, append = TRUE)
}

# K cycle 数据文件 — 与主分析 J/P_ 一一对应
modules <- c(
  "DEMO",   "BMX",   "BPX",  "SMQ",  "ALQ", "DIQ", "BPQ", "MCQ", "PAQ",
  "RXQ_RX",
  "PBCD",                  # Pb / Cd / Hg / Se / Mn (LBXBPB, LBXBCD, LBXTHG, LBXBSE, LBXBMN)
  "LUX",                   # FibroScan (LUXSMED, LUXCAPM)
  "BIOPRO", "CBC",
  "INS", "GLU", "GHB",
  "HDL", "TCHOL", "TRIGLY",
  "DR1TOT", "DR2TOT"
)

tasks <- list()
for (mod in modules) {
  tasks[[length(tasks) + 1]] <- list(
    filename = paste0(mod, "_L"),
    year = 2021, kind = "xpt"
  )
}

log_line(sprintf("总文件数: %d", length(tasks)))

download_task <- function(task, raw_dir) {
  is_real_xpt <- function(path) {
    if (!file.exists(path)) return(FALSE)
    if (file.size(path) < 1024) return(FALSE)
    con <- file(path, "rb"); on.exit(close(con))
    bytes <- readChar(con, 6, useBytes = TRUE)
    identical(bytes, "HEADER")
  }

  destfile <- file.path(raw_dir, paste0(task$filename, ".xpt"))
  if (is_real_xpt(destfile)) {
    return(list(file = task$filename, status = "exists"))
  }
  # K cycle url 路径: /Nchs/Data/Nhanes/Public/2021/DataFiles/<file>.xpt
  urls <- c(
    paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", task$year,
           "/DataFiles/", task$filename, ".xpt"),
    paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", task$year,
           "/DataFiles/", task$filename, ".XPT")
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
  list(file = task$filename, status = "FAIL")
}

library(parallel)
n_workers <- min(6, length(tasks))
log_line(sprintf("启动 %d worker 并行下载...", n_workers))

cl <- makeCluster(n_workers)
on.exit(stopCluster(cl), add = TRUE)

t0 <- Sys.time()
results <- parLapplyLB(cl, tasks, download_task, raw_dir = raw_dir)
elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
stopCluster(cl)

status_tbl <- table(sapply(results, function(r) r$status))
log_line(sprintf("\n并行下载结束 (%.1f 秒):", elapsed))
for (s in names(status_tbl)) log_line(sprintf("  %-10s : %d", s, status_tbl[[s]]))

fails <- sapply(results[sapply(results, function(r) r$status == "FAIL")],
                function(r) r$file)
if (length(fails) > 0) {
  log_line("\n失败文件:")
  for (f in fails) log_line(sprintf("  - %s", f))
}

n_xpt  <- length(list.files(raw_dir, pattern = "\\.xpt$"))
log_line(sprintf("\nK cycle: %d xpt 文件下载完成 → %s/", n_xpt, raw_dir))
