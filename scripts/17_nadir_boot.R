# 006 — RCS Nadir 95% Bootstrap CI for Se → CAP U-shape (OSF v1.5 S9)
# 1000× percentile bootstrap on natural-spline (df=4) lm Se→CAP model
# + sensitivity adjusting for ALT (LBXSATSI)
# 2026-05-18 v2 — avoid rms::datadist (silent failure inside boot env)

suppressPackageStartupMessages({
  library(splines)
})
set.seed(20260516)
cat("==========================================================\n")
cat("006 RCS Nadir 1000x Bootstrap 95% CI (v2 ns-based)\n")
cat("OSF v1.5 S9 (2026-05-18)\n")
cat("==========================================================\n\n")

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE)

log_file <- file.path("logs", paste0("nadir_v2_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

load("data/processed/nhanes_final.RData")
df <- if (exists("nhanes_final")) nhanes_final else nhanes_main
df <- as.data.frame(df)
df <- df[!is.na(df$LBXBSE) & !is.na(df$LUXCAPM), ]
cat("Analytic N =", nrow(df), "\n\n")

# ----------------------------------------------------------
# helper: one bootstrap nadir using splines::ns(df=4) lm
# ----------------------------------------------------------
nadir_one_unadj <- function(d) {
  idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
  d2 <- d[idx, c("LBXBSE", "LUXCAPM"), drop = FALSE]
  d2 <- d2[complete.cases(d2), ]
  if (nrow(d2) < 100) return(NA_real_)
  f <- tryCatch(lm(LUXCAPM ~ ns(LBXBSE, df = 4), data = d2), error = function(e) NULL)
  if (is.null(f)) return(NA_real_)
  q_lo <- as.numeric(quantile(d2$LBXBSE, 0.02))
  q_hi <- as.numeric(quantile(d2$LBXBSE, 0.98))
  grid <- seq(q_lo, q_hi, length.out = 500)
  p <- tryCatch(predict(f, newdata = data.frame(LBXBSE = grid)), error = function(e) NULL)
  if (is.null(p)) return(NA_real_)
  grid[which.min(p)]
}

nadir_one_adj <- function(d) {
  d <- d[!is.na(d$LBXSATSI), ]
  idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
  d2 <- d[idx, c("LBXBSE", "LUXCAPM", "LBXSATSI"), drop = FALSE]
  d2 <- d2[complete.cases(d2), ]
  if (nrow(d2) < 100) return(NA_real_)
  f <- tryCatch(lm(LUXCAPM ~ ns(LBXBSE, df = 4) + LBXSATSI, data = d2),
                error = function(e) NULL)
  if (is.null(f)) return(NA_real_)
  q_lo <- as.numeric(quantile(d2$LBXBSE, 0.02))
  q_hi <- as.numeric(quantile(d2$LBXBSE, 0.98))
  grid <- seq(q_lo, q_hi, length.out = 500)
  p <- tryCatch(predict(f, newdata = data.frame(LBXBSE = grid,
                                                  LBXSATSI = median(d2$LBXSATSI))),
                error = function(e) NULL)
  if (is.null(p)) return(NA_real_)
  grid[which.min(p)]
}

# ----------------------------------------------------------
# 1. Unadjusted
# ----------------------------------------------------------
cat("[1] Unadjusted nadir bootstrap (1000 iter)\n")
ptm <- proc.time()
nadirs1 <- vapply(seq_len(1000), function(i) nadir_one_unadj(df), numeric(1))
elapsed1 <- as.numeric((proc.time() - ptm)[3])
cat("  Bootstrap finished in", round(elapsed1, 1), "sec\n")

valid1 <- nadirs1[!is.na(nadirs1)]
cat("  Valid replicates:", length(valid1), "of 1000\n")
if (length(valid1) >= 10) {
  pe1 <- median(valid1); ci1 <- quantile(valid1, c(0.025, 0.975))
  cat("  Nadir (median):", round(pe1, 1), "µg/L\n")
  cat("  Bootstrap 95% CI:", round(ci1[1], 1), "-", round(ci1[2], 1), "µg/L\n")
  prob_in1 <- mean(valid1 >= 130 & valid1 <= 170)
  cat("  P(nadir in prespecified 130-170 µg/L):", round(prob_in1, 3), "\n")
} else { pe1 <- NA; ci1 <- c(NA, NA); prob_in1 <- NA }

out1 <- data.frame(model = "Unadjusted ns(df=4) lm",
                   estimate = pe1, lower95 = ci1[1], upper95 = ci1[2],
                   prob_in_130_170 = prob_in1,
                   B = 1000, N_valid = length(valid1),
                   elapsed_sec = round(elapsed1, 2))

# ----------------------------------------------------------
# 2. ALT-adjusted
# ----------------------------------------------------------
if ("LBXSATSI" %in% names(df) && sum(!is.na(df$LBXSATSI)) >= 1000) {
  cat("\n[2] ALT-adjusted nadir bootstrap\n")
  ptm <- proc.time()
  nadirs2 <- vapply(seq_len(1000), function(i) nadir_one_adj(df), numeric(1))
  elapsed2 <- as.numeric((proc.time() - ptm)[3])
  cat("  Bootstrap finished in", round(elapsed2, 1), "sec\n")

  valid2 <- nadirs2[!is.na(nadirs2)]
  cat("  Valid replicates:", length(valid2), "of 1000\n")
  if (length(valid2) >= 10) {
    pe2 <- median(valid2); ci2 <- quantile(valid2, c(0.025, 0.975))
    cat("  Nadir (ALT-adj):", round(pe2, 1), "µg/L\n")
    cat("  Bootstrap 95% CI:", round(ci2[1], 1), "-", round(ci2[2], 1), "µg/L\n")
    prob_in2 <- mean(valid2 >= 130 & valid2 <= 170)
  } else { pe2 <- NA; ci2 <- c(NA, NA); prob_in2 <- NA }
  out2 <- data.frame(model = "ALT-adjusted ns(df=4) lm",
                     estimate = pe2, lower95 = ci2[1], upper95 = ci2[2],
                     prob_in_130_170 = prob_in2,
                     B = 1000, N_valid = length(valid2),
                     elapsed_sec = round(elapsed2, 2))
} else {
  cat("\n[2] ALT skipped (no LBXSATSI or N<1000)\n")
  out2 <- data.frame(model = "ALT-adjusted ns(df=4) lm",
                     estimate = NA, lower95 = NA, upper95 = NA,
                     prob_in_130_170 = NA, B = 0, N_valid = 0, elapsed_sec = 0)
}

out_all <- rbind(out1, out2)
write.csv(out_all, "output/tables/nadir_bootstrap.csv", row.names = FALSE)
saveRDS(list(unadj = nadirs1, adj = if(exists("nadirs2")) nadirs2 else NULL),
        "output/tables/nadir_bootstrap_raw.rds")

cat("\n=== Summary ===\n")
print(out_all)
cat("\n==========================================================\n")
cat("Nadir bootstrap (v2 ns-based) complete\n")
cat("Output: output/tables/nadir_bootstrap.csv + .rds\n")
cat("==========================================================\n")
