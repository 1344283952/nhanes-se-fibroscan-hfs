# ============================================================
# 006v2 / 18_mixture_qgcomp.R — qgcomp 5-metal mixture analysis
# Pb / Cd / Hg / Mn / Se → CAP / LSM / HFS  (Keil 2020 qgcomp)
# 2026-05-18 v2 升级
# ============================================================
library(qgcomp); library(dplyr); library(survey)
set.seed(20260518)

cat("==== 18_mixture_qgcomp ====\n")
load("data/processed/nhanes_final.RData")

# 主 cohort 已有 5 metals + outcomes, 100% valid
d <- nhanes_final %>%
  filter(!is.na(pb_ugdl), !is.na(cd_ugl), !is.na(hg_ugl), !is.na(mn_ugl), !is.na(se_ugl))
cat(sprintf("Analytic n: %d\n", nrow(d)))

# Mixture: 5 metals on log scale (qgcomp quantizes internally)
mixture <- c("ln_pb","ln_cd","ln_hg","ln_mn","ln_se_blood")

# Covariates — Pearl-backdoor cov_pre set (manuscript §2.4)
covs <- c("RIDAGEYR","sex_male","race","education","pir","SMQ020","ALQ111","DR1TKCAL")
# 把 NA codes 转 missing
for (cv in covs) if (cv %in% colnames(d)) d[[cv]] <- ifelse(d[[cv]] %in% c(7,9,77,99,777,999), NA, d[[cv]])
d <- d %>% filter(complete.cases(d[, c("RIDAGEYR","sex_male","DR1TKCAL")]))
cat(sprintf("After cov filter: n=%d\n", nrow(d)))

run_qgcomp <- function(outcome_var, family = gaussian()) {
  fml <- as.formula(paste(outcome_var, "~", paste(c(mixture, covs), collapse = "+")))
  fit <- tryCatch(
    qgcomp::qgcomp.noboot(fml, expnms = mixture, data = d, family = family, q = 4),
    error = function(e) { cat("  qgcomp failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)
  ws <- fit$pos.weights; if (is.null(ws)) ws <- numeric()
  wn <- fit$neg.weights; if (is.null(wn)) wn <- numeric()
  list(
    outcome = outcome_var, psi = fit$coef["psi1"],
    se = sqrt(fit$var.psi), p = 2*(1-pnorm(abs(fit$coef["psi1"]/sqrt(fit$var.psi)))),
    ci_lo = fit$ci[1], ci_hi = fit$ci[2],
    pos_metals = paste(names(ws), sprintf("%.2f", ws), sep="=", collapse="; "),
    neg_metals = paste(names(wn), sprintf("%.2f", wn), sep="=", collapse="; ")
  )
}

# 5 outcomes: continuous + binary
results <- list()
for (oc in c("cap","lsm","hfs")) {
  if (oc %in% colnames(d)) {
    r <- run_qgcomp(oc, gaussian())
    if (!is.null(r)) { r$type <- "continuous"; results[[paste0("qg_",oc)]] <- r }
  }
}
for (oc in c("steatosis_cap275","fibrosis_lsm8","hfs_high","ggt_high","hyperuricemia")) {
  if (oc %in% colnames(d) && all(d[[oc]] %in% c(0,1,NA))) {
    r <- run_qgcomp(oc, binomial())
    if (!is.null(r)) { r$type <- "binary"; results[[paste0("qg_",oc)]] <- r }
  }
}

if (length(results) > 0) {
  res_df <- do.call(rbind, lapply(results, function(r)
    data.frame(outcome=r$outcome, type=r$type, psi=r$psi,
               se=r$se, ci_lo=r$ci_lo, ci_hi=r$ci_hi, p=r$p,
               pos_metals=r$pos_metals, neg_metals=r$neg_metals)))
  res_df$p_BH <- p.adjust(res_df$p, "BH")
  write.csv(res_df, "output/tables/qgcomp_mixture.csv", row.names=FALSE)
  cat("\n=== qgcomp 结果 ===\n"); print(res_df, row.names=FALSE)
  cat("\n[OK] saved to output/tables/qgcomp_mixture.csv\n")
} else {
  cat("[ERROR] 0 qgcomp 结果\n")
}
