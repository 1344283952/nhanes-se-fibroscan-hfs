# ============================================================
# 006v2 / 19_mixture_wqs.R — WQS regression pos + neg direction
# Carrico 2015 WQS index for 5-metal mixture → liver outcomes
# 2026-05-18 v2 升级
# ============================================================
library(gWQS); library(dplyr)
set.seed(20260518)

cat("==== 19_mixture_wqs ====\n")
load("data/processed/nhanes_final.RData")

mixture <- c("pb_ugdl","cd_ugl","hg_ugl","mn_ugl","se_ugl")
covs <- c("RIDAGEYR","sex_male","education","pir","SMQ020","ALQ111","DR1TKCAL")

d <- nhanes_final %>%
  filter(if_all(all_of(mixture), ~ !is.na(.))) %>%
  filter(!is.na(RIDAGEYR), !is.na(sex_male), !is.na(DR1TKCAL))
cat(sprintf("Analytic n: %d\n", nrow(d)))

# Need education / SMQ020 / ALQ111 as factors
for (v in c("education","SMQ020","ALQ111")) {
  d[[v]] <- as.factor(d[[v]])
}

run_wqs <- function(outcome_var, family_str = "gaussian", direction = "positive") {
  fml <- as.formula(paste(outcome_var, "~ wqs +",
                          paste(c("RIDAGEYR","sex_male","education","pir","SMQ020","ALQ111","DR1TKCAL"),
                                collapse = "+")))
  fit <- tryCatch(
    gWQS::gwqs(fml, mix_name = mixture, data = d, q = 4, validation = 0.6,
               b = 200, b1_pos = (direction == "positive"), b1_constr = TRUE,
               family = family_str, seed = 20260518),
    error = function(e) { cat("  wqs failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)
  fs <- summary(fit$fit)
  coef <- fs$coefficients["wqs", ]
  list(
    outcome = outcome_var, direction = direction,
    beta = coef[1], se = coef[2], z = coef[3], p = coef[4],
    weights = paste(rownames(fit$final_weights),
                    sprintf("%.3f", fit$final_weights$mean_weight),
                    sep = "=", collapse = "; ")
  )
}

outcomes <- list(
  list(var = "cap", family = "gaussian"),
  list(var = "lsm", family = "gaussian"),
  list(var = "ln_ggt", family = "gaussian"),
  list(var = "steatosis_cap275", family = "binomial"),
  list(var = "fibrosis_lsm8", family = "binomial"),
  list(var = "ggt_high", family = "binomial")
)

results <- list()
for (oc in outcomes) {
  if (!(oc$var %in% colnames(d))) next
  cat(sprintf("\n--- %s (%s) ---\n", oc$var, oc$family))
  for (dir in c("positive","negative")) {
    r <- run_wqs(oc$var, oc$family, dir)
    if (!is.null(r)) results[[paste0(oc$var, "_", dir)]] <- r
  }
}

if (length(results) > 0) {
  res_df <- do.call(rbind, lapply(results, function(r)
    data.frame(outcome = r$outcome, direction = r$direction,
               beta = r$beta, se = r$se, z = r$z, p = r$p,
               weights = r$weights)))
  res_df$p_BH <- p.adjust(res_df$p, "BH")
  write.csv(res_df, "output/tables/wqs_mixture.csv", row.names = FALSE)
  cat("\n=== WQS 结果 ===\n"); print(res_df, row.names = FALSE)
  cat("\n[OK] saved to output/tables/wqs_mixture.csv\n")
} else {
  cat("[ERROR] 0 WQS 结果\n")
}
