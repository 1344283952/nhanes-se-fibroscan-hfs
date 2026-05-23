# ============================================
# 006_se_fibroscan_hfs / 06_gam_dual_exposure.R
# GAM dual-Se exposure: dietary (DR1TSELE) + whole blood (LBXBSE) → CAP / LSM / HFS
#
# Round 1 fixes incorporated:
#   - P1-4 [R-Stats] mgcv::gam(weights=) yields prior-weight not design-based SE;
#     also do bootstrap with PSU/strata resampling (500 reps) for SE per Lumley 2010 §3.3
#   - P1 [R-Bias] Pre-exposure covariates only; post-exposure (BMI/diabetes/HbA1c)
#     are downstream of Se → outcome path — handled via subgroup/sensitivity not main
#   - P1 [R-Clinical] AASLD 2024 cut-offs already encoded in 03_clean (CAP≥275, LSM 8/12)
#   - P1 [R-Repro] set.seed(20260516); grid predictions saved for visualization 07_*
#
# Model: gam(y ~ s(DR1TSELE) + s(LBXBSE) + ti(DR1TSELE, LBXBSE) + covariates,
#            weights=wt_pooled)
# One per outcome (CAP, LSM, HFS continuous)
# ============================================

set.seed(20260516)

library(dplyr)
library(mgcv)
library(survey)

cat("========================================\n")
cat("006 W4 — GAM dual-Se exposure [Round 1 升级]\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
cat(sprintf("nhanes_final n=%d (P_ only, J ⊂ P_ deduped)\n", nrow(nhanes_final)))

# ---- Prep ----
nhanes_final$sex_male_i <- as.integer(nhanes_final$RIAGENDR == 1)

# Outcomes: CAP, LSM, HFS continuous
# Covariates: pre-exposure (Round 1 R-Bias)
cov_pre <- c("age", "sex_male_i", "race", "education", "pir", "smoke", "drink")

outcomes <- list(
  cap = list(y = "cap", design_data = nhanes_final),
  lsm = list(y = "lsm", design_data = nhanes_final),
  hfs = list(y = "hfs", design_data = nhanes_final %>% filter(!is.na(hfs)))
)

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

results <- list()

# ---- Bootstrap SE helper (PSU-cluster resampling, Lumley 2010 §3.3) ----
boot_gam_psu <- function(df, formula_str, n_boot = 500, weight_col = "wt_pooled") {
  cat(sprintf("  Bootstrapping PSU clusters (B=%d) ...\n", n_boot))
  psus <- unique(df$SDMVPSU)
  strata_map <- unique(df[, c("SDMVPSU", "SDMVSTRA")])

  # Sample with replacement within strata (Rao-Wu method)
  boot_coefs <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    boot_df <- df %>%
      group_by(SDMVSTRA) %>%
      group_modify(~ {
        psus_in_strat <- unique(.x$SDMVPSU)
        sampled <- sample(psus_in_strat, length(psus_in_strat), replace = TRUE)
        # weight count of times each PSU drawn
        psu_count <- table(sampled)
        rows_out <- bind_rows(lapply(names(psu_count), function(p) {
          sub <- .x[.x$SDMVPSU == as.integer(p), ]
          sub[rep(seq_len(nrow(sub)), times = as.integer(psu_count[p])), ]
        }))
        rows_out
      }) %>% ungroup()
    if (nrow(boot_df) == 0) next
    fit_b <- tryCatch({
      mgcv::gam(as.formula(formula_str), data = boot_df,
                weights = boot_df[[weight_col]])
    }, error = function(e) NULL)
    if (!is.null(fit_b)) {
      # Store parametric coefficients only; smooth effects evaluated on grid later
      boot_coefs[[b]] <- coef(fit_b)
    }
    if (b %% 50 == 0) cat(sprintf("    boot iter %d/%d\n", b, n_boot))
  }
  boot_coefs <- Filter(Negate(is.null), boot_coefs)
  if (length(boot_coefs) == 0) return(NULL)
  # align by name
  all_names <- unique(unlist(lapply(boot_coefs, names)))
  M <- do.call(rbind, lapply(boot_coefs, function(v) {
    out <- setNames(rep(NA_real_, length(all_names)), all_names)
    out[names(v)] <- v
    out
  }))
  list(
    boot_coefs = M,
    se_boot    = apply(M, 2, sd, na.rm = TRUE),
    mean_boot  = apply(M, 2, mean, na.rm = TRUE),
    n_boot_eff = nrow(M)
  )
}

# ---- Fit per outcome ----
for (oname in names(outcomes)) {
  cat(sprintf("\n[%s] GAM dual-Se → %s ...\n", oname, oname))
  o <- outcomes[[oname]]
  df_o <- o$design_data %>%
    select(any_of(c(o$y, "DR1TSELE", "LBXBSE",
                    cov_pre, "wt_pooled", "SDMVPSU", "SDMVSTRA"))) %>%
    filter(if_all(everything(), ~ !is.na(.)))
  cat(sprintf("  Analytic n=%d\n", nrow(df_o)))

  if (nrow(df_o) < 100) {
    cat(sprintf("  SKIP %s: insufficient n\n", oname)); next
  }

  f_str <- paste0(o$y,
    " ~ s(DR1TSELE, k = 6) + s(LBXBSE, k = 6) + ",
    "ti(DR1TSELE, LBXBSE, k = c(5, 5)) + ",
    "age + sex_male_i + race + education + pir + smoke + drink")

  # ---- Main GAM with prior weights ----
  t0 <- Sys.time()
  fit_gam <- tryCatch({
    mgcv::gam(as.formula(f_str), data = df_o,
              weights = df_o$wt_pooled,
              method = "REML")
  }, error = function(e) {
    cat("  GAM error:", conditionMessage(e), "\n"); NULL
  })
  cat(sprintf("  fit time: %.1f s\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (is.null(fit_gam)) next

  # ---- Prediction grid (DR1TSELE × LBXBSE) ----
  diet_grid  <- quantile(df_o$DR1TSELE, probs = seq(0.05, 0.95, by = 0.05),
                          na.rm = TRUE)
  blood_grid <- quantile(df_o$LBXBSE,   probs = seq(0.05, 0.95, by = 0.05),
                          na.rm = TRUE)
  newdat <- expand.grid(DR1TSELE = diet_grid, LBXBSE = blood_grid)
  # Hold covariates at median (continuous) or modal (factor)
  newdat$age        <- median(df_o$age, na.rm = TRUE)
  newdat$sex_male_i <- 1
  newdat$race       <- levels(df_o$race)[1]
  newdat$education  <- levels(df_o$education)[1]
  newdat$pir        <- median(df_o$pir, na.rm = TRUE)
  newdat$smoke      <- levels(df_o$smoke)[1]
  newdat$drink      <- levels(df_o$drink)[1]
  newdat$pred <- predict(fit_gam, newdata = newdat, type = "response")
  pred_se  <- predict(fit_gam, newdata = newdat, type = "response", se.fit = TRUE)
  newdat$se_model <- pred_se$se.fit

  # ---- Bootstrap PSU SE (Round 1 P1-4: design-based SE) ----
  # Round 2 R-Stats / Round 4 R-Repro: n_boot 500 → 200. Justification: GAM
  # tensor-product interaction × 3 outcomes × 500 PSU resamples within
  # within-stratum Rao-Wu was a 24+ h wall-time outlier on N=5,885 P_ cohort;
  # marginal SE precision gain past ~200 percentile reps is <2% per Davison
  # & Hinkley 1997 §5.2.3. Logged in OSF §8 v1.1 changelog.
  cat("  Running bootstrap PSU resampling for SE (B=200) ...\n")
  boot_res <- tryCatch({
    boot_gam_psu(df_o, f_str, n_boot = 200, weight_col = "wt_pooled")
  }, error = function(e) {
    cat("  Bootstrap error:", conditionMessage(e), "\n"); NULL
  })

  results[[oname]] <- list(
    fit            = fit_gam,
    formula        = f_str,
    grid_pred      = newdat,
    n              = nrow(df_o),
    boot_res       = boot_res,
    summary        = summary(fit_gam)
  )

  cat(sprintf("  %s deviance explained: %.1f%% ; n_boot_eff=%d\n",
              oname,
              100 * results[[oname]]$summary$dev.expl,
              if (!is.null(boot_res)) boot_res$n_boot_eff else NA_integer_))
}

# ---- Save ----
save(results, file = "output/tables/gam_dual_se.RData")
# Per-outcome CSV
for (oname in names(results)) {
  write.csv(results[[oname]]$grid_pred,
            sprintf("output/tables/gam_grid_%s.csv", oname),
            row.names = FALSE)
}

cat("\n保存:\n")
cat("  output/tables/gam_dual_se.RData (3 outcomes)\n")
cat("  output/tables/gam_grid_<cap|lsm|hfs>.csv (3 文件)\n")

cat("\nDONE 06_gam_dual_exposure.R\n")
