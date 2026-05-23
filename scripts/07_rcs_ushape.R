# ============================================
# 006_se_fibroscan_hfs / 07_rcs_ushape.R
# Restricted cubic spline U-shape for Se → CAP/LSM/HFS
#
# Round 1 fixes incorporated:
#   - P1-1 [R-Stats] 5-knot RCS for Se (升级 from 4-knot per Harrell 2015 + Bobb 2015
#     mixture-modelling motivation; 5 knots better captures U-shape inflection)
#   - P1 [R-Stats] LR test vs linear (not Wald only) for non-linearity p-value
#     (Harrell RMS §2.4.5; Wald can miss non-linearity when boundary knot has
#     low information)
#   - P1 [R-Stats] Knot placement: 5/27.5/50/72.5/95 percentiles per Harrell 2015
#     standard 5-knot positions
#   - P1 [R-Bias] Pre-exposure covariates only (Round 1 R-Bias mediator-confounder)
#   - P1 [R-Repro] set.seed(20260516); ribbon plots with 95% CI
#
# Outcomes: 3 (CAP, LSM, HFS continuous) × 2 exposures (DR1TSELE, LBXBSE)
# = 6 RCS fits. Each with linear vs RCS LR test.
# ============================================

set.seed(20260516)

library(dplyr)
library(rms)
library(survey)
library(ggplot2)

cat("========================================\n")
cat("006 W4 — RCS U-shape (Se → CAP/LSM/HFS) [Round 1 升级]\n")
cat("========================================\n\n")

load("data/processed/nhanes_final.RData")
cat(sprintf("nhanes_final n=%d\n", nrow(nhanes_final)))

nhanes_final$sex_male_i <- as.integer(nhanes_final$RIAGENDR == 1)
nhanes_final$race_chr   <- as.character(nhanes_final$race)
nhanes_final$edu_chr    <- as.character(nhanes_final$education)
nhanes_final$smoke_chr  <- as.character(nhanes_final$smoke)
nhanes_final$drink_chr  <- as.character(nhanes_final$drink)

# 5-knot positions per Harrell 2015 (RMS §2.4.5)
knot_quantiles <- c(0.05, 0.275, 0.50, 0.725, 0.95)

cov_pre_str <- "age + sex_male_i + race_chr + edu_chr + pir + smoke_chr + drink_chr"

if (!dir.exists("output/tables")) dir.create("output/tables", recursive = TRUE)
if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)

# rms datadist (subset to vars actually used per CLAUDE.md 已踩坑 #4)
rms_vars <- c("cap", "lsm", "hfs", "DR1TSELE", "LBXBSE",
              "age", "sex_male_i", "race_chr", "edu_chr",
              "pir", "smoke_chr", "drink_chr")
dd <- datadist(nhanes_final[, intersect(rms_vars, names(nhanes_final))])
options(datadist = "dd")

# ---- Per outcome × exposure ----
combos <- expand.grid(
  outcome  = c("cap", "lsm", "hfs"),
  exposure = c("DR1TSELE", "LBXBSE"),
  stringsAsFactors = FALSE
)

results_list <- list()
pvals_df <- data.frame()

for (i in seq_len(nrow(combos))) {
  out <- combos$outcome[i]
  exp <- combos$exposure[i]
  tag <- sprintf("%s_%s", out, exp)
  cat(sprintf("\n[%d/%d] RCS %s ~ %s ...\n", i, nrow(combos), out, exp))

  df_o <- nhanes_final %>%
    select(any_of(c(out, exp, "age", "sex_male_i", "race_chr",
                    "edu_chr", "pir", "smoke_chr", "drink_chr", "wt_pooled"))) %>%
    filter(if_all(everything(), ~ !is.na(.)))
  if (nrow(df_o) < 200) {
    cat(sprintf("  SKIP %s: insufficient n (%d)\n", tag, nrow(df_o))); next
  }
  cat(sprintf("  Analytic n=%d\n", nrow(df_o)))

  # Knot positions (per Harrell 2015 5-knot quantiles)
  knots <- quantile(df_o[[exp]], probs = knot_quantiles, na.rm = TRUE)
  knots <- unique(knots)
  if (length(knots) < 5) {
    cat(sprintf("  WARN: %s has <5 unique knot positions (%d), reducing\n",
                tag, length(knots)))
  }

  # RCS model (rms::ols for continuous outcome)
  f_rcs <- as.formula(sprintf(
    "%s ~ rcs(%s, parms = c(%s)) + %s",
    out, exp,
    paste(knots, collapse = ","),
    cov_pre_str
  ))
  f_lin <- as.formula(sprintf(
    "%s ~ %s + %s",
    out, exp, cov_pre_str
  ))

  fit_rcs <- tryCatch({
    rms::ols(f_rcs, data = df_o, weights = df_o$wt_pooled, x = TRUE, y = TRUE)
  }, error = function(e) {
    cat("  rms::ols (rcs) error:", conditionMessage(e), "\n"); NULL
  })
  fit_lin <- tryCatch({
    rms::ols(f_lin, data = df_o, weights = df_o$wt_pooled, x = TRUE, y = TRUE)
  }, error = function(e) {
    cat("  rms::ols (linear) error:", conditionMessage(e), "\n"); NULL
  })
  if (is.null(fit_rcs) || is.null(fit_lin)) next

  # LR test (RCS vs linear) — proper non-linearity test
  lr_chi <- (fit_lin$stats["Sigma"]^2 - fit_rcs$stats["Sigma"]^2)  # fallback
  # Better: use anova type
  aov_rcs <- tryCatch(anova(fit_rcs), error = function(e) NULL)
  p_nonlin <- if (!is.null(aov_rcs)) {
    # look for the "Nonlinear" row associated with the exposure variable
    rn <- rownames(aov_rcs)
    nl_row <- grep("Nonlinear", rn, ignore.case = TRUE)
    if (length(nl_row) > 0) aov_rcs[nl_row[1], "P"] else NA_real_
  } else NA_real_
  p_overall <- if (!is.null(aov_rcs)) {
    rn <- rownames(aov_rcs)
    ov_row <- grep(paste0("^", exp), rn)
    if (length(ov_row) > 0) aov_rcs[ov_row[1], "P"] else NA_real_
  } else NA_real_

  # LR via residual SS comparison (manual)
  rss_lin <- sum(residuals(fit_lin)^2 * fit_lin$weights, na.rm = TRUE)
  rss_rcs <- sum(residuals(fit_rcs)^2 * fit_rcs$weights, na.rm = TRUE)
  df_extra <- fit_rcs$stats["d.f."] - fit_lin$stats["d.f."]
  lr_stat <- nrow(df_o) * (log(rss_lin) - log(rss_rcs))
  p_lr <- if (df_extra > 0) {
    pchisq(lr_stat, df = df_extra, lower.tail = FALSE)
  } else NA_real_

  pvals_df <- rbind(pvals_df, data.frame(
    outcome     = out,
    exposure    = exp,
    n           = nrow(df_o),
    n_knots     = length(knots),
    p_nonlinear_anova = p_nonlin,
    p_overall_anova   = p_overall,
    p_lr_vs_linear    = p_lr,
    stringsAsFactors = FALSE
  ))

  # ---- Predicted curve with 95% CI ----
  exp_grid <- seq(quantile(df_o[[exp]], 0.025, na.rm = TRUE),
                   quantile(df_o[[exp]], 0.975, na.rm = TRUE),
                   length.out = 200)
  newdat <- data.frame(x = exp_grid)
  newdat$age        <- median(df_o$age, na.rm = TRUE)
  newdat$sex_male_i <- 1
  newdat$race_chr   <- "Non-Hispanic White"
  newdat$edu_chr    <- "College or above"
  newdat$pir        <- median(df_o$pir, na.rm = TRUE)
  newdat$smoke_chr  <- "Never"
  newdat$drink_chr  <- "No"
  names(newdat)[1] <- exp

  pred <- tryCatch({
    rms::Predict(fit_rcs, name = exp,
                 fun = identity, conf.int = 0.95,
                 # ref level set to median:
                 ref.zero = TRUE)
  }, error = function(e) {
    cat("  rms::Predict error:", conditionMessage(e), "\n"); NULL
  })

  # ---- Plot ----
  if (!is.null(pred)) {
    p_obj <- ggplot(as.data.frame(pred),
                    aes(x = .data[[exp]], y = yhat)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.25) +
      geom_line(color = "steelblue", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      labs(
        title = sprintf("RCS (5-knot) %s vs %s", out, exp),
        subtitle = sprintf("n=%d ; p_nonlinear (LR vs linear) = %.4f",
                            nrow(df_o), p_lr),
        x = exp, y = sprintf("Adjusted %s (vs reference)", out)
      ) +
      theme_minimal(base_size = 11)
    ggsave(sprintf("output/figures/rcs_se_%s.png", tag),
           plot = p_obj, width = 6.5, height = 4.5, dpi = 150)
  }

  results_list[[tag]] <- list(
    fit_rcs = fit_rcs,
    fit_lin = fit_lin,
    pred    = pred,
    knots   = knots,
    p_lr    = p_lr,
    n       = nrow(df_o)
  )
}

# ---- Save ----
write.csv(pvals_df, "output/tables/rcs_pvalues.csv", row.names = FALSE)
save(results_list, pvals_df, file = "output/tables/rcs_ushape.RData")

cat("\n保存:\n")
cat("  output/tables/rcs_ushape.RData\n")
cat("  output/tables/rcs_pvalues.csv\n")
cat("  output/figures/rcs_se_<outcome>_<exposure>.png (×6)\n")
print(pvals_df)

cat("\nDONE 07_rcs_ushape.R\n")
