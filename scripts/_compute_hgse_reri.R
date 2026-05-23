# ============================================================
# Compute RERI / AP / SI for Hg × Se interaction on GGT outcome
# R4 P1 fix (2026-05-23)
# Reference: VanderWeele & Knol (2014) — Tutorial on interaction
# RERI = RR11 - RR10 - RR01 + 1
# AP   = RERI / RR11
# SI   = (RR11 - 1) / ((RR10 - 1) + (RR01 - 1))
# 95% CI via delta-method on logit interaction model (Hosmer-Lemeshow 1992)
# ============================================================
library(dplyr); library(survey); library(msm)
set.seed(20260523)

cat("==== Compute RERI/AP/SI for Hg × Se on GGT outcome ====\n")
load("data/processed/nhanes_final.RData")

d <- nhanes_final %>%
  filter(!is.na(hg_ugl), !is.na(se_ugl), !is.na(ggt_high),
         !is.na(RIDAGEYR), !is.na(sex_male))
cat(sprintf("Analytic n: %d\n", nrow(d)))

# Dichotomise Hg and Se at medians (VanderWeele standard practice for additive interaction)
med_hg <- median(d$hg_ugl, na.rm = TRUE)
med_se <- median(d$se_ugl, na.rm = TRUE)
d$hg_hi <- as.integer(d$hg_ugl > med_hg)
d$se_hi <- as.integer(d$se_ugl > med_se)
cat(sprintf("Median Hg: %.3f µg/L; Median Se: %.2f µg/L\n", med_hg, med_se))
cat(sprintf("Cells (Hg_hi, Se_hi): (0,0)=%d  (1,0)=%d  (0,1)=%d  (1,1)=%d\n",
            sum(d$hg_hi==0 & d$se_hi==0), sum(d$hg_hi==1 & d$se_hi==0),
            sum(d$hg_hi==0 & d$se_hi==1), sum(d$hg_hi==1 & d$se_hi==1)))

# Logistic regression: ggt_high ~ hg_hi * se_hi + covariates
covs <- c("RIDAGEYR","sex_male","race","education","pir","SMQ020","ALQ111","DR1TKCAL")
fml <- as.formula(paste("ggt_high ~ hg_hi * se_hi +",
                        paste(covs, collapse = "+")))
fit <- glm(fml, data = d, family = binomial())
cf <- summary(fit)$coefficients
print(cf[c("hg_hi","se_hi","hg_hi:se_hi"), ])

b1 <- coef(fit)["hg_hi"]
b2 <- coef(fit)["se_hi"]
b3 <- coef(fit)["hg_hi:se_hi"]

RR10 <- exp(b1)
RR01 <- exp(b2)
RR11 <- exp(b1 + b2 + b3)
RERI <- RR11 - RR10 - RR01 + 1
AP   <- RERI / RR11
SI   <- (RR11 - 1) / ((RR10 - 1) + (RR01 - 1))

# Delta-method 95% CIs (Hosmer-Lemeshow 1992 / VanderWeele-Knol 2014)
V <- vcov(fit)[c("hg_hi","se_hi","hg_hi:se_hi"),
              c("hg_hi","se_hi","hg_hi:se_hi")]
se_RERI <- deltamethod(~ exp(x1+x2+x3) - exp(x1) - exp(x2) + 1,
                       mean = c(b1,b2,b3), cov = V)
se_AP   <- deltamethod(~ (exp(x1+x2+x3) - exp(x1) - exp(x2) + 1) / exp(x1+x2+x3),
                       mean = c(b1,b2,b3), cov = V)
se_lnSI <- deltamethod(~ log((exp(x1+x2+x3) - 1) / ((exp(x1) - 1) + (exp(x2) - 1))),
                       mean = c(b1,b2,b3), cov = V)
lnSI <- log(SI)

ci_RERI <- c(RERI - 1.96*se_RERI, RERI + 1.96*se_RERI)
ci_AP   <- c(AP   - 1.96*se_AP,   AP   + 1.96*se_AP)
ci_SI   <- c(exp(lnSI - 1.96*se_lnSI), exp(lnSI + 1.96*se_lnSI))

cat(sprintf("\n=== Hg × Se RERI/AP/SI on GGT-high (n=%d) ===\n", nrow(d)))
cat(sprintf("RR10 (Hg high, Se low):  %.3f\n", RR10))
cat(sprintf("RR01 (Hg low, Se high):  %.3f\n", RR01))
cat(sprintf("RR11 (both high):        %.3f\n", RR11))
cat(sprintf("RERI: %.3f (95%% CI %.3f, %.3f)\n", RERI, ci_RERI[1], ci_RERI[2]))
cat(sprintf("AP:   %.3f (95%% CI %.3f, %.3f)\n", AP,   ci_AP[1],   ci_AP[2]))
cat(sprintf("SI:   %.3f (95%% CI %.3f, %.3f)\n", SI,   ci_SI[1],   ci_SI[2]))

# Save
out <- data.frame(
  measure = c("RERI", "AP", "SI"),
  estimate = c(RERI, AP, SI),
  ci_lo = c(ci_RERI[1], ci_AP[1], ci_SI[1]),
  ci_hi = c(ci_RERI[2], ci_AP[2], ci_SI[2]),
  interpretation = c(
    "Relative excess risk due to interaction; >0 = positive additive interaction",
    "Attributable proportion due to interaction; 0-1 = positive additive",
    "Synergy index; >1 = synergy, <1 = antagonism, =1 = additivity"),
  note = c(
    "Hosmer-Lemeshow 1992 delta-method 95% CI",
    "Hosmer-Lemeshow 1992 delta-method 95% CI",
    "Log-SI delta-method 95% CI, exponentiated"),
  stringsAsFactors = FALSE
)
write.csv(out, "output/tables/hgse_reri.csv", row.names = FALSE)
cat("\n[OK] saved output/tables/hgse_reri.csv\n")
cat("\n=== Caveats per VanderWeele-Knol (2014) ===\n")
cat("- Dichotomisation at median is a coarse simplification of continuous Hg/Se exposure;\n")
cat("  the primary metric for §3.9 remains the multiplicative OR 1.58 from the continuous\n")
cat("  Hg × Se interaction model. The additive-scale layer reported here is offered for\n")
cat("  public-health-interpretation completeness rather than as a competing primary inference.\n")
