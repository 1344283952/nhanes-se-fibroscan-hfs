# Verification helper - re-run after any clean step to confirm numerical truth
# Outputs key data integrity numbers needed by manuscript
load("data/processed/nhanes_final.RData")
cat("=== 006 Data integrity verification ===\n")
cat(sprintf("Cohort N = %d\n\n", nrow(nhanes_final)))

cat(sprintf("LBXBSE  : median=%.2f  IQR Q1=%.2f  Q3=%.2f  n=%d\n",
  median(nhanes_final$LBXBSE, na.rm=TRUE),
  as.numeric(quantile(nhanes_final$LBXBSE, 0.25, na.rm=TRUE)),
  as.numeric(quantile(nhanes_final$LBXBSE, 0.75, na.rm=TRUE)),
  sum(!is.na(nhanes_final$LBXBSE))))
cat(sprintf("DR1TSELE: median=%.2f  IQR Q1=%.2f  Q3=%.2f  n=%d\n\n",
  median(nhanes_final$DR1TSELE, na.rm=TRUE),
  as.numeric(quantile(nhanes_final$DR1TSELE, 0.25, na.rm=TRUE)),
  as.numeric(quantile(nhanes_final$DR1TSELE, 0.75, na.rm=TRUE)),
  sum(!is.na(nhanes_final$DR1TSELE))))

cat(sprintf("CAP>=275 n=%d (%.2f%%)\n",
  sum(nhanes_final$steatosis_cap275, na.rm=TRUE),
  100*mean(nhanes_final$steatosis_cap275, na.rm=TRUE)))
cat(sprintf("LSM>=8   n=%d (%.2f%%)\n",
  sum(nhanes_final$fibrosis_lsm8, na.rm=TRUE),
  100*mean(nhanes_final$fibrosis_lsm8, na.rm=TRUE)))
cat(sprintf("LSM>=12  n=%d (%.2f%%)\n\n",
  sum(nhanes_final$fibrosis_lsm12, na.rm=TRUE),
  100*mean(nhanes_final$fibrosis_lsm12, na.rm=TRUE)))

cat("HFS chain (sequential intersections):\n")
cat(sprintf("  step 1: hfs non-NA                                     = %d\n",
  sum(!is.na(nhanes_final$hfs))))
cat(sprintf("  step 2: + fibrosis_lsm8 non-NA                         = %d\n",
  sum(!is.na(nhanes_final$hfs) & !is.na(nhanes_final$fibrosis_lsm8))))
cat(sprintf("  step 3: + wt_saf_pooled non-NA & >0                    = %d\n",
  sum(!is.na(nhanes_final$hfs) & !is.na(nhanes_final$fibrosis_lsm8) &
      !is.na(nhanes_final$wt_saf_pooled) & nhanes_final$wt_saf_pooled>0)))
cat(sprintf("  step 4: + fib4 non-NA (= ROC denom; AUROC analytic N)  = %d\n",
  sum(!is.na(nhanes_final$hfs) & !is.na(nhanes_final$fibrosis_lsm8) &
      !is.na(nhanes_final$wt_saf_pooled) & nhanes_final$wt_saf_pooled>0 &
      !is.na(nhanes_final$fib4))))

hfs_anal_set <- !is.na(nhanes_final$hfs) & !is.na(nhanes_final$fibrosis_lsm8) &
                !is.na(nhanes_final$wt_saf_pooled) & nhanes_final$wt_saf_pooled>0 &
                !is.na(nhanes_final$fib4)
cat(sprintf("\nIn ROC analytic set (HFS+LSM8+FIB4+wt_saf>0): n=%d\n", sum(hfs_anal_set)))
cat(sprintf("  HFS>=0.47 in this set: n=%d\n",
  sum(hfs_anal_set & nhanes_final$hfs >= 0.47)))
cat(sprintf("  LSM>=8 events in this set: n=%d (%.2f%%)\n",
  sum(hfs_anal_set & nhanes_final$fibrosis_lsm8==1),
  100*mean(nhanes_final$fibrosis_lsm8[hfs_anal_set]==1)))

cat("\nMetALD groups (CAP>=275 primary):\n")
print(table(nhanes_final$metald_group, useNA="ifany"))
