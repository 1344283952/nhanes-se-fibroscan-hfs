# 006: Selenium x FibroScan x HFS (NHANES Pre-pandemic 2017-March 2020)

Replication code for: Dual-source selenium (dietary DR1TSELE + serum LBXBSE) x hepatic steatosis/stiffness x Hepamet HFS prediction.

## Key findings
- H2 confirmed: Se U-shape on CAP via mgcv tensor *P* = 1.27e-11
- H5 confirmed: Hepamet HFS AUROC = 0.731; dAUROC vs FIB-4 = 0.057 (95% CI 0.034-0.081)

## Reproducibility
All scripts run from project root with `Rscript scripts/0X_*.R` or `Rscript scripts/run_all.R`.
Raw NHANES data is public domain - re-download via `scripts/01_download_data.R`.

## License: MIT
## AI disclosure: code + draft assisted by Claude Opus 4.7 (1M context) per COPE 2025.
