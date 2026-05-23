# Selenium U-shape, Mercury-Selenium Antagonism, and Five-Metal Mixture in Steatotic Liver Disease (NHANES 2017–2023)

Reproducible analytic pipeline for the manuscript:

> **Selenium U-shape, mercury–selenium antagonism, and a five-metal mixture in steatotic liver disease: NHANES 2017–2023 with multi-cycle validation**
> Submitted to *Environment International* (Elsevier, IF 11.0).

## What this repository contains

The full R pipeline for our analysis of dual-source selenium (whole-blood + dietary), four co-exposed trace metals (Pb, Cd, Hg, Mn), and Hg-Se molar antagonism against FibroScan-defined hepatic steatosis (CAP ≥ 275 dB/m), significant fibrosis (LSM ≥ 8 kPa), and the Hepamet Fibrosis Score (HFS) in **5,885 NHANES Pre-pandemic (2017–March 2020) US adults**, with NHANES cycle K (Aug 2021–Aug 2023, *n* = 3,310) external validation.

- **Primary exposures**: whole-blood Se (LBXBSE), dietary Se (DR1TSELE), Hg/Se molar ratio (Ralston 2018), 5-metal mixture (Pb / Cd / Hg / Mn / Se)
- **Primary outcomes**: CAP ≥ 275 dB/m, LSM ≥ 8 kPa, continuous HFS / FIB-4 / NFS / APRI
- **Methods**: survey-weighted GAM tensor-product + RCS 5-knot non-linearity; qgcomp + WQS mixture triangulation; 4-metal BKMR (drop-Mn) with convergence diagnostics; CMAverse 4-way decomposition for Se → GGT → fibrosis; cycle K external validation under TRIPOD-AI 2024

## Repository layout

```
scripts/             R analytic pipeline (run 00 → 23)
  00_install_packages.R    one-time install
  01_download_data.R       download NHANES + biomarker modules
  02_merge_data.R          merge across two cycles
  03_clean_data.R          exposure + outcome coding
  04_survey_design.R       svydesign setup
  05_table1.R              Table 1 baseline
  06_gam_dual_exposure.R   dual-source GAM tensor + RCS
  07_rcs_ushape.R          5-knot RCS U-shape
  08_metald_stratified.R   MetALD substratum
  09_dag.R                 directed acyclic graph
  09_hfs_predict.R         HFS calibration + ROC
  10_subgroup_forest.R     stratified + interaction
  11_sensitivity.R         9 sensitivity analyses
  12_ratio_analysis.R      Se/Zn dietary ratio
  13_cross_classification.R high/low exposure × outcome
  14_consort.R             CONSORT flow chart
  14_evalue.R              VanderWeele-Ding E-value
  15_mr_se_masld.R         MR (literature-cite reference)
  15_primary_fdr.R         BH-FDR within outcome family
  16_mice_mi.R             multiple imputation m=20
  17_nadir_boot.R          1000× nadir bootstrap
  18_mixture_qgcomp.R      qgcomp mixture estimator
  19_mixture_wqs.R         WQS mixture estimator
  20_hgse_antagonism.R     Hg-Se molar antagonism
  21_cmaverse_redox.R      4-way decomposition
  22_mr_refit.R            MR local refit (token-gated)
  23_bkmr_se_metals_checkpoint.R   4-metal BKMR
  run_all.R                end-to-end orchestrator
data/processed/      intermediate .RData (start here to skip 01-03)
output/
  tables/            primary + supplementary tables (CSV + XLSX)
  figures/           CONSORT, DAG, RCS, GAM heatmap, forest, BKMR
references.bib       BibTeX entries (cited in manuscript)
```

## How to reproduce

Software: **R 4.6** + **Rtools 4.6** (for compiling packages).

```r
# 1. Install packages (one-time, ~30-60 min)
Rscript scripts/00_install_packages.R

# 2. Download raw NHANES data (~15 min)
Rscript scripts/01_download_data.R

# 3. Run the full pipeline (~30-45 min)
Rscript scripts/run_all.R
```

Or, to skip the data download/merge step, load the analytic sample directly:

```r
load("data/processed/nhanes_main.RData")  # → 5,885 main cohort
load("data/processed/nhanes_design_main.RData") # → svydesign object
# then run any scripts from 05_table1.R onward
```

## Reproducibility notes

- All scripts initialise with `set.seed(20260516)`.
- BKMR runs 4-metal drop-Mn (10,000 iter × 2 chains); convergence diagnostics reported in §3.8 + Supplementary §S15.
- CMAverse 4-way decomposition uses 500 bootstrap iterations.
- All survey-weighted analyses use the Pre-pandemic file (P_) 3.2-year sample weight (WTMECPRP).

## Data

- **NHANES 2017–March 2020 Pre-pandemic file (P_)**: public domain. CDC/NCHS. https://wwwn.cdc.gov/nchs/nhanes/
- **NHANES cycle K (Aug 2021–Aug 2023)**: public domain.

This repository does **not** include the raw `.XPT` NHANES files. `scripts/01_download_data.R` downloads them into `data/raw/` (which is `.gitignore`-d).

## Authors

Jie Li (first author), Xiubo Sun, Jing Zhang, Lijie Zhai (Department of Pharmacy / Department of Obstetrics and Gynecology, The Second Hospital of Jilin University, Changchun, Jilin Province, China); **Ling Yu** (corresponding author, Department of Pharmacy, yulingyxb@jlu.edu.cn, ORCID 0000-0001-7362-3581).

### Funding
Jilin Provincial Higher Education Research Project (Grant No. JGJX2021D37), awarded to Ling Yu.

## License: MIT

## Declaration of Generative AI and AI-assisted technologies in the writing process

AI-assisted writing tools were used in two scope-limited ways:

(i) code implementation for the pre-submission sensitivity analyses (Fine-Gray competing-risk model, g-formula causal mediation, qgcomp, WQS, Hg-Se antagonism, CMAverse, and 4-metal BKMR) after the analytic plan was finalised by the authors;

(ii) sentence-level language polishing of the Methods and Results sections only.

The Background, Discussion, and Conclusions sections were authored without AI assistance. All study design choices, statistical model selection, scientific decisions, citation choices, numerical claims, and interpretations were independently made by the authors and verified against the underlying statistical outputs. Disclosure conforms to COPE 2023 guidance on authorship and AI tools (accessed 2026).
