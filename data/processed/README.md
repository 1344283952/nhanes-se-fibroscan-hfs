# Processed Data

Re-assembled wide tables (`nhanes_*.RData`) are excluded from public release due to NHANES re-identification gray-zone considerations. Rebuild via:

```r
Rscript scripts/02_merge_data.R   # downloads + merges
Rscript scripts/03_clean_data.R   # cleans + recodes
```

Raw NHANES `.XPT` files are public domain at <https://wwwn.cdc.gov/nchs/nhanes/>. The merged wide table is bit-for-bit reproducible from these inputs using the included scripts.
