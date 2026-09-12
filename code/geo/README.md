# GEO analysis code

This directory contains the R scripts used to retrieve GEO records, audit sample metadata, preprocess the reported microarray and RNA-seq cohorts, fit cohort-specific models, perform the whole-blood meta-analysis, and analyse the DLPFC validation cohort.

Run the scripts from the repository root. The reported cohorts are GSE98793, GSE251778, GSE19738, GSE39653, and GSE101521. The archived scripts also contain audit definitions for datasets evaluated during study development but not used in the final five-cohort manuscript analysis; those additional datasets must not be interpreted as reported validation cohorts.

The scripts expect downloaded GEO files under `data/raw_geo/` and generate intermediate expression matrices under `data/clean/`. These potentially large, reproducible downloads and matrices are not versioned. The complete gene-level effect tables used for the reported results are archived under `data/geo/full_effect_tables/`.

Core execution order:

```bash
Rscript code/geo/02_download_geo.R
Rscript code/geo/03_build_metadata.R
Rscript code/geo/04_preprocess_microarray.R
Rscript code/geo/05_preprocess_rnaseq.R
Rscript code/geo/06_fit_cohort_models.R
Rscript code/geo/07_meta_analysis.R
Rscript code/geo/09_brain_validation.R
```

The exact numerical random seed is retained in the scripts as a reproducibility parameter. Complete package versions from the analysis environment are archived in `environment/geo/environment.lock.tsv`.
