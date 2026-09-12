# GEO data provenance

Raw expression files are not redistributed. They can be obtained from the NCBI Gene Expression Omnibus using the following accession numbers.

| Accession | Tissue and role | Reported analysis sample |
| --- | --- | --- |
| GSE98793 | Whole blood, primary cohort | 128 MDD and 64 control samples |
| GSE251778 | Whole blood, primary cohort | 80 MDD and 89 control samples; RNA-seq subseries |
| GSE19738 | Whole blood, primary cohort | 33 MDD and 34 control basal samples; stimulated repeats excluded |
| GSE39653 | PBMC, separate validation cohort | 21 MDD and 24 control samples; bipolar samples excluded |
| GSE101521 | DLPFC, separate brain validation cohort | Complete-case model: 28 MDD and 28 control samples |

`full_effect_tables/` contains the complete gene-level cohort results and whole-blood meta-analysis table underlying the reported transcriptomic findings. These derived tables are provided so that reported effect sizes, P values, multiple-testing results, and Figure 3 can be checked without redistributing the original expression matrices.

GSE101521 contains more eligible total-RNA donors before complete-case covariate filtering; the manuscript reports the 56 samples used in the age-, sex-, RIN-, and PMI-adjusted complete-case model.
