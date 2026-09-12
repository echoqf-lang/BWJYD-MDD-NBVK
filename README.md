# BWJYD-MDD-NBVK

Reproducible data and code accompanying the computational evaluation of Bawei Jieyu Decoction (BWJYD) in major depressive disorder (MDD), including network pharmacology, cross-cohort transcriptomic analysis, network-based virtual knockout (NBVK), molecular docking, and ADMET prediction.

## Repository status

This repository is being prepared for manuscript submission. Files and documentation are undergoing provenance, licensing, privacy, and reproducibility review. Results should not be treated as a final archival release until a versioned release and DOI are provided.

## Repository contents

- `code/nbvk_v1/` and `code/nbvk_v2/`: network virtual-knockout analysis code and tests.
- `code/geo/`: GEO retrieval, preprocessing, cohort modelling, meta-analysis, and brain-validation code.
- `data/network_inputs/`: fixed PPI edge tables and node metrics consumed by the NBVK analyses.
- `data/geo/`: GEO accession documentation and complete gene-level derived effect tables.
- `results/nbvk_v1/` and `results/nbvk_v2/`: frozen NBVK outputs, including null distributions.
- `results/geo/`: cohort metadata and compact transcriptomic result tables underlying Figure 3 and Tables S5–S6.
- `environment/`: software versions and reproducible environment specifications.
- `docking_reproducibility/`: audited docking inputs, configurations, validation scripts, and checksums.
- `docs/nbvk/`: frozen analysis plan, deviations, audit notes, evidence boundaries, and rerun instructions.

NBVK-specific execution order and file-to-manuscript mapping are documented in `docs/nbvk/README.md`.

## Quick start

Run the following commands from the repository root to test the core implementation and rebuild the archived NBVK inference tables:

```bash
Rscript code/nbvk_v1/tests/test_nbvk_core.R
python code/nbvk_v2/tests/test_combination_matcher.py
Rscript code/nbvk_v1/02_rebuild_inference_from_frozen_nulls.R
```

The complete execution order, including computationally intensive permutation analyses, is provided in `docs/nbvk/README.md`. Software versions and fixed numerical random seeds are recorded under `environment/` and `docs/nbvk/`.

## Data provenance

Public transcriptomic datasets are referenced by GEO accession number rather than duplicated as raw downloads. Complete author-derived gene-level effect tables are provided under `data/geo/full_effect_tables/`, with compact manuscript-facing summaries under `results/geo/`. The fixed network inputs used by NBVK are provided in `data/network_inputs/`. Database-derived files remain subject to their source databases' terms of use.

## Evidence boundary

The repository supports reproducibility of computational analyses. Network associations, transcriptomic observations, docking poses, and in silico ADMET predictions do not by themselves establish therapeutic efficacy, direct target engagement, or causal pharmacological mechanisms.

## Citation

Machine-readable citation metadata are provided in `CITATION.cff`. The associated article citation and an archival DOI will be added to the versioned public release.

## License

Original analysis code and original project documentation are released under the MIT License. Third-party database records, molecular structures, and software outputs remain subject to their source-specific terms. See `THIRD_PARTY_NOTICES.md` before reuse or redistribution.
