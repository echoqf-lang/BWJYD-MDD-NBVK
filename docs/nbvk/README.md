# NBVK reproducibility package

This package contains the two successive network-based virtual knockout (NBVK) analyses reported in the manuscript. NBVK v1 evaluates single-, pair-, and triple-node deletion against uniform, marginal topology-matched, hub, and sensitivity baselines. NBVK v2 conditions the triple-node comparison on induced structure and pairwise distance categories.

## Fixed inputs

The scripts consume the following files in `data/network_inputs/`:

- `string_edges_score700.csv`: STRING edges at `combined_score >= 0.700`.
- `string_edges_score400.csv`: STRING edges at `combined_score >= 0.400`.
- `full_ppi_node_metrics_score700.csv`: fixed node metrics for the score-0.700 network.
- `string_id_mapping_audit.csv`: retained mapping audit for the original 1,309 candidates.
- `anchor_set_membership.csv`: membership audit for ARG1, LCN2, and LTF.

These files are analysis inputs derived from the audited STRING query. They are not evidence of direct compound-target binding.

The frozen plan retains the original working paths under `ppi_review_20260726/results/`. In this repository, those five frozen inputs map directly to files with the same basenames under `data/network_inputs/`. The historical paths are preserved in the plan as provenance and are not commands to be executed in the public repository.

## Run labels and random seeds

Human-readable run labels (`run01`, `run02`, and `run03`) are separated from numerical random seeds. See `NBVK_SEED_MAP.tsv`. These are three run records, not three unique numerical seeds: NBVK v2 deliberately reuses the fixed primary seed. Renaming a file or run label does not alter the numerical seed used by the algorithm.

## Execution

Run commands from the repository root. Passing the repository root as the last argument is optional because each script resolves it from its own location.

```bash
Rscript code/nbvk_v1/tests/test_nbvk_core.R
python code/nbvk_v2/tests/test_combination_matcher.py

Rscript code/nbvk_v1/01_run_primary_nbvk.R
Rscript code/nbvk_v1/02_rebuild_inference_from_frozen_nulls.R
Rscript code/nbvk_v1/03_run_sensitivity_nbvk.R
Rscript code/nbvk_v1/04_run_variant_topology_inference.R
python code/nbvk_v1/05_posthoc_combination_stress_test.py

Rscript code/nbvk_v2/01_prepare_network_features.R
python code/nbvk_v2/02_match_combinations.py
Rscript code/nbvk_v2/03_run_global_damage.R
```

The full permutation workflows can be computationally expensive. To audit the primary inference without regenerating the 10,000-draw null distributions, run `02_rebuild_inference_from_frozen_nulls.R` against the archived null files.

The recorded Python package versions are listed in `environment/nbvk/python_requirements.txt`; the original R session information is retained under `environment/nbvk/`.

## Manuscript mapping

- `results/nbvk_v1/primary/observed_nbvk.csv`: single-, pair-, and triple-node observed perturbations.
- `results/nbvk_v1/primary/damage_inference.csv`: NBVK v1 random and marginal topology-matched inference.
- `results/nbvk_v1/exploratory_posthoc/`: explicitly post hoc structural-context analyses.
- `results/nbvk_v2/global_damage_inference.csv`: strict and relaxed NBVK v2 triplet inference reported in Figure 4 and Table S7.
- `results/nbvk_v2/global_damage_null_relax2.csv`: exact relaxed comparison distribution.

The strict NBVK v2 comparison retained 599 triplets, below the prespecified minimum of 1,000. The adequately supported relaxed comparison did not show a significant combination-specific global network effect. These analyses do not establish pharmacological synergy.

## Provenance and deviations

The frozen design is in `preregistered_analysis_plan.md`. Implementation corrections and deviations are retained in `nbvk_v1_deviations.md`, `nbvk_v1_implementation_audit.md`, `nbvk_v2_deviations.md`, and `nbvk_v2_structural_checkpoint.md`. The original numerical seeds are recorded rather than concealed or changed.
