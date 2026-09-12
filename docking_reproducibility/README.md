# Supplementary Data S1

## Docking inputs, configurations, validation outputs, and per-seed poses

This archive accompanies the manuscript on the in-silico evaluation of Bawei Jieyu Decoction in major depressive disorder. It contains the audited docking materials used for Figure 5 and Tables S8–S9.

## Reproducibility scope

The archive preserves the corrected three-target workflow finalized on 27 July 2026:

- ARG1: PDB 3KV2; co-crystallized NNH (nor-NOHA) redocking validation; puerarin (PubChem CID 5281807) candidate docking.
- LCN2: PDB 1QQS; DKA (decanoic acid) redocking audit; linoleate (PubChem CID 5280450) exploratory pose generation.
- LTF: PDB 1LFG; frozen fpocket N- and C-lobe boxes; linoleate with palmitate and stearate fatty-acid class controls.

Public filenames use neutral run identifiers rather than date-like seed strings. The exact AutoDock Vina seed values remain fixed and must be supplied unchanged to reproduce the stochastic searches:

| Run identifier | AutoDock Vina seed |
| --- | ---: |
| `run01` | `20260727` |
| `run02` | `20260728` |
| `run03` | `20260729` |

These integers are algorithm parameters, not dates, molecular identifiers, or pose identifiers. The machine-readable mapping is provided in `RUN_SEED_MAP.tsv`.

Stearate required a separate three-dimensional coordinate-generation step. RDKit ETKDGv3 used `randomSeed=20260727`, followed by MMFF optimization. This RDKit seed controls initial conformer generation and is distinct from the three Vina docking seeds.

## Main software and fixed docking parameters

- AutoDock Vina 1.2.7.
- ARG1 candidate docking: 22 × 22 × 22 Å co-crystal-centred box; exhaustiveness 64; 20 modes; energy range 6 kcal/mol.
- LCN2 candidate docking: 26 × 22 × 22 Å DKA-centred box; exhaustiveness 64; 20 modes; energy range 6 kcal/mol.
- LTF candidate/control docking: frozen N-lobe and C-lobe boxes; exhaustiveness 64; 20 modes; energy range 6 kcal/mol.
- Vina runs: `run01`, `run02`, and `run03`, corresponding to seeds `20260727`, `20260728`, and `20260729`, respectively.
- Representative poses were selected using cross-seed recurrence and geometry review, not docking score alone.

The target-specific box centres, receptor-preparation rules, scripts, RMSD calculations, and sensitivity summaries are retained within the corresponding target directories.

## Directory structure

### `ARG1_3KV2/`

- `00_raw_inputs/`: downloaded PDB structure.
- `01_receptor_clean/` and `02_receptor_pdbqt/`: prepared receptor files.
- `03_NNH_reference/` and `04_NNH_pdbqt/`: crystallographic validation ligand.
- `05_redocking/`: complete NNH outputs for all three seeds.
- `06_RMSD_validation/`: symmetry-corrected RMSD results and validation status.
- `07_puerarin_input/`: audited puerarin structure and docking input.
- `08_puerarin_docking/`: complete puerarin outputs for all three seeds.
- `09_PLIP_PyMOL/`: representative complex and interaction reports used for Figure 5.
- `configs/`: box definition, preparation, RMSD, comparison, and rendering scripts.

### `LCN2_1QQS/`

- `00_raw_inputs/` to `04_DKA_pdbqt/`: receptor and DKA validation inputs.
- `05_redocking/` and `06_RMSD_validation/`: unbiased/focused validation attempts and the recorded validation failure.
- `07_linoleic_input/`: audited linoleate input and explicit exploratory-status record.
- `08_linoleic_docking/`: complete linoleate outputs for all three seeds.
- `09_PLIP_PyMOL/`: representative exploratory complex and interaction reports.
- `configs/`: receptor-repair rules, box definitions, preparation, RMSD, and comparison scripts.

### `LTF_1LFG/`

- `00_raw_inputs/` to `04_ligand_pdbqt/`: receptor and three audited fatty-acid inputs.
- `05_Nlobe_docking/` and `06_Clobe_docking/`: complete linoleate, palmitate, and stearate outputs for all three seeds in both frozen boxes.
- `07_sensitivity_summary/`: score ranges, cross-seed convergence, geometry checks, and interpretation boundary.
- `08_PLIP_PyMOL/`: representative N- and C-lobe linoleate complexes and interaction reports.
- `configs/`: frozen protocol, preparation, comparison, quality-control, and rendering scripts.
- `pocket_detection/`: fpocket output retained to document prospective pocket selection.

### `SUMMARY/`

- `ARG1_LCN2_LTF_docking_summary.txt`: audited three-target narrative and numerical summary.
- `Table_S8_ligand_registry.csv`: ligand identities and structure provenance.
- `Table_S9_docking_validation_and_results.csv`: validation status and target-specific docking results.

### Root files

- `FILE_MANIFEST.tsv`: relative path and file size for every archived file.
- `SHA256SUMS.txt`: SHA-256 checksum for every archived file except the checksum file itself.

## Interpretation boundaries

- The ARG1 3KV2/NNH protocol passed the predefined redocking criterion in all three seeds; the puerarin poses support a reproducible structural hypothesis only.
- The LCN2 unbiased DKA redocking protocol failed. LCN2–linoleate files are retained as exploratory pose-generation outputs, not validated binding evidence.
- LTF had no fatty-acid co-crystal reference, and the class-control score ranges overlapped. LTF files are exploratory and do not establish linoleate selectivity.
- Vina scores must not be compared across different receptors.
- PLIP contacts and docking poses are computational annotations, not experimental target engagement.
- The archive excludes the superseded seed-42 screening runs and the obsolete ligand-free molecular-dynamics workflow.

## Relationship to the manuscript

- Figure 5A summarizes the target-specific evidence hierarchy.
- Figure 5B–E uses representative complexes stored in the target-specific `09_PLIP_PyMOL/` or `08_PLIP_PyMOL/` directories.
- Table S8 provides the ligand registry.
- Table S9 summarizes validation and candidate docking results.

When redistributing raw PDB or PubChem-derived files, users should retain the original database identifiers and cite the corresponding source databases and structure records.
