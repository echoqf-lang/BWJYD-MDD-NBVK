from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolAlign


ROOT = Path("./ARG1_3KV2")
SEEDS = ["20260727", "20260728", "20260729"]

poses = {}

for seed in SEEDS:
    path = (
        ROOT
        / "08_puerarin_docking"
        / f"puerarin_seed{seed}_poses.sdf"
    )
    supplier = Chem.SDMolSupplier(
        str(path), removeHs=False, sanitize=True
    )
    pose = supplier[0]
    if pose is None:
        raise RuntimeError(f"Could not read top pose for seed {seed}.")
    poses[seed] = Chem.RemoveHs(pose)

for i, seed_a in enumerate(SEEDS):
    for seed_b in SEEDS[i + 1 :]:
        rmsd = rdMolAlign.CalcRMS(
            poses[seed_a],
            poses[seed_b],
            maxMatches=100000,
            symmetrizeConjugatedTerminalGroups=True,
        )
        print(f"{seed_a} vs {seed_b}: {rmsd:.4f} Angstrom")
