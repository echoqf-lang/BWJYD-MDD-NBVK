import csv
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign


if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: python calculate_redocking_rmsd.py "
        "REFERENCE.sdf POSES.sdf OUTPUT.csv"
    )

reference_path = Path(sys.argv[1])
poses_path = Path(sys.argv[2])
output_path = Path(sys.argv[3])

reference = Chem.MolFromMolFile(
    str(reference_path), removeHs=False, sanitize=True
)
if reference is None:
    raise RuntimeError("Could not read reference SDF.")

reference = Chem.RemoveHs(reference)
reference_symbols = [atom.GetSymbol() for atom in reference.GetAtoms()]
reference_xyz = np.array(reference.GetConformer().GetPositions())

supplier = Chem.SDMolSupplier(
    str(poses_path), removeHs=False, sanitize=True
)

rows = []

for mode, pose in enumerate(supplier, start=1):
    if pose is None:
        raise RuntimeError(f"Could not read pose {mode}.")

    pose = Chem.RemoveHs(pose)
    pose_symbols = [atom.GetSymbol() for atom in pose.GetAtoms()]

    if pose_symbols != reference_symbols:
        raise RuntimeError(
            f"Atom-order mismatch in mode {mode}:\n"
            f"reference={reference_symbols}\npose={pose_symbols}"
        )

    pose_xyz = np.array(pose.GetConformer().GetPositions())
    mapped_rmsd = np.sqrt(
        np.mean(np.sum((pose_xyz - reference_xyz) ** 2, axis=1))
    )

    symmetry_rmsd = rdMolAlign.CalcRMS(
        pose,
        reference,
        maxMatches=100000,
        symmetrizeConjugatedTerminalGroups=True,
    )

    affinity = ""
    for property_name in (
        "free_energy",
        "vina_affinity",
        "Affinity",
        "affinity",
    ):
        if pose.HasProp(property_name):
            affinity = pose.GetProp(property_name)
            break

    rows.append(
        {
            "mode": mode,
            "affinity_kcal_mol": affinity,
            "mapped_RMSD_A": f"{mapped_rmsd:.4f}",
            "symmetry_RMSD_A": f"{symmetry_rmsd:.4f}",
            "passes_2A": symmetry_rmsd < 2.0,
        }
    )

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote: {output_path}")
print("mode mapped_RMSD symmetry_RMSD pass")
for row in rows:
    print(
        row["mode"],
        row["mapped_RMSD_A"],
        row["symmetry_RMSD_A"],
        row["passes_2A"],
    )
