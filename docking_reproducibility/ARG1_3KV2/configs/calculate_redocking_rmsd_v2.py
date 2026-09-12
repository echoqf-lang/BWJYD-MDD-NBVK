import csv
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign


if len(sys.argv) != 4:
    raise SystemExit(
        "Usage: python calculate_redocking_rmsd_v2.py "
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
reference_xyz = np.array(reference.GetConformer().GetPositions())

supplier = Chem.SDMolSupplier(
    str(poses_path), removeHs=False, sanitize=True
)

rows = []

for mode, pose in enumerate(supplier, start=1):
    if pose is None:
        raise RuntimeError(f"Could not read pose {mode}.")

    pose = Chem.RemoveHs(pose)
    pose_xyz = np.array(pose.GetConformer().GetPositions())

    matches = pose.GetSubstructMatches(
        reference,
        uniquify=False,
        useChirality=True,
        maxMatches=100000,
    )

    if not matches:
        raise RuntimeError(
            f"No topology-based atom mapping found for mode {mode}."
        )

    mapped_rmsds = []

    for match in matches:
        matched_pose_xyz = pose_xyz[list(match)]
        rmsd = np.sqrt(
            np.mean(
                np.sum(
                    (matched_pose_xyz - reference_xyz) ** 2,
                    axis=1,
                )
            )
        )
        mapped_rmsds.append(rmsd)

    topology_rmsd = min(mapped_rmsds)

    symmetry_rmsd = rdMolAlign.CalcRMS(
        pose,
        reference,
        maxMatches=100000,
        symmetrizeConjugatedTerminalGroups=True,
    )

    affinity = ""
    for property_name in pose.GetPropNames():
        if property_name.lower() in {
            "free_energy",
            "vina_affinity",
            "affinity",
        }:
            affinity = pose.GetProp(property_name)
            break

    rows.append(
        {
            "mode": mode,
            "affinity_kcal_mol": affinity,
            "topology_RMSD_A": f"{topology_rmsd:.4f}",
            "symmetry_RMSD_A": f"{symmetry_rmsd:.4f}",
            "passes_2A": symmetry_rmsd < 2.0,
            "mapping_count": len(matches),
        }
    )

with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote: {output_path}")
print("mode topology_RMSD symmetry_RMSD pass mappings")
for row in rows:
    print(
        row["mode"],
        row["topology_RMSD_A"],
        row["symmetry_RMSD_A"],
        row["passes_2A"],
        row["mapping_count"],
    )
