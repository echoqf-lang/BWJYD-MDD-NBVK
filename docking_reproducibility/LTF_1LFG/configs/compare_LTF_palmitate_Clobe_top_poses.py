import math
from pathlib import Path
from rdkit import Chem

root = Path(
    "./"
    "06_Clobe_docking"
)

files = {
    "20260727": root / "palmitate_Clobe_run01_poses.sdf",
    "20260728": root / "palmitate_Clobe_run02_poses.sdf",
    "20260729": root / "palmitate_Clobe_run03_poses.sdf",
}


def first_heavy_pose(path):
    supplier = Chem.SDMolSupplier(str(path), removeHs=True)
    mol = supplier[0]
    if mol is None:
        raise RuntimeError(f"Could not read first pose from {path}")

    conformer = mol.GetConformer()
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coordinates = [
        (
            conformer.GetAtomPosition(i).x,
            conformer.GetAtomPosition(i).y,
            conformer.GetAtomPosition(i).z,
        )
        for i in range(mol.GetNumAtoms())
    ]
    return elements, coordinates


poses = {seed: first_heavy_pose(path) for seed, path in files.items()}

pairs = [
    ("20260727", "20260728"),
    ("20260727", "20260729"),
    ("20260728", "20260729"),
]

for first, second in pairs:
    elements_a, coordinates_a = poses[first]
    elements_b, coordinates_b = poses[second]

    if elements_a != elements_b:
        raise RuntimeError(f"Atom-order mismatch: {first} vs {second}")

    squared = [
        sum((a[i] - b[i]) ** 2 for i in range(3))
        for a, b in zip(coordinates_a, coordinates_b)
    ]
    rmsd = math.sqrt(sum(squared) / len(squared))
    maximum = math.sqrt(max(squared))

    print(
        f"{first} vs {second}: "
        f"RMSD={rmsd:.4f} Angstrom, "
        f"maximum_atom_deviation={maximum:.4f} Angstrom"
    )
