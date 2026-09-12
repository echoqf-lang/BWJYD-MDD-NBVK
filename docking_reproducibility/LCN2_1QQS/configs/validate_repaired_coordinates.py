import math
from pathlib import Path

old_path = Path(
    "LCN2_1QQS/"
    "01_receptor_clean/LCN2_1QQS_chainA_clean.pdb"
)
new_path = Path(
    "LCN2_1QQS/"
    "01_receptor_clean/LCN2_1QQS_chainA_sidechains_repaired.pdb"
)


def read_atoms(path):
    atoms = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            key = (
                line[21:22],
                line[22:26].strip(),
                line[12:16].strip(),
            )
            atoms[key] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
    return atoms


old_atoms = read_atoms(old_path)
new_atoms = read_atoms(new_path)
common = set(old_atoms) & set(new_atoms)

squared_distances = []
for key in common:
    old_xyz = old_atoms[key]
    new_xyz = new_atoms[key]
    squared_distances.append(
        sum((old_xyz[i] - new_xyz[i]) ** 2 for i in range(3))
    )

rmsd = math.sqrt(sum(squared_distances) / len(squared_distances))
maximum = math.sqrt(max(squared_distances))

print(f"Original atoms compared: {len(common)}")
print(f"Coordinate RMSD: {rmsd:.6f} Angstrom")
print(f"Maximum deviation: {maximum:.6f} Angstrom")
print(f"New atoms added: {len(new_atoms) - len(old_atoms)}")
