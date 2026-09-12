import math
from pathlib import Path

pdb_path = Path(
    "LCN2_1QQS/"
    "00_raw_inputs/1QQS_raw.pdb"
)

protein = []
ligand = []

with pdb_path.open() as handle:
    for line in handle:
        record = line[0:6].strip()
        chain = line[21:22]
        residue = line[17:20]
        resid = line[22:26].strip()

        xyz = (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        ) if record in {"ATOM", "HETATM"} else None

        if record == "ATOM" and chain == "A":
            protein.append(xyz)

        if (
            record == "HETATM"
            and chain == "A"
            and residue == "DKA"
            and resid == "181"
        ):
            ligand.append(xyz)

partner = [
    (
        y + 27.145,
        -x + 81.435,
        z - 30.495,
    )
    for x, y, z in protein
]


def minimum_distance(first, second):
    return min(
        math.dist(atom1, atom2)
        for atom1 in first
        for atom2 in second
    )


print(f"DKA atoms: {len(ligand)}")
print(f"Original protein atoms: {len(protein)}")
print(
    "Minimum DKA-to-original-chain distance:",
    f"{minimum_distance(ligand, protein):.3f} Angstrom",
)
print(
    "Minimum DKA-to-BIOMT-partner distance:",
    f"{minimum_distance(ligand, partner):.3f} Angstrom",
)
print(
    "Partner atoms within 4.0 Angstrom:",
    sum(
        any(math.dist(atom, lig) <= 4.0 for lig in ligand)
        for atom in partner
    ),
)
