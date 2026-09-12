import csv
import math
from pathlib import Path

from rdkit import Chem


RECEPTOR = Path(
    "./01_receptor_clean/LTF_1LFG_chainA_2Fe_2CO3_repaired.pdb"
)
SEARCH_DIRS = [
    Path("./05_Nlobe_docking"),
    Path("./06_Clobe_docking"),
]
OUTPUT = Path(
    "./07_sensitivity_summary/LTF_top_pose_geometry_QC.csv"
)


def distance(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def read_receptor(path):
    protein = []
    receptor_nonhydrogen = []
    iron = []

    with path.open() as handle:
        for line in handle:
            record = line[0:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue

            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain = line[21:22].strip()
            residue_id = line[22:26].strip()
            element = line[76:78].strip()

            if not element:
                element = "".join(
                    character for character in atom_name if character.isalpha()
                )[:1]

            if element.upper() == "H":
                continue

            coordinates = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
            atom = {
                "name": atom_name,
                "residue": residue_name,
                "chain": chain,
                "residue_id": residue_id,
                "coordinates": coordinates,
            }

            receptor_nonhydrogen.append(atom)

            if record == "ATOM":
                protein.append(atom)

            if residue_name == "FE":
                iron.append(atom)

    return protein, receptor_nonhydrogen, iron


def read_first_pose(path):
    supplier = Chem.SDMolSupplier(str(path), removeHs=True)
    molecule = supplier[0]

    if molecule is None:
        raise RuntimeError(f"Could not read first pose from {path}")

    conformer = molecule.GetConformer()
    atoms = []

    for index, atom in enumerate(molecule.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            continue

        position = conformer.GetAtomPosition(index)
        atoms.append(
            {
                "name": f"{atom.GetSymbol()}{index + 1}",
                "coordinates": (position.x, position.y, position.z),
            }
        )

    return atoms


protein_atoms, receptor_atoms, iron_atoms = read_receptor(RECEPTOR)

if len(iron_atoms) != 2:
    raise RuntimeError(
        f"Expected 2 Fe atoms in receptor, found {len(iron_atoms)}"
    )

rows = []

for directory in SEARCH_DIRS:
    for path in sorted(directory.glob("*_poses.sdf")):
        ligand_atoms = read_first_pose(path)

        minimum_protein = float("inf")
        nearest_ligand = ""
        nearest_protein = ""
        clashes_below_1_5 = 0
        contacts_below_2_0 = 0

        for ligand_atom in ligand_atoms:
            for protein_atom in protein_atoms:
                current = distance(
                    ligand_atom["coordinates"],
                    protein_atom["coordinates"],
                )

                if current < minimum_protein:
                    minimum_protein = current
                    nearest_ligand = ligand_atom["name"]
                    nearest_protein = (
                        f'{protein_atom["name"]} '
                        f'{protein_atom["chain"]}:'
                        f'{protein_atom["residue"]}:'
                        f'{protein_atom["residue_id"]}'
                    )

                if current < 1.5:
                    clashes_below_1_5 += 1

                if current < 2.0:
                    contacts_below_2_0 += 1

        minimum_receptor = min(
            distance(ligand_atom["coordinates"], receptor_atom["coordinates"])
            for ligand_atom in ligand_atoms
            for receptor_atom in receptor_atoms
        )

        fe_distances = []
        for fe_atom in iron_atoms:
            minimum_fe = min(
                distance(
                    ligand_atom["coordinates"],
                    fe_atom["coordinates"],
                )
                for ligand_atom in ligand_atoms
            )
            fe_distances.append(minimum_fe)

        name = path.stem.replace("_poses", "")
        parts = name.split("_")
        ligand = parts[0]
        lobe = parts[1]
        seed = parts[2].replace("seed", "")

        metal_flag = (
            "YES"
            if min(fe_distances) <= 3.0
            else "NEAR"
            if min(fe_distances) <= 4.0
            else "NO"
        )

        rows.append(
            {
                "ligand": ligand,
                "lobe": lobe,
                "seed": seed,
                "minimum_protein_distance_A": f"{minimum_protein:.3f}",
                "nearest_pair": f"{nearest_ligand}--{nearest_protein}",
                "protein_contacts_below_1.5A": clashes_below_1_5,
                "protein_contacts_below_2.0A": contacts_below_2_0,
                "minimum_full_receptor_distance_A": f"{minimum_receptor:.3f}",
                "minimum_distance_to_FE701_A": f"{fe_distances[0]:.3f}",
                "minimum_distance_to_FE702_A": f"{fe_distances[1]:.3f}",
                "possible_metal_artifact": metal_flag,
                "source_file": str(path),
            }
        )

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote: {OUTPUT}")
print(f"Top poses checked: {len(rows)}")

for row in rows:
    print(
        row["ligand"],
        row["lobe"],
        row["seed"],
        f'protein_min={row["minimum_protein_distance_A"]}',
        f'below_1.5A={row["protein_contacts_below_1.5A"]}',
        f'below_2.0A={row["protein_contacts_below_2.0A"]}',
        f'FE701={row["minimum_distance_to_FE701_A"]}',
        f'FE702={row["minimum_distance_to_FE702_A"]}',
        f'metal_flag={row["possible_metal_artifact"]}',
    )
