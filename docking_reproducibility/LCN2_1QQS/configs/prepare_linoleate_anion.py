from pathlib import Path
from rdkit import Chem

source = Path(
    "LCN2_1QQS/"
    "07_linoleic_input/linoleic_acid_CID5280450_pubchem_3D.sdf"
)
output = Path(
    "LCN2_1QQS/"
    "07_linoleic_input/linoleic_acid_CID5280450_anion_H.sdf"
)

mol = Chem.MolFromMolFile(str(source), removeHs=False)
if mol is None:
    raise RuntimeError("Could not read PubChem linoleic acid SDF")

mol = Chem.RemoveHs(mol)
heavy_before = Chem.RemoveHs(mol)
conf_before = heavy_before.GetConformer()
coords_before = [
    (
        conf_before.GetAtomPosition(i).x,
        conf_before.GetAtomPosition(i).y,
        conf_before.GetAtomPosition(i).z,
    )
    for i in range(heavy_before.GetNumAtoms())
]

matches = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[O;H1]"))
if len(matches) != 1:
    raise RuntimeError(f"Expected one carboxylic acid group, found {len(matches)}")

acid_oxygen = mol.GetAtomWithIdx(matches[0][2])
acid_oxygen.SetNumExplicitHs(0)
acid_oxygen.SetNoImplicit(True)
acid_oxygen.SetFormalCharge(-1)

Chem.SanitizeMol(mol)
mol = Chem.AddHs(mol, addCoords=True)

heavy_after = Chem.RemoveHs(mol)
conf_after = heavy_after.GetConformer()
sum_sq = 0.0

for i, before in enumerate(coords_before):
    after = conf_after.GetAtomPosition(i)
    sum_sq += (
        (before[0] - after.x) ** 2
        + (before[1] - after.y) ** 2
        + (before[2] - after.z) ** 2
    )

rmsd = (sum_sq / len(coords_before)) ** 0.5

mol.SetProp("_Name", "Linoleate_CID5280450")
Chem.MolToMolFile(mol, str(output))

print(f"Wrote: {output}")
print(f"Heavy atoms: {mol.GetNumHeavyAtoms()}")
print(f"Hydrogens: {mol.GetNumAtoms() - mol.GetNumHeavyAtoms()}")
print(f"Formal charge: {Chem.GetFormalCharge(mol)}")
print(f"Heavy-atom coordinate RMSD: {rmsd:.6f} Angstrom")
print(f"Isomeric SMILES: {Chem.MolToSmiles(mol, isomericSmiles=True)}")
