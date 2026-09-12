from pathlib import Path
from rdkit import Chem

source = Path(
    "LCN2_1QQS/"
    "03_DKA_reference/DKA_A181_RCSB_instance.sdf"
)
output = Path(
    "LCN2_1QQS/"
    "03_DKA_reference/DKA_A181_crystal_anion_H.sdf"
)

mol = Chem.MolFromMolFile(str(source), removeHs=False)
if mol is None:
    raise RuntimeError("Could not read the DKA reference SDF")

heavy_before = Chem.RemoveHs(mol).GetConformer()
coords_before = [
    (
        heavy_before.GetAtomPosition(i).x,
        heavy_before.GetAtomPosition(i).y,
        heavy_before.GetAtomPosition(i).z,
    )
    for i in range(heavy_before.GetNumAtoms())
]

matches = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[O;H1]"))
if len(matches) != 1:
    raise RuntimeError(f"Expected one carboxylic acid group, found {len(matches)}")

acid_oxygen = mol.GetAtomWithIdx(matches[0][2])
acid_oxygen.SetNoImplicit(True)
acid_oxygen.SetNumExplicitHs(0)
acid_oxygen.SetFormalCharge(-1)

Chem.SanitizeMol(mol)
mol_h = Chem.AddHs(mol, addCoords=True)

heavy_after = Chem.RemoveHs(mol_h)
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
mol_h.SetProp("_Name", "DKA_A181_crystal_carboxylate")
Chem.MolToMolFile(mol_h, str(output))

print(f"Wrote: {output}")
print(f"Heavy atoms: {mol_h.GetNumHeavyAtoms()}")
print(f"Hydrogens added: {mol_h.GetNumAtoms() - mol_h.GetNumHeavyAtoms()}")
print(f"Formal charge: {Chem.GetFormalCharge(mol_h)}")
print(f"Heavy-atom coordinate RMSD: {rmsd:.6f} Angstrom")
