import sys
import math
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

source = Path(sys.argv[1])
output = Path(sys.argv[2])
generate_3d = len(sys.argv) > 3 and sys.argv[3] == "--generate-3d"

mol = Chem.MolFromMolFile(str(source), removeHs=False)
if mol is None:
    raise RuntimeError(f"Could not read {source}")

mol = Chem.RemoveHs(mol)

if generate_3d:
    mol_h = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 20260727
    params.useRandomCoords = True

    result = AllChem.EmbedMolecule(mol_h, params)
    if result != 0:
        raise RuntimeError(f"ETKDG embedding failed with code {result}")

    optimization = AllChem.MMFFOptimizeMolecule(
        mol_h,
        maxIters=2000,
    )
    mol = Chem.RemoveHs(mol_h)
else:
    optimization = None

conf_before = mol.GetConformer()
coords_before = [
    (
        conf_before.GetAtomPosition(i).x,
        conf_before.GetAtomPosition(i).y,
        conf_before.GetAtomPosition(i).z,
    )
    for i in range(mol.GetNumAtoms())
]

matches = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[O;H1]"))
if len(matches) != 1:
    raise RuntimeError(
        f"Expected one carboxylic acid group, found {len(matches)}"
    )

acid_oxygen = mol.GetAtomWithIdx(matches[0][2])
acid_oxygen.SetNumExplicitHs(0)
acid_oxygen.SetNoImplicit(True)
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

rmsd = math.sqrt(sum_sq / len(coords_before))

mol_h.SetProp("_Name", output.stem)
Chem.MolToMolFile(mol_h, str(output))

print(f"Wrote: {output}")
print(f"Heavy atoms: {mol_h.GetNumHeavyAtoms()}")
print(f"Hydrogens: {mol_h.GetNumAtoms() - mol_h.GetNumHeavyAtoms()}")
print(f"Formal charge: {Chem.GetFormalCharge(mol_h)}")
print(f"Post-preparation heavy-atom RMSD: {rmsd:.6f} Angstrom")
print(f"Generated 3D: {generate_3d}")
print(f"ETKDG random seed: {20260727 if generate_3d else 'not used'}")
print(f"MMFF optimization status: {optimization}")
print(f"Isomeric SMILES: {Chem.MolToSmiles(mol_h, isomericSmiles=True)}")
