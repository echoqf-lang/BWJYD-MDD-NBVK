import sys
from pathlib import Path

import numpy as np
from rdkit import Chem


if len(sys.argv) != 3:
    raise SystemExit("Usage: python add_explicit_hydrogens.py INPUT.sdf OUTPUT.sdf")

input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])

mol = Chem.MolFromMolFile(str(input_path), removeHs=False, sanitize=True)
if mol is None:
    raise RuntimeError(f"Could not read: {input_path}")

heavy_count = mol.GetNumAtoms()
heavy_before = np.array(mol.GetConformer().GetPositions())

mol_h = Chem.AddHs(mol, addCoords=True)

heavy_after = np.array(mol_h.GetConformer().GetPositions())[:heavy_count]
heavy_rmsd = np.sqrt(
    np.mean(np.sum((heavy_before - heavy_after) ** 2, axis=1))
)

mol_h.SetProp("hydrogen_method", "RDKit Chem.AddHs(addCoords=True)")
mol_h.SetProp("source_file", input_path.name)

writer = Chem.SDWriter(str(output_path))
writer.write(mol_h)
writer.close()

print(f"Wrote: {output_path}")
print(f"Heavy atoms: {heavy_count}")
print(f"Hydrogens added: {mol_h.GetNumAtoms() - heavy_count}")
print(f"Total atoms: {mol_h.GetNumAtoms()}")
print(f"Formal charge: {Chem.GetFormalCharge(mol_h)}")
print(f"Heavy-atom coordinate RMSD: {heavy_rmsd:.6f} Angstrom")

