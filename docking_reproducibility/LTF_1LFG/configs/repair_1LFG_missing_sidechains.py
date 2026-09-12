from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile

source = Path(
    "LTF_1LFG/"
    "00_raw_inputs/1LFG_raw.pdb"
)
output = Path(
    "LTF_1LFG/"
    "01_receptor_clean/LTF_1LFG_chainA_protein_repaired.pdb"
)

fixer = PDBFixer(filename=str(source))
fixer.removeHeterogens(False)
fixer.findMissingResidues()

if fixer.missingResidues:
    raise RuntimeError(
        f"Unexpected missing residue blocks: {fixer.missingResidues}"
    )

fixer.findMissingAtoms()

expected = {
    ("A", "86", "ARG"): {"CG", "CD", "NE", "CZ", "NH1", "NH2"},
    ("A", "637", "GLU"): {"CG", "CD", "OE1", "OE2"},
}

observed = {
    (res.chain.id, res.id, res.name): {atom.name for atom in atoms}
    for res, atoms in fixer.missingAtoms.items()
}

if observed != expected:
    raise RuntimeError(
        f"Unexpected missing atoms.\nExpected: {expected}\nObserved: {observed}"
    )

fixer.addMissingAtoms()

with output.open("w") as handle:
    PDBFile.writeFile(
        fixer.topology,
        fixer.positions,
        handle,
        keepIds=True,
    )

print(f"Wrote: {output}")
print("Repaired residues:", len(expected))
print("Added heavy atoms:", sum(len(atoms) for atoms in expected.values()))
