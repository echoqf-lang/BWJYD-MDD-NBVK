from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile

source = Path(
    "LCN2_1QQS/"
    "01_receptor_clean/LCN2_1QQS_chainA_clean.pdb"
)
output = Path(
    "LCN2_1QQS/"
    "01_receptor_clean/LCN2_1QQS_chainA_sidechains_repaired.pdb"
)

fixer = PDBFixer(filename=str(source))
fixer.findMissingResidues()

if fixer.missingResidues:
    raise RuntimeError(
        f"Unexpected missing residue blocks: {fixer.missingResidues}"
    )

fixer.findMissingAtoms()

expected = {
    ("A", "4", "THR"): {"CG2", "OG1"},
    ("A", "72", "ARG"): {"CG", "CD", "NE", "CZ", "NH1", "NH2"},
    ("A", "73", "LYS"): {"CG", "CD", "CE", "NZ"},
    ("A", "74", "LYS"): {"CG", "CD", "CE", "NZ"},
    ("A", "98", "LYS"): {"CG", "CD", "CE", "NZ"},
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
