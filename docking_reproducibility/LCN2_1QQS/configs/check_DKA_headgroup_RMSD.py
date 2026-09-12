import sys
import math
from rdkit import Chem

reference_path = sys.argv[1]
poses_path = sys.argv[2]

reference = Chem.MolFromMolFile(reference_path, removeHs=True)
poses = Chem.SDMolSupplier(poses_path, removeHs=True)

if reference is None:
    raise RuntimeError("Could not read reference molecule")

query = Chem.MolFromSmarts("C(=O)O")
reference_match = reference.GetSubstructMatch(query)

if len(reference_match) != 3:
    raise RuntimeError("Reference carboxyl group not found")

ref_conf = reference.GetConformer()


def xyz(conf, atom_index):
    point = conf.GetAtomPosition(atom_index)
    return (point.x, point.y, point.z)


def rmsd(first, second):
    return math.sqrt(
        sum(
            sum((a[i] - b[i]) ** 2 for i in range(3))
            for a, b in zip(first, second)
        )
        / len(first)
    )


ref_head = [xyz(ref_conf, i) for i in reference_match]

for mode, pose in enumerate(poses, start=1):
    if pose is None:
        continue

    pose_match = pose.GetSubstructMatch(query)
    if len(pose_match) != 3:
        raise RuntimeError(f"Mode {mode}: carboxyl group not found")

    pose_conf = pose.GetConformer()
    pose_head = [xyz(pose_conf, i) for i in pose_match]

    direct = rmsd(ref_head, pose_head)
    swapped = rmsd(
        ref_head,
        [pose_head[0], pose_head[2], pose_head[1]],
    )

    print(
        f"mode={mode} "
        f"carboxyl_RMSD={min(direct, swapped):.4f} Angstrom"
    )
