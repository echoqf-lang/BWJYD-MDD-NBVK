import math
from pathlib import Path

root = Path(
    "pocket_detection/"
    "LTF_1LFG_chainA_2Fe_2CO3_repaired_out/pockets"
)

iron = {
    "FE701": (-27.958, 23.026, 2.254),
    "FE702": (13.172, 26.254, 16.064),
}

for number in range(1, 10):
    vertex_path = root / f"pocket{number}_vert.pqr"
    atom_path = root / f"pocket{number}_atm.pdb"

    vertices = []
    with vertex_path.open() as handle:
        for line in handle:
            if line.startswith("ATOM"):
                vertices.append(
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                )

    residues = set()
    with atom_path.open() as handle:
        for line in handle:
            if line.startswith("ATOM"):
                chain = line[21:22]
                resid = int(line[22:26])
                resname = line[17:20].strip()
                residues.add((chain, resid, resname))

    xs = [xyz[0] for xyz in vertices]
    ys = [xyz[1] for xyz in vertices]
    zs = [xyz[2] for xyz in vertices]

    center = (
        (min(xs) + max(xs)) / 2,
        (min(ys) + max(ys)) / 2,
        (min(zs) + max(zs)) / 2,
    )
    spans = (
        max(xs) - min(xs),
        max(ys) - min(ys),
        max(zs) - min(zs),
    )

    n_count = sum(1 for _, resid, _ in residues if resid <= 333)
    c_count = sum(1 for _, resid, _ in residues if resid >= 345)
    lobe = "N-lobe" if n_count > c_count else "C-lobe"

    fe_distances = {
        name: math.dist(center, xyz)
        for name, xyz in iron.items()
    }

    residue_numbers = sorted(resid for _, resid, _ in residues)

    print(f"Pocket {number}: {lobe}")
    print(
        "  center:",
        f"{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}",
    )
    print(
        "  vertex spans:",
        f"{spans[0]:.3f} {spans[1]:.3f} {spans[2]:.3f}",
    )
    print(
        "  residues:",
        len(residues),
        f"range={min(residue_numbers)}-{max(residue_numbers)}",
        f"N_count={n_count}",
        f"C_count={c_count}",
    )
    print(
        "  center-to-Fe:",
        " ".join(
            f"{name}={distance:.3f}"
            for name, distance in fe_distances.items()
        ),
    )
