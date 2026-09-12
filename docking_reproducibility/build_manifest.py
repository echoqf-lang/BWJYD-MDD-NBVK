#!/usr/bin/env python3
"""Regenerate the Supplementary Data S1 file manifest and SHA-256 checksums."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "FILE_MANIFEST.tsv"
CHECKSUMS = ROOT / "SHA256SUMS.txt"
EXCLUDED = {MANIFEST.name, CHECKSUMS.name}


def archived_files() -> list[Path]:
    return sorted(
        path
        for path in ROOT.rglob("*")
        if path.is_file() and path.name not in EXCLUDED
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    files = archived_files()
    manifest_lines = ["relative_path\tsize_bytes"]
    checksum_lines = []
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        manifest_lines.append(f"{relative}\t{path.stat().st_size}")
        checksum_lines.append(f"{sha256(path)}  {relative}")

    MANIFEST.write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    CHECKSUMS.write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
