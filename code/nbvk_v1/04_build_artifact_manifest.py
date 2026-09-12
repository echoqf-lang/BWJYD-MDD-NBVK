from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "NBVK_FILE_MANIFEST.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


included = []
for folder in (
    "code/geo",
    "code/nbvk_v1",
    "code/nbvk_v2",
    "data/geo",
    "data/network_inputs",
    "results/geo",
    "results/nbvk_v1",
    "results/nbvk_v2",
    "docs/nbvk",
    "environment/geo",
    "environment/nbvk",
):
    for path in sorted((ROOT / folder).rglob("*")):
        if not path.is_file() or path == OUTPUT:
            continue
        included.append(
            {
                "relative_path": str(path.relative_to(ROOT)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )

manifest = {
    "created_at": datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(),
    "platform": platform.platform(),
    "artifact_count": len(included),
    "artifacts": included,
}
OUTPUT.write_text(
    json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
    encoding="utf-8",
)
print(f"Wrote {OUTPUT} with {len(included)} artifacts")
