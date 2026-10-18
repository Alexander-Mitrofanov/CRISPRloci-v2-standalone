"""Run the byte-verified archived Python detector against the pinned native profile.

Invoke with PYTHONPATH pointing at the extracted baseline src. No archived source
is edited; its existing fasta36 setting selects the new dependency binary.
"""

from __future__ import annotations

import json
import platform
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path

import crisprtracrrna
from crisprtracrrna.adapters import sha256_file
from crisprtracrrna.config import load_settings
from crisprtracrrna.pipeline import atomic_json, doctor, run

ROOT = Path(__file__).resolve().parents[1]
ROUND = ROOT / "benchmarks/council_round2"


def main():
    output = ROUND / "archived_detector_pinned_native"
    output.mkdir(exist_ok=False)
    config = ROOT / "local.reproducible.toml"
    settings = load_settings(config, {"threads": 1, "mode": "complete"})
    source_root = Path(crisprtracrrna.__file__).parent
    expected = json.loads((ROUND / "baseline/source_sha256.json").read_text())
    for path in source_root.glob("*.py"):
        relative = "src/crisprtracrrna/" + path.name
        assert sha256_file(path) == expected[relative], relative
    source_data = source_root / "data"
    for relative, digest in json.loads(
        (ROOT / "src/crisprtracrrna/data/manifest.json").read_text()
    ).items():
        assert sha256_file(source_data / relative) == digest
    entries = [
        {"label": e["label"], "fasta": e["fasta"]}
        for e in json.loads((ROOT / "benchmarks/literature/loci.json").read_text())
    ]
    manifest = {
        "created_at": datetime.now(UTC).isoformat(),
        "profile": "archived_default_pinned_native",
        "configuration_file": str(config),
        "configuration_sha256": sha256_file(config),
        "settings": settings.public(),
        "workers": 2,
        "threads_per_locus": 1,
        "protocol_sha256": sha256_file(ROUND / "protocol.json"),
        "protocol_amendment_sha256": sha256_file(ROUND / "protocol_amendment_v2.json"),
        "native_execution_manifest_sha256": sha256_file(ROUND / "native_execution_manifest.json"),
        "baseline_source_archive_sha256": sha256_file(ROUND / "baseline/source.zip"),
        "source_sha256": {
            "src/crisprtracrrna/" + p.name: sha256_file(p) for p in source_root.glob("*.py")
        },
        "source_location": str(source_root),
        "resource_relocation": str(source_data.resolve()),
        "inputs": [{**e, "fasta_sha256": sha256_file(ROOT / e["fasta"])} for e in entries],
        "models": {str(p): sha256_file(p) for p in settings.model_paths()},
        "executables": {
            k: {"path": p, "sha256": sha256_file(Path(p))}
            for k, p in doctor(settings)["executables"].items()
        },
        "upstream_lock_sha256": sha256_file(ROOT / "upstream-lock.json"),
        "python": platform.python_version(),
        "platform": platform.platform(),
    }
    atomic_json(output / "run_manifest.json", manifest)
    started = time.monotonic()

    def execute(entry):
        before = time.monotonic()
        try:
            report = run(ROOT / entry["fasta"], output / entry["label"], settings)
            return {
                "label": entry["label"],
                "status": "complete",
                "elapsed_seconds": time.monotonic() - before,
                "prediction_count": report["summary"]["prediction_count"],
            }
        except Exception as exc:
            return {
                "label": entry["label"],
                "status": "failed",
                "error": str(exc),
                "elapsed_seconds": time.monotonic() - before,
            }

    rows = []
    with ThreadPoolExecutor(2) as pool, (output / "results.jsonl").open("w") as handle:
        for row in pool.map(execute, entries):
            rows.append(row)
            handle.write(json.dumps(row) + "\n")
            handle.flush()
            print(json.dumps(row), flush=True)
    manifest.update(
        {
            "batch_wall_seconds": time.monotonic() - started,
            "completed_at": datetime.now(UTC).isoformat(),
            "results": rows,
        }
    )
    atomic_json(output / "run_manifest.json", manifest)
    assert len(rows) == 79 and all(r["status"] == "complete" for r in rows)


if __name__ == "__main__":
    main()
