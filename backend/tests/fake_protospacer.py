"""Subprocess fixture for the real queue/containment/publication lifecycle."""

import hashlib
import json
from pathlib import Path
import sys


def main():
    args = sys.argv[1:]
    output = Path(args[2])
    native = output / "protospacer"
    native.mkdir(parents=True)
    fixture = Path(__file__).parent / "fixtures/protospacer"
    report = json.loads((fixture / "analysis.json").read_text())
    expected = args[args.index("--spacer-database-sha256") + 1]
    report["reference"]["database_sha256"] = expected
    retained = json.loads((fixture / "inputs.json").read_text())
    retained["reference"] = report["reference"]
    (native / "analysis.json").write_text(json.dumps(report))
    (native / "inputs.json").write_text(json.dumps(retained))
    (native / "matches.tsv").write_text("private native table\n")
    files = {}
    for name in ("analysis.json", "inputs.json", "matches.tsv"):
        data = (native / name).read_bytes()
        files[name] = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    (native / "manifest.json").write_text(json.dumps({"schema_version": "1.0.0", "kind": "crisprspacer.database-search", "files": files}))
    document = {"schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"}, "status": "completed", "workflow": "protospacer", "stages": {"crisprspacer": {"status": "completed", "counts": report["counts"], "database_sha256": expected}}}
    (output / "integration_result.json").write_text(json.dumps(document))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
