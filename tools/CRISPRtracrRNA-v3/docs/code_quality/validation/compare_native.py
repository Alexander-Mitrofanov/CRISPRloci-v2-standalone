"""Compare baseline and final native test artifacts without ignoring prediction fields."""

import csv
import hashlib
import json
import sys
from pathlib import Path

baseline, final, destination = map(Path, sys.argv[1:])
results = []
for entry in json.loads((baseline.parent / "native-report-index.json").read_text()):
    relative = Path(entry["path"])
    before = json.loads((baseline / relative).read_text())
    after = json.loads((final / relative).read_text())
    checks = {"summary": before["summary"] == after["summary"]}
    for index, (old, new) in enumerate(zip(before["records"], after["records"], strict=True)):
        for key in (
            "predictions", "raw_model_hits", "raw_genomic_model_hits", "arrays",
            "short_rescue_proposals",
        ):
            checks[f"record_{index}_{key}"] = old.get(key) == new.get(key)
    for name in ("predictions.fasta", "predictions.gff3"):
        checks[name] = (
            (baseline / relative.parent / name).read_bytes()
            == (final / relative.parent / name).read_bytes()
        )
    def csv_rows(root):
        with (root / relative.parent / "summary.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            # Temporary fixture files live beneath different pytest basetemps.
            row["source_file"] = row["source_file"].replace(str(root), "<TEST_ROOT>")
        return rows
    checks["summary_csv_normalized_test_root"] = csv_rows(baseline) == csv_rows(final)
    manifest = json.loads((final / relative.parent / "manifest.json").read_text())
    checks["manifest_hashes"] = all(
        hashlib.sha256((final / relative.parent / name).read_bytes()).hexdigest() == digest
        for name, digest in manifest["files"].items()
    )
    checks["record_json_consistency"] = all(
        json.loads((final / relative.parent / "work" / row["key"] / "record.json").read_text())
        == row for row in after["records"]
    )
    results.append({"path": str(relative), "checks": checks, "summary": after["summary"]})
result = {"passed": all(all(row["checks"].values()) for row in results), "runs": results}
destination.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
sys.exit(0 if result["passed"] else 1)
