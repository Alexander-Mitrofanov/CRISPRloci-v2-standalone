"""Explicit occurrence grouping and offline association transactions."""

from __future__ import annotations

import csv
import shutil
from contextlib import closing
from pathlib import Path

from . import __version__
from .engine import run
from .io import InputError, digest, fasta_bytes, groups_value, json_value, read, transaction, write_json
from .reference import PUBLIC_FIELDS, load
from .report import parse


def associate(
    input_path: Path,
    group_path: Path,
    reference: Path,
    output: Path,
    fdr: float = 0.05,
    threads: int = 1,
    memory_mb: int = 2048,
) -> dict:
    if (
        type(fdr) not in (int, float)
        or not 0.001 <= fdr <= 0.5
        or type(threads) is not int
        or not 1 <= threads <= 8
        or type(memory_mb) is not int
        or not 1024 <= memory_mb <= 8192
    ):
        raise InputError("Invalid FDR or resource limits")
    input_raw, group_raw = read(input_path), read(group_path)
    records = fasta_bytes(input_raw)
    assignments = groups_value(json_value(group_raw), records)
    manifest, targets = load(reference)
    by_id = {row["id"]: row for row in records}
    query: dict[str, dict] = {}
    native_groups = {}
    group_rows = []
    with closing(targets), transaction(output) as temp:
        work = temp / "work"
        work.mkdir()
        files = []
        for index, group in enumerate(assignments):
            native = f"group{index:04d}.fasta"
            native_groups[native] = group["id"]
            unique: dict[str, list[str]] = {}
            unsupported = []
            for ident in group["members"]:
                row = by_id[ident]
                if len(row["sequence"]) < 27 or set(row["sequence"]) - set("ACGT"):
                    unsupported.append(ident)
                    continue
                unique.setdefault(row["sequence"], []).append(ident)
            group_rows.append(
                {
                    **group,
                    "supported_occurrence_count": len(group["members"]) - len(unsupported),
                    "unique_sequence_count": len(unique),
                    "unsupported_occurrence_ids": unsupported,
                }
            )
            if not unique:
                continue
            path = work / native
            chunks = []
            for seq, ids in unique.items():
                native_id = f"q{len(query):06d}"
                query[native_id] = {"sequence": seq, "occurrence_ids": ids, "group_id": group["id"]}
                chunks.append(f">{native_id}\n{seq}\n")
            path.write_text("".join(chunks))
            files.append(str(path))
        rows = []
        calibration_receipt = None
        calibration = "not_run"
        warnings = []
        if files:
            run(
                [
                    "createsetdb",
                    *files,
                    str(work / "query"),
                    str(work / "query_tmp"),
                    "--extractorf-spacer",
                    "1",
                    "--threads",
                    str(threads),
                ],
                temp / "engine.log",
            )
            run(
                [
                    "predictmatch",
                    str(work / "query"),
                    str(reference / "db/target"),
                    str(reference / "db/control"),
                    str(temp / "native.tsv"),
                    str(work / "predict_tmp"),
                    "--fmt",
                    "2",
                    "--report-fdr",
                    "1",
                    "--report-pam",
                    "0",
                    "--fdr",
                    str(fdr),
                    "--threads",
                    str(threads),
                    "--split-memory-limit",
                    f"{memory_mb}M",
                ],
                temp / "engine.log",
            )
            receipts = [
                path
                for path in (work / "predict_tmp").glob("*/match.calibration.json")
                if not path.parent.is_symlink()
            ]
            if len(receipts) != 1:
                raise InputError("Native calibration receipt missing or ambiguous")
            calibration_receipt = json_value(read(receipts[0], 4096))
            if (
                set(calibration_receipt)
                != {"schema_version", "status", "target_score_count", "control_score_count"}
                or calibration_receipt["schema_version"] != "spacepharer.calibration/1.0.0"
                or calibration_receipt["status"] not in {"estimated", "unavailable"}
            ):
                raise InputError("Invalid native calibration receipt")
            for count in ("target_score_count", "control_score_count"):
                if (
                    type(calibration_receipt[count]) is not int
                    or not 0 <= calibration_receipt[count] <= 10_000_000
                ):
                    raise InputError("Invalid native calibration count")
            calibrated = calibration_receipt["status"] == "estimated"
            if (
                calibrated
                and min(calibration_receipt["target_score_count"], calibration_receipt["control_score_count"])
                < 2
            ):
                raise InputError("Invalid native calibration support")
            calibration = calibration_receipt["status"]
            if not calibrated:
                warnings.append(
                    "Control evidence cannot estimate FDR; candidates are not accepted associations."
                )
            rows = parse(
                read(temp / "native.tsv", 16 * 1024 * 1024).decode(),
                query,
                targets,
                native_groups,
                calibrated,
                fdr,
            )
        report: dict = {
            "schema_version": "spacepharer.association/2.0.0",
            "tool_version": __version__,
            "input_sha256": digest(input_raw),
            "groups_sha256": digest(group_raw),
            "requested_fdr": fdr,
            "calibration_status": calibration,
            "calibration": calibration_receipt,
            "warnings": warnings,
            "reference": {key: manifest[key] for key in (*PUBLIC_FIELDS, "manifest_sha256")},
            "records": records,
            "groups": group_rows,
            "associations": rows,
            "interpretation": "Spacer similarity supports a group–target association hypothesis, not active infection, cleavage or a PAM prediction.",
        }
        write_json(temp / "report.json", report)
        with (temp / "associations.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(
                ["group_id", "target_id", "status", "combined_score", "native_fdr", "native_hit_count"]
            )
            for row in rows:
                writer.writerow(
                    [
                        row[k]
                        for k in [
                            "group_id",
                            "target_id",
                            "status",
                            "combined_score",
                            "native_fdr",
                            "native_hit_count",
                        ]
                    ]
                )
        shutil.rmtree(work)
    return report
