"""Versioned file contracts for unmodified CRISPRidentify v2 and CasAndra."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, cast

from .config import Settings
from .errors import ContractError, TracrError
from .models import Array, Interval, Sequence
from .process import Runner
from .sequence import normalize


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _decode_json_object(text: str) -> dict[str, Any]:
    """Decode an upstream record without ambiguous fields or nonfinite numbers."""

    def unique_fields(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key {key!r}")
            result[key] = value
        return result

    def invalid_constant(value: str) -> None:
        raise ValueError(f"Non-finite JSON number {value}")

    def finite_float(value: str) -> float:
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"Non-finite JSON number {value}")
        return number

    data = json.loads(
        text,
        object_pairs_hook=unique_fields,
        parse_constant=invalid_constant,
        parse_float=finite_float,
    )
    if not isinstance(data, dict):
        raise ValueError("Expected JSON object")
    return cast(dict[str, Any], data)


def load_json(path: Path) -> dict[str, Any]:
    try:
        return _decode_json_object(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ContractError(f"Cannot read structured report {path}: {exc}") from exc


def parse_identify(path: Path, source: Sequence) -> tuple[list[Array], dict[str, Any]]:
    report = load_json(path)
    try:
        schema = report["schema"]
        if schema["name"] != "CRISPRidentify-report" or schema["version"].split(".")[0] != "1":
            raise ValueError(f"Unsupported CRISPRidentify schema {schema}")
        coords = report["coordinate_system"]
        if (
            coords["indexing"] != 1
            or coords["end_inclusive"] is not True
            or coords["reference"] != "input_forward_strand"
        ):
            raise ValueError("Unsupported CRISPRidentify coordinates")
        reported_source = report["source"]
        if (
            reported_source["sha256"] != source.sha256
            or reported_source["id"] != source.id
            or reported_source["length"] != len(source.bases)
        ):
            raise ValueError("CRISPRidentify source identity/hash does not match the input")
        if report["validation"]["status"] != "valid":
            raise ValueError("CRISPRidentify reports invalid source reconstruction")
        arrays = []
        seen = set()
        for row in report["arrays"]:
            interval = Interval.from_public(
                row["source_interval"]["start"], row["source_interval"]["end"]
            )
            if interval.end > len(source.bases):
                raise ValueError("Array exceeds contig boundary")
            if row["validation"]["status"] != "valid":
                raise ValueError("Invalid array reconstruction")
            if (
                row.get("source_sequence", source.bases[interval.start : interval.end])
                != source.bases[interval.start : interval.end]
            ):
                raise ValueError("Array sequence does not reconstruct from source")
            strand = row["orientation"].get("strand")
            if strand not in {"+", "-", None}:
                raise ValueError(f"Invalid array orientation {strand}")
            score = row["certainty_score"]
            if type(score) not in {int, float} or not 0 <= score <= 1 or not math.isfinite(score):
                raise ValueError("Invalid array score")
            identifier = row["id"]
            if not isinstance(identifier, str) or not identifier:
                raise ValueError("Array id must be a nonempty string")
            try:
                repeat = normalize(row["consensus_repeat"])
            except TracrError as exc:
                raise ContractError(
                    f"Incompatible CRISPRidentify report {path}: invalid consensus repeat: {exc}"
                ) from exc
            if row["category"] not in {"Bona-fide", "Possible"}:
                continue
            if identifier in seen:
                raise ValueError("Duplicate array id")
            seen.add(identifier)
            arrays.append(
                Array(
                    identifier,
                    interval,
                    repeat,
                    strand,
                    row["category"],
                    float(score),
                )
            )
        return arrays, report
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        raise ContractError(f"Incompatible CRISPRidentify report {path}: {exc}") from exc


def run_identify(
    fasta: Path, source: Sequence, work: Path, settings: Settings, runner: Runner
) -> tuple[list[Array], dict[str, Any]]:
    output = work / "identify"
    runner.run(
        [
            *settings.identify_command,
            "--file",
            str(fasta),
            "--result_folder",
            str(output),
            "--preset",
            settings.identify_preset,
            "--cas",
            "False",
            "--cpu",
            str(settings.threads),
        ]
    )
    reports = sorted(output.rglob("report.json"))
    if len(reports) != 1:
        raise ContractError(f"Expected one CRISPRidentify report, got {len(reports)} in {output}")
    return parse_identify(reports[0], source)


def parse_casandra(directory: Path, fasta: Path, source: Sequence) -> dict[str, Any]:
    run = load_json(directory / "run.json")
    manifest = load_json(directory / "manifest.json")
    try:
        if (
            run["program"] != "CasAndra"
            or run["schema_version"] != 5
            or manifest["schema_version"] != 5
        ):
            raise ValueError("Unsupported CasAndra genome report schema (expected 5)")
        input_digest = sha256_file(fasta)
        if not any(
            row["sha256"] == input_digest and row["kind"] == "genome_fasta" for row in run["inputs"]
        ):
            raise ValueError("CasAndra input digest does not match")
        for count in ("genes", "cas_cassettes"):
            if type(run[count]) is not int or run[count] < 0:
                raise ValueError(f"CasAndra {count} count must be a nonnegative integer")
        result: dict[str, Any] = {"run": run, "status": "computed"}
        for filename, key in (("proteins.jsonl", "genes"), ("cassettes.jsonl", "cassettes")):
            path = directory / filename
            if sha256_file(path) != manifest["files"][filename]["sha256"]:
                raise ValueError(f"CasAndra checksum mismatch: {filename}")
            rows = []
            seen = set()
            with path.open(encoding="utf-8") as handle:
                for number, line in enumerate(handle, 1):
                    if not line.strip():
                        continue
                    try:
                        row = _decode_json_object(line)
                        interval = Interval.from_public(
                            row["start_1based"], row["end_1based_inclusive"]
                        )
                        if interval.end > len(source.bases) or row["contig_id"] != source.id:
                            raise ValueError(
                                "CasAndra coordinates or contig do not match the source"
                            )
                        if key == "genes":
                            if type(row["prediction"]["is_cas"]) is not bool:
                                raise ValueError("CasAndra is_cas must be a boolean")
                        else:
                            if not isinstance(row["final_classification"], dict):
                                raise ValueError("Missing CasAndra cassette classification")
                            identifier = row["cassette_id"]
                            if not isinstance(identifier, str) or not identifier:
                                raise ValueError("CasAndra cassette id must be a nonempty string")
                            if identifier in seen:
                                raise ValueError("Duplicate CasAndra cassette id")
                            seen.add(identifier)
                    except (KeyError, TypeError, ValueError) as exc:
                        raise ContractError(
                            f"Incompatible CasAndra record {path}:{number}: {exc}"
                        ) from exc
                    rows.append(row)
            result[key] = rows
        if len(result["genes"]) != run["genes"] or len(result["cassettes"]) != run["cas_cassettes"]:
            raise ValueError("CasAndra record count mismatch")
        return result
    except (OSError, KeyError, TypeError, ValueError) as exc:
        raise ContractError(f"Incompatible CasAndra output in {directory}: {exc}") from exc


def run_casandra(
    fasta: Path, source: Sequence, work: Path, settings: Settings, runner: Runner
) -> dict[str, Any]:
    output = work / "casandra"
    runner.run(
        [
            *settings.casandra_command,
            "predict-genome",
            "--genome",
            str(fasta),
            "--output",
            str(output),
            "--threads",
            str(settings.threads),
            "--gene-mode",
            settings.gene_mode,
        ]
    )
    return parse_casandra(output, fasta, source)


def nearest_cas(interval: Interval, cas: dict[str, Any], model_type: str) -> dict[str, Any] | None:
    rows = [
        row for row in cas["cassettes"] if row["final_classification"].get("type") == model_type
    ]
    if not rows:
        return None
    row = min(
        rows,
        key=lambda r: (
            interval.distance(Interval.from_public(r["start_1based"], r["end_1based_inclusive"])),
            r["cassette_id"],
        ),
    )
    location = Interval.from_public(row["start_1based"], row["end_1based_inclusive"])
    return {
        "id": row["cassette_id"],
        "interval": location.public(),
        "distance_bases": interval.distance(location),
        "classification": row["final_classification"],
    }
