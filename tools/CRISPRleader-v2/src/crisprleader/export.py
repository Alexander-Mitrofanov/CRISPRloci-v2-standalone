"""Atomic data-only export with explicit circular segments and stable identifiers."""

from __future__ import annotations

import csv
import ctypes
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any
from urllib.parse import quote

from .errors import InputError


def _identifier(*parts: str) -> str:
    return "leader-" + hashlib.sha256("\0".join(parts).encode()).hexdigest()[:20]


def _cell(value: Any) -> Any:
    """Prevent formula evaluation when TSV is opened in spreadsheet software."""
    if value is None:
        return ""
    if isinstance(value, str) and value[:1] in {"=", "+", "-", "@"}:
        return "'" + value
    return value


def _features(report: dict[str, Any]) -> list[tuple[str, dict[str, Any], str]]:
    features = []
    for candidate in report["candidates"]:
        features.append((_identifier(candidate["source_id"], candidate["id"], "flank"),
                         candidate, "extracted_flank"))
    for row in report["predictions"]:
        if row["core"] is not None:
            features.append((_identifier(row["source_id"], row["array_id"], "core"),
                             row["core"], "predicted_conserved_core"))
        for index, motif in enumerate(row["motifs"]):
            features.append((_identifier(row["source_id"], row["array_id"], "motif", str(index)),
                             motif["feature"], "motif_occurrence"))
    return features


def _publish_no_replace(source: Path, target: Path) -> None:
    """Linux atomic no-clobber rename; fail closed when unsupported."""
    library = ctypes.CDLL(None, use_errno=True)
    rename = getattr(library, "renameat2", None)
    if rename is None:
        raise InputError("Atomic directory publication requires Linux renameat2 support")
    rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    rename.restype = ctypes.c_int
    if rename(-100, os.fsencode(source), -100, os.fsencode(target), 1) != 0:
        code = ctypes.get_errno()
        raise OSError(code, os.strerror(code), str(target))


def write_json_report(document: dict[str, Any], output: str | Path) -> Path:
    """Atomically publish one data-only report without clobbering any destination."""
    target = Path(output).absolute()
    if target.exists() or target.is_symlink():
        raise InputError("Output path already exists; choose a new report filename")
    if not target.parent.is_dir():
        raise InputError("Output parent directory does not exist")
    descriptor, name = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(document, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
        _publish_no_replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    return target


def write_bundle(report: dict[str, Any], output: str | Path) -> Path:
    """Export a validated report; atomically publish without replacing any target."""
    target = Path(output).absolute()
    if target.exists() or target.is_symlink():
        raise InputError("Output path already exists; choose a new directory")
    if not target.parent.is_dir():
        raise InputError("Output parent directory does not exist")
    temp = Path(tempfile.mkdtemp(prefix=f".{target.name}.", dir=target.parent))
    try:
        (temp / "result.json").write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        with (temp / "leaders.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["source_id", "array_id", "status", "candidate_id", "strand",
                             "segments_0based_halfopen", "core_sequence", "reasons"])
            for row in report["predictions"]:
                core = row.get("core") or {}
                values = [row["source_id"], row["array_id"], row["status"],
                          row["selected_candidate_id"], core.get("strand"),
                          json.dumps(core.get("segments")), core.get("sequence"),
                          ";".join(row["reasons"])]
                writer.writerow([_cell(value) for value in values])
        with (temp / "motifs.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["source_id", "array_id", "candidate_id", "model_id", "domain",
                             "local_start", "local_end", "strand", "score", "pvalue", "qvalue",
                             "source_segments_0based_halfopen"])
            for row in report["predictions"]:
                for hit in row["motifs"]:
                    values = [row["source_id"], row["array_id"], hit["candidate_id"], hit["model_id"],
                              hit["domain"], hit["start"], hit["end"], hit["feature"]["strand"],
                              hit["score"], hit["pvalue"], hit["qvalue"],
                              json.dumps(hit["feature"]["segments"])]
                    writer.writerow([_cell(value) for value in values])
        features = _features(report)
        with (temp / "sequences.fasta").open("w") as fasta, \
             (temp / "features.gff3").open("w") as gff, \
             (temp / "features.bed").open("w") as bed:
            gff.write("##gff-version 3\n")
            gff.write("# Linked segment IDs share feature_id; sequences.fasta stores oriented full arcs.\n")
            for source in report["sources"]:
                gff.write(f"##sequence-region {quote(source['id'], safe='')} 1 {source['length']}\n")
            for identifier, feature, kind in features:
                if not feature["sequence"]:
                    continue
                fasta.write(f">{identifier} kind={kind} source={quote(feature['source_id'], safe='')}\n")
                for offset in range(0, len(feature["sequence"]), 80):
                    fasta.write(feature["sequence"][offset:offset + 80] + "\n")
                segments = feature["segments"]
                for number, segment in enumerate(segments, 1):
                    sid = f"{identifier}.part{number}"
                    attributes = (f"ID={sid};feature_id={identifier};part={number}/{len(segments)};"
                                  f"evidence_kind={kind};array_id={quote(feature['array_id'], safe='')}")
                    feature_type = "sequence_motif" if kind == "motif_occurrence" else "sequence_feature"
                    gff.write("\t".join([quote(feature["source_id"], safe=""), "CRisprleader-v2",
                                         feature_type, str(segment["start"] + 1), str(segment["end"]),
                                         ".", feature["strand"], ".", attributes]) + "\n")
                    bed.write("\t".join([feature["source_id"], str(segment["start"]),
                                         str(segment["end"]), sid, "0", feature["strand"]]) + "\n")
        (temp / "README.txt").write_text(
            "CRisprleader v2 results\n"
            f"Analysis mode: {report['analysis_mode']}\n"
            "result.json is the complete scientific evidence record.\n"
            "Native JSON/BED: zero-based half-open. GFF3: one-based inclusive.\n"
            "Origin-crossing features use linked parts, never min-max intervals.\n"
            "FASTA identifiers join GFF3 feature_id; sequences follow feature strand.\n"
            "Extracted flanks are context windows, not inferred biological leader lengths.\n"
            "Predicted conserved cores/motifs do not establish activity or full leader boundaries.\n"
            "TSV cells beginning with spreadsheet formula markers are escaped with an apostrophe.\n"
        )
        manifest = {path.name: {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                               "size_bytes": path.stat().st_size}
                    for path in sorted(temp.iterdir()) if path.is_file()}
        (temp / "manifest.json").write_text(json.dumps({"schema_version": "1.0.0", "files": manifest},
                                                       indent=2, sort_keys=True) + "\n")
        if target.exists() or target.is_symlink():
            raise InputError("Output path appeared during analysis; refusing to overwrite")
        _publish_no_replace(temp, target)
    except BaseException:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    return target
