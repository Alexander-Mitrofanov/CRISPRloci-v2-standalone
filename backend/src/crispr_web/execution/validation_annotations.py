"""Independent validation of CasAndra and CRISPRtracrRNA against immutable input.

This module deliberately does not import Integration or scientific tool packages.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import stat

from .validation_identify import _validate_identify_reports

_LIMIT = 20_000_000
_KEY = re.compile(r"f[0-9]{6}\.s[0-9]{6}")
_DNA = frozenset("ACGTRYSWKMBDHVN")
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def _payload(root, path, checkpoint=None):
    relative = path.relative_to(root)
    current = root
    for part in relative.parts:
        current /= part
        if stat.S_ISLNK(current.lstat().st_mode):
            raise ValueError("Symlink in scientific output")
    metadata = path.stat()
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_size > _LIMIT:
        raise ValueError("Invalid scientific output size/type")
    data = bytearray()
    with path.open("rb") as stream:
        while block := stream.read(65536):
            if checkpoint:
                checkpoint()
            data.extend(block)
            if len(data) > _LIMIT:
                raise ValueError("Scientific output grew past the limit")
    return bytes(data)


def _strict_json(payload):
    def object_pairs(pairs):
        value = {}
        for key, item in pairs:
            if key in value:
                raise ValueError("Duplicate JSON key")
            value[key] = item
        return value

    def number(value):
        parsed = float(value)
        if not math.isfinite(parsed):
            raise ValueError("Nonfinite JSON number")
        return parsed

    def invalid(value):
        raise ValueError("Invalid JSON constant: " + value)

    return json.loads(
        payload,
        object_pairs_hook=object_pairs,
        parse_float=number,
        parse_constant=invalid,
    )


def _json(root, path, checkpoint=None):
    document = _strict_json(_payload(root, path, checkpoint))
    if not isinstance(document, dict):
        raise ValueError("Scientific report must be an object")
    return document


def _fasta(payload):
    lines = payload.decode("ascii").splitlines()
    if (
        not lines
        or not lines[0].startswith(">")
        or any(line.startswith(">") for line in lines[1:])
    ):
        raise ValueError("Expected one normalized FASTA record")
    source_id = lines[0][1:].split()[0]
    sequence = "".join(lines[1:]).upper()
    if not sequence or set(sequence) - _DNA:
        raise ValueError("Invalid normalized FASTA sequence")
    return source_id, sequence


def _interval(row, length, start_key="start", end_key="end"):
    start, end = row.get(start_key), row.get(end_key)
    if (
        type(start) is not int
        or type(end) is not int
        or not 1 <= start <= end <= length
    ):
        raise ValueError("Scientific coordinates do not match input")
    return start, end


def _count(value):
    if type(value) is not int or value < 0:
        raise ValueError("Invalid scientific count")
    return value


def _metric(value, *, nonnegative=False):
    if value is not None and (
        type(value) not in {int, float}
        or not math.isfinite(value)
        or (nonnegative and value < 0)
    ):
        raise ValueError("Invalid scientific numeric field")


def _cas(root, directory, payload, source_id, sequence, checkpoint):
    run = _json(root, directory / "run.json", checkpoint)
    manifest = _json(root, directory / "manifest.json", checkpoint)
    if (
        run.get("program") != "CasAndra"
        or type(run.get("schema_version")) is not int
        or run["schema_version"] != 5
    ):
        raise ValueError("Unsupported CasAndra report")
    if manifest.get("schema_version") != 5 or not any(
        item.get("kind") == "genome_fasta"
        and item.get("sha256") == hashlib.sha256(payload).hexdigest()
        for item in run["inputs"]
    ):
        raise ValueError("CasAndra input/manifest mismatch")
    counts = {}
    cas_genes = 0
    for filename, key, id_key in [
        ("proteins.jsonl", "genes", "protein_id"),
        ("cassettes.jsonl", "cas_cassettes", "cassette_id"),
    ]:
        raw = _payload(root, directory / filename, checkpoint)
        if hashlib.sha256(raw).hexdigest() != manifest["files"][filename]["sha256"]:
            raise ValueError("CasAndra output checksum mismatch")
        seen = set()
        for line in raw.splitlines():
            if not line.strip():
                continue
            if checkpoint:
                checkpoint()
            row = _strict_json(line)
            if row["contig_id"] != source_id:
                raise ValueError("CasAndra output belongs to a different contig")
            _interval(row, len(sequence), "start_1based", "end_1based_inclusive")
            identity = row[id_key]
            if not isinstance(identity, str) or not identity or identity in seen:
                raise ValueError("Invalid or duplicate CasAndra record id")
            seen.add(identity)
            if key == "genes":
                prediction = row["prediction"]
                if type(prediction.get("is_cas")) is not bool or row["strand"] not in {
                    "+",
                    "-",
                }:
                    raise ValueError("Invalid CasAndra protein call/strand")
                if not isinstance(prediction.get("classification"), dict):
                    raise ValueError("Missing CasAndra protein classification")
                start, end = row["start_1based"], row["end_1based_inclusive"]
                nucleotides = sequence[start - 1 : end]
                if row["strand"] == "-":
                    nucleotides = nucleotides.translate(_COMPLEMENT)[::-1]
                # Pyrodigal reports caller-normalized CDS DNA: each submitted
                # ambiguity symbol becomes N. The run SHA above still binds
                # the exact submitted FASTA; normalization never changes it.
                if (row.get("caller") or {}).get("program") == "Pyrodigal":
                    nucleotides = "".join(base if base in "ACGT" else "N" for base in nucleotides)
                if (
                    row.get("nucleotide_sequence") != nucleotides
                    or row.get("nucleotide_sha256")
                    != hashlib.sha256(nucleotides.encode()).hexdigest()
                ):
                    raise ValueError(
                        "CasAndra gene sequence differs from submitted input"
                    )
                _metric(prediction.get("score_margin"))
                cas_genes += prediction["is_cas"]
            else:
                classification = row["final_classification"]
                if not isinstance(classification, dict):
                    raise ValueError("Missing cassette classification")
                confidence = classification.get("confidence")
                if confidence is not None and (
                    type(confidence) not in {int, float} or not 0 <= confidence <= 1
                ):
                    raise ValueError("Invalid cassette confidence")
                _count(row["cas_gene_count"])
        if len(seen) != _count(run[key]):
            raise ValueError("CasAndra output count mismatch")
        counts[key] = len(seen)
    if cas_genes != _count(run["cas_proteins"]):
        raise ValueError("CasAndra Cas-protein count mismatch")
    return {
        "gene_count": counts["genes"],
        "cas_gene_count": cas_genes,
        "cassette_count": counts["cas_cassettes"],
    }


def _tracr(claimed, root, sources, stages, checkpoint):
    output = root / "tracrrna"
    report = _json(root, output / "report.json", checkpoint)
    manifest = _json(root, output / "manifest.json", checkpoint)
    mode = "complete" if claimed.mode == "loci" else "model"
    if (
        report.get("schema") != {"name": "CRISPRtracrRNA-report", "version": "3.0.0"}
        or report.get("status") != "complete"
    ):
        raise ValueError("Invalid tracrRNA completion report")
    if manifest.get("algorithm") != "sha256":
        raise ValueError("Unsupported tracrRNA manifest")
    config = report["configuration"]
    if (
        config.get("mode") != mode
        or config.get("model_type") != claimed.tracr_model_type
    ):
        raise ValueError("tracrRNA analysis policy mismatch")
    for flag, value in claimed.tool_options.get("crispr_tracrrna", {}).items():
        if config.get(flag[2:].replace("-", "_")) != value:
            raise ValueError("tracrRNA scientific option mismatch")
    if config.get("short_anti_repeat_rescue") is not False:
        raise ValueError("Experimental tracrRNA policy was unexpectedly enabled")
    coordinates = report["coordinate_system"]
    if (
        coordinates.get("indexing") != 1
        or coordinates.get("end_inclusive") is not True
        or coordinates.get("reference") != "input_forward_strand"
    ):
        raise ValueError("Unsupported tracrRNA coordinates")
    seen, keys, prediction_ids = set(), set(), set()
    total = 0
    for record in report["records"]:
        source = record["source"]
        source_id = source["id"]
        key = record["key"]
        if (
            not isinstance(key, str)
            or not _KEY.fullmatch(key)
            or key in keys
            or source_id in seen
            or source_id not in sources
        ):
            raise ValueError("Invalid tracrRNA record identity")
        sequence, original = sources[source_id]
        if (source.get("length"), source.get("sha256"), source.get("file_sha256")) != (
            len(sequence),
            hashlib.sha256(sequence.encode()).hexdigest(),
            hashlib.sha256(original).hexdigest(),
        ):
            raise ValueError("tracrRNA source binding mismatch")
        seen.add(source_id)
        keys.add(key)
        relative = f"work/{key}/record.json"
        raw = _payload(root, output / relative, checkpoint)
        if (
            hashlib.sha256(raw).hexdigest() != manifest["files"][relative]
            or _strict_json(raw) != record
        ):
            raise ValueError("tracrRNA retained record mismatch")
        for prediction in record["predictions"]:
            start, end = _interval(prediction["interval"], len(sequence))
            strand = prediction.get("strand")
            if strand not in {"+", "-"}:
                raise ValueError("Invalid tracrRNA strand")
            expected = sequence[start - 1 : end]
            if strand == "-":
                expected = expected.translate(_COMPLEMENT)[::-1]
            if prediction.get("sequence") != expected:
                raise ValueError("tracrRNA sequence does not reconstruct from source")
            _metric((prediction.get("ranking") or {}).get("score"))
            for name in ("identity", "coverage"):
                value = (prediction.get("anti_repeat") or {}).get(name)
                _metric(value, nonnegative=True)
                if value is not None and value > 1:
                    raise ValueError("Invalid tracrRNA anti-repeat fraction")
            _metric((prediction.get("interaction") or {}).get("energy_kcal_mol"))
            _metric(prediction.get("best_score"))
            _metric(prediction.get("best_evalue"), nonnegative=True)
            _metric((prediction.get("model_hit") or {}).get("evalue"), nonnegative=True)
            identity = prediction["export_id"]
            if (
                not isinstance(identity, str)
                or not re.fullmatch(re.escape(key) + r"\.p[0-9]{6}", identity)
                or identity in prediction_ids
            ):
                raise ValueError("Invalid tracrRNA export identity")
            prediction_ids.add(identity)
            total += 1
        if mode == "complete":
            native = output / "work" / key
            internal_payload = _payload(root, native / "input.fasta", checkpoint)
            if _fasta(internal_payload) != (key, sequence):
                raise ValueError(
                    "tracrRNA internal input differs from submitted sequence"
                )
            bindings = {
                key: (len(sequence), hashlib.sha256(sequence.encode()).hexdigest())
            }
            if not _validate_identify_reports(
                root, bindings, checkpoint, reports_root=native / "identify"
            ):
                raise ValueError("Nested CRISPRidentify output is invalid")
            counts = _cas(
                root, native / "casandra", internal_payload, key, sequence, checkpoint
            )
            if (
                record["cas"].get("status") != "computed"
                or record["cas"].get("cas_gene_count") != counts["cas_gene_count"]
            ):
                raise ValueError(
                    "Nested CasAndra report disagrees with tracrRNA context"
                )
            cassette_rows = [
                _strict_json(line)
                for line in _payload(
                    root, native / "casandra/cassettes.jsonl", checkpoint
                ).splitlines()
                if line.strip()
            ]
            if record["cas"].get("cassettes") != cassette_rows or record["cas"].get(
                "run"
            ) != _json(root, native / "casandra/run.json", checkpoint):
                raise ValueError(
                    "tracrRNA context differs from retained CasAndra evidence"
                )
            for name in ("crispridentify", "casandra"):
                if (
                    stages[name].get("status") != "completed"
                    or stages[name].get("provided_by") != "tracrrna_complete"
                ):
                    raise ValueError("Missing complete-mode provenance")
    expected_counts = {"sequence_count": len(sources), "prediction_count": total}
    if seen != set(sources) or report.get("summary") != expected_counts:
        raise ValueError("tracrRNA output counts do not cover the input")
    stage = stages["tracrrna"]
    if (
        stage.get("status") != "completed"
        or stage.get("record_count") != len(sources)
        or stage.get("prediction_count") != total
    ):
        raise ValueError("tracrRNA integration counts disagree with native output")
    if stage.get("mode") != mode or stage.get("model_type") != claimed.tracr_model_type:
        raise ValueError("tracrRNA integration policy mismatch")


def validate_annotation_pipeline(
    claimed, root, expected, integration, input_root, checkpoint=None
):
    if input_root is None or integration.get("status") != "completed":
        return False
    try:
        selection = {"cas": "casandra", "tracrrna": "tracrrna", "loci": "loci"}[
            claimed.mode
        ]
        if integration.get("workflow") != selection:
            return False
        sources = {}
        for path in sorted(input_root.iterdir()):
            raw = _payload(input_root, path, checkpoint)
            identity, sequence = _fasta(raw)
            if identity in sources or expected.get(identity) != (
                len(sequence),
                hashlib.sha256(sequence.encode()).hexdigest(),
            ):
                return False
            sources[identity] = (sequence, raw)
        if set(sources) != set(expected):
            return False
        stages = integration["stages"]
        if claimed.mode != "cas":
            _tracr(claimed, root, sources, stages, checkpoint)
            return True
        totals = {"gene_count": 0, "cas_gene_count": 0, "cassette_count": 0}
        directories = sorted((root / "casandra").iterdir())
        if len(directories) != len(sources):
            return False
        consumed = set()
        by_hash = {
            hashlib.sha256(raw).hexdigest(): (identity, sequence, raw)
            for identity, (sequence, raw) in sources.items()
        }
        for directory in directories:
            run = _json(root, directory / "run.json", checkpoint)
            hashes = [
                item["sha256"]
                for item in run["inputs"]
                if item.get("kind") == "genome_fasta"
            ]
            if len(hashes) != 1 or hashes[0] not in by_hash or hashes[0] in consumed:
                return False
            consumed.add(hashes[0])
            identity, sequence, raw = by_hash[hashes[0]]
            for key, value in _cas(
                root, directory, raw, identity, sequence, checkpoint
            ).items():
                totals[key] += value
        stage = stages["casandra"]
        return stage.get("status") == "completed" and all(
            stage.get(key) == value for key, value in totals.items()
        )
    except (OSError, ValueError, KeyError, TypeError, IndexError, AttributeError):
        return False
