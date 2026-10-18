"""Integration's source and numeric gate, independent of the scientific package."""

from __future__ import annotations

import hashlib
import math

from .annotations import read_object

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def digest(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def source_sequences(path, molecule):
    sources = {}
    for file in sorted(path.iterdir()) if path.is_dir() else [path]:
        if file.is_symlink() or file.stat().st_size > 100_000_000:
            raise ValueError("Invalid repeat source file")
        identifier = None
        for line in file.read_text(encoding="ascii").splitlines():
            if line.startswith(">"):
                identifier = line[1:].split()[0]
                if identifier in sources:
                    raise ValueError("Duplicate repeat source")
                sources[identifier] = ""
            elif line.strip():
                if identifier is None:
                    raise ValueError("Invalid FASTA")
                sources[identifier] += line.strip().upper()
    alphabet = set("ACGURYSWKMBDHVN" if molecule == "RNA" else "ACGTRYSWKMBDHVN")
    if not sources or any(not seq or set(seq) - alphabet for seq in sources.values()):
        raise ValueError("Invalid repeat molecule alphabet")
    return sources


def _pairs(rows, sequence):
    length = len(sequence)
    found = {}
    if not isinstance(rows, list) or len(rows) > 250000:
        raise ValueError("Invalid pair table")
    for i, j, value in rows:
        if (
            type(i) is not int
            or type(j) is not int
            or not 0 <= i < j < length
            or type(value) not in {float, int}
            or not math.isfinite(value)
            or not 0 < value <= 1
            or j - i < 4
            or sequence[i] + sequence[j] not in {"AU", "UA", "GC", "CG", "GU", "UG"}
            or (i, j) in found
        ):
            raise ValueError("Invalid pairing probability")
        found[i, j] = value
    return found


def _finite(value):
    return type(value) in {int, float} and math.isfinite(value)


def _unpaired(values, length, pairs=None):
    if not isinstance(values, list) or len(values) != length:
        raise ValueError("Incomplete unpaired probabilities")
    paired = [0.0] * length
    if pairs is not None:
        for (i, j), value in pairs.items():
            paired[i] += value
            paired[j] += value
    for index, value in enumerate(values):
        if not _finite(value) or not 0 <= value <= 1:
            raise ValueError("Invalid unpaired probability")
        if pairs is not None and (paired[index] > 1 + 1e-8 or abs(value - max(0, 1 - paired[index])) > 1e-8):
            raise ValueError("Inconsistent global pair marginals")


def _isolated(result, sequence):
    if set(sequence) - set("ACGU"):
        if result != {"status": "unsupported", "reason": "ambiguous_bases"}:
            raise ValueError("Ambiguous folding must remain unavailable")
        return None, []
    if result.get("status") != "completed" or result.get("method") != "global_partition_function":
        raise ValueError("Missing isolated fold")
    if any(not _finite(result.get(key)) for key in ("mfe_kcal_mol", "ensemble_kcal_mol")):
        raise ValueError("Invalid folding energy")
    pairs = _pairs(result["pairs"], sequence)
    _unpaired(result["unpaired"], len(sequence), pairs)
    structure, stack, motif = result["structure"], [], []
    if not isinstance(structure, str) or len(structure) != len(sequence):
        raise ValueError("Structure length mismatch")
    for index, char in enumerate(structure):
        if char == "(":
            stack.append(index)
        elif char == ")" and stack:
            left = stack.pop()
            if index - left < 4 or sequence[left] + sequence[index] not in {"AU", "UA", "GC", "CG", "GU", "UG"}:
                raise ValueError("Invalid MFE pairing")
            motif.append([left, index])
        elif char != ".":
            raise ValueError("Invalid structure")
    motif.sort()
    if stack or motif != result["motif_pairs"]:
        raise ValueError("Invalid motif")
    return pairs, motif


def validate_repeat_report(report, fasta, reports, molecule, categories, *, tool_options=None):
    options = tool_options or {}
    requested_model = {
        "temperature": options.get("--temperature", 37),
        "window": options.get("--window", 150),
        "span": options.get("--span", 80),
        "no_lonely_pairs": options.get("--no-lonely-pairs", False),
    }
    if (
        report.get("schema") != {"name": "CRISPRrepeat-report", "version": "1.0.0"}
        or report.get("status") != "completed"
        or report.get("molecule") != molecule
        or report.get("input_kind") != ("arrays" if reports else "repeats")
        or report.get("tool") not in ({"name": "CRISPRrepeat", "version": "1.0.0"},
                                      {"name": "CRISPRrepeat", "version": "1.0.1"})
    ):
        raise ValueError("Unsupported repeat output contract")
    if any(report.get("model", {}).get(key) != value for key, value in requested_model.items()):
        raise ValueError("Repeat effective model differs from requested settings")
    sources = source_sequences(fasta, molecule)
    if report["sources"] != [{"id": k, "length": len(v), "sha256": digest(v)} for k, v in sources.items()]:
        raise ValueError("Repeat source identities disagree")
    orientations = ["forward"] if molecule == "RNA" else ["forward", "reverse"]
    if report["orientation_policy"] != ("forward" if molecule == "RNA" else "both"):
        raise ValueError("Repeat orientation policy mismatch")
    expected = []
    array_keys = {}
    if reports:
        native = {}
        for path in sorted(reports.rglob("report.json")):
            doc = read_object(reports, path)
            source_id = doc["source"]["id"]
            for array in doc["arrays"]:
                if array["category"] in categories:
                    native[source_id, array["id"]] = array
        if len(native) != len(report["arrays"]):
            raise ValueError("Repeat array coverage mismatch")
        for array in report["arrays"]:
            raw = native.pop((array["source_id"], array["native_id"]))
            start, end = raw["source_interval"]["start"] - 1, raw["source_interval"]["end"]
            seq = sources[array["source_id"]]
            if (array["start"], array["end"], array["source_sequence"]) != (start, end, seq[start:end]):
                raise ValueError("Repeat array reconstruction failed")
            if len(array["repeats"]) != raw["repeat_count"]:
                raise ValueError("Repeat count mismatch")
            array_keys[array["id"]] = array
            for strand in orientations:
                for repeat, original in zip(array["repeats"], raw["repeats"], strict=True):
                    rs, re = original["source_interval"]["start"] - 1, original["source_interval"]["end"]
                    if (repeat["start"], repeat["end"], repeat["sequence"]) != (rs, re, seq[rs:re]):
                        raise ValueError("Repeat reconstruction failed")
                    expected.append((array["source_id"], array["id"], repeat["index"], strand, rs, re, seq[rs:re]))
    else:
        if report["arrays"] or report["contexts"]:
            raise ValueError("Repeat-only result invents context")
        expected = [
            (key, None, None, strand, None, None, sequence)
            for key, sequence in sources.items()
            for strand in orientations
        ]
    contexts, context_keys = {}, set()
    if len(report["contexts"]) != len(array_keys) * len(orientations):
        raise ValueError("Incomplete repeat context coverage")
    for number, context in enumerate(report["contexts"], 1):
        key = context["array_id"], context["orientation"]
        if key in context_keys or key[1] not in orientations or context["id"] != f"context_{number:06d}":
            raise ValueError("Invalid or duplicate context identity")
        context_keys.add(key)
        array = array_keys[context["array_id"]]
        seq = array["source_sequence"]
        if context["orientation"] == "reverse":
            seq = seq.translate(_COMPLEMENT)[::-1]
        if context["rna_sequence"] != seq.replace("T", "U") or context["rna_sha256"] != digest(context["rna_sequence"]):
            raise ValueError("Context reconstruction failed")
        window = min(len(seq), requested_model["window"])
        span = min(window, requested_model["span"])
        if (context.get("requested_window"), context.get("requested_span"), context.get("window"), context.get("span")) != (
            requested_model["window"], requested_model["span"], window, span
        ):
            raise ValueError("Repeat local folding settings differ from request")
        rna = context["rna_sequence"]
        reason = "ambiguous_bases" if set(rna) - set("ACGU") else "context_too_short" if len(rna) < 5 else None
        if context.get("method") != "local_window_pair_averages" or context["status"] != ("unsupported" if reason else "completed"):
            raise ValueError("Invalid context folding method or status")
        if reason:
            if context.get("reason") != reason or "pairs" in context or "unpaired" in context:
                raise ValueError("Unsupported context invents folding evidence")
            pairs = {}
        else:
            pairs = _pairs(context["pairs"], rna)
            if any(j - i + 1 > span for i, j in pairs):
                raise ValueError("Context pair exceeds configured span")
            _unpaired(context["unpaired"], len(seq))
        contexts[context["id"]] = (context, pairs)
    consensus_keys = set()
    for consensus in report["consensus_folds"]:
        key = consensus["array_id"], consensus["orientation"]
        if key not in context_keys or key in consensus_keys:
            raise ValueError("Invalid consensus fold identity")
        consensus_keys.add(key)
        sequence = array_keys[key[0]]["consensus_repeat"]
        if key[1] == "reverse":
            sequence = sequence.translate(_COMPLEMENT)[::-1]
        rna = sequence.replace("T", "U")
        if consensus["rna_sequence"] != rna:
            raise ValueError("Consensus sequence mismatch")
        _isolated(consensus["isolated"], rna)
    if consensus_keys != context_keys:
        raise ValueError("Incomplete consensus fold coverage")
    if len(expected) != len(report["instances"]):
        raise ValueError("Repeat instance coverage mismatch")
    for number, (row, bound) in enumerate(zip(report["instances"], expected, strict=True), 1):
        if any(value is not None and type(value) is not int for value in (
            row["repeat_index"], row["start"], row["end"], row["context_start"]
        )):
            raise ValueError("Repeat coordinates must be integers")
        if (
            tuple(
                row[key]
                for key in ("source_id", "array_id", "repeat_index", "orientation", "start", "end", "source_sequence")
            )
            != bound
        ):
            raise ValueError("Repeat instance provenance mismatch")
        sequence = bound[-1]
        if row["orientation"] == "reverse":
            sequence = sequence.translate(_COMPLEMENT)[::-1]
        if row["rna_sequence"] != sequence.replace("T", "U") or row["rna_sha256"] != digest(row["rna_sequence"]):
            raise ValueError("Repeat RNA mapping mismatch")
        array = array_keys.get(bound[1])
        offset = (bound[4] - array["start"] if bound[3] == "forward" else array["end"] - bound[5]) if array else None
        terminal = bool(array and bound[2] in {1, len(array["repeats"])})
        if row["id"] != f"repeat_{number:06d}" or row["context_start"] != offset or type(row["terminal"]) is not bool or row["terminal"] != terminal:
            raise ValueError("Repeat context offset or identity mismatch")
        context, q = contexts[row["context_id"]] if array else (None, {})
        if array and (context["array_id"], context["orientation"]) != (array["id"], bound[3]):
            raise ValueError("Repeat bound to a different context")
        if not array and row["context_id"] is not None:
            raise ValueError("Isolated repeat invents a context")
        values, motif = _isolated(row["isolated"], row["rna_sequence"])
        reason, iso, ctx = None, None, None
        if values is None:
            reason = "ambiguous_bases"
        elif not motif:
            reason = "no_reference_pairs"
        else:
            iso = sum(values.get(tuple(pair), 0) for pair in motif) / len(motif)
            if context is None:
                reason = "context_not_supplied"
            elif context["status"] != "completed":
                reason = context["reason"]
            elif any(j - i + 1 > context["span"] for i, j in motif):
                reason = "motif_exceeds_context_span"
            else:
                ctx = sum(q.get((i + offset, j + offset), 0) for i, j in motif) / len(motif)
        comp = row["comparison"]
        if comp["motif_pairs"] != motif or comp["motif_origin"] != "isolated_mfe_hypothesis" or comp["reason"] != reason or comp["status"] != ("unsupported" if reason else "completed"):
            raise ValueError("Invalid comparison evidence status")
        for key, value in {"isolated_pair_support": iso, "context_pair_support": ctx,
                           "support_change": ctx - iso if ctx is not None else None}.items():
            if (value is None and comp.get(key) is not None) or (value is not None and (not _finite(comp.get(key)) or abs(comp[key] - value) > 1e-8)):
                raise ValueError("Comparison arithmetic mismatch")
    if any(type(value) is not int for value in report["counts"].values()) or report["counts"] != {
        "sources": len(sources),
        "arrays": len(array_keys),
        "contexts": len(contexts),
        "instances": len(expected),
    }:
        raise ValueError("Repeat summary counts mismatch")
