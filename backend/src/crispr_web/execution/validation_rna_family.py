"""Independent immutable-input, alignment and native-fold gate for RNA families.

Scientific engines and Integration validators are intentionally not imported.
"""
from collections import Counter
import hashlib
import math
import re

from .validation_annotations import _payload, _strict_json

_POLICY = {"homology": "user_asserted_not_tested", "orientation": "as_submitted_5prime_to_3prime",
           "boundaries": "user_supplied_not_experimentally_verified", "deduplication": "exact_rna_sequence",
           "weights": "equal_per_distinct_sequence", "temperature_celsius": 37,
           "energy_parameters": "RNA_Turner2004", "covariance_factor": 1.0, "noncanonical_factor": 1.0,
           "probability_threshold": 0.01, "pair_coordinates": "one_based_alignment_columns",
           "covariation": "descriptive_not_statistical_significance", "activity": "not_assessed"}
_LIMITS = {"min_records": 2, "max_records": 50, "min_length": 30, "max_length": 500,
           "max_alignment_columns": 1000, "max_input_bytes": 100000}
_CANONICAL = {"AU", "UA", "GC", "CG", "GU", "UG"}


def _sha(value):
    return hashlib.sha256(value.encode("ascii") if isinstance(value, str) else value).hexdigest()


def _finite(value, lower=None, upper=None):
    return (type(value) in (int, float) and math.isfinite(value)
            and (lower is None or value >= lower) and (upper is None or value <= upper))


def _fasta(data, alphabet):
    result, key = {}, None
    for line in data.decode("ascii").splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            if not re.fullmatch(r"[A-Za-z0-9_.:-]{1,200}", key) or key in result:
                raise ValueError("Invalid family identifier")
            result[key] = ""
        elif line.strip():
            if key is None:
                raise ValueError("Missing family header")
            result[key] += line.strip().upper()
    if not result or any(not seq or set(seq) - alphabet for seq in result.values()):
        raise ValueError("Invalid family sequence")
    return result


def validate_rna_family_pipeline(claimed, root, expected, integration, input_root, checkpoint=None):
    try:
        return _validate(claimed, root, expected, integration, input_root, checkpoint)
    except (OSError, ValueError, TypeError, KeyError, IndexError, AttributeError, UnicodeError, RecursionError):
        return False


def _validate(claimed, root, expected, integration, input_root, checkpoint):
    if input_root is None or claimed.molecule not in ("DNA", "RNA") or integration.get("workflow") != "rna_family":
        return False
    records = {}
    for path in sorted(input_root.iterdir()):
        incoming = _fasta(_payload(input_root, path, checkpoint), set("ACGU" if claimed.molecule == "RNA" else "ACGT"))
        if len(incoming) != 1:
            return False
        for key, seq in incoming.items():
            if key in records or expected.get(key) != (len(seq), _sha(seq)):
                return False
            records[key] = seq
    if set(records) != set(expected) or not 2 <= len(records) <= 50 or any(not 30 <= len(seq) <= 500 for seq in records.values()):
        return False
    doc = _strict_json(_payload(root, root / "rna-family/report.json", checkpoint))
    if (doc["schema"] != {"name": "rnaalifold.family", "version": "1.0.0"}
            or doc["tool"] not in ({"name": "RNAalifold-v2", "version": "1.0.0"}, {"name": "RNAalifold-v2", "version": "1.0.1"}) or doc["status"] != "completed"
            or doc["policy"] != _POLICY or doc["limits"] != _LIMITS
            or doc["input"] != {"sha256": _sha("".join(f">{key}\n{seq}\n" for key, seq in records.items())), "molecule": claimed.molecule}
            or doc["engine"] != {"name": "ViennaRNA", "version": "2.7.0"}
            or doc["aligner"]["name"] != "MAFFT" or doc["aligner"]["strategy"] != "FFT-NS-2"
            or not re.fullmatch(r"7\.[0-9]+", doc["aligner"]["version"])
            or type(doc["aligner"]["threads"]) is not int or not 1 <= doc["aligner"]["threads"] <= 8):
        return False
    group_sequences, groups, sources = [], [], []
    for key, seq in records.items():
        sequence = seq.replace("T", "U")
        if sequence not in group_sequences:
            group_sequences.append(sequence)
            groups.append({"id": f"family{len(groups)+1:06d}", "member_ids": []})
        index = group_sequences.index(sequence)
        groups[index]["member_ids"].append(key)
        sources.append({"id": key, "length": len(seq), "sha256": _sha(seq), "rna_sha256": _sha(sequence),
                        "representative_id": groups[index]["id"]})
    if not 2 <= len(groups) <= 50 or doc["sources"] != sources or len(doc["alignment"]) != len(groups):
        return False
    aligned = []
    for index, row in enumerate(doc["alignment"]):
        sequence = row["sequence"]
        if (row["id"] != groups[index]["id"] or row["member_ids"] != groups[index]["member_ids"]
                or sequence.replace("-", "") != group_sequences[index] or set(sequence)-set("ACGU-")
                or not 30 <= len(sequence) <= 1000):
            return False
        aligned.append(sequence)
    if len({len(seq) for seq in aligned}) != 1:
        return False
    length = len(aligned[0])
    columns, consensus = [], ""
    for index in range(length):
        counts = {base: sum(seq[index] == base for seq in aligned) for base in "ACGU-"}
        maximum = max(counts[base] for base in "ACGU")
        if not maximum:
            return False
        winners = [base for base in "ACGU" if counts[base] == maximum]
        consensus += winners[0] if len(winners) == 1 else "N"
        columns.append({"column": index+1, "counts": counts})
    fold = doc["fold"]
    if (doc["columns"] != columns or fold["consensus"] != consensus or fold["temperature_celsius"] != 37
            or not _finite(fold["mfe_score_kcal_mol"]) or not _finite(fold["ensemble_free_energy_kcal_mol"])
            or any(type(row["column"]) is not int or any(type(v) is not int for v in row["counts"].values())
                   for row in doc["columns"])):
        return False
    structure = fold["structure"]
    if len(structure) != length:
        return False
    opening, pairs = [], []
    for index, mark in enumerate(structure, 1):
        if mark == "(":
            opening.append(index)
        elif mark == ")":
            if not opening:
                return False
            left = opening.pop()
            if index-left <= 3:
                return False
            pairs.append((left, index))
        elif mark != ".":
            return False
    if opening:
        return False
    pairs.sort()
    probabilities, mass = {}, [0.0]*(length+1)
    if len(doc["probabilities"]) > length*(length-1)//2:
        return False
    previous = (0, 0)
    for row in doc["probabilities"]:
        if checkpoint:
            checkpoint()
        left, right, probability = row["left"], row["right"], row["probability"]
        if (type(left) is not int or type(right) is not int or not 1 <= left < right <= length
                or (left, right) <= previous or not _finite(probability, 0, 1)
                or probability < 0.01 and (left, right) not in pairs):
            return False
        probabilities[left, right] = probability
        mass[left] += probability
        mass[right] += probability
        previous = left, right
    if max(mass) > 1.00001 or any(pair not in probabilities for pair in pairs):
        return False
    evidence = []
    for left, right in pairs:
        types = Counter(seq[left-1]+seq[right-1] for seq in aligned)
        canonical = sum(types[pair] for pair in _CANONICAL)
        gaps = sum(count for pair, count in types.items() if "-" in pair)
        covarying = 0
        for a_index, a in enumerate(sorted(types)):
            for b in sorted(types)[a_index+1:]:
                if a in _CANONICAL and b in _CANONICAL and a[0] != b[0] and a[1] != b[1]:
                    covarying += types[a]*types[b]
        evidence.append({"left": left, "right": right, "probability": probabilities[left, right],
                         "canonical_sequences": canonical, "noncanonical_sequences": len(aligned)-canonical-gaps,
                         "gapped_sequences": gaps, "pair_types": dict(sorted(types.items())),
                         "compensatory_sequence_pairs": covarying})
    counts = {"sources": len(sources), "unique_sequences": len(groups), "alignment_columns": length,
              "mfe_pairs": len(pairs), "retained_probabilities": len(probabilities)}
    if (doc["pairs"] != evidence or doc["counts"] != counts or any(type(v) is not int for v in doc["counts"].values())
            or any(any(type(row[key]) is not int for key in ("left", "right", "canonical_sequences", "noncanonical_sequences",
                       "gapped_sequences", "compensatory_sequence_pairs")) or any(type(v) is not int for v in row["pair_types"].values())
                   for row in doc["pairs"])):
        return False
    native_alignment = _fasta(_payload(root, root / "rna-family/work/aligned.fasta", checkpoint), set("ACGT-"))
    if list(native_alignment) != [row["id"] for row in doc["alignment"]] or any(
            native_alignment[row["id"]].replace("T", "U") != row["sequence"] for row in doc["alignment"]):
        return False
    native_fold = _strict_json(_payload(root, root / "rna-family/work/fold.json", checkpoint))
    if (native_fold["engine"] != "ViennaRNA" or native_fold["version"] != "2.7.0"
            or native_fold["alignment_sha256"] != _sha("".join(f'>{row["id"]}\n{row["sequence"]}\n' for row in doc["alignment"]))
            or any(native_fold[key] != fold[key] for key in ("structure", "mfe_score_kcal_mol", "ensemble_free_energy_kcal_mol"))
            or native_fold["probabilities"] != doc["probabilities"]):
        return False
    stage = integration["stages"]["rnaalifold"]
    return stage["status"] == "completed" and stage["counts"] == counts
