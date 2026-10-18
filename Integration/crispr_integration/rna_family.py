"""Homologous RNA-family execution with source-bound alignment evidence."""
from collections import Counter
import hashlib
import math
from pathlib import Path
import re

from .protospacer import _json, _payload

POLICY = {"homology": "user_asserted_not_tested", "orientation": "as_submitted_5prime_to_3prime",
          "boundaries": "user_supplied_not_experimentally_verified", "deduplication": "exact_rna_sequence",
          "weights": "equal_per_distinct_sequence", "temperature_celsius": 37,
          "energy_parameters": "RNA_Turner2004", "covariance_factor": 1.0, "noncanonical_factor": 1.0,
          "probability_threshold": 0.01, "pair_coordinates": "one_based_alignment_columns",
          "covariation": "descriptive_not_statistical_significance", "activity": "not_assessed"}
LIMITS = {"min_records": 2, "max_records": 50, "min_length": 30, "max_length": 500,
          "max_alignment_columns": 1000, "max_input_bytes": 100000}
CANONICAL = {"AU", "UA", "GC", "CG", "GU", "UG"}


def sha(value):
    return hashlib.sha256(value.encode("ascii") if isinstance(value, str) else value).hexdigest()


def records_from_input(source, molecule):
    records = {}
    for path in sorted(source.iterdir()) if source.is_dir() else [source]:
        data = _payload(source if source.is_dir() else source.parent, path)
        if len(data) > 100000:
            raise ValueError("RNA input byte limit exceeded")
        key = None
        for line in data.decode("ascii").splitlines():
            if line.startswith(">"):
                key = line[1:].split()[0]
                if not re.fullmatch(r"[A-Za-z0-9_.:-]{1,200}", key) or key in records:
                    raise ValueError("Invalid RNA source identifier")
                records[key] = ""
            elif line.strip():
                if key is None:
                    raise ValueError("RNA sequence before header")
                records[key] += line.strip().upper()
    alphabet = set("ACGU" if molecule == "RNA" else "ACGT")
    if not 2 <= len(records) <= 50 or any(not 30 <= len(s) <= 500 or set(s) - alphabet for s in records.values()):
        raise ValueError("Supply 2–50 homologous sequences of 30–500 unambiguous nucleotides")
    if len({s.replace("T", "U") for s in records.values()}) < 2:
        raise ValueError("Supply at least two distinct RNA sequences")
    return records


def validate_report(doc, records, molecule, input_sha256):
    if (doc["schema"] != {"name": "rnaalifold.family", "version": "1.0.0"}
            or doc["tool"] not in ({"name": "RNAalifold-v2", "version": "1.0.0"}, {"name": "RNAalifold-v2", "version": "1.0.1"}) or doc["status"] != "completed"
            or doc["policy"] != POLICY or doc["limits"] != LIMITS
            or doc["input"] != {"sha256": input_sha256, "molecule": molecule}
            or doc["engine"] != {"name": "ViennaRNA", "version": "2.7.0"}
            or doc["aligner"]["name"] != "MAFFT" or doc["aligner"]["strategy"] != "FFT-NS-2"
            or not re.fullmatch(r"7\.[0-9]+", doc["aligner"]["version"])
            or type(doc["aligner"]["threads"]) is not int or not 1 <= doc["aligner"]["threads"] <= 8):
        raise ValueError("RNA family identity mismatch")
    groups, sources = {}, []
    for key, seq in records.items():
        rna = seq.replace("T", "U")
        if rna not in groups:
            groups[rna] = {"id": f"family{len(groups)+1:06d}", "member_ids": []}
        groups[rna]["member_ids"].append(key)
        sources.append({"id": key, "length": len(seq), "sha256": sha(seq), "rna_sha256": sha(rna),
                        "representative_id": groups[rna]["id"]})
    if doc["sources"] != sources or len(doc["alignment"]) != len(groups) or not 2 <= len(groups) <= 50:
        raise ValueError("RNA source or duplicate membership mismatch")
    aligned = []
    for row, (rna, group) in zip(doc["alignment"], groups.items(), strict=True):
        if row["id"] != group["id"] or row["member_ids"] != group["member_ids"] or row["sequence"].replace("-", "") != rna:
            raise ValueError("RNA alignment does not reconstruct input")
        aligned.append(row["sequence"])
    lengths = {len(s) for s in aligned}
    if len(lengths) != 1 or not 30 <= len(aligned[0]) <= 1000 or any(set(s)-set("ACGU-") for s in aligned):
        raise ValueError("Invalid RNA alignment shape")
    consensus, columns = "", []
    for i, values in enumerate(zip(*aligned, strict=True), 1):
        counts = {key: values.count(key) for key in "ACGU-"}
        maximum = max(counts[k] for k in "ACGU")
        if not maximum:
            raise ValueError("All-gap RNA alignment column")
        winners = [key for key in "ACGU" if counts[key] == maximum]
        consensus += winners[0] if len(winners) == 1 else "N"
        columns.append({"column": i, "counts": counts})
    fold = doc["fold"]
    if (doc["columns"] != columns or fold["consensus"] != consensus or fold["temperature_celsius"] != 37
            or any(type(row["column"]) is not int or any(type(v) is not int for v in row["counts"].values())
                   for row in doc["columns"])):
        raise ValueError("Invalid RNA column evidence")
    structure = fold["structure"]
    if len(structure) != len(consensus) or any(type(fold[k]) not in (int, float) or not math.isfinite(fold[k])
            for k in ("mfe_score_kcal_mol", "ensemble_free_energy_kcal_mol")):
        raise ValueError("Invalid comparative structure or energy")
    stack, pairs = [], []
    for i, char in enumerate(structure, 1):
        if char == "(":
            stack.append(i)
        elif char == ")" and stack:
            pairs.append((stack.pop(), i))
        elif char != ".":
            raise ValueError("Unbalanced RNA structure")
    if stack:
        raise ValueError("Unbalanced RNA structure")
    pair_set = set(pairs)
    if len(doc["probabilities"]) > len(consensus)*(len(consensus)-1)//2:
        raise ValueError("Too many RNA pair probabilities")
    probabilities, masses, last = {}, Counter(), (0, 0)
    for row in doc["probabilities"]:
        left, right, value = row["left"], row["right"], row["probability"]
        key = (left, right)
        if (type(left) is not int or type(right) is not int or not 1 <= left < right <= len(consensus)
                or key <= last or type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1
                or value < 0.01 and key not in pair_set):
            raise ValueError("Invalid RNA pair probability")
        probabilities[key] = value
        masses[left] += value
        masses[right] += value
        last = key
    if any(v > 1.00001 for v in masses.values()) or not pair_set <= probabilities.keys():
        raise ValueError("Invalid RNA probability mass or missing MFE pair")
    evidence = []
    for left, right in sorted(pairs):
        types = Counter(row[left-1]+row[right-1] for row in aligned)
        canonical = sum(v for k, v in types.items() if k in CANONICAL)
        gaps = sum(v for k, v in types.items() if "-" in k)
        compensatory = sum(types[a]*types[b] for a in types for b in types if a < b and a in CANONICAL
                           and b in CANONICAL and a[0] != b[0] and a[1] != b[1])
        evidence.append({"left": left, "right": right, "probability": probabilities[left, right],
                         "canonical_sequences": canonical, "gapped_sequences": gaps,
                         "noncanonical_sequences": len(aligned)-canonical-gaps, "pair_types": dict(sorted(types.items())),
                         "compensatory_sequence_pairs": compensatory})
    counts = {"sources": len(sources), "unique_sequences": len(aligned), "alignment_columns": len(consensus),
              "mfe_pairs": len(pairs), "retained_probabilities": len(probabilities)}
    if (doc["pairs"] != evidence or doc["counts"] != counts or any(type(v) is not int for v in doc["counts"].values())
            or any(any(type(row[key]) is not int for key in ("left", "right", "canonical_sequences", "noncanonical_sequences",
                       "gapped_sequences", "compensatory_sequence_pairs")) or any(type(v) is not int for v in row["pair_types"].values())
                   for row in doc["pairs"])):
        raise ValueError("RNA pair evidence or counts differ")


def run_rna_family_pipeline(input_path, output_path, *, runner_config, cpu=1, timeout_seconds=3600, molecule="DNA"):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process
    if molecule not in ("DNA", "RNA"):
        raise ValueError("Invalid RNA-family input molecule")
    records = records_from_input(Path(input_path).absolute(), molecule)
    job = _empty_job_root(output_path)
    (job / "logs").mkdir()
    (job / "rna-family-input").mkdir()
    fasta = job / "rna-family-input/sequences.fasta"
    fasta.write_text("".join(f">{key}\n{seq}\n" for key, seq in records.items()), encoding="ascii")
    result = {"schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"}, "workflow": "rna_family",
              "status": "running", "stages": {}, "tool_options": runner_config.tool_options}
    _atomic_json(job / "integration_result.json", result)
    try:
        stage = _run_process("rnaalifold", [*runner_config.rnaalifold, "--input", str(fasta), "--output", str(job / "rna-family"),
            "--molecule", molecule, "--threads", str(min(cpu, 8)), "--timeout-seconds", str(max(1, min(math.ceil(timeout_seconds), 14400)))],
            job / "logs", timeout_seconds=timeout_seconds)
        doc = _json(_payload(job, job / "rna-family/report.json"))
        validate_report(doc, records, molecule, sha(fasta.read_bytes()))
        stage["counts"] = doc["counts"]
        result["stages"]["rnaalifold"] = stage
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("rnaalifold", str(exc)) from exc
    _atomic_json(job / "integration_result.json", result)
    return result
