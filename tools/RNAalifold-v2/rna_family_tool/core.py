"""Source-bound, redundancy-aware comparative RNA evidence."""
from collections import Counter
from itertools import combinations
import hashlib
import math
import re
from typing import TypedDict

PROBABILITY_THRESHOLD = 0.01
POLICY = {"homology": "user_asserted_not_tested", "orientation": "as_submitted_5prime_to_3prime",
          "boundaries": "user_supplied_not_experimentally_verified", "deduplication": "exact_rna_sequence",
          "weights": "equal_per_distinct_sequence", "temperature_celsius": 37,
          "energy_parameters": "RNA_Turner2004", "covariance_factor": 1.0, "noncanonical_factor": 1.0,
          "probability_threshold": PROBABILITY_THRESHOLD, "pair_coordinates": "one_based_alignment_columns",
          "covariation": "descriptive_not_statistical_significance", "activity": "not_assessed"}
LIMITS = {"min_records": 2, "max_records": 50, "min_length": 30, "max_length": 500,
          "max_alignment_columns": 1000, "max_input_bytes": 100000}
CANONICAL = frozenset(("AU", "UA", "GC", "CG", "GU", "UG"))
ID = re.compile(r"[A-Za-z0-9_.:-]{1,200}\Z")


class Representative(TypedDict):
    id: str
    member_ids: list[str]
    sequence: str


def sha(value: str | bytes) -> str:
    return hashlib.sha256(value.encode("ascii") if isinstance(value, str) else value).hexdigest()


def fasta(payload: bytes, *, molecule: str = "RNA", aligned: bool = False) -> dict[str, str]:
    if molecule not in {"RNA", "DNA"}:
        raise ValueError("Unsupported molecule; choose RNA or DNA")
    records: dict[str, str] = {}
    key: str | None = None
    sequence: list[str] = []
    alphabet = set("ACGU" if molecule == "RNA" else "ACGT") | ({"-"} if aligned else set())
    if len(payload) > LIMITS["max_input_bytes"]:
        raise ValueError("Input exceeds byte limit")
    lines = payload.decode("ascii").splitlines()
    if any(line.startswith(">") and not line[1:].strip() for line in lines):
        raise ValueError("Empty FASTA header")
    for line in lines + [">"]:
        if line.startswith(">"):
            if key is not None:
                seq = "".join(sequence).upper()
                if not seq or set(seq) - alphabet or key in records:
                    raise ValueError("Invalid or duplicate FASTA sequence")
                records[key] = seq
                if len(records) > LIMITS["max_records"]:
                    raise ValueError("Too many RNA family records")
            if line == ">":
                key = None
                continue
            key, sequence = line[1:].split()[0], []
            if not ID.fullmatch(key):
                raise ValueError("Invalid FASTA identifier")
        elif line.strip():
            if key is None:
                raise ValueError("Sequence before header")
            sequence.append(line.strip())
    if not aligned and (not LIMITS["min_records"] <= len(records) <= LIMITS["max_records"] or
                        any(not LIMITS["min_length"] <= len(s) <= LIMITS["max_length"] for s in records.values())):
        raise ValueError("Supply 2–50 sequences of 30–500 nucleotides")
    return records


def representatives(records: dict[str, str]):
    groups: dict[str, Representative] = {}
    sources = []
    for key, sequence in records.items():
        rna = sequence.replace("T", "U")
        if rna not in groups:
            groups[rna] = {"id": f"family{len(groups)+1:06d}", "member_ids": [], "sequence": rna}
        groups[rna]["member_ids"].append(key)
        sources.append({"id": key, "length": len(sequence), "sha256": sha(sequence),
                        "rna_sha256": sha(rna), "representative_id": groups[rna]["id"]})
    if len(groups) < 2:
        raise ValueError("At least two distinct RNA sequences are required after deduplication")
    return sources, list(groups.values())


def structure_pairs(structure: str) -> list[tuple[int, int]]:
    stack, result = [], []
    for i, char in enumerate(structure, 1):
        if char == "(":
            stack.append(i)
        elif char == ")":
            if not stack:
                raise ValueError("Unbalanced consensus structure")
            result.append((stack.pop(), i))
        elif char != ".":
            raise ValueError("Unsupported structure notation")
    if stack:
        raise ValueError("Unbalanced consensus structure")
    return sorted(result)


def columns(alignment):
    result, consensus = [], []
    for index, values in enumerate(zip(*alignment, strict=True), 1):
        observed = Counter(values)
        counts = {key: observed.get(key, 0) for key in "ACGU-"}
        maximum = max(counts[key] for key in "ACGU")
        bases = [key for key in "ACGU" if counts[key] == maximum and maximum]
        consensus.append(bases[0] if len(bases) == 1 else "N")
        result.append({"column": index, "counts": counts})
    return "".join(consensus), result


def pair_evidence(alignment, left, right, probability):
    types = Counter(row[left-1] + row[right-1] for row in alignment)
    canonical = sum(count for pair, count in types.items() if pair in CANONICAL)
    gaps = sum(count for pair, count in types.items() if "-" in pair)
    canonical_types = [(pair, count) for pair, count in types.items() if pair in CANONICAL]
    compensatory = sum(a_count * b_count for (a, a_count), (b, b_count) in combinations(canonical_types, 2)
                       if a[0] != b[0] and a[1] != b[1])
    return {"left": left, "right": right, "probability": probability, "canonical_sequences": canonical,
            "noncanonical_sequences": len(alignment)-canonical-gaps, "gapped_sequences": gaps,
            "pair_types": dict(sorted(types.items())), "compensatory_sequence_pairs": compensatory}


def alignment_rows(groups: list[Representative], aligned: dict[str, str]) -> list[Representative]:
    """Bind native alignment order and every ungapped RNA to its representative."""
    if list(aligned) != [row["id"] for row in groups] or len({len(s) for s in aligned.values()}) != 1:
        raise ValueError("MAFFT changed sequence identity or alignment shape")
    result: list[Representative] = []
    for row in groups:
        value = aligned[row["id"]].replace("T", "U")
        if value.replace("-", "") != row["sequence"] or len(value) > LIMITS["max_alignment_columns"]:
            raise ValueError("MAFFT changed sequence or exceeded alignment bounds")
        result.append({**row, "sequence": value})
    if any(all(base == "-" for base in column) for column in zip(*(row["sequence"] for row in result), strict=True)):
        raise ValueError("All-gap alignment column")
    return result


def fold(alignment: list[str]):
    import RNA
    if RNA.__version__ != "2.7.0":
        raise ValueError("ViennaRNA 2.7.0 is required")
    if (not LIMITS["min_records"] <= len(alignment) <= LIMITS["max_records"] or
            len({len(row) for row in alignment}) != 1 or
            not 1 <= len(alignment[0]) <= LIMITS["max_alignment_columns"]):
        raise ValueError("Invalid comparative alignment shape")
    RNA.params_load_RNA_Turner2004()
    # Explicit constructor settings also rebuild native pair tables. RNA.md()
    # alone inherits mutable process-wide defaults, including energy_set.
    model = RNA.md(
        temperature=37.0, cv_fact=1.0, nc_fact=1.0, compute_bpp=1,
        dangles=2, betaScale=1.0, pf_smooth=1, special_hp=1,
        noLP=0, noGU=0, noGUclosure=0, logML=0, circ=0, circ_penalty=1,
        gquad=0, uniq_ML=0, energy_set=0, backtrack=1, backtrack_type="F",
        max_bp_span=-1, min_loop_size=3, window_size=-1, oldAliEn=0, ribo=0,
        sfact=1.07, salt=1.021, saltMLLower=6, saltMLUpper=24,
        saltDPXInit=99999, saltDPXInitFact=-45.324, helical_rise=2.8, backbone_length=6.0,
    )
    # The native comparative Boltzmann table includes family-size scaling.
    # Snapshot it with an evaluation-only compound (no DP matrices), then copy
    # both tables so later process-global parameter replacement cannot leak in.
    snapshot = RNA.fold_compound(alignment, model, RNA.OPTION_PF | RNA.OPTION_EVAL_ONLY)
    compound = RNA.fold_compound(alignment, model, RNA.OPTION_MFE | RNA.OPTION_PF)
    compound.params_subst(snapshot.params)
    compound.exp_params_subst(snapshot.exp_params)
    del snapshot
    structure, mfe = compound.mfe()
    if len(structure) != len(alignment[0]):
        raise ValueError("Consensus structure length differs from alignment")
    mfe_pairs = set(structure_pairs(structure))
    if not math.isfinite(mfe):
        raise ValueError("Nonfinite comparative folding energy")
    compound.exp_params_rescale(mfe)
    _, ensemble = compound.pf()
    if not math.isfinite(ensemble):
        raise ValueError("Nonfinite comparative folding energy")
    matrix = compound.bpp()
    if len(matrix) != len(structure) + 1 or any(len(row) != len(matrix) for row in matrix):
        raise ValueError("Invalid base-pair probability matrix shape")
    probabilities = []
    for left in range(1, len(structure)+1):
        for right in range(left+1, len(structure)+1):
            value = float(matrix[left][right])
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError("Invalid base-pair probability")
            if value >= PROBABILITY_THRESHOLD or (left, right) in mfe_pairs:
                probabilities.append({"left": left, "right": right, "probability": value})
    return {"engine": "ViennaRNA", "version": RNA.__version__, "structure": structure,
            "mfe_score_kcal_mol": float(mfe), "ensemble_free_energy_kcal_mol": float(ensemble),
            "probabilities": probabilities}
