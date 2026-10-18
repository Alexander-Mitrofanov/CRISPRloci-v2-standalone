"""Independent backend gate: source CDSs, genetic codes and alignment evidence.

No Integration or scientific-runtime imports are permitted here. Static NCBI
code assignments are checked locally; the scientific engine uses Biopython.
"""

import hashlib
from itertools import product, combinations
import math
import re

BASES = "TCAG"
STANDARD = dict(
    zip(
        ("".join(c) for c in product(BASES, repeat=3)),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=True,
    )
)
# NCBI genetic-code assignments; no initiator-M substitution (CasAndra policy).
CHANGES = {
    1: {},
    2: {"AGA": "*", "AGG": "*", "TGA": "W", "ATA": "M"},
    3: {"TGA": "W", "CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T", "ATA": "M"},
    4: {"TGA": "W"},
    5: {"TGA": "W", "ATA": "M", "AGA": "S", "AGG": "S"},
    6: {"TAA": "Q", "TAG": "Q"},
    9: {"TGA": "W", "AAA": "N", "AGA": "S", "AGG": "S"},
    10: {"TGA": "C"},
    11: {},
    12: {"CTG": "S"},
    13: {"TGA": "W", "ATA": "M", "AGA": "G", "AGG": "G"},
    14: {"TAA": "Y", "TGA": "W", "AAA": "N", "AGA": "S", "AGG": "S"},
    15: {"TAG": "Q"},
    16: {"TAG": "L"},
    21: {"TGA": "W", "ATA": "M", "AAA": "N", "AGA": "S", "AGG": "S"},
    22: {"TCA": "*", "TAG": "L"},
    23: {"TTA": "*"},
    24: {"TGA": "W", "AGA": "S", "AGG": "K"},
    25: {"TGA": "G"},
    26: {"CTG": "A"},
    27: {"TGA": "W", "TAA": "Q", "TAG": "Q"},
    28: {"TAA": "Q", "TAG": "Q", "TGA": "W"},
    29: {"TAA": "Y", "TAG": "Y"},
    30: {"TAA": "E", "TAG": "E"},
    31: {"TAA": "E", "TAG": "E", "TGA": "W"},
    32: {"TAG": "W"},
    33: {"TAA": "Y", "TGA": "W", "AGA": "S", "AGG": "K"},
}
IUPAC = dict(
    zip(
        "ACGTRYSWKMBDHVN",
        (
            "A",
            "C",
            "G",
            "T",
            "AG",
            "CT",
            "GC",
            "AT",
            "GT",
            "AC",
            "CGT",
            "AGT",
            "ACT",
            "ACG",
            "ACGT",
        ),
        strict=True,
    )
)
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
GENE_FIELDS = {
    "id",
    "start",
    "end",
    "strand",
    "phase",
    "translation_table",
    "partial_5prime",
    "partial_3prime",
    "translation",
}
LIMITS = {
    "loci": 8,
    "window_nt": 50000,
    "genes": 120,
    "protein_aa": 3000,
    "protein_pairs": 6000,
    "alignment_cells": 250000000,
}
ALIGNER = {
    "mode": "global",
    "substitution_matrix": "BLOSUM62",
    "open_gap_score": -10.0,
    "extend_gap_score": -0.5,
}
INTERPRETATION = "Protein similarity and observed gene order; not orthology, function, or a phylogeny."


def sha(value):
    return hashlib.sha256(
        value.encode() if isinstance(value, str) else value
    ).hexdigest()


def translate(dna, table):
    code = {**STANDARD, **CHANGES[table]}
    result = []
    for offset in range(0, len(dna), 3):
        codon = dna[offset : offset + 3]
        residues = {
            code["".join(bases)] for bases in product(*(IUPAC[x] for x in codon))
        }
        result.append(
            next(iter(residues))
            if len(residues) == 1
            else "B"
            if residues == {"D", "N"}
            else "Z"
            if residues == {"E", "Q"}
            else "J"
            if residues == {"I", "L"}
            else "X"
        )
    return "".join(result).removesuffix("*")


def _id(value):
    if not isinstance(value, str) or not re.fullmatch(
        r"[A-Za-z0-9_.:|+%=\-]{1,512}", value
    ):
        raise ValueError("Invalid locus or gene identity")


def expected_loci(features, records):
    if (
        set(features) != {"schema", "coordinate_system", "loci"}
        or features["schema"] != {"name": "clinker.features", "version": "1.0.0"}
        or features["coordinate_system"] != "1-based-inclusive-source-forward"
    ):
        raise ValueError("Invalid locus input contract")
    if not isinstance(features["loci"], list) or not 2 <= len(features["loci"]) <= 8:
        raise ValueError("Locus count outside bounds")
    loci = []
    seen = set()
    count = 0
    for locus in features["loci"]:
        if set(locus) != {"id", "source_id", "start", "end", "genes"}:
            raise ValueError("Invalid locus fields")
        _id(locus["id"])
        _id(locus["source_id"])
        if locus["id"] in seen:
            raise ValueError("Duplicate locus identity")
        seen.add(locus["id"])
        seq = records[locus["source_id"]]
        start, end = locus["start"], locus["end"]
        if (
            type(start) is not int
            or type(end) is not int
            or not 1 <= start <= end <= len(seq)
            or end - start + 1 > 50000
        ):
            raise ValueError("Invalid locus boundaries")
        genes = []
        ids = set()
        previous = None
        if not isinstance(locus["genes"], list) or len(locus["genes"]) > 120:
            raise ValueError("Invalid gene count")
        for gene in locus["genes"]:
            if set(gene) != GENE_FIELDS:
                raise ValueError("Invalid CDS fields")
            _id(gene["id"])
            gs, ge = gene["start"], gene["end"]
            phase = gene["phase"]
            table = gene["translation_table"]
            if (
                gene["id"] in ids
                or type(gs) is not int
                or type(ge) is not int
                or not start <= gs <= ge <= end
                or type(phase) is not int
                or phase not in (0, 1, 2)
                or type(table) is not int
                or table not in CHANGES
                or gene["strand"] not in ("+", "-")
                or any(
                    type(gene[k]) is not bool
                    for k in ("partial_5prime", "partial_3prime")
                )
            ):
                raise ValueError(
                    "Invalid CDS coordinates, code, phase or partial evidence"
                )
            ids.add(gene["id"])
            order = (gs, ge, gene["id"])
            if previous is not None and order <= previous:
                raise ValueError("CDS order differs")
            previous = order
            nucleotide = seq[gs - 1 : ge]
            if gene["strand"] == "-":
                nucleotide = nucleotide.translate(COMP)[::-1]
            coding = nucleotide[phase:]
            if (
                phase
                and not gene["partial_5prime"]
                or len(coding) % 3
                and not gene["partial_3prime"]
            ):
                raise ValueError("Incomplete CDS lacks partial evidence")
            coding = coding[: len(coding) // 3 * 3]
            protein = gene["translation"]
            if (
                not isinstance(protein, str)
                or not 1 <= len(protein) <= 3000
                or set(protein) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
                or translate(coding, table) != protein
            ):
                raise ValueError("Protein differs from source translation")
            genes.append(
                {
                    **gene,
                    "nucleotide_sha256": sha(nucleotide),
                    "protein_sha256": sha(protein),
                }
            )
        count += len(genes)
        loci.append(
            {
                **locus,
                "genes": genes,
                "sequence_sha256": sha(seq[start - 1 : end]),
                "source_length": len(seq),
                "at_source_start": start == 1,
                "at_source_end": end == len(seq),
            }
        )
    if count > 120:
        raise ValueError("Total CDS limit exceeded")
    return loci


def validate_report(report, features, records, input_digest, features_digest, identity):
    if (
        set(report)
        != {
            "schema",
            "tool",
            "coordinate_system",
            "input_sha256",
            "features_sha256",
            "options",
            "limits",
            "sources",
            "loci",
            "links",
            "evaluated_protein_pairs",
            "alignment_cells",
            "interpretation",
        }
        or report["schema"] != {"name": "clinker.comparison", "version": "1.0.0"}
        or report["tool"]
        not in ({"name": "clinker-v2", "version": "2.0.0", "biopython": "1.85"}, {"name": "clinker-v2", "version": "2.0.1", "biopython": "1.85"})
        or report["coordinate_system"] != "1-based-inclusive-source-forward"
        or report["input_sha256"] != input_digest
        or report["features_sha256"] != features_digest
        or report["options"] != {"identity": identity, "jobs": 1, "aligner": ALIGNER}
        or report["limits"] != LIMITS
        or report["interpretation"] != INTERPRETATION
    ):
        raise ValueError("Unsupported clinker output or run identity")
    loci = expected_loci(features, records)
    if report["loci"] != loci or report["sources"] != [
        {"id": k, "length": len(v), "sha256": sha(v)} for k, v in records.items()
    ]:
        raise ValueError("Source-bound features differ")
    allowed = {}
    cells = 0
    for a, b in combinations(loci, 2):
        for x, y in product(a["genes"], b["genes"]):
            allowed[a["id"], x["id"], b["id"], y["id"]] = (
                x["translation"],
                y["translation"],
            )
            cells += len(x["translation"]) * len(y["translation"])
    if (
        len(allowed) > 6000
        or cells > 250000000
        or type(report["evaluated_protein_pairs"]) is not int
        or report["evaluated_protein_pairs"] != len(allowed)
        or type(report["alignment_cells"]) is not int
        or report["alignment_cells"] != cells
    ):
        raise ValueError("Comparison work or counts differ")
    seen = set()
    groups = [set(x) for x in ("GAVLI", "FYW", "CM", "ST", "KRH", "DENQ", "P")]
    if not isinstance(report["links"], list) or len(report["links"]) > len(allowed):
        raise ValueError("Alignment count differs")
    for link in report["links"]:
        if set(link) != {
            "query_locus",
            "query_gene",
            "target_locus",
            "target_gene",
            "identity",
            "similarity",
            "score",
            "aligned_query",
            "aligned_target",
        }:
            raise ValueError("Unexpected alignment fields")
        key = tuple(
            link[k]
            for k in ("query_locus", "query_gene", "target_locus", "target_gene")
        )
        if key not in allowed or key in seen:
            raise ValueError("Unbound or duplicate protein alignment")
        seen.add(key)
        a, b = link["aligned_query"], link["aligned_target"]
        if (
            not isinstance(a, str)
            or not isinstance(b, str)
            or not 1 <= len(a) == len(b) <= 6000
            or (a.replace("-", ""), b.replace("-", "")) != allowed[key]
            or any(x == y == "-" for x, y in zip(a, b, strict=True))
        ):
            raise ValueError("Gapped alignment differs from verified proteins")
        matches = sum(x == y for x, y in zip(a, b, strict=True))
        similar = sum(
            x != y and any(x in group and y in group for group in groups)
            for x, y in zip(a, b, strict=True)
        )
        if any(
            type(link[k]) not in (int, float) or not math.isfinite(link[k])
            for k in ("identity", "similarity", "score")
        ):
            raise ValueError("Invalid alignment number")
        if (
            link["score"] != alignment_score(a, b)
            or link["identity"] != matches / len(a)
            or link["similarity"] != (matches + similar) / len(a)
            or link["identity"] < identity
        ):
            raise ValueError("Alignment identity/similarity differs")
    return {
        "loci": len(loci),
        "genes": sum(len(x["genes"]) for x in loci),
        "links": len(seen),
        "evaluated_protein_pairs": len(allowed),
    }


def validate_locus_compare_pipeline(
    claimed, root, expected, integration, input_root, checkpoint=None
):
    from .validation_annotations import _payload, _strict_json, _fasta, _cas
    import json

    try:
        if (
            input_root is None
            or claimed.molecule != "DNA"
            or integration.get("workflow") != "locus_compare"
        ):
            return False
        records = {}
        for path in sorted(input_root.iterdir()):
            key, sequence = _fasta(_payload(input_root, path, checkpoint))
            if key in records:
                return False
            records[key] = sequence
        if (
            not 2 <= len(records) <= 8
            or any(len(v) > 50000 for v in records.values())
            or {k: (len(v), sha(v)) for k, v in records.items()} != expected
        ):
            return False
        directory_names = {f"record-{index:06d}" for index in range(1, len(records) + 1)}
        if ({path.name for path in (root / "casandra").iterdir()} != directory_names
                or {path.name for path in (root / "locus-compare-input/sources").iterdir()}
                != {name + ".fasta" for name in directory_names}):
            return False
        loci = []
        cas_counts = {"gene_count": 0, "cas_gene_count": 0, "cassette_count": 0}
        for index, (source, sequence) in enumerate(records.items(), 1):
            if checkpoint:
                checkpoint()
            dna = f">{source}\n{sequence}\n".encode()
            path = root / "locus-compare-input/sources" / f"record-{index:06d}.fasta"
            if _payload(root, path, checkpoint) != dna:
                return False
            directory = root / "casandra" / f"record-{index:06d}"
            counts = _cas(root, directory, dna, source, sequence, checkpoint)
            for key, value in counts.items():
                cas_counts[key] += value
            rows = [
                _strict_json(line)
                for line in _payload(
                    root, directory / "proteins.jsonl", checkpoint
                ).splitlines()
                if line.strip()
            ]
            genes = []
            for row in rows:
                if row["contig_id"] != source:
                    return False
                genes.append(
                    {
                        "id": row["protein_id"],
                        "start": row["start_1based"],
                        "end": row["end_1based_inclusive"],
                        "strand": row["strand"],
                        "phase": row["phase"],
                        "translation_table": row["translation_table"],
                        "partial_5prime": row["partial_5prime"],
                        "partial_3prime": row["partial_3prime"],
                        "translation": row["protein_sequence"],
                    }
                )
            genes.sort(key=lambda g: (g["start"], g["end"], g["id"]))
            loci.append(
                {
                    "id": f"locus{index:03d}",
                    "source_id": source,
                    "start": 1,
                    "end": len(sequence),
                    "genes": genes,
                }
            )
        features = {
            "schema": {"name": "clinker.features", "version": "1.0.0"},
            "coordinate_system": "1-based-inclusive-source-forward",
            "loci": loci,
        }
        serialized = (
            json.dumps(features, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode()
        native = b"".join(f">{k}\n{v}\n".encode() for k, v in records.items())
        if (
            _payload(root, root / "locus-compare-input/features.json", checkpoint)
            != serialized
            or _payload(root, root / "locus-compare-input/input.fasta", checkpoint)
            != native
        ):
            return False
        report = _strict_json(
            _payload(root, root / "locus-compare/comparison.json", checkpoint)
        )
        identity = claimed.tool_options.get("clinker", {}).get("--identity", 0.3)
        counts = validate_report(
            report, features, records, sha(native), sha(serialized), identity
        )
        stage = integration["stages"]["clinker"]
        return (
            stage["status"] == "completed"
            and stage["counts"] == counts
            and integration["stages"]["casandra"]
            == {"status": "completed", **cas_counts}
        )
    except (
        OSError,
        ValueError,
        TypeError,
        KeyError,
        AttributeError,
        IndexError,
        UnicodeError,
        RecursionError,
    ):
        return False


# BLOSUM62 standard published substitution scores (NCBI distribution).
MATRIX_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZX*"
MATRIX_ROWS = (
    (
        4,
        -1,
        -2,
        -2,
        0,
        -1,
        -1,
        0,
        -2,
        -1,
        -1,
        -1,
        -1,
        -2,
        -1,
        1,
        0,
        -3,
        -2,
        0,
        -2,
        -1,
        0,
        -4,
    ),
    (
        -1,
        5,
        0,
        -2,
        -3,
        1,
        0,
        -2,
        0,
        -3,
        -2,
        2,
        -1,
        -3,
        -2,
        -1,
        -1,
        -3,
        -2,
        -3,
        -1,
        0,
        -1,
        -4,
    ),
    (
        -2,
        0,
        6,
        1,
        -3,
        0,
        0,
        0,
        1,
        -3,
        -3,
        0,
        -2,
        -3,
        -2,
        1,
        0,
        -4,
        -2,
        -3,
        3,
        0,
        -1,
        -4,
    ),
    (
        -2,
        -2,
        1,
        6,
        -3,
        0,
        2,
        -1,
        -1,
        -3,
        -4,
        -1,
        -3,
        -3,
        -1,
        0,
        -1,
        -4,
        -3,
        -3,
        4,
        1,
        -1,
        -4,
    ),
    (
        0,
        -3,
        -3,
        -3,
        9,
        -3,
        -4,
        -3,
        -3,
        -1,
        -1,
        -3,
        -1,
        -2,
        -3,
        -1,
        -1,
        -2,
        -2,
        -1,
        -3,
        -3,
        -2,
        -4,
    ),
    (
        -1,
        1,
        0,
        0,
        -3,
        5,
        2,
        -2,
        0,
        -3,
        -2,
        1,
        0,
        -3,
        -1,
        0,
        -1,
        -2,
        -1,
        -2,
        0,
        3,
        -1,
        -4,
    ),
    (
        -1,
        0,
        0,
        2,
        -4,
        2,
        5,
        -2,
        0,
        -3,
        -3,
        1,
        -2,
        -3,
        -1,
        0,
        -1,
        -3,
        -2,
        -2,
        1,
        4,
        -1,
        -4,
    ),
    (
        0,
        -2,
        0,
        -1,
        -3,
        -2,
        -2,
        6,
        -2,
        -4,
        -4,
        -2,
        -3,
        -3,
        -2,
        0,
        -2,
        -2,
        -3,
        -3,
        -1,
        -2,
        -1,
        -4,
    ),
    (
        -2,
        0,
        1,
        -1,
        -3,
        0,
        0,
        -2,
        8,
        -3,
        -3,
        -1,
        -2,
        -1,
        -2,
        -1,
        -2,
        -2,
        2,
        -3,
        0,
        0,
        -1,
        -4,
    ),
    (
        -1,
        -3,
        -3,
        -3,
        -1,
        -3,
        -3,
        -4,
        -3,
        4,
        2,
        -3,
        1,
        0,
        -3,
        -2,
        -1,
        -3,
        -1,
        3,
        -3,
        -3,
        -1,
        -4,
    ),
    (
        -1,
        -2,
        -3,
        -4,
        -1,
        -2,
        -3,
        -4,
        -3,
        2,
        4,
        -2,
        2,
        0,
        -3,
        -2,
        -1,
        -2,
        -1,
        1,
        -4,
        -3,
        -1,
        -4,
    ),
    (
        -1,
        2,
        0,
        -1,
        -3,
        1,
        1,
        -2,
        -1,
        -3,
        -2,
        5,
        -1,
        -3,
        -1,
        0,
        -1,
        -3,
        -2,
        -2,
        0,
        1,
        -1,
        -4,
    ),
    (
        -1,
        -1,
        -2,
        -3,
        -1,
        0,
        -2,
        -3,
        -2,
        1,
        2,
        -1,
        5,
        0,
        -2,
        -1,
        -1,
        -1,
        -1,
        1,
        -3,
        -1,
        -1,
        -4,
    ),
    (
        -2,
        -3,
        -3,
        -3,
        -2,
        -3,
        -3,
        -3,
        -1,
        0,
        0,
        -3,
        0,
        6,
        -4,
        -2,
        -2,
        1,
        3,
        -1,
        -3,
        -3,
        -1,
        -4,
    ),
    (
        -1,
        -2,
        -2,
        -1,
        -3,
        -1,
        -1,
        -2,
        -2,
        -3,
        -3,
        -1,
        -2,
        -4,
        7,
        -1,
        -1,
        -4,
        -3,
        -2,
        -2,
        -1,
        -2,
        -4,
    ),
    (
        1,
        -1,
        1,
        0,
        -1,
        0,
        0,
        0,
        -1,
        -2,
        -2,
        0,
        -1,
        -2,
        -1,
        4,
        1,
        -3,
        -2,
        -2,
        0,
        0,
        0,
        -4,
    ),
    (
        0,
        -1,
        0,
        -1,
        -1,
        -1,
        -1,
        -2,
        -2,
        -1,
        -1,
        -1,
        -1,
        -2,
        -1,
        1,
        5,
        -2,
        -2,
        0,
        -1,
        -1,
        0,
        -4,
    ),
    (
        -3,
        -3,
        -4,
        -4,
        -2,
        -2,
        -3,
        -2,
        -2,
        -3,
        -2,
        -3,
        -1,
        1,
        -4,
        -3,
        -2,
        11,
        2,
        -3,
        -4,
        -3,
        -2,
        -4,
    ),
    (
        -2,
        -2,
        -2,
        -3,
        -2,
        -1,
        -2,
        -3,
        2,
        -1,
        -1,
        -2,
        -1,
        3,
        -3,
        -2,
        -2,
        2,
        7,
        -1,
        -3,
        -2,
        -1,
        -4,
    ),
    (
        0,
        -3,
        -3,
        -3,
        -1,
        -2,
        -2,
        -3,
        -3,
        3,
        1,
        -2,
        1,
        -1,
        -2,
        -2,
        0,
        -3,
        -1,
        4,
        -3,
        -2,
        -1,
        -4,
    ),
    (
        -2,
        -1,
        3,
        4,
        -3,
        0,
        1,
        -1,
        0,
        -3,
        -4,
        0,
        -3,
        -3,
        -2,
        0,
        -1,
        -4,
        -3,
        -3,
        4,
        1,
        -1,
        -4,
    ),
    (
        -1,
        0,
        0,
        1,
        -3,
        3,
        4,
        -2,
        0,
        -3,
        -3,
        1,
        -1,
        -3,
        -1,
        0,
        -1,
        -3,
        -2,
        -2,
        1,
        4,
        -1,
        -4,
    ),
    (
        0,
        -1,
        -1,
        -1,
        -2,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -1,
        -2,
        0,
        0,
        -2,
        -1,
        -1,
        -1,
        -1,
        -1,
        -4,
    ),
    (
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        -4,
        1,
    ),
)


def alignment_score(one, two):
    total = 0.0
    previous_gap = None
    for a, b in zip(one, two, strict=True):
        gap = "query" if a == "-" else "target" if b == "-" else None
        if gap:
            total += -0.5 if gap == previous_gap else -10.0
        elif a in MATRIX_ALPHABET and b in MATRIX_ALPHABET:
            total += MATRIX_ROWS[MATRIX_ALPHABET.index(a)][MATRIX_ALPHABET.index(b)]
        previous_gap = gap
    return total
