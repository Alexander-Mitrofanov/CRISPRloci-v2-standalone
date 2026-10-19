"""Offline, source-verified comparison for bounded genomic neighborhoods.

The scientific alignment is clinker global BLOSUM62 with its identity/similarity
definition. This entry point does not call the HTML renderer or HTTP server.
"""
from __future__ import annotations

import argparse
import ctypes
import errno
import hashlib
import io
import json
import math
import os
import re
import shutil
import sys
import tempfile
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path

from Bio import Align, SeqIO
from Bio import __version__ as bio_version
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from clinker import __version__
from clinker.align import compute_identity, extend_matrix_alphabet

SCHEMA = {"name": "clinker.features", "version": "1.0.0"}
INTERPRETATION = "Protein similarity and observed gene order; not orthology, function, or a phylogeny."
LIMITS = {"loci": 8, "window_nt": 50000, "genes": 120, "protein_aa": 3000,
          "protein_pairs": 6000, "alignment_cells": 250000000}
ALIGNER = {"mode": "global", "substitution_matrix": "BLOSUM62",
           "open_gap_score": -10.0, "extend_gap_score": -0.5}
GENE_KEYS = {"id", "start", "end", "strand", "phase", "translation_table",
             "partial_5prime", "partial_3prime", "translation"}


def digest(value):
    return hashlib.sha256(value.encode() if isinstance(value, str) else value).hexdigest()


def read_bytes(path, maximum):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_size > maximum:
        raise ValueError("Input must be a bounded regular file")
    with path.open("rb") as handle:
        raw = handle.read(maximum + 1)
    if len(raw) > maximum:
        raise ValueError("Input exceeds size limit")
    return raw


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate JSON key")
        result[key] = value
    return result


def identifier(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9_.:|+%=\-]{1,512}", value):
        raise ValueError("Invalid evidence identifier")
    return value


def integer(value, minimum, maximum):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError("Integer is outside permitted range")
    return value


def validate_features(document, records):
    if (not isinstance(document, dict) or set(document) != {"schema", "coordinate_system", "loci"}
            or document["schema"] != SCHEMA
            or document["coordinate_system"] != "1-based-inclusive-source-forward"):
        raise ValueError("Invalid feature document schema")
    loci = document["loci"]
    if not isinstance(loci, list) or not 2 <= len(loci) <= LIMITS["loci"]:
        raise ValueError("Comparison requires 2–8 bounded loci")
    output = []
    seen = set()
    gene_total = 0
    for locus in loci:
        if not isinstance(locus, dict) or set(locus) != {"id", "source_id", "start", "end", "genes"}:
            raise ValueError("Invalid locus fields")
        lid = identifier(locus["id"])
        source = identifier(locus["source_id"])
        if lid in seen or source not in records:
            raise ValueError("Duplicate locus or absent source")
        seen.add(lid)
        sequence = records[source]
        start = integer(locus["start"], 1, len(sequence))
        end = integer(locus["end"], start, min(len(sequence), start + LIMITS["window_nt"] - 1))
        genes = locus["genes"]
        if not isinstance(genes, list) or len(genes) > LIMITS["genes"]:
            raise ValueError("Invalid gene collection")
        checked = []
        ids = set()
        previous = None
        for gene in genes:
            if not isinstance(gene, dict) or set(gene) != GENE_KEYS:
                raise ValueError("Invalid CDS fields")
            gid = identifier(gene["id"])
            if gid in ids:
                raise ValueError("Duplicate gene identity within locus")
            ids.add(gid)
            gstart = integer(gene["start"], start, end)
            gend = integer(gene["end"], gstart, end)
            key = (gstart, gend, gid)
            if previous is not None and key <= previous:
                raise ValueError("Genes must be in source-coordinate order")
            previous = key
            phase = integer(gene["phase"], 0, 2)
            table = integer(gene["translation_table"], 1, 33)
            if gene["strand"] not in ("+", "-") or any(type(gene[k]) is not bool for k in ("partial_5prime", "partial_3prime")):
                raise ValueError("Invalid strand or partial-CDS evidence")
            if phase and not gene["partial_5prime"]:
                raise ValueError("Nonzero phase requires an explicitly partial 5-prime CDS")
            nt = sequence[gstart - 1:gend]
            if gene["strand"] == "-":
                nt = str(Seq(nt).reverse_complement())
            coding = nt[phase:]
            if len(coding) % 3:
                if not gene["partial_3prime"]:
                    raise ValueError("Incomplete codon requires an explicitly partial 3-prime CDS")
                coding = coding[:len(coding) // 3 * 3]
            try:
                translated = str(Seq(coding).translate(table=table, cds=False)).removesuffix("*")
            except (ValueError, KeyError) as exc:
                raise ValueError("Unsupported genetic code or invalid coding DNA") from exc
            protein = gene["translation"]
            if (not isinstance(protein, str) or not 1 <= len(protein) <= LIMITS["protein_aa"]
                    or set(protein) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO") or protein != translated):
                raise ValueError("CDS translation differs from source DNA or is unsupported")
            checked.append({**gene, "nucleotide_sha256": digest(nt), "protein_sha256": digest(protein)})
        gene_total += len(checked)
        output.append({**locus, "genes": checked, "sequence_sha256": digest(sequence[start-1:end]),
                       "source_length": len(sequence), "at_source_start": start == 1,
                       "at_source_end": end == len(sequence)})
    if gene_total > LIMITS["genes"]:
        raise ValueError("Total CDS count exceeds bounded-comparison limit")
    pairs = cells = 0
    for one, two in combinations(output, 2):
        pairs += len(one["genes"]) * len(two["genes"])
        cells += (sum(len(gene["translation"]) for gene in one["genes"])
                  * sum(len(gene["translation"]) for gene in two["genes"]))
    if pairs > LIMITS["protein_pairs"] or cells > LIMITS["alignment_cells"]:
        raise ValueError("Protein comparison exceeds work limit; supply smaller neighborhoods")
    return output, pairs, cells


def compare(fasta_raw, features_raw, identity=0.3):
    if type(identity) not in (int, float) or not math.isfinite(identity) or not 0 <= identity <= 1:
        raise ValueError("Identity cutoff must be finite and in [0,1]")
    records = {}
    for record in SeqIO.parse(io.StringIO(fasta_raw.decode("utf-8")), "fasta"):
        rid = identifier(record.id)
        sequence = str(record.seq).upper()
        if rid in records or not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError("Duplicate, empty or invalid nucleotide source")
        records[rid] = sequence
    if not records or len(records) > 20 or sum(map(len, records.values())) > 2000000:
        raise ValueError("Source sequence limit exceeded")
    document = json.loads(features_raw, object_pairs_hook=unique_object,
                          parse_constant=lambda _: (_ for _ in ()).throw(ValueError("Nonfinite JSON value")))
    loci, pairs, cells = validate_features(document, records)
    aligner = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    aligner.substitution_matrix = extend_matrix_alphabet(substitution_matrices.load("BLOSUM62"))

    @lru_cache(maxsize=128)
    def align_proteins(query, target):
        # Ordered keys retain Biopython's first optimum and gap tie placement.
        # Per-run cache is bounded; evidence counts still cover every gene pair.
        alignment = aligner.align(query, target)[0]
        aligned = (alignment[0], alignment[1])
        ident, similar = compute_identity(aligned)
        return ident, similar, float(alignment.score), *aligned

    links = []
    for one, two in combinations(loci, 2):
        for a, b in product(one["genes"], two["genes"]):
            ident, similar, score, query, target = align_proteins(a["translation"], b["translation"])
            if ident < identity:
                continue
            links.append({"query_locus": one["id"], "query_gene": a["id"],
                          "target_locus": two["id"], "target_gene": b["id"],
                          "identity": ident, "similarity": similar,
                          "score": score, "aligned_query": query,
                          "aligned_target": target})
    report = {"schema": {"name": "clinker.comparison", "version": "1.0.0"},
              "tool": {"name": "clinker-v2", "version": __version__, "biopython": bio_version},
              "coordinate_system": "1-based-inclusive-source-forward",
              "input_sha256": digest(fasta_raw), "features_sha256": digest(features_raw),
              "options": {"identity": identity, "jobs": 1, "aligner": ALIGNER.copy()}, "limits": LIMITS.copy(),
              "sources": [{"id": k, "length": len(v), "sha256": digest(v)} for k, v in records.items()],
              "loci": loci, "links": links, "evaluated_protein_pairs": pairs,
              "alignment_cells": cells, "interpretation": INTERPRETATION}
    return report, records


def write_exports(directory, report, records):
    (directory / "comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    keys = ("query_locus", "query_gene", "target_locus", "target_gene", "identity", "similarity", "score", "aligned_query", "aligned_target")
    with (directory / "alignments.tsv").open("w") as handle:
        handle.write("\t".join(keys) + "\n")
        for link in report["links"]:
            cells = []
            for key in keys:
                value = link[key]
                cell = str(value)
                numeric = key in {"identity", "similarity", "score"} and type(value) in (int, float) and math.isfinite(value)
                residues = key in {"aligned_query", "aligned_target"} and re.fullmatch(r"[ABCDEFGHIKLMNPQRSTVWXYZJUO-]+", cell)
                if not (numeric or residues) and cell.lstrip().startswith(("=", "+", "-", "@")):
                    cell = "'" + cell
                cells.append(cell)
            handle.write("\t".join(cells) + "\n")
    gbk = []
    proteins = []
    for i, locus in enumerate(report["loci"], 1):
        record = SeqRecord(Seq(records[locus["source_id"]][locus["start"]-1:locus["end"]]),
                           id=f"locus{i:03d}", name=f"locus{i:03d}", description=locus["id"])
        record.annotations = {"molecule_type": "DNA", "topology": "linear",
                              "comment": f"Source {locus['source_id']}:{locus['start']}-{locus['end']}; exported feature coordinates are local to this window."}
        for j, gene in enumerate(locus["genes"], 1):
            record.features.append(SeqFeature(SimpleLocation(gene["start"]-locus["start"], gene["end"]-locus["start"]+1,
                                                             strand=1 if gene["strand"] == "+" else -1), type="CDS",
                qualifiers={"locus_tag": [gene["id"]], "translation": [gene["translation"]],
                            "transl_table": [str(gene["translation_table"])], "codon_start": [str(gene["phase"] + 1)],
                            "note": [f"partial_5prime={gene['partial_5prime']};partial_3prime={gene['partial_3prime']}"]}))
            proteins.append(SeqRecord(Seq(gene["translation"]), id=f"locus{i:03d}_gene{j:03d}",
                                      description=f"locus={locus['id']} gene={gene['id']} source={locus['source_id']}:{gene['start']}-{gene['end']} strand={gene['strand']}"))
        gbk.append(record)
    SeqIO.write(gbk, directory / "loci.gbk", "genbank")
    SeqIO.write(proteins, directory / "proteins.faa", "fasta")


def _rename_new_directory(source: Path, destination: Path) -> None:
    """Atomically publish without replacing a concurrently created directory."""
    if os.name == "nt":
        # Windows rename already refuses to replace an existing destination.
        source.rename(destination)
        return
    library = ctypes.CDLL(None, use_errno=True)
    if sys.platform.startswith("linux") and hasattr(library, "renameat2"):
        rename = library.renameat2
        rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
        rename.restype = ctypes.c_int
        # AT_FDCWD, RENAME_NOREPLACE: the existence check occurs inside rename.
        result = rename(-100, os.fsencode(source), -100, os.fsencode(destination), 1)
    elif sys.platform == "darwin" and hasattr(library, "renamex_np"):
        rename = library.renamex_np
        rename.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        rename.restype = ctypes.c_int
        result = rename(os.fsencode(source), os.fsencode(destination), 4)  # RENAME_EXCL
    else:
        raise ValueError("This platform does not support atomic directory publication without replacement")
    if result:
        error = ctypes.get_errno()
        if error in (errno.EEXIST, errno.ENOTEMPTY):
            raise ValueError("Output must be a new non-symlink directory")
        raise OSError(error, os.strerror(error), str(destination))



def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("compare")
    command.add_argument("--input", required=True, type=Path)
    command.add_argument("--features", required=True, type=Path)
    command.add_argument("--outdir", required=True, type=Path)
    command.add_argument("--identity", type=float, default=0.3)
    args = parser.parse_args(argv)
    staging = None
    try:
        if args.outdir.exists() or args.outdir.is_symlink():
            raise ValueError("Refusing to overwrite output directory")
        report, records = compare(read_bytes(args.input, 4000000), read_bytes(args.features, 4000000), args.identity)
        args.outdir.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=".clinker-", dir=args.outdir.parent))
        write_exports(staging, report, records)
        _rename_new_directory(staging, args.outdir)
        staging = None
    except (ValueError, OSError, UnicodeError) as exc:
        parser.exit(2, f"clinker-v2: {exc}\n")
    finally:
        if staging is not None:
            shutil.rmtree(staging)


if __name__ == "__main__":
    main()
