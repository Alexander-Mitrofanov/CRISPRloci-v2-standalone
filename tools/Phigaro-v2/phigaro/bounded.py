"""Bounded offline worker entry point; original Phigaro scoring, explicit evidence."""

from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import io
import json
import math
import os
import re
import signal
import subprocess
import sys
import threading
from collections import defaultdict
from pathlib import Path

from Bio.Data import CodonTable
from Bio.Seq import Seq

from . import const
from .finder.v2 import calc_scores, scan_phages

VERSION = "2.4.1"
PARAMETERS = {
    "window": 32,
    "threshold_min": 45.39,
    "threshold_max": 46.0,
    "max_evalue": 0.00445,
    "penalty_black": 2.2,
    "penalty_white": 0.7,
    "mode": "basic",
    "gene_mode": "single",
    "minimum_record_bases": 20000,
    "gc_source": "prodigal_reported_gene_gc",
    "protein_translation": "prodigal_strict_ambiguous_X",
}
BLACK_PROFILES = frozenset(const.DEFAULT_BLACK_LIST)
WHITE_PROFILES = frozenset(const.DEFAULT_WHITE_LIST)
CODE = dict(CodonTable.unambiguous_dna_by_id[11].forward_table)
CODE.update({codon: "*" for codon in CodonTable.unambiguous_dna_by_id[11].stop_codons})
INTERPRETATION = (
    "Candidate prophage intervals in assessed sequences; boundaries are gene-based. "
    "Absence of a prediction does not establish absence of prophages. "
    "Neither inducibility nor active CRISPR targeting is established."
)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def file_digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def read_fasta_snapshot(path):
    path = Path(path)
    if (
        path.is_symlink()
        or not path.is_file()
        or not 0 < path.stat().st_size <= 5_000_000
    ):
        raise ValueError("Input must be a bounded regular FASTA file")
    chunks = {}
    ident = None
    bases = 0
    with path.open("rb") as stream:
        payload = stream.read(5_000_001)
    if not payload or len(payload) > 5_000_000:
        raise ValueError("Input must be a bounded regular FASTA file")
    for line in io.StringIO(payload.decode("ascii"), newline=None):
        if line.startswith(">"):
            tokens = line[1:].split()
            ident = tokens[0] if tokens else ""
            if (
                not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:-]{0,199}", ident)
                or ident in chunks
            ):
                raise ValueError("Invalid or duplicate FASTA identifier")
            chunks[ident] = []
            if len(chunks) > 20:
                raise ValueError(
                    "Input is limited to 20 records and 2 million bases"
                )
        elif line.strip():
            if ident is None:
                raise ValueError("Sequence appears before a FASTA header")
            sequence = line.strip().upper()
            chunks[ident].append(sequence)
            bases += len(sequence)
            if bases > 2_000_000:
                raise ValueError(
                    "Input is limited to 20 records and 2 million bases"
                )
    if not chunks:
        raise ValueError("Input is limited to 20 records and 2 million bases")
    records = {key: "".join(parts) for key, parts in chunks.items()}
    if any(not seq or set(seq) - set("ACGTRYSWKMBDHVN") for seq in records.values()):
        raise ValueError("Empty or invalid DNA sequence")
    return records, sha256(payload)


def read_fasta(path):
    return read_fasta_snapshot(path)[0]


def _stop_native(signum, frame):
    # Raising lets the command's finally block reap its isolated group.
    raise SystemExit(128 + signum)


def _run_native(command, *, timeout, check=True, capture_output=False, **kwargs):
    """Contain the native command and every descendant on all exit paths."""
    if capture_output:
        kwargs.update(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    child = subprocess.Popen(command, start_new_session=True, **kwargs)
    previous = None
    if threading.current_thread() is threading.main_thread():
        previous = signal.signal(signal.SIGTERM, _stop_native)
    try:
        stdout, stderr = child.communicate(timeout=timeout)
    finally:
        with contextlib.suppress(ProcessLookupError):
            os.killpg(child.pid, signal.SIGKILL)
        child.wait()
        for stream in (child.stdin, child.stdout, child.stderr):
            if stream is not None:
                stream.close()
        if previous is not None:
            signal.signal(signal.SIGTERM, previous)
    if check and child.returncode:
        raise subprocess.CalledProcessError(child.returncode, command, stdout, stderr)
    return subprocess.CompletedProcess(command, child.returncode, stdout, stderr)


def process(command, work, name, timeout):
    with (
        (work / f"{name}.stdout").open("wb") as out,
        (work / f"{name}.stderr").open("wb") as err,
    ):
        _run_native(command, stdout=out, stderr=err, check=True, timeout=timeout)


def parse_genes(path, sequences, native_ids):
    genes = []
    seen = set()
    header, protein = None, []

    def append():
        if header is None:
            return
        parts = header.split(" # ")
        if len(parts) != 5:
            raise ValueError("Malformed Prodigal protein header")
        native_id, start, end, strand, metadata = parts
        if len(genes) >= 10000 or native_id in seen:
            raise ValueError("Too many or duplicate predicted genes")
        key = native_id.rsplit("_", 1)[0]
        source_id = native_ids.get(key)
        if source_id is None:
            raise ValueError("Prodigal returned an unknown source")
        start, end, strand = int(start), int(end), int(strand)
        seq = sequences[source_id]
        if not 1 <= start <= end <= len(seq) or strand not in (-1, 1):
            raise ValueError("Prodigal coordinates are outside source DNA")
        fields = dict(item.split("=", 1) for item in metadata.split(";") if "=" in item)
        partial = fields.get("partial")
        if partial not in {"00", "01", "10", "11"}:
            raise ValueError("Prodigal partial-gene evidence is unavailable")
        dna = seq[start - 1 : end]
        oriented = dna if strand == 1 else str(Seq(dna).reverse_complement())
        if len(oriented) % 3:
            raise ValueError("Prodigal CDS length is not a multiple of three")
        predicted = "".join(protein).removesuffix("*")
        # Native Prodigal emits X for every non-ACGT codon, including GCN/AAR.
        translated = "".join(
            CODE.get(oriented[i : i + 3], "X") for i in range(0, len(oriented), 3)
        ).removesuffix("*")
        # Prodigal uses methionine for complete alternative initiation codons.
        partial_5prime = partial[0 if strand == 1 else 1] == "1"
        partial_3prime = partial[1 if strand == 1 else 0] == "1"
        if not partial_5prime and oriented[:3] not in {
            "TTG",
            "CTG",
            "ATT",
            "ATC",
            "ATA",
            "ATG",
            "GTG",
        }:
            raise ValueError("Complete CDS has no supported initiation codon")
        if not partial_3prime and oriented[-3:] not in {"TAA", "TAG", "TGA"}:
            raise ValueError("Complete CDS has no terminal stop codon")
        if not partial_5prime and translated:
            translated = "M" + translated[1:]
        if not predicted or predicted != translated or "*" in predicted:
            raise ValueError("Protein is inconsistent with its source CDS")
        gc = float(fields["gc_cont"])
        if not math.isfinite(gc) or not 0 <= gc <= 1:
            raise ValueError("Invalid gene GC fraction")
        observed_gc = (dna.count("G") + dna.count("C")) / len(dna)
        # Preserve the native GC covariate used to calibrate Phigaro thresholds.
        # Prodigal 2.6.3 can differ from exact CDS GC on reverse-strand genes.
        seen.add(native_id)
        genes.append(
            {
                "id": native_id,
                "source_id": source_id,
                "start": start,
                "end": end,
                "strand": "+" if strand == 1 else "-",
                "translation_table": 11,
                "partial_5prime": partial_5prime,
                "partial_3prime": partial_3prime,
                "sequence_sha256": sha256(dna.encode()),
                "protein_sha256": sha256(predicted.encode()),
                "translation": predicted,
                "gc": gc,
                "source_gc": observed_gc,
                "best_hit": None,
            }
        )

    with Path(path).open(encoding="ascii") as stream:
        for raw in stream:
            line = raw.strip()
            if line.startswith(">"):
                append()
                header, protein = line[1:], []
            elif line:
                if header is None:
                    raise ValueError("Protein appears before a Prodigal header")
                protein.append(line)
        append()
    return genes


def read_hits(path, genes, profiles):
    by_id = {gene["id"]: gene for gene in genes}
    if Path(path).stat().st_size > 50_000_000:
        raise ValueError("HMMER evidence exceeds byte limit")
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(maxsplit=18)
            if len(fields) < 18 or fields[0] not in by_id or fields[2] not in profiles:
                raise ValueError("Unknown or malformed HMMER evidence")
            evalue, score = float(fields[4]), float(fields[5])
            if not math.isfinite(evalue) or evalue < 0 or not math.isfinite(score):
                raise ValueError("Nonfinite HMMER evidence")
            hit = {"profile": fields[2], "evalue": evalue, "score": score}
            gene = by_id[fields[0]]
            previous = gene["best_hit"]
            if previous is None or (evalue, -score, fields[2]) < (
                previous["evalue"],
                -previous["score"],
                previous["profile"],
            ):
                gene["best_hit"] = hit


def supported_hit(gene):
    hit = gene["best_hit"]
    return hit is not None and hit["evalue"] <= PARAMETERS["max_evalue"]


def profile_marker(profile):
    if profile in BLACK_PROFILES:
        return 1 - PARAMETERS["penalty_black"]
    if profile in WHITE_PROFILES:
        return 1 + PARAMETERS["penalty_white"]
    return 1.0


def regions_for(sequences, genes):
    """Score each source independently; retain inclusive terminal genes."""
    by_source = defaultdict(list)
    for gene in genes:
        by_source[gene["source_id"]].append(gene)
    regions = []
    for source_id, seq in sequences.items():
        group = sorted(
            by_source[source_id], key=lambda g: (g["start"], g["end"], g["id"])
        )
        supported = [supported_hit(gene) for gene in group]
        marker = [
            profile_marker(gene["best_hit"]["profile"]) if hit else 0.0
            for gene, hit in zip(group, supported)
        ]
        gc = [gene["gc"] if hit else 0.0 for gene, hit in zip(group, supported)]
        window = PARAMETERS["window"]
        scores = [
            a * b for a, b in zip(calc_scores(marker, window), calc_scores(gc, window))
        ]
        for begin, end in scan_phages(
            scores, PARAMETERS["threshold_min"], PARAMETERS["threshold_max"]
        ):
            selected = group[begin : end + 1]
            start, stop = selected[0]["start"], max(g["end"] for g in selected)
            sequence = seq[start - 1 : stop]
            regions.append(
                {
                    "id": f"region{len(regions) + 1:06d}",
                    "source_id": source_id,
                    "start": start,
                    "end": stop,
                    "sequence": sequence,
                    "sequence_sha256": sha256(sequence.encode()),
                    "gene_ids": [g["id"] for g in selected],
                    "supporting_gene_ids": [
                        g["id"]
                        for g, hit in zip(selected, supported[begin : end + 1])
                        if hit
                    ],
                    "first_gene_at_record_edge": begin == 0,
                    "last_gene_at_record_edge": end == len(group) - 1,
                    "peak_score": max(scores[begin : end + 1]),
                }
            )
    return regions


def predict(args):
    records, input_digest = read_fasta_snapshot(args.input)
    models = Path(args.models_hmm)
    if (
        models.is_symlink()
        or not models.is_file()
        or not 0 < models.stat().st_size <= 2_000_000_000
    ):
        raise ValueError("Model bundle must be an operator-owned bounded regular file")
    if (
        not re.fullmatch(r"[0-9a-f]{64}", args.models_sha256)
        or file_digest(models) != args.models_sha256
    ):
        raise ValueError("Model bundle differs from its pinned SHA-256")
    profiles = set()
    with models.open(encoding="ascii") as stream:
        for line in stream:
            if line.startswith("NAME "):
                name = line.split()[1]
                if not re.fullmatch(r"VOG\d{4,6}[a-z]?", name) or name in profiles:
                    raise ValueError("Unexpected or duplicate pVOG profile")
                profiles.add(name)
    if not 1 <= len(profiles) <= 20000:
        raise ValueError("Invalid number of pVOG profiles")
    output = Path(args.output)
    if output.exists() and (
        output.is_symlink() or not output.is_dir() or any(output.iterdir())
    ):
        raise ValueError("Output directory must be new or empty")
    output.mkdir(parents=True, exist_ok=True)
    work = output / "work"
    work.mkdir()
    eligible = {
        key: seq
        for key, seq in records.items()
        if len(seq) >= PARAMETERS["minimum_record_bases"]
    }
    native_ids = {
        f"record{i:06d}": source_id for i, source_id in enumerate(eligible, 1)
    }
    genes = []
    if eligible:
        fasta = work / "input.fasta"
        fasta.write_text(
            "".join(
                f">{key}\n{eligible[source_id]}\n"
                for key, source_id in native_ids.items()
            )
        )
        proteins = work / "proteins.faa"
        process(
            [
                args.prodigal,
                "-i",
                str(fasta),
                "-a",
                str(proteins),
                "-p",
                "single",
                "-g",
                "11",
                "-q",
            ],
            work,
            "prodigal",
            args.timeout_seconds,
        )
        genes = parse_genes(proteins, records, native_ids)
        if genes:
            hits = work / "hits.tsv"
            process(
                [
                    args.hmmsearch,
                    "--cpu",
                    str(args.threads),
                    "--notextw",
                    "--noali",
                    "--tblout",
                    str(hits),
                    str(models),
                    str(proteins),
                ],
                work,
                "hmmsearch",
                args.timeout_seconds,
            )
            read_hits(hits, genes, profiles)
    report = {
        "schema": {"name": "phigaro.regions", "version": "1.0.0"},
        "tool": {"name": "Phigaro-v2", "version": VERSION},
        "coordinate_system": "1-based-inclusive-source-forward",
        "input_sha256": input_digest,
        "reference": {
            "name": "pVOGs",
            "sha256": args.models_sha256,
            "profile_count": len(profiles),
        },
        "parameters": dict(PARAMETERS),
        "sources": [
            {
                "id": key,
                "length": len(seq),
                "sha256": sha256(seq.encode()),
                "status": "assessed"
                if key in eligible
                else "ineligible_short_sequence",
            }
            for key, seq in records.items()
        ],
        "genes": genes,
        "regions": regions_for(records, genes),
        "interpretation": INTERPRETATION,
    }
    (output / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    with (output / "regions.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(
            ["id", "source_id", "start", "end", "gene_count", "supporting_gene_count"]
        )
        for region in report["regions"]:
            writer.writerow(
                [region[k] for k in ["id", "source_id", "start", "end"]]
                + [len(region["gene_ids"]), len(region["supporting_gene_ids"])]
            )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=VERSION)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("predict")
    command.add_argument("--input", required=True)
    command.add_argument("--output", required=True)
    command.add_argument("--models-hmm", required=True)
    command.add_argument("--models-sha256", required=True)
    command.add_argument("--prodigal", default="prodigal")
    command.add_argument("--hmmsearch", default="hmmsearch")
    command.add_argument("--threads", type=int, choices=range(1, 9), default=1)
    command.add_argument(
        "--timeout-seconds",
        type=int,
        choices=range(1, 14401),
        default=3600,
        metavar="SECONDS",
        help="Time limit in seconds (1-14400; default: 3600).",
    )
    args = parser.parse_args(argv)
    try:
        result = predict(args)
    except (ValueError, OSError, KeyError, subprocess.SubprocessError) as exc:
        print(f"Phigaro analysis failed: {type(exc).__name__}", file=sys.stderr)
        return 2
    print(
        json.dumps({"regions": len(result["regions"]), "genes": len(result["genes"])})
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
