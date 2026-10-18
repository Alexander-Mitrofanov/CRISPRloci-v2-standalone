"""Offline HMM/system-rule execution with original-source coordinates.

Modified 2026-10-12 for CRISPRloci. GPL-3.0; see LICENSE.
No model discovery, downloads, shell invocation, embeddings or model updates.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from importlib.metadata import version
import logging
import math
import os
from pathlib import Path
import re
import subprocess

import pyrodigal
from macsypy.scripts import macsyfinder

from . import __version__ as VERSION
MODEL = {"id": "defense-finder-models-1.2.3", "version": "1.2.3",
         "commit": "5bdb41257fc7234fb8c3b70e3d362c52541860dd", "license": "GPL-3.0",
         "scope": "DefenseFinder_1-5_and_RM;CasFinder_and_newer_models_not_assessed"}
POLICY = {"coordinates": "one_based_inclusive", "topology": "linear",
          "gene_caller": "pyrodigal-3.7.1-meta", "coverage_profile": 0.4,
          "scope": "submitted_sequences_only", "functional_activity": "not_assessed"}
GROUPS = [f"DefenseFinder_{n}" for n in range(1, 6)] + ["RM"]
LIMIT = 20_000_000
INPUT_LIMIT = 4_000_000
MAX_RECORDS = 20
MAX_BASES = 2_000_000
DNA = frozenset("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def sha(value):
    return hashlib.sha256(value.encode("ascii") if isinstance(value, str) else value).hexdigest()


def bounded_bytes(path: Path, limit: int, label: str) -> bytes:
    """Read one regular file with a byte bound, including growth after stat()."""
    if path.is_symlink() or not path.is_file() or path.stat().st_size > limit:
        raise ValueError(f"Invalid or oversized {label}")
    with path.open("rb") as stream:
        data = stream.read(limit + 1)
    if len(data) > limit:
        raise ValueError(f"Invalid or oversized {label}")
    return data


def parse_fasta(data: bytes) -> dict[str, str]:
    """Join sequence chunks once; source order defines stable native replicons."""
    chunks: dict[str, list[str]] = {}
    key = None
    bases = 0
    for line in data.decode("ascii").splitlines():
        if line.startswith(">"):
            header = line[1:].split()
            if not header or not re.fullmatch(r"[A-Za-z0-9_.:-]{1,100}", header[0]):
                raise ValueError("Invalid source ID")
            key = header[0]
            if key in chunks or len(chunks) >= MAX_RECORDS:
                raise ValueError("Duplicate source ID or source count limit exceeded")
            chunks[key] = []
        elif line.strip():
            if key is None:
                raise ValueError("Missing FASTA header")
            sequence = line.strip().upper()
            bases += len(sequence)
            if bases > MAX_BASES or set(sequence) - DNA:
                raise ValueError("Expected DNA records and at most 2 Mb")
            chunks[key].append(sequence)
    if not chunks or any(not sequence for sequence in chunks.values()):
        raise ValueError("Expected 1-20 nonempty DNA records")
    return {key: "".join(sequence) for key, sequence in chunks.items()}


def read_fasta(path: Path) -> dict[str, str]:
    return parse_fasta(bounded_bytes(path, INPUT_LIMIT, "FASTA"))


def verify_models(directory):
    """Require every model byte to match the source-controlled release manifest."""
    manifest = Path(__file__).with_name("models.lock.json").read_bytes()
    lock = json.loads(manifest)
    package = directory / "defense-finder-models"
    for name, digest in lock["files"].items():
        path = package / name
        if path.is_symlink() or not path.is_file() or sha(path.read_bytes()) != digest:
            raise ValueError("DefenseFinder model bundle failed pinned identity check")
    actual = {p.relative_to(package).as_posix() for part in ("profiles", "definitions")
              for p in (package / part).rglob("*") if p.is_file()}
    wanted = {p for p in lock["files"] if p.startswith(("profiles/", "definitions/"))}
    if actual != wanted:
        raise ValueError("Unexpected model files")
    return {**MODEL, "manifest_sha256": sha(manifest)}


def call_genes(records):
    if pyrodigal.__version__ != "3.7.1":
        raise ValueError("Pyrodigal release mismatch")
    caller = pyrodigal.GeneFinder(meta=True)
    genes = []
    for source_number, (source_id, sequence) in enumerate(records.items(), 1):
        native = sorted(caller.find_genes(sequence.encode("ascii")), key=lambda g: (g.begin, g.end, g.strand))
        for position, gene in enumerate(native, 1):
            protein = gene.translate(include_stop=False)
            # Pyrodigal's gene.sequence() replaces IUPAC ambiguity with N.
            # Bind CDS evidence to the exact oriented submitted sequence instead.
            nucleotides = sequence[gene.begin - 1:gene.end]
            if gene.strand == -1:
                nucleotides = nucleotides.translate(COMPLEMENT)[::-1]
            genes.append({"id": f"r{source_number:06d}_{position:06d}", "source_id": source_id,
                          "replicon": f"r{source_number:06d}", "position": position,
                          "start": gene.begin, "end": gene.end, "strand": gene.strand,
                          "translation_table": gene.translation_table,
                          "partial_begin": bool(gene.partial_begin), "partial_end": bool(gene.partial_end),
                          "initiator_methionine": gene.start_type != "Edge",
                          "cds_sha256": sha(nucleotides), "protein": protein, "protein_sha256": sha(protein)})
            if len(genes) > 10_000:
                raise ValueError("Protein count limit exceeded")
    return genes


def native_lines(path: Path):
    """Stream bounded native tables, including headerless HMM extracts."""
    if path.is_symlink() or not path.is_file() or path.stat().st_size > LIMIT:
        raise ValueError("Invalid native table")
    with path.open(encoding="ascii", newline="") as stream:
        size = 0
        while line := stream.readline(LIMIT - size + 1):
            size += len(line)
            if size > LIMIT:
                raise ValueError("Native table byte budget exceeded")
            if not line.startswith("#"):
                yield line


def table(path: Path):
    yield from csv.DictReader(native_lines(path), delimiter="\t")


def number(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError("Nonfinite native score")
    return value


def bound_gene(lookup, ranks, identity, replicon, position, length):
    gene = lookup.get(identity)
    if (gene is None or replicon != gene["replicon"] or int(position) != ranks[identity]
            or int(length) != len(gene["protein"])):
        raise ValueError("Native gene identity mismatch")
    return gene


def profile_hit(row, lookup, ranks):
    if len(row) != 11:
        raise ValueError("Malformed HMM hit")
    gene = bound_gene(lookup, ranks, *row[:4])
    hit = {"gene_id": row[0], "profile": row[4], "i_evalue": number(row, 5),
           "score": number(row, 6), "profile_coverage": number(row, 7),
           "sequence_coverage": number(row, 8), "protein_start": int(row[9]), "protein_end": int(row[10])}
    if (not hit["profile"] or hit["i_evalue"] < 0
            or not 1 <= hit["protein_start"] <= hit["protein_end"] <= len(gene["protein"])):
        raise ValueError("Invalid HMM hit evidence")
    return hit


def collect(work: Path, genes):
    """Keep all system assignments and the first highest-scoring hit per profile.

    Systems can overlap or share a protein. Membership is distinct within each
    system; equal-score profile hits keep the original family/path traversal tie.
    """
    lookup = {gene["id"]: gene for gene in genes}
    if len(lookup) != len(genes):
        raise ValueError("Duplicate native gene ID")
    # MacSyFinder positions are global FASTA ranks; exported positions stay local.
    ranks = {gene["id"]: rank for rank, gene in enumerate(genes, 1)}
    components = []
    hit_map: dict[tuple[str, str], dict] = {}
    system_map = {}
    members: dict[str, set[str]] = {}
    for family in GROUPS:
        directory = work / family
        for row in table(directory / "best_solution.tsv"):
            gene = bound_gene(lookup, ranks, row["hit_id"], row["replicon"], row["hit_pos"], row["hit_seq_len"])
            system_id, model = row["sys_id"], row["model_fqn"]
            parts = model.split("/")
            if not system_id or len(parts) < 4 or parts[:2] != ["defense-finder-models", family]:
                raise ValueError("Native system model mismatch")
            wholeness, score = number(row, "sys_wholeness"), number(row, "sys_score")
            if system_id not in system_map:
                system_map[system_id] = {"id": system_id, "source_id": gene["source_id"],
                    "type": parts[-2], "subtype": parts[-1], "model": model,
                    "gene_ids": [], "wholeness": wholeness, "score": score}
                members[system_id] = set()
            system = system_map[system_id]
            identity = system["source_id"], system["model"], system["wholeness"], system["score"]
            if identity != (gene["source_id"], model, wholeness, score):
                raise ValueError("Inconsistent native system evidence")
            if gene["id"] not in members[system_id]:
                system["gene_ids"].append(gene["id"])
                members[system_id].add(gene["id"])
            components.append({"system_id": system_id, "gene_id": gene["id"], "profile": row["gene_name"],
                "status": row["hit_status"], "score": number(row, "hit_score"),
                "i_evalue": number(row, "hit_i_eval"), "profile_coverage": number(row, "hit_profile_cov")})
        for path in sorted((directory / "hmmer_results").glob("*.res_hmm_extract")):
            for row in csv.reader(native_lines(path), delimiter="\t"):
                if not row:
                    continue
                hit = profile_hit(row, lookup, ranks)
                key = hit["gene_id"], hit["profile"]
                if key not in hit_map or hit["score"] > hit_map[key]["score"]:
                    hit_map[key] = hit
    systems = list(system_map.values())
    for system in systems:
        system["gene_ids"].sort(key=lambda key: lookup[key]["position"])
        system["start"] = min(lookup[key]["start"] for key in system["gene_ids"])
        system["end"] = max(lookup[key]["end"] for key in system["gene_ids"])
    systems.sort(key=lambda s: (s["source_id"], s["start"], s["id"]))
    components.sort(key=lambda c: (c["system_id"], c["gene_id"], c["profile"]))
    return systems, components, [hit_map[key] for key in sorted(hit_map)]


def reject_ambient_configuration(models_dir):
    """Reject undeclared files consumed by pinned MacSyFinder 2.1.4.

    Do not read their contents or change the caller's environment. VIRTUAL_ENV
    takes precedence over MACSY_CONF in MacSyFinder's own resolver.
    """
    virtual = os.environ.get("VIRTUAL_ENV")
    system = Path(virtual) / "etc/macsyfinder/macsyfinder.conf" if virtual else Path(os.environ.get("MACSY_CONF") or "/etc/macsyfinder/macsyfinder.conf")
    candidates = [system, Path.home() / ".macsyfinder/macsyfinder.conf", Path.cwd() / "macsyfinder.conf"]
    candidates.extend(Path(models_dir) / "defense-finder-models" / family / "model_conf.xml" for family in GROUPS)
    if any(path.exists() or path.is_symlink() for path in candidates):
        raise ValueError("Undeclared MacSyFinder configuration is not supported by defense-context")


def run(input_path, output_path, models_dir, hmmer, workers):
    if type(workers) is not int or not 1 <= workers <= 8:
        raise ValueError("Worker count must be in 1..8")
    if version("macsyfinder") != "2.1.4":
        raise ValueError("MacSyFinder release mismatch")
    reject_ambient_configuration(models_dir)
    engine = subprocess.run([str(hmmer), "-h"], capture_output=True, text=True, timeout=10, check=True)
    if "HMMER 3.4 " not in engine.stdout:
        raise ValueError("HMMER release mismatch")
    input_data = bounded_bytes(input_path, INPUT_LIMIT, "FASTA")
    records = parse_fasta(input_data)
    model = verify_models(models_dir)
    output_path.mkdir(parents=True, exist_ok=False)
    work = output_path / "work"
    work.mkdir()
    genes = call_genes(records)
    protein_path = work / "proteins.faa"
    protein_path.write_text("".join(f'>{g["id"]}\n{g["protein"]}\n' for g in genes))
    systems, components, hits = [], [], []
    if genes:
        for family in GROUPS:
            try:
                status = macsyfinder.main(args=["--db-type", "gembase", "--replicon-topology", "linear",
                    "--sequence-db", str(protein_path), "--models-dir", str(models_dir),
                    "--models", f"defense-finder-models/{family}", "all", "--out-dir", str(work / family),
                    "--worker", str(workers), "--coverage-profile", "0.4", "--exchangeable-weight", "1",
                    "--hmmer", str(hmmer), "--mute"])
                if status not in (None, 0):
                    raise RuntimeError("MacSyFinder did not complete")
            finally:
                logger = logging.getLogger("macsypy")
                for handler in logger.handlers[:]:
                    handler.close()
                    logger.removeHandler(handler)
        systems, components, hits = collect(work, genes)
    doc = {"schema": {"name": "defensefinder.context", "version": "1.0.0"},
           "tool": {"name": "DefenseFinder-v2", "version": VERSION}, "execution_completed": True,
           "model": model, "policy": POLICY,
           "input": {"sha256": sha(input_data), "record_count": len(records)},
           "sources": [{"id": key, "length": len(seq), "sha256": sha(seq)} for key, seq in records.items()],
           "genes": genes, "systems": systems, "components": components, "hits": hits,
           "counts": {"sources": len(records), "genes": len(genes), "systems": len(systems),
                      "components": len(components), "hits": len(hits)}}
    data = json.dumps(doc, sort_keys=True, allow_nan=False) + "\n"
    if len(data.encode()) > LIMIT:
        raise ValueError("Report byte budget exceeded")
    (output_path / "report.json").write_text(data)
    return doc


def main():
    parser = argparse.ArgumentParser(description="Offline source-bound DefenseFinder HMM context")
    parser.add_argument("--version", action="version", version=f"DefenseFinder-v2 {VERSION}")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--models-dir", required=True, type=Path)
    parser.add_argument("--hmmer", required=True, type=Path)
    parser.add_argument("--workers", type=int, choices=range(1, 9), default=1)
    args = parser.parse_args()
    run(args.input.absolute(), args.output.absolute(), args.models_dir.absolute(), args.hmmer.absolute(), args.workers)
