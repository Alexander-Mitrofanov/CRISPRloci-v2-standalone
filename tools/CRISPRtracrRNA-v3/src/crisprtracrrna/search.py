"""Sequence, covariance-model, interaction, and terminator tool adapters."""

from __future__ import annotations

import csv
import json
import math
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, cast

from .config import Settings
from .errors import ContractError
from .models import AntiRepeat, Array, Candidate, Hit, Interval, Sequence, Strand
from .process import Runner
from .sequence import reverse_complement, subsequence, write_fasta
from .short_rescue_policy import SHORT_RESCUE_POLICY


def finite(value: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Non-finite scientific result")
    return number


def parse_cm_table(path: Path, lengths: dict[str, int], model: str) -> list[Hit]:
    """Parse Infernal 1.x cmscan --fmt 1; seqfrom/seqto encode strand direction."""
    hits = []
    try:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.split(maxsplit=17)
                if len(cols) < 17:
                    raise ValueError("Truncated cmscan table row")
                query, strand = cols[2], cols[9]
                start, end = int(cols[7]), int(cols[8])
                if strand not in {"+", "-"} or query not in lengths:
                    raise ValueError(f"Unknown cmscan sequence or strand: {query}, {strand}")
                if (strand == "+" and start > end) or (strand == "-" and start < end):
                    raise ValueError("cmscan strand disagrees with coordinate direction")
                interval = Interval.from_public(min(start, end), max(start, end))
                if interval.end > lengths[query]:
                    raise ValueError("cmscan hit exceeds sequence")
                evalue = finite(cols[15])
                if evalue < 0:
                    raise ValueError("Negative E-value")
                hits.append(
                    Hit(
                        query,
                        interval,
                        cast(Strand, strand),
                        finite(cols[14]),
                        evalue,
                        model + ":" + cols[0],
                    )
                )
    except (OSError, IndexError, ValueError) as exc:
        raise ContractError(f"Invalid Infernal table {path}: {exc}") from exc
    return hits


def scan_models(
    records: list[Sequence], work: Path, settings: Settings, runner: Runner
) -> list[Hit]:
    if not records:
        return []
    fasta = work / "cm-input.fasta"
    write_fasta(fasta, records)
    lengths = {r.id: len(r.bases) for r in records}
    hits = []
    for index, model in enumerate(settings.model_paths()):
        table = work / f"cm-{index}.tbl"
        runner.run(
            [
                settings.cmscan,
                "--cpu",
                str(settings.threads),
                "--noali",
                "--fmt",
                "1",
                "--tblout",
                str(table),
                str(model),
                str(fasta),
            ]
        )
        hits.extend(parse_cm_table(table, lengths, model.name))
    return hits


def parse_alignment_table(
    path: Path,
    queries: dict[str, tuple[Array, str, Strand]],
    source: Sequence,
    arrays: list[Array],
    settings: Settings,
    genome_is_query: bool = False,
    short: bool = False,
) -> list[AntiRepeat]:
    results = []
    try:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.split()
                if len(cols) != 12:
                    raise ValueError("Expected 12-column BLAST-format alignment")
                qid, sid = cols[0:2]
                length, mismatches, gap_opens = map(int, cols[3:6])
                if length <= 0 or min(mismatches, gap_opens) < 0 or mismatches + gap_opens > length:
                    raise ValueError("Invalid alignment counts")
                if finite(cols[10]) < 0:
                    raise ValueError("Negative alignment E-value")
                finite(cols[11])
                qstart, qend, sstart, send = map(int, cols[6:10])
                if genome_is_query:
                    qid, sid = sid, qid
                    qstart, qend, sstart, send = sstart, send, qstart, qend
                if qid not in queries or sid != source.id:
                    raise ValueError("Unknown sequence in alignment")
                array, repeat, repeat_strand = queries[qid]
                interval = Interval.from_public(min(sstart, send), max(sstart, send))
                q_interval = Interval.from_public(min(qstart, qend), max(qstart, qend))
                if interval.end > len(source.bases) or q_interval.end > len(repeat):
                    raise ValueError("Alignment exceeds sequence length")
                if max(interval.length, q_interval.length) > length:
                    raise ValueError("Alignment span exceeds alignment length")
                # Query span, not gapped alignment length, keeps coverage within [0, 1].
                coverage = q_interval.length / len(repeat)
                identity = finite(cols[2]) / 100
                if not 0 <= identity <= 1:
                    raise ValueError("Invalid alignment identity")
                if short:
                    if (
                        coverage >= settings.coverage
                        or q_interval.length < SHORT_RESCUE_POLICY.min_query_span
                        or identity < SHORT_RESCUE_POLICY.min_identity
                    ):
                        continue
                elif coverage < settings.coverage or identity < settings.similarity:
                    continue
                if any(interval.overlap(a.interval) for a in arrays):
                    continue
                relative_forward = (qend - qstart) * (send - sstart) >= 0
                complementary: Strand = "-" if relative_forward else "+"
                strands: tuple[Strand, ...] = (
                    (complementary,)
                    if short
                    else (complementary, "+" if complementary == "-" else "-")
                )
                for strand in strands:
                    results.append(
                        AntiRepeat(
                            array,
                            repeat,
                            repeat_strand,
                            interval,
                            strand,
                            identity,
                            coverage,
                            strand == complementary,
                        )
                    )
    except (OSError, ValueError, KeyError) as exc:
        raise ContractError(f"Invalid anti-repeat alignment table {path}: {exc}") from exc
    return results


def search_anti_repeats(
    source: Sequence,
    fasta: Path,
    arrays: list[Array],
    work: Path,
    settings: Settings,
    runner: Runner,
    short_proposals: list[AntiRepeat] | None = None,
) -> list[AntiRepeat]:
    if not arrays:
        return []
    search_fasta = fasta
    if settings.mask_arrays:
        # Perfect repeat copies inside arrays otherwise dominate FASTA36's
        # secondary-hit search. Post-filtering those alignments is too late:
        # valid, weaker anti-repeats may never have been reported.
        bases = list(source.bases)
        for array in arrays:
            bases[array.interval.start : array.interval.end] = "N" * array.interval.length
        search_fasta = work / "anti-search.fasta"
        write_fasta(search_fasta, [Sequence(source.id, source.description, "".join(bases))])
    queries: dict[str, tuple[Array, str, Strand]] = {}
    for index, array in enumerate(arrays):
        for strand in ("+", "-"):
            repeat = array.repeat if strand == "+" else reverse_complement(array.repeat)
            queries[f"repeat{index:06d}{'f' if strand == '+' else 'r'}"] = (array, repeat, strand)
    repeats = work / "repeats.fasta"
    write_fasta(repeats, [Sequence(k, k, v[1]) for k, v in queries.items()])
    output = work / "anti-repeat.tsv"
    if settings.search_backend == "blastn":
        runner.run(
            [
                settings.blastn,
                "-task",
                "blastn-short",
                "-query",
                str(repeats),
                "-subject",
                str(search_fasta),
                "-strand",
                "both",
                "-dust",
                "no",
                "-soft_masking",
                "false",
                "-evalue",
                "100",
                "-outfmt",
                "6",
                "-out",
                str(output),
            ]
        )
    else:
        # FASTA36's legacy library-name parser has a fixed-size path buffer.
        # Absolute paths of 120 characters can be corrupted (".fastax").
        # Our per-record directory already owns both inputs, so use short names.
        runner.run(
            [
                settings.fasta36,
                "-T",
                str(settings.threads),
                "-m",
                "8",
                repeats.name,
                search_fasta.name,
            ],
            cwd=work,
            stdout=output,
        )
    hits = parse_alignment_table(output, queries, source, arrays, settings)
    # Both repeat orientations are already searched. The old additional pass
    # with an entire genome as the FASTA36 query is redundant for strand coverage
    # and triggers native band-alignment corruption on some real loci (Kki).
    # Keep the short repeat as query; malformed output still fails closed.
    if short_proposals is not None:
        short_proposals.extend(
            deduplicate_anti_repeats(
                parse_alignment_table(output, queries, source, arrays, settings, short=True)
            )
        )
    return deduplicate_anti_repeats(hits)


def deduplicate_anti_repeats(hits: list[AntiRepeat]) -> list[AntiRepeat]:
    unique: dict[tuple[str, int, int, str, str], AntiRepeat] = {}
    for hit in hits:
        key = (hit.array.id, hit.interval.start, hit.interval.end, hit.repeat_strand, hit.strand)
        incumbent = unique.get(key)
        if incumbent is None or (hit.identity * hit.coverage, hit.identity) > (
            incumbent.identity * incumbent.coverage,
            incumbent.identity,
        ):
            unique[key] = hit
    return [unique[k] for k in sorted(unique)]


def make_candidates(
    source: Sequence, anti_repeats: list[AntiRepeat], settings: Settings, prefix: str = "candidate"
) -> list[Candidate]:
    results = []
    for index, anti in enumerate(anti_repeats):
        if anti.strand == "+":
            interval = Interval(
                max(0, anti.interval.start - settings.upstream_length),
                min(len(source.bases), anti.interval.end + settings.tail_length),
            )
            local = Interval(
                anti.interval.start - interval.start, anti.interval.end - interval.start
            )
        else:
            interval = Interval(
                max(0, anti.interval.start - settings.tail_length),
                min(len(source.bases), anti.interval.end + settings.upstream_length),
            )
            local = Interval(interval.end - anti.interval.end, interval.end - anti.interval.start)
        results.append(
            Candidate(
                f"{prefix}{index:06d}",
                anti,
                interval,
                subsequence(source.bases, interval, anti.strand),
                local,
            )
        )
    return results


def parse_erpin(path: Path, lengths: dict[str, int]) -> dict[str, list[Hit]]:
    results: dict[str, list[Hit]] = defaultdict(list)
    query = None
    try:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    query = line[1:].split()[0]
                    if query not in lengths:
                        raise ValueError(f"Unknown ERPIN sequence {query}")
                elif re.match(r"^(FW|RC)\s", line):
                    if query is None:
                        raise ValueError("ERPIN hit before sequence header")
                    cols = line.split()
                    start, end = map(int, cols[2].split(".."))
                    interval = Interval.from_public(min(start, end), max(start, end))
                    if interval.end > lengths[query]:
                        raise ValueError("ERPIN hit exceeds sequence")
                    score, evalue = finite(cols[3]), finite(cols[4])
                    if evalue < 0:
                        raise ValueError("Negative ERPIN E-value")
                    results[query].append(
                        Hit(
                            query,
                            interval,
                            "+" if cols[0] == "FW" else "-",
                            score,
                            evalue,
                            "rho-indep",
                        )
                    )
    except (OSError, ValueError, IndexError) as exc:
        raise ContractError(f"Invalid ERPIN output {path}: {exc}") from exc
    return dict(results)


def terminators(
    candidates: list[Candidate], work: Path, settings: Settings, runner: Runner
) -> dict[str, list[Hit]]:
    if not candidates:
        return {}
    fasta = work / "terminator-input.fasta"
    write_fasta(fasta, [Sequence(c.id, c.id, c.sequence) for c in candidates])
    output = work / "erpin.txt"
    runner.run(
        [
            settings.erpin,
            settings.erpin_model,
            str(fasta),
            "-1,4",
            "-add",
            "2",
            "4",
            "2",
            "-pcw",
            "3.0",
            "-cutoff",
            "100%",
        ],
        stdout=output,
    )
    return parse_erpin(output, {c.id: len(c.sequence) for c in candidates})


def _interaction_row(row: dict[str, str], target_length: int, query_length: int) -> dict[str, Any]:
    target = Interval.from_public(int(row["start1"]), int(row["end1"]))
    query = Interval.from_public(int(row["start2"]), int(row["end2"]))
    if target.end > target_length or query.end > query_length:
        raise ValueError("IntaRNA interaction exceeds sequence")
    structure = row["hybridDP"]
    if structure.count("&") != 1:
        raise ValueError("IntaRNA hybrid must have exactly one ampersand")
    left, right = structure.split("&")
    if set(left) - {".", "("} or set(right) - {".", ")"}:
        raise ValueError("Invalid IntaRNA hybrid characters")
    if len(left) != target.length or len(right) != query.length:
        raise ValueError("IntaRNA hybrid lengths do not match coordinates")
    paired_target = [target.start + i + 1 for i, value in enumerate(left) if value == "("]
    paired_repeat = [query.start + i + 1 for i, value in enumerate(right) if value == ")"]
    if not paired_target or len(paired_target) != len(paired_repeat):
        raise ValueError("IntaRNA hybrid must have equal positive paired counts")
    return {
        "target_interval": target.public(),
        "repeat_interval": query.public(),
        "energy_kcal_mol": finite(row["E"]),
        "hybrid_dot_bracket": structure,
        "target_paired_positions": paired_target,
        "repeat_paired_positions": paired_repeat,
    }


def _interaction_key(row: dict[str, Any]) -> tuple[float, int, int, int, int, str]:
    """A total key makes selection independent of native CSV row order."""
    return (
        row["energy_kcal_mol"],
        row["target_interval"]["start"],
        row["target_interval"]["end"],
        row["repeat_interval"]["start"],
        row["repeat_interval"]["end"],
        row["hybrid_dot_bracket"],
    )


def _parse_interactions(
    path: Path,
    targets: dict[str, int],
    query_id: str,
    query_length: int,
    with_ids: bool,
) -> dict[str, dict[str, Any] | None]:
    results: dict[str, dict[str, Any] | None] = dict.fromkeys(targets)
    seen: dict[tuple[str, int, int, int, int, str], float] = {}
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter=";")
            required = {"start1", "end1", "start2", "end2", "E", "hybridDP"}
            if with_ids:
                required |= {"id1", "id2"}
            columns = reader.fieldnames or []
            if not required <= set(columns) or len(columns) != len(set(columns)):
                raise ValueError("Missing or duplicate IntaRNA output columns")
            for row in reader:
                if None in row or any(row.get(name) in (None, "") for name in required):
                    raise ValueError("Incomplete IntaRNA output row")
                identity = row["id1"] if with_ids else next(iter(targets))
                if identity not in targets or (with_ids and row["id2"] != query_id):
                    raise ValueError("Unknown IntaRNA target or query ID")
                value = _interaction_row(row, targets[identity], query_length)
                key = _interaction_key(value)
                site = (identity, *key[1:])
                if site in seen:
                    raise ValueError("Duplicate or inconsistent IntaRNA interaction record")
                seen[site] = key[0]
                current = results[identity]
                if current is None or key < _interaction_key(current):
                    results[identity] = value
    except (OSError, KeyError, TypeError, ValueError) as exc:
        raise ContractError(f"Invalid IntaRNA output {path}: {exc}") from exc
    return results


def parse_intarna(path: Path, target_length: int, query_length: int) -> dict[str, Any] | None:
    """Compatibility parser; validates all rows, including nonwinning rows."""
    return _parse_interactions(path, {"target": target_length}, "query", query_length, False)[
        "target"
    ]


def parse_intarna_batch(
    path: Path,
    targets: dict[str, int],
    query_id: str,
    query_length: int,
) -> dict[str, dict[str, Any] | None]:
    return _parse_interactions(path, targets, query_id, query_length, True)


def interactions(
    candidates: list[Candidate],
    work: Path,
    settings: Settings,
    runner: Runner,
    anchored: bool = False,
) -> dict[str, dict[str, Any] | None]:
    """Evaluate exactly one oriented repeat per target, never a query cross product.

    IntaRNA 3.2.0 restricts seedTRange to a single target. Anchored calls therefore
    remain individual even if their repeat and interval constraints are identical.
    Size 1 is an explicit serial path; failures never silently retry.
    """
    if len({c.id for c in candidates}) != len(candidates):
        raise ContractError("Duplicate IntaRNA candidate IDs")
    groups: dict[str, list[Candidate]] = defaultdict(list)
    for candidate in candidates:
        groups[candidate.anti.repeat].append(candidate)
    results: dict[str, dict[str, Any] | None] = {}
    batch_size = 1 if anchored else settings.intarna_batch_size
    for group_index, repeat in enumerate(sorted(groups)):
        members = sorted(groups[repeat], key=lambda c: c.id)
        for chunk_index, offset in enumerate(range(0, len(members), batch_size)):
            chunk = members[offset : offset + batch_size]
            stem = f"intarna-{group_index:06d}-{chunk_index:06d}"
            target_ids = {f"t{i:06d}": c for i, c in enumerate(chunk)}
            query_id = "q000000"
            target_path, query_path = work / f"{stem}.targets.fasta", work / f"{stem}.query.fasta"
            write_fasta(target_path, [Sequence(k, k, c.sequence) for k, c in target_ids.items()])
            write_fasta(query_path, [Sequence(query_id, query_id, repeat)])
            (work / f"{stem}.ids.json").write_text(
                json.dumps({k: c.id for k, c in target_ids.items()}, indent=2) + "\n"
            )
            output = work / f"{stem}.csv"
            command = [
                settings.intarna,
                "-t",
                str(target_path),
                "-q",
                str(query_path),
                "--outMode",
                "C",
                "--outCsvCols",
                "id1,id2,start1,end1,start2,end2,E,hybridDP",
                "--outNumber",
                "1",
                "--out",
                str(output),
                "--threads",
                str(settings.threads),
            ]
            if anchored:
                local = chunk[0].anti_local.public()
                command.extend(["--seedTRange", f"{local['start']}-{local['end']}"])
            runner.run(command)
            parsed = parse_intarna_batch(
                output, {k: len(c.sequence) for k, c in target_ids.items()}, query_id, len(repeat)
            )
            results.update({target_ids[k].id: value for k, value in parsed.items()})
    return results


def interaction(
    candidate: Candidate, work: Path, settings: Settings, runner: Runner
) -> dict[str, Any] | None:
    """Compatibility entry point for a single candidate."""
    single_work = Path(tempfile.mkdtemp(prefix="intarna-single-", dir=work))
    return interactions([candidate], single_work, settings, runner)[candidate.id]
