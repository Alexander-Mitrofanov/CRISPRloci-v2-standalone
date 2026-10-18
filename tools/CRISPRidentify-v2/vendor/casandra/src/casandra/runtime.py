"""Dependency-light PyHMMER runtime used by the final predictor."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import pyhmmer


def parse_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split(maxsplit=1)[0]
                chunks = []
            elif name is None:
                raise ValueError("FASTA sequence occurs before a header")
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def scan_profiles(
    queries: list[dict[str, str]],
    hmms: list[pyhmmer.plan7.HMM],
    metadata: dict[str, dict[str, object]],
    threads: int,
) -> list[dict[str, object]]:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = [
        pyhmmer.easel.TextSequence(
            name=query["sequence_id"].encode(), sequence=query["sequence"]
        ).digitize(alphabet)
        for query in queries
    ]
    output: list[dict[str, object]] = []
    hits_iterator = pyhmmer.hmmer.hmmscan(
        sequences, hmms, cpus=threads, E=1e6, domE=1e6, incE=1e6, incdomE=1e6
    )
    for query, hits in zip(queries, hits_iterator, strict=True):
        best_by_use = {"primary_search": -200.0, "hard_negative": -200.0}
        best_profile: str | None = None
        best_positive = -200.0
        label_scores: dict[str, dict[str, float]] = {
            "class": {},
            "type": {},
            "subtype": {},
        }
        hit_count = 0
        for hit in hits:
            profile_id = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
            profile = metadata.get(profile_id)
            if profile is None:
                continue
            hit_count += 1
            score = float(hit.score)
            use = str(profile["intended_use"])
            best_by_use[use] = max(best_by_use[use], score)
            if use != "primary_search":
                continue
            if score > best_positive:
                best_positive = score
                best_profile = profile_id
            for head in label_scores:
                labels = profile["labels"][head]
                if len(labels) == 1:
                    label = str(labels[0])
                    label_scores[head][label] = max(label_scores[head].get(label, -200.0), score)
        output.append(
            {
                "sequence_id": query["sequence_id"],
                "positive_score": best_by_use["primary_search"],
                "negative_score": best_by_use["hard_negative"],
                "margin": best_by_use["primary_search"] - best_by_use["hard_negative"],
                "best_positive_profile": best_profile,
                "classification": {
                    head: max(scores, key=scores.get) if scores else None
                    for head, scores in label_scores.items()
                },
                "hit_count": hit_count,
            }
        )
    return output
