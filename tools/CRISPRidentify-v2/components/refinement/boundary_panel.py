"""Curated repeat phase for newly extended, source-supported paths."""

from collections import Counter

from components.paths import TOOLS_ROOT
from components.sequence.utils import reverse_complement


def oriented_consensuses(value):
    entries = value.upper().split(",")
    if len(entries) > 256 or any(
        not 20 <= len(s) <= 60 or set(s) - set("ACGT") for s in entries
    ):
        raise ValueError(
            "provide up to 256 comma-separated ACGT repeat consensuses of 20–60 bases"
        )
    return sorted(set(entries) | {reverse_complement(s) for s in entries})


def configured_consensuses(value):
    if value is None:
        panel = TOOLS_ROOT / "blasting" / "curated_boundary_consensuses.txt"
        value = ",".join(panel.read_text(encoding="ascii").splitlines())
    if value.lower() == "none":
        return []
    return oriented_consensuses(value)


def resolve_phase(dna, path, motifs):
    """Require one alternative phase and at least two exact source exemplars."""
    if not motifs:
        return None
    widths = {end - start for start, end in path}
    if len(widths) != 1:
        return None
    (width,) = widths
    choices = {}
    for left in (-1, 0, 1):
        for right in (-1, 0, 1):
            seq_width = width + right - left
            panel = {s for s in motifs if len(s) == seq_width}
            if not panel or path[0][0] + left < 0 or path[-1][1] + right > len(dna):
                continue
            projected = tuple((a + left, b + right) for a, b in path)
            if any(b >= c for (_, b), (c, _) in zip(projected, projected[1:])):
                continue
            counts = Counter(dna[a:b] for a, b in projected)
            matches = []
            for motif in panel.intersection(counts):
                if counts[motif] < 2:
                    continue
                support = sum(
                    n
                    for s, n in counts.items()
                    if sum(x != y for x, y in zip(s, motif)) <= seq_width // 5
                )
                if 4 * support >= 3 * len(path):
                    matches.append(motif)
            if matches:
                choices[projected] = sorted(matches)
    if path in choices or len(choices) != 1:
        return None
    return next(iter(choices.items()))
