"""CPU-efficient PyHMMER scan with a validation-selected reporting cutoff."""

from __future__ import annotations

import pyhmmer


def scan_profiles_fast(queries, hmms, metadata, threads, report_evalue):
    """Scan all profiles while reporting only hits below ``report_evalue``.

    The original runtime requested every hit up to E=1e6 and then discarded
    nearly all of them in Python.  Detection only needs the best reported
    positive and hard-negative hits.  The cutoff is selected on validation
    data and frozen in the model card.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = [
        pyhmmer.easel.TextSequence(name=q["sequence_id"].encode(), sequence=q["sequence"]).digitize(alphabet)
        for q in queries
    ]
    output = []
    iterator = pyhmmer.hmmer.hmmscan(
        sequences, hmms, cpus=threads,
        E=report_evalue, domE=report_evalue,
        incE=report_evalue, incdomE=report_evalue,
    )
    for query, hits in zip(queries, iterator, strict=True):
        best_by_use = {"primary_search": -200.0, "hard_negative": -200.0}
        best_profile = None
        best_positive = -200.0
        label_scores = {"class": {}, "type": {}, "subtype": {}}
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
                best_positive, best_profile = score, profile_id
            for head in label_scores:
                labels = profile["labels"][head]
                if len(labels) == 1:
                    label = str(labels[0])
                    label_scores[head][label] = max(label_scores[head].get(label, -200.0), score)
        output.append({
            "sequence_id": query["sequence_id"],
            "positive_score": best_by_use["primary_search"],
            "negative_score": best_by_use["hard_negative"],
            "margin": best_by_use["primary_search"] - best_by_use["hard_negative"],
            "best_positive_profile": best_profile,
            "classification": {head: max(scores, key=scores.get) if scores else None for head, scores in label_scores.items()},
            "hit_count": hit_count,
        })
    return output

