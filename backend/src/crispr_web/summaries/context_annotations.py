"""Bounded coordinate joins; proximity does not establish biological function."""

from .annotations import _rows

MAX_NEAREST_CASSETTES = 32


def cas_neighbors(root, intervals, checkpoint=None):
    if len(intervals) > 20_000:
        raise ValueError("Context interval budget exceeded")
    by_source = {}
    cassette_count = 0
    for directory in sorted((root / "casandra").glob("record-*")):
        for row in _rows(root, directory / "cassettes.jsonl", checkpoint):
            cassette_count += 1
            if cassette_count > 10_000:
                raise ValueError("Cas context cassette budget exceeded")
            by_source.setdefault(row["contig_id"], []).append({
                "id": row["cassette_id"], "start": row["start_1based"], "end": row["end_1based_inclusive"]})
    result = []
    for interval in intervals:
        nearest, count, candidates = None, 0, []
        for cassette in by_source.get(interval["source_id"], []):
            if checkpoint:
                checkpoint()
            overlap = max(0, min(interval["end"], cassette["end"]) - max(interval["start"], cassette["start"]) + 1)
            distance = 0 if overlap else max(cassette["start"] - interval["end"] - 1, interval["start"] - cassette["end"] - 1)
            if nearest is None or distance < nearest:
                nearest, count, candidates = distance, 0, []
            if distance == nearest:
                count += 1
                if len(candidates) < MAX_NEAREST_CASSETTES:
                    candidates.append({"cassette_id": cassette["id"], "start": cassette["start"], "end": cassette["end"],
                                       "overlap_bases": overlap, "intervening_bases": distance})
        result.append({"interval_id": interval["id"], "source_id": interval["source_id"],
                       "nearest_cassettes": candidates, "nearest_cassette_count": count,
                       "nearest_cassettes_truncated": count > len(candidates),
                       "interpretation": "Coordinate proximity in submitted sequence only; function and disruption not assessed."})
    return result
