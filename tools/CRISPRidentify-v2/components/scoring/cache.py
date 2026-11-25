"""Reuse scores only for unchanged candidates with the same provenance."""

from collections import deque


class _CandidateScoreSnapshot:
    __slots__ = ("candidate", "entry", "row")

    def __init__(self, entry):
        self.entry = entry
        self.candidate = entry[1]
        self.row = tuple(entry)


class CandidateScoreCache:
    """Reuse scores only while the exact scored entry remains unchanged.

    Enhancement replaces the candidate at index one of an existing score row.
    The row object's identity is therefore a stable provenance token even when
    structural filtering moves the row to another category. Equality with the
    original candidate is an independent fence: provenance alone never permits
    reuse after sequence, alignment, spacer, or coordinate changes.
    """

    def __init__(self, categories):
        self._snapshots = {}
        self.reused_count = 0
        self.recomputed_count = 0
        for _category_index, _interval, _entry_index, entry in self._entries(
                categories):
            self._snapshots.setdefault(
                id(entry), _CandidateScoreSnapshot(entry)
            )

    @staticmethod
    def _entries(categories):
        for category_index, category in enumerate(categories):
            for interval in sorted(category):
                for entry_index, entry in enumerate(category[interval]):
                    if not isinstance(entry, (tuple, list)) or len(entry) < 2:
                        raise RuntimeError(
                            "cannot cache a malformed scored candidate entry"
                        )
                    yield category_index, interval, entry_index, entry

    @staticmethod
    def _candidate_is_unchanged(snapshot, entry):
        return (
            snapshot is not None
            and snapshot.entry is entry
            and snapshot.candidate == entry[1]
        )

    def merge_after_enhancement(self, categories, score_changed_candidates):
        """Return all score rows in legacy traversal order.

        ``score_changed_candidates`` receives only candidates that failed the
        provenance-plus-structural-equality fence. It must return the ordinary
        ``{interval: [score rows...]}`` representation in candidate order.
        """

        slots = []
        changed_candidates = {}
        reused_count = 0
        recomputed_count = 0

        for _category_index, interval, _entry_index, entry in self._entries(
                categories):
            candidate = entry[1]
            snapshot = self._snapshots.get(id(entry))
            if self._candidate_is_unchanged(snapshot, entry):
                cached_row = list(snapshot.row)
                # Preserve the final post-enhancement object even when it is a
                # newly constructed but structurally equal candidate.
                cached_row[1] = candidate
                slots.append((interval, candidate, cached_row))
                reused_count += 1
            else:
                changed_candidates.setdefault(interval, []).append(candidate)
                slots.append((interval, candidate, None))
                recomputed_count += 1

        if changed_candidates:
            fresh_rows = score_changed_candidates(changed_candidates)
        else:
            fresh_rows = {}

        unexpected_intervals = set(fresh_rows).difference(changed_candidates)
        if unexpected_intervals:
            raise RuntimeError(
                "fresh evaluation returned unexpected intervals: {}".format(
                    sorted(unexpected_intervals)
                )
            )
        fresh_queues = {
            interval: deque(fresh_rows.get(interval, ()))
            for interval in changed_candidates
        }

        merged_rows = {}
        for interval, candidate, cached_row in slots:
            if cached_row is None:
                queue = fresh_queues[interval]
                if not queue:
                    raise RuntimeError(
                        "fresh evaluation omitted a changed candidate"
                    )
                row = list(queue.popleft())
                if len(row) < 2 or row[1] is not candidate:
                    raise RuntimeError(
                        "fresh evaluation changed candidate order or identity"
                    )
                row[1] = candidate
            else:
                row = cached_row
            merged_rows.setdefault(interval, []).append(row)

        if any(fresh_queues.values()):
            raise RuntimeError("fresh evaluation returned surplus score rows")

        self.reused_count = reused_count
        self.recomputed_count = recomputed_count
        return merged_rows
