"""Reconcile supplied-template candidates before structural fallback and scoring."""

from components.evidence.discovery_origins import inherit_discovery_origins, representation
from components.refinement.candidates import SameStartEndFilter


def reconcile_template_candidates(mapping):
    # Import locally: this helper uses exactly the historical locus partition.
    from components.refinement.stage import DetectionRefinement, _SourceLocusKey

    is_parent = DetectionRefinement._is_local_chain_family_parent
    interval = DetectionRefinement._fuzzy_source_interval

    def signature(candidate):
        normalized = DetectionRefinement._normalize_fuzzy_candidate(candidate)
        return (representation(normalized), is_parent(candidate),
                getattr(candidate, "source_path_classification_eligible", True))

    def template(candidate):
        return getattr(candidate, "discovery_route", None) == "template"

    if not any(template(c) for rows in mapping.values() for c in rows):
        return mapping

    # Self representatives always win exact cross-route equality. Preserve
    # pre-existing self grouping even when it has equal geometry in two keys.
    self_by_value = {}
    for candidates in mapping.values():
        for candidate in candidates:
            if not template(candidate):
                self_by_value.setdefault(signature(candidate), []).append(candidate)
    seen_templates, deduplicated, novel = {}, {}, set()
    for key, candidates in mapping.items():
        kept = []
        for candidate in candidates:
            if not template(candidate):
                kept.append(candidate)
                continue
            value = signature(candidate)
            existing = self_by_value.get(value) or seen_templates.get(value)
            if existing:
                for representative in existing:
                    inherit_discovery_origins(representative, (representative, candidate))
                continue
            seen_templates[value] = [candidate]
            kept.append(candidate)
            if not is_parent(candidate) and getattr(candidate, "source_path_classification_eligible", True) is True:
                novel.add(id(candidate))
        if kept or not candidates:
            deduplicated[key] = kept
    if not novel:
        return deduplicated

    # Normalize each historical key first, and keep diagnostic parents apart.
    deduplicated = SameStartEndFilter(deduplicated, preserve_roles=True).output()
    components = []
    for key, candidates in deduplicated.items():
        local = [c for c in candidates if not is_parent(c)]
        parents = [c for c in candidates if is_parent(c)]
        parts = DetectionRefinement._source_locus_components(local)
        parts += DetectionRefinement._source_locus_components(parents)
        if not parts:
            components.append((key, [], False, False))
        for part in parts:
            spans = [interval(c) for c in part]
            part_key = key if len(parts) == 1 else _SourceLocusKey(
                min(s[0] for s in spans), max(s[1] for s in spans))
            components.append((part_key, part, any(id(c) in novel for c in part), is_parent(part[0])))

    parent = list(range(len(components)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    # Sweep actual local candidate spans, never padded keys or parent envelopes.
    spans = sorted((interval(c), i) for i, (_, rows, _, restricted) in enumerate(components)
                   if not restricted for c in rows)
    maximum_end, representative = -1, None
    for (start, end), index in spans:
        if start < maximum_end:
            parent[root(index)] = root(representative)
        if start >= maximum_end or end > maximum_end:
            maximum_end, representative = end, index
    groups = {}
    for i in range(len(components)):
        groups.setdefault(root(i), []).append(i)
    affected = {i: group for group in groups.values()
                if any(components[j][2] for j in group) for i in group}
    output, emitted = {}, set()
    for i, (key, rows, _, _) in enumerate(components):
        if i not in affected:
            output[key] = rows
            continue
        group = affected[i]
        if group[0] in emitted:
            continue
        emitted.add(group[0])
        combined = [c for j in group for c in components[j][1]]
        bounds = [interval(c) for c in combined]
        new_key = _SourceLocusKey(min(s[0] for s in bounds), max(s[1] for s in bounds))
        output[new_key] = SameStartEndFilter({new_key: combined}, preserve_roles=True).output()[new_key]
    return output
