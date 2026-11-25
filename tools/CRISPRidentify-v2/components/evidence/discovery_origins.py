"""Observational discovery ancestry. This never grants classification authority."""

from dataclasses import dataclass

from components.detection.repeat_templates import digest


def representation(candidate):
    starts = getattr(candidate, "list_absolute_start", None)
    if starts is None:
        starts = candidate.list_repeat_starts
    gapped = getattr(candidate, "list_gaped_repeats", None)
    if gapped is None:
        gapped = candidate.list_repeats_gaped
    return (tuple(starts), tuple(candidate.list_repeats), tuple(gapped), tuple(candidate.list_spacers))


@dataclass(frozen=True)
class DiscoveryBirth:
    source_sha256: str
    panel_sha256: str
    pattern: str
    task_id: str
    representation: tuple


def mark_discovery_birth(candidate, source_hash, panel_hash, pattern, task_id):
    birth = DiscoveryBirth(source_hash, panel_hash, pattern, task_id, representation(candidate))
    candidate.discovery_route = "template"
    candidate.discovery_origins = (representation(candidate), (birth,))
    return {"id": digest(birth.representation), "pattern": pattern, "task_id": task_id,
            "representation": birth.representation}


def inherit_discovery_origins(target, parents):
    births = set()
    for parent in parents:
        origins = getattr(parent, "discovery_origins", None)
        if origins is None:
            continue
        if origins[0] != representation(parent):
            raise ValueError("stale discovery ancestry representation")
        births.update(origins[1])
    if births:
        target.discovery_origins = (representation(target), tuple(sorted(births, key=repr)))
