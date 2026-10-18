"""Current restriction evidence shared by internal completion routes."""

from copy import deepcopy
import hashlib

from components.internal_repeat_completion import (
    _json, _restrictions, _source_valid, _validate_lineage, parent_authorization,
)
from components.repeat_boundary_repair import source_geometry
from components.terminal_repair_completion import _overlap


AUTHORITY_FIELDS = (
    "source_path_classification_eligible", "source_path_classification_reason",
    "source_path_restriction_reasons", "source_path_provenance",
)
CURRENT_INVENTORY_POLICY = "pretransaction-current-authority-inventory-v1"


def authority_state(candidate):
    return {name: getattr(candidate, name) for name in AUTHORITY_FIELDS if hasattr(candidate, name)}


def unrestricted(source, candidate):
    """An absent legacy marker is distinct from a present malformed marker."""
    try:
        _source_valid(source, candidate)
        if (hasattr(candidate, "source_path_classification_eligible")
                and candidate.source_path_classification_eligible is not True):
            return False
        if _restrictions(candidate):
            return False
        lineage = getattr(candidate, "source_path_provenance", None)
        if hasattr(candidate, "source_path_provenance"):
            _validate_lineage(source, candidate, lineage)
            if any(item.derivation == "local_chain_family_parent" for item in lineage):
                return False
        return True
    except (ValueError, TypeError, KeyError, AttributeError):
        return False


def current_by_geometry(categories):
    result = {}
    for mapping in categories:
        for rows in mapping.values():
            for row in rows:
                result.setdefault(source_geometry(row[1]), []).append(row[1])
    return result


def capture_current_authority_inventory(source, categories, pending):
    """Snapshot all current rows overlapping any pending target, independently.

    No eligibility or phase-family filter is applied. Duplicate derivations and
    absent versus present authority attributes must survive report replay.
    """
    targets = sorted({source_geometry(candidate) for candidate in pending})
    rows = [{"geometry": geometry, "authority": deepcopy(authority_state(candidate))}
            for geometry, candidates in current_by_geometry(categories).items()
            if any(_overlap(geometry, target) for target in targets)
            for candidate in candidates]
    return {"policy": CURRENT_INVENTORY_POLICY, "complete": True,
            "source_sha256": hashlib.sha256(source.encode("ascii")).hexdigest(),
            "target_geometries": targets, "rows": sorted(rows, key=_json)}


def parent_current_allowed(source, geometry, current, structural, source_hash):
    for candidate in current.get(geometry, ()):
        try:
            parent_authorization(source, candidate, structural, source_hash)
        except (ValueError, TypeError, KeyError, AttributeError):
            # Ordinary source-only current rows have no lineage to discharge.
            if not unrestricted(source, candidate):
                return False
    return True
