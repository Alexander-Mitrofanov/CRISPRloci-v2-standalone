"""One compact scientific result for every successful public analysis mode."""

from ..exchange import MAX_EXCHANGE_BYTES
from .compact.annotations import cas_data, leader_data, loci_data, tracr_data
from .compact.common import Reader, array_data
from .compact.evolution import evolution_data
from .compact.repeats import repeat_data
from .compact.repeat_type import repeat_type_data
from .compact.association import association_data
from .compact.array_compare import array_compare_data
from .compact.locus_compare import locus_compare_data
from .compact.prophage_context import prophage_context_data
from .compact.insertion_context import insertion_context_data
from .compact.defense_context import defense_context_data
from .compact.rna_family import rna_family_data
from .compact.search import protospacer_data, repeat_map_data, viral_data
from .safety import _contains_private_path, _private_markers
from .sequence_ports import write_document


PROJECTORS = {
    "locus_compare": locus_compare_data,
    "prophage_context": prophage_context_data,
    "insertion_context": insertion_context_data,
    "defense_context": defense_context_data,
    "rna_family": rna_family_data,
    "array_compare": array_compare_data,
    "detection": array_data,
    "cas": cas_data,
    "tracrrna": tracr_data,
    "loci": loci_data,
    "leader": leader_data,
    "repeats": repeat_data,
    "repeat_context": repeat_data,
    "protospacer": protospacer_data,
    "viral_search": viral_data,
    "repeat_map": repeat_map_data,
    "repeat_type": repeat_type_data,
    "spacer_association": association_data,
}
MODULES = frozenset({*PROJECTORS, "orientation", "reconstruction"})
SCHEMA = {"name": "crisprloci.module-result", "version": "1.0.0"}


def module_result(root, mode, status, *, checkpoint=None):
    if mode not in MODULES or status not in {
        "completed",
        "completed_no_eligible_groups",
    }:
        raise ValueError("compact results require a successful supported analysis")
    read = Reader(root, checkpoint)
    data = (
        evolution_data(read, mode, status)
        if mode in {"orientation", "reconstruction"}
        else PROJECTORS[mode](read)
    )
    return {
        "schema": dict(SCHEMA),
        "module": mode,
        "status": status,
        "complete": True,
        "data": data,
    }


def write_module_result(
    root, mode, status, *, max_bytes, private_roots=(), checkpoint=None
):
    """Never turn an unavailable export into apparently empty scientific data."""
    if status not in {"completed", "completed_no_eligible_groups"}:
        return []
    limit = min(max_bytes, MAX_EXCHANGE_BYTES)
    try:
        document = module_result(root, mode, status, checkpoint=checkpoint)
        if _contains_private_path(document, _private_markers((*private_roots, root))):
            raise ValueError("compact output contains private data")
        write_document(
            root, "result.json", document, max_bytes=limit, checkpoint=checkpoint
        )
        return []
    except (ValueError, KeyError, TypeError, OSError, UnicodeError, RecursionError):
        # No partial rows or empty-results substitution. A small status document
        # keeps the per-mode artifact contract explicit when bounded export fails.
        fallback = {
            "schema": dict(SCHEMA),
            "module": mode,
            "status": status,
            "complete": False,
            "reason": "compact_export_unavailable",
        }
        try:
            write_document(
                root, "result.json", fallback, max_bytes=limit, checkpoint=checkpoint
            )
        except (ValueError, OSError):
            pass
        return [
            "Compact result data could not be exported within validation and size limits; consult the detailed artifacts."
        ]
