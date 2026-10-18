"""Connect CRISPRidentify output to SpacerPlacer and CRISPR-evOr."""

from .adapter import adapt_reports
from .parser import (
    ArrayRecord,
    ReportFormatError,
    SourceBinding,
    load_fasta_source_bindings,
    load_identified_arrays,
)

__all__ = [
    "ArrayRecord",
    "ReportFormatError",
    "SourceBinding",
    "adapt_reports",
    "load_fasta_source_bindings",
    "load_identified_arrays",
]
__version__ = "0.15.0"
