"""Compatibility imports for scientific-output validation modules."""

from .validation_adapter import _adapter_group_names
from .validation_common import _supported_schema_major
from .validation_pipeline import _validated_pipeline_status

__all__ = [
    "_adapter_group_names",
    "_supported_schema_major",
    "_validated_pipeline_status",
]
