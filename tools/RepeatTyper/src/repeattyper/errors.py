"""Typed failures for stable CLI exit codes."""


class InputError(ValueError):
    """Invalid or oversized input (exit 2)."""


class ModelError(ValueError):
    """Missing, changed or incompatible model (exit 3)."""


class OutputError(ValueError):
    """Output cannot be published atomically (exit 4)."""
