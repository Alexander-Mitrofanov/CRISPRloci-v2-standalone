"""Errors exposed by the CLI without hiding external-tool failures."""


class TracrError(Exception):
    """Invalid input, incompatible data, or failed scientific computation."""


class ContractError(TracrError):
    """An upstream result violates the supported data contract."""


class ToolError(TracrError):
    """An external tool failed or exceeded its timeout."""
