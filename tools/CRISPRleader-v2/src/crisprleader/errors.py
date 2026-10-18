"""Errors distinguish broken execution/input from scientific no-call results."""


class LeaderError(Exception):
    """Base actionable error."""


class InputError(LeaderError):
    """Invalid or unbound scientific input."""


class ReferenceError(LeaderError):
    """Missing, incompatible, corrupt, or invalid reference bundle."""


class ExecutionError(LeaderError):
    """External execution failed, exceeded a bound, or was cancelled."""
