"""Private exec helper: resource limits without preexec_fn in threaded callers."""

from __future__ import annotations
import os
import resource
import sys


def main() -> None:
    # This fresh single-threaded interpreter is in the process group established
    # by process.run. exec preserves that identity, file handles and the limit.
    from .io import MAX_BYTES

    resource.setrlimit(resource.RLIMIT_FSIZE, (MAX_BYTES, MAX_BYTES))
    try:
        os.execvp(sys.argv[1], sys.argv[1:])
    except OSError:
        raise SystemExit(127) from None


if __name__ == "__main__":
    main()
