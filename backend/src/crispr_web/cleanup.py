"""One-shot retention cleanup for a systemd timer."""

from __future__ import annotations

from .config import Settings
import os

from .service import JobService, require_unprivileged_process


def main() -> None:
    require_unprivileged_process()
    os.umask(0o077)
    removed = JobService(Settings.from_env()).cleanup_expired()
    print(f"removed {removed} expired CRISPR web job(s)")


if __name__ == "__main__":
    main()
