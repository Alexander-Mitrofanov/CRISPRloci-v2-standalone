#!/usr/bin/env python3
"""Run the selected immutable CRISPRloci release (loopback by default)."""

from __future__ import annotations

import argparse
import fcntl
import os
from pathlib import Path
import signal
import subprocess
import time

from manage_release import environment, require_gate, verify
from runtime.common import ROOT, RUNTIME


def serve(args):
    RUNTIME.mkdir(exist_ok=True)
    # Hold the lock before resolving 'current': selection and process lifetime
    # are atomic with respect to activation. Each process keeps it after exec.
    with (RUNTIME / "active.lock").open("r") as lock:
        fcntl.flock(lock, fcntl.LOCK_SH)
        release = args.release.resolve(strict=True)
        document = verify(release)
        env = environment(
            release, args.data_root, args.frontend_origin, verified=document
        )
        if getattr(args, "catalog", None) is not None:
            # Reference data is operator-selected and independent of the job queue
            # and immutable scientific tool release. Never fetch on a request.
            env["CRISPR_WEB_CATALOG_DATABASE"] = str(args.catalog.absolute())
        require_gate(release)
        python = release / "environment/bin/python"
        commands = {
            "api": [
                python,
                "-m",
                "uvicorn",
                "crispr_web.api:create_app",
                "--factory",
                "--host",
                args.host,
                "--port",
                str(args.port),
                "--no-proxy-headers",
                "--no-access-log",
                "--loop",
                "uvloop",
                "--http",
                "httptools",
                "--ws",
                "none",
                "--limit-concurrency",
                "20",
                "--timeout-keep-alive",
                "5",
            ],
            "worker": [python, "-m", "crispr_web.worker"],
            "cleanup": [python, "-m", "crispr_web.cleanup"],
        }
        if args.service != "all":
            os.set_inheritable(lock.fileno(), True)
            os.execve(python, [str(part) for part in commands[args.service]], env)
        children = []
        stopping = False

        def stop(signum, _frame):
            nonlocal stopping
            stopping = True
            for process in children:
                if process.poll() is None:
                    process.send_signal(signum)

        signal.signal(signal.SIGINT, stop)
        signal.signal(signal.SIGTERM, stop)
        try:
            for name in ("api", "worker"):
                children.append(
                    subprocess.Popen(
                        commands[name], cwd=ROOT, env=env, pass_fds=(lock.fileno(),)
                    )
                )
            while not stopping and all(process.poll() is None for process in children):
                time.sleep(0.25)
        finally:
            stop(signal.SIGTERM, None)
            for process in children:
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
        return max(abs(process.returncode) for process in children)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release", type=Path, default=RUNTIME / "current")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(os.getenv("CRISPR_WEB_DATA_ROOT", ROOT / ".runtime-data")),
    )
    parser.add_argument(
        "--frontend-origin",
        default=os.getenv("CRISPR_WEB_CORS_ORIGINS", "http://127.0.0.1:5173"),
    )
    parser.add_argument(
        "--service", choices=["all", "api", "worker", "cleanup"], default="all"
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument(
        "--catalog", type=Path,
        help="immutable catalogue registry (also CRISPR_WEB_CATALOG_DATABASE)",
    )
    raise SystemExit(serve(parser.parse_args()))
