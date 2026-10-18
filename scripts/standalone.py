#!/usr/bin/env python3
"""Start the validated standalone API, worker and local browser interface."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time
from urllib.error import URLError
from urllib.request import urlopen
import webbrowser

ROOT = Path(__file__).resolve().parents[1]


def port_available(port):
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", port))


def start(no_browser=False, api_port=8000, frontend_port=5173):
    from installer.engine import local_environment

    settings = ROOT / ".installer/installation.json"
    if not settings.is_file():
        raise ValueError("Run bash install.sh and finish installation first")
    installed = json.loads(settings.read_text())
    if installed.get("schema_version") != 1 or installed.get("project") != str(ROOT):
        raise ValueError("The installation moved or has incompatible settings; rerun the installer in its final location")
    if api_port == frontend_port or any(not 1024 <= p <= 65535 for p in (api_port, frontend_port)):
        raise ValueError("Choose two distinct local ports between 1024 and 65535")
    for port in (api_port, frontend_port):
        port_available(port)
    env = local_environment(ROOT)
    env.pop("CRISPR_INSTALL_PROGRESS", None)
    env.update(installed["environment"])
    env.update(CRISPR_WEB_WORKER_CPU=str(installed["choices"]["cpu"]), CRISPR_WEB_HSTS="false",
               VITE_API_BASE_URL=f"http://127.0.0.1:{api_port}")
    bootstrap = ROOT / ".installer/bootstrap/bin"
    origin = f"http://127.0.0.1:{frontend_port}"
    commands = [
        [bootstrap / "python", ROOT / "scripts/serve.py", "--release", installed["release"],
         "--data-root", installed["data_root"], "--frontend-origin", origin, "--port", str(api_port)],
        [bootstrap / "node", ROOT / "frontend/node_modules/vite/bin/vite.js", "--host", "127.0.0.1",
         "--port", str(frontend_port), "--strictPort"],
    ]
    children = []
    stopping = False

    def stop(_signum=None, _frame=None):
        nonlocal stopping
        stopping = True
        for process in children:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)

    signal.signal(signal.SIGINT, stop)
    signal.signal(signal.SIGTERM, stop)
    try:
        for command, cwd in zip(commands, [ROOT, ROOT / "frontend"], strict=True):
            children.append(subprocess.Popen([str(arg) for arg in command], cwd=cwd, env=env, start_new_session=True))
        deadline = time.monotonic() + 300
        opened = False
        while not stopping and all(process.poll() is None for process in children):
            if not opened:
                try:
                    with urlopen(f"http://127.0.0.1:{api_port}/api/v1/health", timeout=2) as response:
                        healthy = response.status == 200
                    with urlopen(origin, timeout=2) as response:
                        healthy = healthy and response.status == 200
                    if healthy:
                        print("CRISPRloci is ready at " + origin + ". Press Ctrl+C to stop.", flush=True)
                        if not no_browser:
                            webbrowser.open(origin)
                        opened = True
                except (URLError, OSError):
                    if time.monotonic() > deadline:
                        raise ValueError("The local services did not become ready; check the startup messages")
            time.sleep(0.25)
    finally:
        stop()
        for process in children:
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
    return max((abs(p.returncode) for p in children), default=1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--api-port", type=int, default=8000)
    parser.add_argument("--frontend-port", type=int, default=5173)
    args = parser.parse_args()
    try:
        raise SystemExit(start(args.no_browser, args.api_port, args.frontend_port))
    except (ValueError, OSError) as error:
        sys.exit(str(error))
