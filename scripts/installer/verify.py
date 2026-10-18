"""Check the candidate with selected references before making it the default."""
from __future__ import annotations

import json
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time
from urllib.error import URLError
from urllib.request import urlopen

from .engine import local_environment
from runtime.progress import emit


def verify_local(root, work, release, value):
    env = local_environment(root)
    for key in value["references"]:
        receipt = json.loads((work / (key + ".reference.json")).read_text())
        env.update(receipt["environment"])
    env.update(CRISPR_WEB_HSTS="false", CRISPR_WEB_WORKER_CPU=str(value["cpu"]))
    with socket.socket() as socket_:
        socket_.bind(("127.0.0.1", 0))
        port = socket_.getsockname()[1]
    active_lock = root / ".runtime/active.lock"
    active_lock.parent.mkdir(exist_ok=True)
    if active_lock.is_symlink():
        raise ValueError("The runtime service lock must not be a symlink")
    with active_lock.open("a"):
        pass
    data = work / ("local-verification-" + str(time.time_ns()))
    process = subprocess.Popen([sys.executable, str(root / "scripts/serve.py"), "--release", str(release),
                                "--data-root", str(data), "--port", str(port)], cwd=root, env=env)

    def request(route):
        with urlopen(f"http://127.0.0.1:{port}" + route, timeout=5) as response:
            payload = response.read(1024**2 + 1)
            if len(payload) > 1024**2:
                raise ValueError("Local verification response is oversized")
            return json.loads(payload)

    try:
        emit("Starting the candidate with your selected databases")
        deadline = time.monotonic() + 300
        while process.poll() is None and time.monotonic() < deadline:
            try:
                health = request("/api/v1/health")
                if health.get("release_id") == Path(release).name and health.get("worker") == "ok":
                    break
            except (URLError, OSError):
                pass
            time.sleep(0.5)
        else:
            raise ValueError("The candidate did not become ready with the selected reference configuration")
        emit("Checking analysis modes and catalogue", 0, 2, "checks")
        modes = request("/api/v1/config")["modes"]
        expected_modes = json.loads((root / "frontend/src/contracts/module-ports-v1.json").read_text())
        if set(modes) != set(expected_modes):
            raise ValueError("The selected reference setup does not enable every analysis mode")
        emit("Checking analysis modes and catalogue", 1, 2, "checks")
        if request("/api/v1/catalog").get("available") is not True:
            raise ValueError("The selected catalogue is not available through the API")
        emit("Checking analysis modes and catalogue", 2, 2, "checks")
    finally:
        if process.poll() is None:
            process.send_signal(signal.SIGTERM)
        try:
            process.wait(timeout=45)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
