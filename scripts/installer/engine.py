"""Run one explicit, resumable local installation with bounded live state."""
from __future__ import annotations

import copy
import json
import math
import os
from pathlib import Path
import selectors
import signal
import subprocess
import sys
import threading
import time

try:
    import fcntl
except ImportError:  # The GUI can explain unsupported native Windows installs.
    fcntl = None

from .download import download
from .plan import ROOT, choices, stages, system, work_directory

MAMBA_URL = "https://github.com/mamba-org/micromamba-releases/releases/download/2.3.3-0/micromamba-linux-64"
MAMBA_SHA256 = "9496f94a8b78c536573c93d946ec9bba74bd9ff79ee55aaa4b546e30db8f511b"
PREFIX = "CRISPR_PROGRESS "


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def local_environment(root):
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(("CRISPR_WEB_", "CRISPR_RELEASE_", "VITE_"))
           and key not in {"PYTHON", "PYTHONPATH", "PYTHONHOME", "CONDA_PREFIX", "CONDA_DEFAULT_ENV", "MAFFT_BINARIES"}}
    base = root / ".installer"
    env.update(PATH=str(base / "bin") + os.pathsep + str(base / "bootstrap/bin") + os.pathsep + env.get("PATH", ""),
               MAMBA_ROOT_PREFIX=str(base / "mamba"), PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1",
               XDG_CACHE_HOME=str(base / "cache"), PIP_CACHE_DIR=str(base / "cache/pip"),
               npm_config_cache=str(base / "cache/npm"),
               PYTHONUNBUFFERED="1", CRISPR_INSTALL_PROGRESS="1", PIP_DISABLE_PIP_VERSION_CHECK="1")
    return env


class Engine:
    def __init__(self, root=ROOT):
        self.root = Path(root)
        self.directory = self.root / ".installer"
        self.lock = threading.RLock()
        self.cancelled = threading.Event()
        self.thread = None
        self.process = None
        self.state = {"status": "idle", "stages": [], "logs": [], "error": None,
                      "completed": 0, "total": 0, "progress": None}

    def snapshot(self):
        with self.lock:
            return copy.deepcopy(self.state)

    def change(self, **fields):
        with self.lock:
            self.state.update(fields)

    def log(self, line):
        line = line.strip()[-4096:]
        if not line:
            return
        if line.startswith(PREFIX):
            try:
                progress = json.loads(line[len(PREFIX):])
                self.progress(progress["label"], progress["completed"], progress["total"], progress["unit"])
                return
            except (ValueError, KeyError, TypeError):
                pass
        with self.lock:
            self.state["logs"] = [*self.state["logs"][-149:], line]
            if line.startswith("PASS "):
                self.state["progress"] = {"label": "Validated scientific jobs", "completed":
                                          self.state.get("gate_cases", 0) + 1, "total": None, "unit": "jobs"}
                self.state["gate_cases"] = self.state["progress"]["completed"]
            if line.startswith("Acquired or verified "):
                import re

                match = re.match(r"Acquired or verified (\d+)/(\d+) records", line)
                if match:
                    self.progress("Downloading phage records", int(match[1]), int(match[2]), "records")

    def progress(self, label, completed=0, total=None, unit="items"):
        if not isinstance(label, str) or type(completed) not in (int, float) or not math.isfinite(completed) or completed < 0:
            return
        if total is not None and (type(total) not in (int, float) or not math.isfinite(total) or total <= 0 or completed > total):
            total = None
        self.change(progress={"label": label[:300], "completed": completed, "total": total, "unit": str(unit)[:30]})
        if unit == "tools":
            self.change(tool_progress={"completed": completed, "total": total})

    def start(self, raw):
        value = choices(raw)
        if self.directory.is_symlink():
            raise ValueError("The installer state directory must not be a symlink")
        self.directory.mkdir(mode=0o700, exist_ok=True)
        with self.lock:
            if self.thread is not None and self.thread.is_alive():
                raise ValueError("An installation is already running")
            self.cancelled.clear()
            self.state = {"status": "running", "stages": stages(value), "completed": 0,
                          "total": len(stages(value)), "logs": [], "error": None, "progress": None}
            self.thread = threading.Thread(target=self.install, args=(value,), daemon=True)
            self.thread.start()

    def cancel(self):
        self.cancelled.set()

    def command(self, arguments, *, env=None):
        if self.cancelled.is_set():
            raise InterruptedError("Installation cancelled")
        self.log("Running: " + " ".join(str(part) for part in arguments))
        process = subprocess.Popen([str(part) for part in arguments], cwd=self.root,
                                   env=env or local_environment(self.root), stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        self.process = process
        stopped_at = None
        exited_at = None
        buffer = b""
        selector = selectors.DefaultSelector()
        selector.register(process.stdout, selectors.EVENT_READ)
        try:
            while selector.get_map():
                if process.poll() is not None:
                    exited_at = exited_at or time.monotonic()
                    if time.monotonic() - exited_at > 2:
                        # A background descendant must not retain our output pipe
                        # or outlive this fixed installer command.
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                if self.cancelled.is_set() and process.poll() is None:
                    if stopped_at is None:
                        os.killpg(process.pid, signal.SIGTERM)
                        stopped_at = time.monotonic()
                    elif time.monotonic() - stopped_at > 5:
                        os.killpg(process.pid, signal.SIGKILL)
                for key, _ in selector.select(timeout=0.2):
                    chunk = os.read(key.fileobj.fileno(), 16384)
                    if not chunk:
                        selector.unregister(key.fileobj)
                        break
                    buffer += chunk.replace(b"\r", b"\n")
                    while b"\n" in buffer:
                        line, buffer = buffer.split(b"\n", 1)
                        self.log(line.decode("utf-8", errors="replace"))
                    if len(buffer) > 8192:
                        self.log(buffer[:8192].decode("utf-8", errors="replace"))
                        buffer = buffer[8192:]
            if buffer:
                self.log(buffer.decode("utf-8", errors="replace"))
            code = process.wait()
            if self.cancelled.is_set():
                raise InterruptedError("Installation cancelled; completed verified work will be reused on retry")
            if code:
                raise ValueError(f"Command failed (exit {code}). See the installation log for the failing dependency or input.")
        finally:
            selector.close()
            process.stdout.close()
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            if process.poll() is None:
                process.wait()
            self.process = None

    def bootstrap(self):
        binary = self.directory / "bin/micromamba"
        download(MAMBA_URL, binary, MAMBA_SHA256, self.progress, self.cancelled.is_set)
        binary.chmod(0o755)
        # Existing provisioners use the Conda CLI; micromamba supports the same
        # create/list interface. This alias lives only in this installation.
        alias = binary.with_name("conda")
        if not alias.exists():
            alias.symlink_to("micromamba")
        if not alias.is_symlink() or alias.readlink() != Path("micromamba"):
            raise ValueError("The managed Conda alias changed; inspect it before resuming")
        prefix = self.directory / "bootstrap"
        marker = self.directory / "bootstrap.json"
        packages = ["python=3.12", "pip", "nodejs=24", "git", "cxx-compiler", "make", "cmake", "ninja", "pkg-config", "coreutils", "tar", "gzip"]
        if not marker.exists():
            self.progress("Resolving and installing Python 3.12, Node 24 and native build tools")
            action = "install" if (prefix / "conda-meta/history").is_file() else "create"
            self.command([binary, action, "--yes", "--prefix", prefix, "--override-channels",
                          "-c", "conda-forge", *packages])
            from .download import file_digest

            atomic_json(marker, {"packages": packages, "python": file_digest(prefix / "bin/python"),
                                 "node": file_digest(prefix / "bin/node")})
        else:
            from .download import file_digest

            receipt = json.loads(marker.read_text())
            if receipt != {"packages": packages, "python": file_digest(prefix / "bin/python"),
                           "node": file_digest(prefix / "bin/node")}:
                raise ValueError("The managed bootstrap environment changed; inspect it before resuming")
        self.command([prefix / "bin/python", "-c", "import sys; assert sys.version_info[:2] == (3, 12)"])

    def install(self, value):
        if fcntl is None:
            self.change(status="failed", error="Scientific installation requires Linux x86-64")
            return
        if self.directory.is_symlink():
            self.change(status="failed", error="The installer state directory must not be a symlink")
            return
        self.directory.mkdir(mode=0o700, exist_ok=True)
        with (self.directory / "install.lock").open("a") as install_lock:
            try:
                fcntl.flock(install_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                work = work_directory(self.root, value)
                work.mkdir(parents=True, exist_ok=True)
                config = work / "choices.json"
                atomic_json(config, value)
                checkpoint = work / "completed.json"
                completed = json.loads(checkpoint.read_text()) if checkpoint.exists() else []
                if not isinstance(completed, list) or any(item not in {s["id"] for s in stages(value)} for item in completed):
                    raise ValueError("Installer checkpoint is invalid")
                for position, stage in enumerate(stages(value)):
                    if self.cancelled.is_set():
                        raise InterruptedError("Installation cancelled")
                    self.change(current=stage["id"], completed=position, progress=None)
                    self.log(stage["title"])
                    if stage["id"] == "preflight":
                        issues = system(self.root, value)["errors"]
                        if issues:
                            raise ValueError(" ".join(issues))
                        active = self.root / ".runtime/active.lock"
                        if active.exists():
                            with active.open("r") as active_lock:
                                try:
                                    fcntl.flock(active_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                                except BlockingIOError as error:
                                    raise ValueError("Stop the local API and worker before installing or changing this standalone runtime") from error
                        self.command([sys.executable, self.root / "scripts/sync_tools.py", "--check"])
                    elif stage["id"] == "bootstrap":
                        self.bootstrap()
                    else:
                        python = self.directory / "bootstrap/bin/python"
                        self.command([self.directory / "bin/micromamba", "run", "--prefix", self.directory / "bootstrap",
                                      python, self.root / "scripts/installer_task.py", stage["id"],
                                      "--work", work, "--config", config])
                    self.change(completed=position + 1)
                    if stage["id"] not in completed:
                        completed.append(stage["id"])
                    atomic_json(checkpoint, completed)
                self.change(status="complete", progress=None)
            except InterruptedError as error:
                self.change(status="cancelled", error=str(error))
            except Exception as error:
                self.log(str(error))
                self.change(status="failed", error=str(error))
            finally:
                atomic_json(self.directory / "last-install.json", self.snapshot())
