"""Authenticated loopback-only browser wizard; no scientific job API endpoints."""
from __future__ import annotations

import argparse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import secrets
import signal
import subprocess
import sys
from urllib.error import URLError
from urllib.request import urlopen
import webbrowser

from .engine import Engine
from .plan import REFERENCES, ROOT, choices, defaults, stages, system

ASSETS = Path(__file__).with_name("assets")
MAX_BODY = 64 * 1024


class InstallerServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, port=0, root=ROOT):
        self.root = Path(root)
        self.token = secrets.token_urlsafe(32)
        self.engine = Engine(self.root)
        self.application = None
        self.application_log = None
        super().__init__(("127.0.0.1", port), Handler)
        self.origin = "http://127.0.0.1:" + str(self.server_port)

    def app_status(self):
        active = self.application is not None and self.application.poll() is None
        ready = False
        if active:
            try:
                with urlopen("http://127.0.0.1:8000/api/v1/health", timeout=0.5) as response:
                    ready = response.status == 200
                with urlopen("http://127.0.0.1:5173", timeout=0.5) as response:
                    ready = ready and response.status == 200
            except (URLError, OSError):
                pass
        return {"active": active, "ready": ready, "url": "http://127.0.0.1:5173",
                "installed": (self.root / ".installer/installation.json").is_file()}

    def start_application(self):
        if self.engine.snapshot()["status"] == "running":
            raise ValueError("Wait for installation to finish before starting CRISPRloci")
        if not self.app_status()["installed"]:
            raise ValueError("Finish installation before starting CRISPRloci")
        if self.application is not None and self.application.poll() is None:
            return
        if self.application_log is not None:
            self.application_log.close()
        self.application_log = (self.root / ".installer/application.log").open("ab")
        self.application = subprocess.Popen([sys.executable, self.root / "scripts/standalone.py", "--no-browser"],
                                            cwd=self.root, stdout=self.application_log, stderr=subprocess.STDOUT,
                                            start_new_session=True)

    def stop_application(self):
        if self.application is not None and self.application.poll() is None:
            os.killpg(self.application.pid, signal.SIGTERM)
            try:
                self.application.wait(timeout=40)
            except subprocess.TimeoutExpired:
                os.killpg(self.application.pid, signal.SIGKILL)
                self.application.wait()
        if self.application_log is not None:
            self.application_log.close()
            self.application_log = None


class Handler(BaseHTTPRequestHandler):
    server_version = "CRISPRlociInstaller"

    def setup(self):
        super().setup()
        self.connection.settimeout(10)

    def log_message(self, *_args):
        pass

    def send(self, status, value, content_type="application/json"):
        payload = json.dumps(value).encode() if content_type == "application/json" else value
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy", "default-src 'none'; script-src 'self'; style-src 'self'; connect-src 'self'; img-src 'self'; base-uri 'none'; frame-ancestors 'none'; form-action 'none'")
        self.end_headers()
        self.wfile.write(payload)

    def authorized(self):
        if self.headers.get("Host") != self.server.origin.removeprefix("http://"):
            self.send(403, {"error": "Open the exact local installer address printed in your terminal"})
            return False
        origin = self.headers.get("Origin")
        if origin is not None and origin != self.server.origin:
            self.send(403, {"error": "Installer requests must come from its own local page"})
            return False
        if self.path.startswith("/api/") and not secrets.compare_digest(
            self.headers.get("Authorization", ""), "Bearer " + self.server.token
        ):
            self.send(401, {"error": "Use the installer link printed in your terminal, including its session fragment"})
            return False
        return True

    def do_GET(self):
        if not self.authorized():
            return
        if self.path == "/api/options":
            self.send(200, {"defaults": defaults(), "references": REFERENCES, "system": system(self.server.root, defaults())})
        elif self.path == "/api/state":
            self.send(200, self.server.engine.snapshot())
        elif self.path == "/api/app":
            self.send(200, self.server.app_status())
        elif self.path in {"/", "/app.js", "/style.css"}:
            name, mime = {"/": ("index.html", "text/html; charset=utf-8"),
                          "/app.js": ("app.js", "text/javascript; charset=utf-8"),
                          "/style.css": ("style.css", "text/css; charset=utf-8")}[self.path]
            self.send(200, (ASSETS / name).read_bytes(), mime)
        else:
            self.send(404, {"error": "Not found"})

    def do_POST(self):
        if not self.authorized():
            return
        try:
            raw_length = self.headers.get("Content-Length", "")
            if len(self.headers.get_all("Content-Length", [])) != 1 or not raw_length.isdecimal() or not 0 < int(raw_length) <= MAX_BODY or self.headers.get("Transfer-Encoding"):
                self.send(413, {"error": "Installer request is empty, oversized or has unsupported framing"})
                return
            if self.headers.get_content_type() != "application/json":
                self.send(415, {"error": "Use a JSON installer request"})
                return
            data = json.loads(self.rfile.read(int(raw_length)))
            if self.path == "/api/plan":
                value = choices(data)
                self.send(200, {"choices": value, "stages": stages(value), "system": system(self.server.root, value)})
            elif self.path == "/api/install":
                value = choices(data)
                issues = system(self.server.root, value)["errors"]
                if issues:
                    raise ValueError(" ".join(issues))
                self.server.engine.start(value)
                self.send(202, self.server.engine.snapshot())
            elif self.path == "/api/cancel":
                self.server.engine.cancel()
                self.send(202, {"status": "cancelling"})
            elif self.path == "/api/app/start":
                self.server.start_application()
                self.send(202, self.server.app_status())
            elif self.path == "/api/app/stop":
                self.server.stop_application()
                self.send(200, self.server.app_status())
            else:
                self.send(404, {"error": "Not found"})
        except (ValueError, TypeError, KeyError, OSError) as error:
            self.send(400, {"error": str(error)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--headless", type=Path, metavar="CHOICES_JSON", help="Run the same installation from explicit saved choices")
    parser.add_argument("--check", action="store_true", help="Print prerequisites without installing anything")
    args = parser.parse_args()

    def terminate(_signum, _frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, terminate)
    if args.check:
        info = system()
        print(json.dumps(info, indent=2))
        return int(bool(info["errors"]))
    if args.headless:
        engine = Engine()
        engine.start(json.loads(args.headless.read_text()))
        try:
            previous = None
            while engine.thread.is_alive():
                state = engine.snapshot()
                update = {"stage": state.get("current"), "completed_stages": state["completed"],
                          "total_stages": state["total"], "progress": state["progress"]}
                if update != previous:
                    print(json.dumps(update), flush=True)
                    previous = update
                engine.thread.join(timeout=1)
        except KeyboardInterrupt:
            engine.cancel()
            engine.thread.join(timeout=65)
            return 130
        state = engine.snapshot()
        print(json.dumps(state, indent=2))
        return int(state["status"] != "complete")
    server = InstallerServer(args.port)
    url = server.origin + "/#" + server.token
    print("Open the local CRISPRloci installer: " + url, flush=True)
    print("The installer runs only on this computer. Press Ctrl+C to stop it and any app started here.", flush=True)
    if not args.no_browser:
        webbrowser.open(url)
    try:
        server.serve_forever(poll_interval=0.25)
    except KeyboardInterrupt:
        pass
    finally:
        server.engine.cancel()
        if server.engine.thread is not None:
            server.engine.thread.join(timeout=65)
        server.stop_application()
        server.server_close()
    return 0
