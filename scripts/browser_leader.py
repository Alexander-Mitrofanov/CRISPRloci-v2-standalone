#!/usr/bin/env python3
"""Check the compiled leader UI against an isolated immutable release.

Starts only temporary local API/worker/preview processes. The selected server
release and its job database are never changed. Build frontend/ first.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import socket
import shutil
import subprocess
import time
from urllib.parse import urlsplit
from urllib.request import urlopen
import zipfile

from playwright.sync_api import expect, sync_playwright

from manage_release import environment, verify
from runtime.common import ROOT, atomic_json, source_inventory


def available_port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def check(release, output, *, exercise_fn=None, label="Leader"):
    release = release.resolve(strict=True)
    document = verify(release)
    output.mkdir(parents=True, exist_ok=False)
    frontend_identity = source_inventory(ROOT / "frontend/dist")
    shutil.copytree(ROOT / "frontend/dist", output / "frontend")
    if source_inventory(output / "frontend") != frontend_identity:
        raise RuntimeError(
            "Frontend build changed while making the browser-test snapshot"
        )
    api_port, ui_port = available_port(), available_port()
    origin = f"http://127.0.0.1:{ui_port}"
    api_origin = f"http://127.0.0.1:{api_port}"
    env = environment(release, output / "data", api_origin=origin, verified=document)
    env.update(
        CRISPR_WEB_WORKER_CPU="2",
        CRISPR_WEB_HSTS="false",
        CRISPR_WEB_RATE_LIMIT_SUBMISSIONS="20",
    )
    python = release / "environment/bin/python"
    commands = [
        (
            "api",
            [
                str(python),
                "-m",
                "uvicorn",
                "crispr_web.api:create_app",
                "--factory",
                "--host",
                "127.0.0.1",
                "--port",
                str(api_port),
                "--no-access-log",
                "--no-proxy-headers",
            ],
            output,
        ),
        ("worker", [str(python), "-m", "crispr_web.worker"], output),
        (
            "preview",
            [
                str(python),
                "-m",
                "http.server",
                str(ui_port),
                "--bind",
                "127.0.0.1",
                "--directory",
                str(output / "frontend"),
            ],
            output,
        ),
    ]
    children, streams = [], []
    try:
        for name, command, cwd in commands:
            stream = (output / f"{name}.log").open("w")
            streams.append(stream)
            children.append(
                subprocess.Popen(
                    command,
                    cwd=cwd,
                    env=env,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
            )
        deadline = time.monotonic() + 60
        while True:
            try:
                with urlopen(
                    api_origin + "/api/v1/health/fresh", timeout=2
                ) as response:
                    if json.load(response)["status"] == "ok":
                        break
            except OSError:
                if time.monotonic() > deadline or any(
                    child.poll() is not None for child in children
                ):
                    raise
                time.sleep(0.5)
        observations = (exercise_fn or exercise)(origin, api_origin, output)
        atomic_json(
            output / "acceptance.json",
            {
                "status": "passed",
                "release_id": release.name,
                "frontend_files": frontend_identity,
                **observations,
            },
        )
        print(f"{label} browser acceptance passed:", output, flush=True)
    finally:
        import os
        import signal

        for child in children:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
        for child in children:
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        for stream in streams:
            stream.close()


def exercise(origin, api_origin, output):
    errors = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            context = browser.new_context(
                viewport={"width": 1440, "height": 1050}, accept_downloads=True
            )
            page = context.new_page()
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.add_init_script(
                "window.cspViolations=[]; document.addEventListener('securitypolicyviolation', e=>window.cspViolations.push(e.violatedDirective));"
            )

            def proxy(route):
                parsed = urlsplit(route.request.url)
                response = route.fetch(
                    url=api_origin
                    + parsed.path
                    + ("?" + parsed.query if parsed.query else "")
                )
                route.fulfill(response=response)

            page.route("https://analysis.example.org/**", proxy)
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            expect(page.get_by_text("Service ready", exact=True)).to_be_visible(
                timeout=30000
            )
            page.get_by_role("link", name="Extract leader context", exact=True).click()
            page.locator("input[type=file]").set_input_files(
                ROOT / "backend/tests/fixtures/leader/source.fasta"
            )
            page.get_by_text("Analysis options", exact=True).click()
            page.get_by_label("Leader context window (nt per side)").fill("5000")
            page.screenshot(path=str(output / "leader-input.png"), full_page=True)
            with page.expect_response(
                lambda response: (
                    response.request.method == "POST"
                    and response.url.endswith("/api/v1/jobs")
                )
            ) as submitted:
                page.get_by_role("button", name="Compute", exact=True).click()
            assert submitted.value.status == 202
            created = submitted.value.json()
            assert created["job"]["options"]["leader_flank_length"] == 5000
            expect(
                page.get_by_role("heading", name="Leader context", exact=True)
            ).to_be_visible(timeout=180000)
            page.get_by_role("tab", name="Leader context", exact=True).click()
            panel = page.get_by_role("tabpanel", name="Leader context", exact=True)
            expect(
                panel.get_by_role("region", name="Leader context table").get_by_role(
                    "row"
                )
            ).to_have_count(3)
            expect(
                panel.get_by_text("Clipped at contig boundary", exact=True)
            ).to_have_count(2)
            expect(
                panel.get_by_text(re.compile("Leader prediction is unavailable"))
            ).to_be_visible()
            panel.get_by_text("View sequences", exact=True).first.click()
            page.screenshot(path=str(output / "leader-desktop.png"), full_page=True)
            page.get_by_role("tab", name="Files & methods", exact=True).click()
            page.locator("summary").filter(has_text="Methods & provenance").click()
            expect(
                page.get_by_role("link", name=re.compile("CRISPRleader-v2 revision"))
            ).to_have_attribute("href", re.compile("/commit/[0-9a-f]{40}$"))
            with page.expect_download() as archive:
                page.get_by_role(
                    "button", name=re.compile("Download complete result bundle")
                ).click()
            bundle_path = output / "leader-result.zip"
            archive.value.save_as(bundle_path)
            with zipfile.ZipFile(bundle_path) as bundle:
                data = json.loads(bundle.read("results/leader-contexts.json"))
                assert (
                    data["complete"] is True
                    and data["context_count"] == len(data["contexts"]) == 2
                )
                assert (
                    data["prediction_performed"] is False
                    and data["flank_length"] == 5000
                )
            recovery = page.url
            page.reload()
            expect(
                page.get_by_role("heading", name="Leader context", exact=True)
            ).to_be_visible(timeout=30000)
            assert page.url == recovery and "#job=" in recovery
            page.get_by_role("tab", name="Leader context", exact=True).click()
            page.set_viewport_size({"width": 390, "height": 844})
            page.get_by_role(
                "heading", name="Leader context", exact=True
            ).scroll_into_view_if_needed()
            assert page.evaluate(
                "document.documentElement.scrollWidth <= window.innerWidth"
            )
            page.screenshot(path=str(output / "leader-mobile.png"), full_page=True)
            assert not page.evaluate("window.cspViolations"), page.evaluate(
                "window.cspViolations"
            )
            assert not errors, errors
            return {
                "job_id": created["job"]["job_id"],
                "checks": [
                    "real leader submission",
                    "5000-nt clipping",
                    "both context sides",
                    "explicit unavailable prediction",
                    "pinned repository provenance",
                    "authenticated complete download",
                    "reload recovery",
                    "390px layout",
                    "compiled CSP",
                    "no JavaScript exceptions",
                ],
            }
        finally:
            browser.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("release", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / ".runtime/browser-leader" / str(time.time_ns()),
    )
    args = parser.parse_args()
    check(args.release, args.output.resolve())
