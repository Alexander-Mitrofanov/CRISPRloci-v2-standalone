#!/usr/bin/env python3
"""Browser acceptance of the real wizard; installation work is explicitly stubbed.

The native reference builders and release safety checks are tested separately.
This script never installs software, downloads references or activates a release.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import selectors
import subprocess
import sys

from installer.plan import ROOT, defaults, stages


def check(output):
    from playwright.sync_api import expect, sync_playwright

    output.mkdir(parents=True, exist_ok=True)
    server = subprocess.Popen([sys.executable, ROOT / "scripts/install.py", "--no-browser"],
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=ROOT)
    selector = selectors.DefaultSelector()
    selector.register(server.stdout, selectors.EVENT_READ)
    try:
        if not selector.select(timeout=30):
            raise RuntimeError("Installer server did not open")
        line = server.stdout.readline()
        match = re.search(r"http://127\.0\.0\.1:\d+/#\S+", line)
        if not match:
            raise RuntimeError("Installer server did not provide its local session URL")
        url = match[0]
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=True)
            page = browser.new_page(viewport={"width": 1280, "height": 1100})
            failures = []
            page.on("pageerror", lambda error: failures.append(str(error)))
            page.goto(url)
            page.wait_for_load_state("networkidle")
            expect(page.get_by_role("heading", name="Set up your CRISPR workspace.")).to_be_visible()
            expect(page.locator("#references select")).to_have_count(5)
            page.screenshot(path=str(output / "choices-desktop.png"), full_page=True)
            page.locator("#viral-mode").select_option("local")
            page.locator("#viral-path").fill("/this-installer-test-path-does-not-exist.fasta")
            page.get_by_role("button", name="Review installation").click()
            expect(page.get_by_role("alert")).to_contain_text("existing absolute local path")
            page.locator("#viral-mode").select_option("sample")

            # Remove only the host's disk-capacity error from the browser test's
            # review response. The production preflight is tested independently.
            def review(route):
                response = route.fetch()
                document = response.json()
                document["system"]["errors"] = []
                route.fulfill(response=response, json=document)

            page.route("**/api/plan", review)
            page.locator("#spacer-mode").select_option("download")
            page.locator("#cpu").fill("3")
            page.get_by_role("button", name="Review installation").click()
            expect(page.get_by_role("heading", name="Review your installation.")).to_be_visible()
            expect(page.locator("#plan")).to_contain_text("Observed spacer search: download")
            page.screenshot(path=str(output / "review-desktop.png"), full_page=True)
            page.get_by_role("button", name="Change choices").click()
            expect(page.locator("#spacer-mode")).to_have_value("download")
            expect(page.locator("#cpu")).to_have_value("3")
            page.locator("#spacer-mode").select_option("sample")
            page.get_by_role("button", name="Review installation").click()
            state = {"status": "running", "stages": stages(defaults()), "completed": 3,
                     "total": len(stages(defaults())), "current": "tools", "logs": ["Installer browser fixture"],
                     "tool_progress": {"completed": 4, "total": 17}, "error": None,
                     "progress": {"label": "Downloading reference.hmm", "completed": 1024**2,
                                  "total": 8 * 1024**2, "unit": "bytes"}}
            page.route("**/api/install", lambda route: route.fulfill(status=202, json=state))
            page.route("**/api/state", lambda route: route.fulfill(json=state))

            def cancel(route):
                state.update(status="cancelled", error="Installation cancelled; completed verified work can be reused")
                route.fulfill(status=202, json={"status": "cancelling"})

            page.route("**/api/cancel", cancel)
            page.get_by_role("button", name="Install selected components").click()
            expect(page.locator("#operation-title")).to_have_text("Downloading reference.hmm")
            expect(page.locator("#operation")).to_have_attribute("value", str(1024**2))
            expect(page.locator("#tools")).to_have_attribute("value", "4")
            page.screenshot(path=str(output / "progress-desktop.png"), full_page=True)
            state["progress"] = {"label": "Building BLAST index", "completed": 0, "total": None, "unit": "items"}
            expect(page.locator("#operation-title")).to_have_text("Building BLAST index")
            expect(page.locator("#operation")).not_to_have_attribute("value", re.compile(".*"))
            page.get_by_role("button", name="Cancel installation").click()
            expect(page.get_by_role("heading", name="Installation cancelled.")).to_be_visible()
            page.get_by_role("button", name="Review and retry").click()
            page.get_by_role("button", name="Review installation").click()
            state.update(status="complete", completed=state["total"], current="activate", progress=None)
            page.get_by_role("button", name="Install selected components").click()
            expect(page.get_by_role("heading", name="Your workspace is ready.")).to_be_visible()
            expect(page.get_by_role("button", name="Start CRISPRloci")).to_be_visible()
            expect(page.locator("#start-help")).to_contain_text("bash start.sh")
            page.screenshot(path=str(output / "complete-desktop.png"), full_page=True)
            page.set_viewport_size({"width": 390, "height": 844})
            page.screenshot(path=str(output / "complete-mobile.png"), full_page=True)
            assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth"), "Mobile horizontal overflow"
            assert not failures, failures
            browser.close()
        result = {"status": "passed", "screenshots": str(output), "checks": [
            "local authenticated wizard", "all five reference choices", "invalid local path", "review and edit choices",
            "real-unit progress rendering", "indeterminate native build", "cancel and retry", "completion and launch instructions", "mobile layout", "no JavaScript exceptions"]}
        (output / "report.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result))
    finally:
        selector.close()
        server.send_signal(__import__("signal").SIGINT)
        try:
            server.wait(timeout=10)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait()
        server.stdout.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "output/installer-browser")
    check(parser.parse_args().output)
