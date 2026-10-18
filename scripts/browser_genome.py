#!/usr/bin/env python3
"""Exercise the compiled genome UI with synthetic capability-protected artifacts.

Build frontend first, then serve frontend/dist locally. This test never starts
scientific tools, accesses production, or reads a real job. Screenshots and
downloaded figures stay in the explicitly selected local output directory.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import time
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]
JOB_ID = "0123456789abcdef0123456789abcdef"
TOKEN = "a" * 43
API_ORIGIN = "https://analysis.example.org"


def exercise(origin: str, output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    fixture = json.loads(subprocess.check_output([
        "node", "--input-type=module", "-e",
        "import {genomeFixture} from './tests/support/genomeFixture.js'; "
        "console.log(JSON.stringify(genomeFixture()));",
    ], cwd=ROOT / "frontend", text=True))
    artifacts = {
        "result-artifact": fixture["result"],
        "input-artifact": fixture["input"],
    }
    job = {
        "job_id": JOB_ID, "mode": "loci", "status": "completed",
        "created_at": "2026-01-01T00:00:00Z", "expires_at": "2099-01-01T00:00:00Z",
        "options": {}, "summary": {},
        "artifacts": [
            {"artifact_id": key, "name": name, "size_bytes": len(json.dumps(artifacts[key]))}
            for key, name in [("result-artifact", "result.json"), ("input-artifact", "input-sequences.json")]
        ],
    }
    reads, requests, errors = [], [], []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            context = browser.new_context(viewport={"width": 1440, "height": 1100}, accept_downloads=True)
            page = context.new_page()
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.on("request", lambda request: requests.append(request.url))
            page.add_init_script("window.cspViolations=[]; document.addEventListener('securitypolicyviolation',"
                                 "e=>window.cspViolations.push({directive:e.violatedDirective,uri:e.blockedURI}));")

            def api(route):
                request = route.request
                path = urlsplit(request.url).path
                assert TOKEN not in request.url
                if "/jobs/" in path:
                    assert request.headers.get("authorization") == f"Bearer {TOKEN}"
                if "/artifacts/" in path:
                    key = path.rsplit("/", 1)[1]
                    reads.append(key)
                    response = artifacts[key]
                elif path.endswith(JOB_ID):
                    response = job
                elif path.endswith("/config"):
                    response = {"api_version": "1.0.0", "modes": ["loci"], "max_total_bases": 2000000}
                else:
                    response = {"status": "ok", "version": "test"}
                route.fulfill(status=200, content_type="application/json", body=json.dumps(response))

            page.route(f"{API_ORIGIN}/**", api)
            page.goto(f"{origin}/#job={JOB_ID}.{TOKEN}")
            page.wait_for_load_state("networkidle")
            tab = page.get_by_role("tab", name="Genome & locus", exact=True)
            expect(tab).to_be_visible()
            assert not reads, "The viewer must not read DNA before it is opened"
            assert not any("/igv" in url for url in requests), "IGV must remain lazy loaded"
            tab.click()
            panel = page.get_by_role("tabpanel", name="Genome & locus", exact=True)
            expect(panel.get_by_role("button", name="Zoom in", exact=True)).to_be_enabled(timeout=30000)
            page.wait_for_load_state("networkidle")
            assert reads == ["result-artifact", "input-artifact"]
            assert panel.locator(".igv-viewport canvas").count() > 0

            array_viewport = panel.locator(".igv-viewport").nth(2)
            extent = array_viewport.bounding_box()
            array_viewport.click(position={"x": extent["width"] / 4, "y": 15})
            expect(panel.get_by_role("region", name="Selected genomic feature")).to_contain_text("array-1")
            array_viewport.click(button="right", position={"x": extent["width"] / 4, "y": 15})
            expect(panel.get_by_text("BLAT reference sequence", exact=True)).not_to_be_visible()

            panel.get_by_label("Feature", exact=True).select_option("feature_0")
            expect(panel.get_by_role("region", name="Selected genomic feature")).to_contain_text("422–580 (159 nt)")
            expect(panel.get_by_role("region", name="Selected genomic feature")).to_contain_text("Not established")
            panel.get_by_role("button", name="Zoom in", exact=True).click()
            panel.get_by_role("button", name="Zoom out", exact=True).click()
            panel.get_by_label("Start", exact=True).fill("450")
            panel.get_by_label("End", exact=True).fill("490")
            panel.get_by_role("button", name="Go to interval", exact=True).click()
            with page.expect_download() as download:
                panel.get_by_role("button", name="Export SVG figure").click()
            figure = output / "genome.svg"
            download.value.save_as(figure)
            svg = figure.read_text()
            assert "CRISPRloci genomic context: alpha" in svg
            assert TOKEN not in svg and "https://analysis.example.org" not in svg
            assert "<script" not in svg and "<foreignObject" not in svg
            panel.get_by_role("button", name="Open evidence tables").click()
            expect(page.get_by_role("tab", name="Arrays", exact=True)).to_have_attribute("aria-selected", "true")
            tab.click()
            assert reads == ["result-artifact", "input-artifact"], "Returning must reuse the viewer"
            expect(panel.get_by_role("region", name="Selected genomic feature")).to_contain_text("array-1")
            page.screenshot(path=str(output / "genome-desktop.png"), full_page=True)

            # Change records, select the same local array ID, and retain the correct source.
            panel.get_by_label("Sequence record").select_option("beta")
            panel.get_by_label("Feature", exact=True).select_option("feature_4")
            expect(panel.locator(".genome-caption")).to_contain_text("beta")
            page.set_viewport_size({"width": 390, "height": 844})
            expect(panel.get_by_role("button", name="Whole record")).to_be_visible()
            panel.get_by_role("button", name="Whole record").focus()
            page.keyboard.press("Enter")
            assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
            page.screenshot(path=str(output / "genome-mobile.png"), full_page=True)
            assert not errors, errors
            assert not page.evaluate("window.cspViolations"), page.evaluate("window.cspViolations")
            allowed = {urlsplit(origin).netloc, urlsplit(API_ORIGIN).netloc}
            assert all(urlsplit(url).netloc in allowed for url in requests), requests

            # Complete:false remains unavailable, not an apparently empty genome.
            artifacts["result-artifact"] = {**fixture["result"], "complete": False}
            page.reload()
            page.get_by_role("tab", name="Genome & locus", exact=True).click()
            expect(page.get_by_role("alert")).to_contain_text("Complete genome evidence is unavailable")
            assert page.locator(".igv-container").count() == 0

            # Exercise the documented client limits, with no server or scientific run.
            large = copy.deepcopy(fixture)
            sequence = "ACGT" * 2_500_000
            large["input"]["records"] = [{"id": "alpha", "sequence": sequence}]
            data = large["result"]["data"]
            data["sources"] = [{"id": "alpha", "length": len(sequence),
                                "sha256": hashlib.sha256(sequence.encode()).hexdigest()}]
            data["detection"]["arrays"] = data["detection"]["arrays"][:1]
            data["cas"]["genes"] = [
                {"id": f"gene_{index}", "source_id": "alpha", "start": index * 200 + 1,
                 "end": index * 200 + 90, "strand": "+", "profile": "synthetic"}
                for index in range(49_993)
            ]
            artifacts["result-artifact"] = large["result"]
            artifacts["input-artifact"] = large["input"]
            page.set_viewport_size({"width": 1440, "height": 1100})
            page.reload()
            page.get_by_role("tab", name="Genome & locus", exact=True).click()
            started = time.monotonic()
            expect(panel.get_by_role("button", name="Zoom in", exact=True)).to_be_enabled(timeout=30000)
            elapsed = time.monotonic() - started
            expect(panel.locator(".genome-caption")).to_contain_text("10,000,000 nt")
            expect(panel.get_by_text("100 of 50000 matching features listed.", exact=False)).to_be_visible()
            assert not page.evaluate("window.cspViolations")
            assert not errors, errors
            return {"status": "passed", "artifacts_authenticated": True, "external_requests": 0,
                    "csp_violations": 0, "desktop_mobile_keyboard_export": True,
                    "max_bases": 10_000_000, "max_features": 50_000,
                    "limit_load_seconds": round(elapsed, 2)}
        finally:
            browser.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = exercise(args.origin.rstrip("/"), args.output)
    (args.output / "acceptance.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
