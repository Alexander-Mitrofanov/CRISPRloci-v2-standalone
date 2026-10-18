#!/usr/bin/env python3
"""Exercise repeat UI against a staged release with temporary local services."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import time
from urllib.parse import urlsplit
import zipfile

from playwright.sync_api import expect, sync_playwright

from browser_leader import check
from runtime.common import ROOT


def exercise(origin, api_origin, output):
    errors, jobs = [], []
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
            for mode, title, fixture in [
                ("repeats", "Analyze repeat sequences", "repeats.fasta"),
                ("rna", "Analyze repeat sequences", "rna.fasta"),
                ("repeat_context", "Compare repeats in array context", "source.fasta"),
            ]:
                page.goto(origin)
                page.wait_for_load_state("networkidle")
                expect(page.get_by_text("Service ready", exact=True)).to_be_visible(
                    timeout=30000
                )
                page.get_by_role("link", name=title, exact=True).click()
                if mode == "rna":
                    page.get_by_role("radio", name=re.compile("RNA repeats")).check()
                page.locator("input[type=file]").set_input_files(
                    ROOT / "backend/tests/fixtures/repeats" / fixture
                )
                with page.expect_response(
                    lambda response: (
                        response.request.method == "POST"
                        and response.url.endswith("/api/v1/jobs")
                    )
                ) as submitted:
                    page.get_by_role("button", name="Compute", exact=True).click()
                assert submitted.value.status == 202
                created = submitted.value.json()
                jobs.append(created["job"]["job_id"])
                assert created["job"]["mode"] == ("repeats" if mode == "rna" else mode)
                assert created["job"]["options"]["molecule"] == (
                    "RNA" if mode == "rna" else "DNA"
                )
                expect(
                    page.get_by_role("tab", name="Repeat evidence", exact=True)
                ).to_be_visible(timeout=180000)
                page.get_by_role("tab", name="Repeat evidence", exact=True).click()
                panel = page.get_by_role("tabpanel", name="Repeat evidence", exact=True)
                table = panel.get_by_role("region", name="Repeat instances table")
                expect(table.get_by_role("row")).to_have_count(
                    2 if mode == "rna" else 7
                )
                expect(
                    panel.get_by_role("region", name="Selected repeat detail")
                ).to_be_visible()
                if mode == "repeats":
                    table.get_by_role("radio").nth(2).check()
                    expect(
                        panel.get_by_text(re.compile("No reference pairs"))
                    ).to_be_visible()
                    table.get_by_role("radio").first.check()
                if mode == "repeat_context":
                    expect(
                        panel.get_by_text(
                            re.compile("Local context: effective window 150")
                        )
                    ).to_be_visible()
                page.screenshot(
                    path=str(output / f"{mode}-desktop.png"), full_page=True
                )
                page.get_by_role("tab", name="Files & methods", exact=True).click()
                page.locator("summary").filter(has_text="Methods & provenance").click()
                expect(
                    page.get_by_role("link", name=re.compile("CRISPRrepeat revision"))
                ).to_have_attribute("href", re.compile("/commit/[0-9a-f]{40}$"))
                with page.expect_download() as downloaded:
                    page.get_by_role(
                        "button", name=re.compile("Download complete result bundle")
                    ).click()
                path = output / f"{mode}-results.zip"
                downloaded.value.save_as(path)
                with zipfile.ZipFile(path) as archive:
                    doc = json.loads(archive.read("results/repeat-analysis.json"))
                    assert doc["counts"]["instances"] == (1 if mode == "rna" else 6)
                    assert doc["orientation_policy"] == (
                        "forward" if mode == "rna" else "both"
                    )
                    assert doc["classification"]["subtype"]["status"] == "unsupported"
                    if mode == "repeat_context":
                        assert doc["counts"]["contexts"] == 2
                    else:
                        assert not doc["contexts"]
                    assert not any(
                        name.startswith("results/repeat/")
                        for name in archive.namelist()
                    )
                recovery = page.url
                page.reload()
                expect(
                    page.get_by_role("tab", name="Repeat evidence", exact=True)
                ).to_be_visible(timeout=30000)
                assert page.url == recovery
                page.get_by_role("tab", name="Repeat evidence", exact=True).click()
                page.set_viewport_size({"width": 390, "height": 844})
                assert page.evaluate(
                    "document.documentElement.scrollWidth <= window.innerWidth"
                )
                page.screenshot(path=str(output / f"{mode}-mobile.png"), full_page=True)
                assert not page.evaluate("window.cspViolations"), page.evaluate(
                    "window.cspViolations"
                )
                page.set_viewport_size({"width": 1440, "height": 1050})
            assert not errors, errors
            return {
                "job_ids": jobs,
                "checks": [
                    "real DNA repeat submission",
                    "real transcribed RNA submission",
                    "real CRISPRidentify array-context submission",
                    "explicit missing motif",
                    "both DNA alternatives",
                    "source-bound context",
                    "pinned repository provenance",
                    "authenticated complete downloads",
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
        default=ROOT / ".runtime/browser-repeats" / str(time.time_ns()),
    )
    args = parser.parse_args()
    check(args.release, args.output.resolve(), exercise_fn=exercise, label="Repeat")
