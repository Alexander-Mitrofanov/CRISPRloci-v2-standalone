#!/usr/bin/env python3
"""Local UI acceptance with synthetic fixtures and an intercepted API; no live jobs."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import re
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright
from offline_example import check_offline_example

ROOT = Path(__file__).resolve().parents[1]
FASTA = ">target_A\nAAACGTTGCAACGATTCGAGTCTTGACTCGAATCGTTGCAACGTCC\n"
TOKEN = "a" * 43
TABLE = b"target_id\ttarget_start\ttarget_end\ntarget_A\t3\t22\n"


def check(origin: str, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    fixture = json.loads((ROOT / "tests/fixtures/protospacer.json").read_text())
    fixture["artifacts"] = [{
        "artifact_id": "matches", "filename": "protospacer-matches.tsv",
        "label": "Complete protospacer matches", "media_type": "text/tab-separated-values",
        "size_bytes": len(TABLE),
    }]
    state = {"advertised": False, "job": fixture, "requests": [], "submissions": []}
    errors = []

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            context = browser.new_context(viewport={"width": 1440, "height": 1100}, accept_downloads=True)
            page = context.new_page()
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.add_init_script("window.cspViolations=[]; document.addEventListener('securitypolicyviolation', e=>window.cspViolations.push(e.violatedDirective));")

            def api(route):
                request = route.request
                path = urlsplit(request.url).path
                state["requests"].append((request.method, path))
                headers = {"Access-Control-Allow-Origin": origin, "Access-Control-Allow-Headers": "authorization,content-type"}
                if request.method == "OPTIONS":
                    route.fulfill(status=204, headers=headers)
                    return
                if path.endswith("/health"):
                    data = {"version": "2.2.0"}
                elif path.endswith("/config"):
                    data = {"modes": ["protospacer"] if state["advertised"] else ["detection"],
                            "max_total_bases": 1000, "max_record_bases": 1000,
                            "max_records": 10, "max_request_bytes": 10000}
                elif path.endswith("/jobs") and request.method == "POST":
                    payload = request.post_data_json
                    assert payload["mode"] == "protospacer"
                    assert payload["sequence"].strip() == FASTA.strip()
                    assert not {"database", "reference", "reference_path", "molecule"} & payload.keys()
                    state["submissions"].append(payload)
                    data = {"job_id": fixture["job_id"], "access_token": TOKEN, "job": state["job"]}
                elif "/jobs/" in path:
                    assert request.headers.get("authorization") == f"Bearer {TOKEN}"
                    assert TOKEN not in request.url
                    if path.endswith("/artifacts/matches"):
                        route.fulfill(status=200, headers=headers, content_type="text/tab-separated-values", body=TABLE)
                        return
                    data = state["job"]
                else:
                    raise AssertionError(f"Unexpected API request: {request.method} {path}")
                route.fulfill(status=202 if request.method == "POST" else 200,
                              headers=headers, content_type="application/json", body=json.dumps(data))

            page.route("https://analysis.example.org/**", api)
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            check_offline_example(page, state, link='Spacer searches', input_label='Contigs or small genomes', sequence=FASTA, submit='Search references')
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            page.goto(origin + "/?method=protospacer")
            page.wait_for_load_state("networkidle")
            page.get_by_label("Contigs or small genomes", exact=True).fill(FASTA)
            expect(page.get_by_role("button", name="Search references", exact=True)).to_be_disabled()
            page.get_by_text("Review input requirements", exact=True).click()
            expect(page.get_by_role("button", name=re.compile("^Find spacers in my sequence")).get_by_text("Reference unavailable on this server; precomputed example available", exact=True)).to_be_visible()
            expect(page.get_by_text("Selected analysis unavailable on this service", exact=False)).to_be_visible()

            state["advertised"] = True
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            page.get_by_role("link", name="Spacer searches", exact=True).click()
            expect(page.get_by_text("Analysis options", exact=True)).to_be_visible()
            expect(page.get_by_label("Reference path", exact=False)).to_have_count(0)
            page.get_by_label("Upload FASTA file", exact=True).set_input_files({
                "name": "synthetic-target.fasta", "mimeType": "text/plain", "buffer": FASTA.encode(),
            })
            page.get_by_role("button", name="Search references", exact=True).click()
            tab = page.get_by_role("tab", name="Protospacer matches", exact=True)
            tab.click()
            panel = page.get_by_role("tabpanel", name="Protospacer matches", exact=True)
            expect(panel.get_by_text("3–22", exact=True)).to_be_visible()
            expect(panel.get_by_text("25–44", exact=True)).to_be_visible()
            expect(panel.get_by_text("Reverse (−)", exact=True)).to_be_visible()
            expect(panel.get_by_text("d" * 64, exact=True)).to_be_visible()
            expect(page.get_by_role("tab", name="Arrays", exact=True)).to_have_count(0)
            page.screenshot(path=output / "desktop.png", full_page=True)
            page.reload()
            page.wait_for_load_state("networkidle")
            page.get_by_role("tab", name="Protospacer matches", exact=True).click()
            expect(panel.get_by_text("Forward (+)", exact=True)).to_be_visible()

            page.get_by_role("tab", name="Files & methods", exact=True).click()
            with page.expect_download() as download:
                page.get_by_role("button", name=re.compile("Complete protospacer matches")).click()
            saved = output / "matches.tsv"
            download.value.save_as(saved)
            assert saved.read_bytes() == TABLE

            page.set_viewport_size({"width": 390, "height": 844})
            page.get_by_role("tab", name="Protospacer matches", exact=True).click()
            assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth + 1")
            expect(panel.get_by_role("region", name="Protospacer match table")).to_be_visible()
            page.screenshot(path=output / "mobile.png", full_page=True)
            page.set_viewport_size({"width": 1440, "height": 1100})

            variants = [
                (None, "Protospacer evidence was not reported"),
                ({"available": True, "counts": {"matches": 0}}, "Search completion was not confirmed"),
                ({"available": True, "search_complete": True, "outcome": "no_assessable_windows", "counts": {"matches": 0}}, "No target windows could be assessed"),
                ({"available": True, "search_complete": True, "outcome": "no_matches", "counts": {"matches": 0}}, "No matches to eligible reference spacers"),
            ]
            for evidence, expected in variants:
                state["job"] = copy.deepcopy(fixture)
                state["job"]["summary"]["protospacer"] = evidence
                page.reload()
                page.wait_for_load_state("networkidle")
                overview = page.get_by_role("tabpanel", name="Overview", exact=True)
                expect(overview.get_by_text(re.compile(expected))).to_be_visible()
                if not evidence or evidence.get("outcome") != "no_matches":
                    expect(overview.get_by_text(re.compile("No matches to eligible"))).to_have_count(0)

            state["job"] = copy.deepcopy(fixture)
            evidence = state["job"]["summary"]["protospacer"]
            evidence["counts"]["matches"] = 120
            evidence["matches"] = [copy.deepcopy(evidence["matches"][0]) for _ in range(110)]
            page.reload()
            page.wait_for_load_state("networkidle")
            page.get_by_role("tab", name="Protospacer matches", exact=True).click()
            expect(panel.get_by_text(re.compile("preview shows 100 of 120"))).to_be_visible()
            expect(panel.get_by_role("row")).to_have_count(101)
            assert len(state["submissions"]) == 1
            assert not errors, errors
            assert not page.evaluate("window.cspViolations"), "CSP violations"
            (output / "acceptance.json").write_text(json.dumps({
                "status": "passed", "api": "intercepted synthetic fixtures", "deployment": False,
                "checks": ["advertised_capability", "direct_link_guard", "fasta_upload", "submission_payload",
                           "relative_strands", "inclusive_coordinates", "reference_hash", "recovery_reload",
                           "header_only_download", "missing_incomplete_no_assessable_no_match", "preview_cap", "mobile_layout", "no_console_or_csp_errors"],
            }, indent=2) + "\n")
            print("Protospacer frontend browser acceptance passed:", output)
        finally:
            browser.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:4187")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check(args.origin.rstrip("/"), args.output)
