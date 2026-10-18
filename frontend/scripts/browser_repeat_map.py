#!/usr/bin/env python3
"""Local browser acceptance for repeat mapping; synthetic API, no live jobs."""

import argparse
import copy
import json
from pathlib import Path
import re
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright
from offline_example import check_offline_example

ROOT = Path(__file__).resolve().parents[1]
TOKEN = "a" * 43
FASTA = ">repeat\nACGTCGATGCTAGCTAGTCGATCG\n"
TABLE = b"query_id\tstatus\nrepeat\texact\n"


def check(origin, output):
    output.mkdir(parents=True, exist_ok=True)
    fixture = json.loads((ROOT / "tests/fixtures/repeat-map.json").read_text())
    fixture["artifacts"] = [{"artifact_id": "mapping", "name": "repeat-mapping.tsv", "media_type": "text/tab-separated-values", "size_bytes": len(TABLE)}]
    state = {"advertised": False, "job": fixture, "submissions": []}
    errors = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            context = browser.new_context(viewport={"width": 1440, "height": 1050}, accept_downloads=True)
            page = context.new_page()
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.add_init_script("window.cspViolations=[];document.addEventListener('securitypolicyviolation',e=>window.cspViolations.push(e.violatedDirective));")

            def api(route):
                request = route.request
                path = urlsplit(request.url).path
                headers = {"Access-Control-Allow-Origin": origin, "Access-Control-Allow-Headers": "authorization,content-type"}
                if request.method == "OPTIONS":
                    route.fulfill(status=204, headers=headers)
                    return
                if path.endswith("/health"):
                    data = {"version": "2.4.0"}
                elif path.endswith("/config"):
                    data = {"modes": ["repeat_map"] if state["advertised"] else ["detection"], "max_total_bases": 20000,
                            "max_record_bases": 20000, "max_records": 1000, "max_request_bytes": 100000}
                elif path.endswith("/jobs") and request.method == "POST":
                    payload = request.post_data_json
                    assert payload["mode"] == "repeat_map" and payload["sequence"].strip() == FASTA.strip()
                    assert "reference" not in payload
                    state["submissions"].append(payload)
                    data = {"job_id": fixture["job_id"], "access_token": TOKEN, "job": state["job"]}
                elif "/jobs/" in path:
                    assert request.headers.get("authorization") == "Bearer " + TOKEN
                    assert TOKEN not in request.url
                    if path.endswith("/artifacts/mapping"):
                        route.fulfill(status=200, headers=headers, content_type="text/tab-separated-values", body=TABLE)
                        return
                    data = state["job"]
                else:
                    raise AssertionError(f"Unexpected request: {path}")
                route.fulfill(status=202 if request.method == "POST" else 200, headers=headers, content_type="application/json", body=json.dumps(data))

            page.route("https://analysis.example.org/**", api)
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            check_offline_example(page, state, link='Map repeats to reference families', input_label='Repeat sequences', sequence=FASTA, submit='Compute')
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            state["advertised"] = True
            page.reload()
            page.wait_for_load_state("networkidle")
            page.get_by_role("link", name="Map repeats to reference families", exact=True).click()
            page.get_by_label("Repeat sequences", exact=True).fill(">long\n" + "A" * 201)
            expect(page.get_by_role("button", name="Compute", exact=True)).to_be_disabled()
            page.get_by_label("Repeat sequences", exact=True).fill(FASTA)
            page.get_by_role("button", name="Compute", exact=True).click()
            page.get_by_role("tab", name="Repeat matches", exact=True).click()
            panel = page.get_by_role("tabpanel", name="Repeat matches", exact=True)
            expect(panel.get_by_text("exact · exact", exact=True)).to_be_visible()
            expect(panel.get_by_text("negative · no_hit", exact=True)).to_be_visible()
            panel.get_by_label("Filter displayed repeats").fill("ambiguous")
            expect(panel.get_by_text("ambiguous · unsupported", exact=True)).to_be_visible()
            expect(panel.get_by_text("exact · exact", exact=True)).to_have_count(0)
            panel.get_by_label("Filter displayed repeats").fill("")
            page.screenshot(path=output / "desktop.png", full_page=True)
            page.reload()
            page.wait_for_load_state("networkidle")
            page.get_by_role("tab", name="Files & methods", exact=True).click()
            with page.expect_download() as download:
                page.get_by_role("button", name=re.compile("repeat-mapping.tsv")).click()
            saved = output / "repeat-mapping.tsv"
            download.value.save_as(saved)
            assert saved.read_bytes() == TABLE
            page.set_viewport_size({"width": 390, "height": 844})
            page.get_by_role("tab", name="Repeat matches", exact=True).click()
            panel.get_by_text("Reference alignments", exact=True).first.click()
            assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth + 1")
            page.screenshot(path=output / "mobile.png", full_page=True)
            state["job"] = copy.deepcopy(fixture)
            state["job"]["summary"]["repeat_map"] = None
            page.reload()
            page.wait_for_load_state("networkidle")
            expect(page.get_by_role("tabpanel", name="Overview", exact=True).get_by_text("A completed repeat mapping was not reported.")).to_be_visible()
            assert not errors, errors
            assert not page.evaluate("window.cspViolations")
            assert len(state["submissions"]) == 1
            (output / "acceptance.json").write_text(json.dumps({"status": "passed", "api": "synthetic intercepted fixtures", "checks": ["reference advertisement", "input bounds", "submission", "exact/near/no_hit/unsupported", "filter", "reload recovery", "authenticated download", "390px layout with alignment", "missing evidence", "no JS/CSP errors"]}, indent=2) + "\n")
            print("Repeat mapping browser acceptance passed")
        finally:
            browser.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:4198")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check(args.origin.rstrip("/"), args.output)
