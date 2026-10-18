#!/usr/bin/env python3
"""Desktop/mobile acceptance using real CCTK public evidence and an intercepted API."""

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
TABLE = b"source_a\tarray_a\tsource_b\tarray_b\tshared_spacers\tjaccard\ncontrol-A\tsame-array-name\tcontrol-B\tsame-array-name\t1\t0.3333333333333333\n"


def check(origin, output):
    output.mkdir(parents=True, exist_ok=True)
    fixture = json.loads((ROOT / "tests/fixtures/array-compare.json").read_text())
    packet = json.loads((ROOT / "tests/fixtures/array-compare-input.json").read_text())
    fixture["artifacts"] = [
        {
            "artifact_id": "sharing",
            "name": "array-sharing.tsv",
            "media_type": "text/tab-separated-values",
            "size_bytes": len(TABLE),
        }
    ]
    state = {"advertised": False, "job": fixture, "submissions": []}
    errors = []
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        try:
            context = browser.new_context(
                viewport={"width": 1440, "height": 1050}, accept_downloads=True
            )
            page = context.new_page()
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.add_init_script(
                "window.cspViolations=[];document.addEventListener('securitypolicyviolation',e=>window.cspViolations.push(e.violatedDirective));"
            )

            def api(route):
                request = route.request
                path = urlsplit(request.url).path
                headers = {
                    "Access-Control-Allow-Origin": origin,
                    "Access-Control-Allow-Headers": "authorization,content-type",
                }
                if request.method == "OPTIONS":
                    route.fulfill(status=204, headers=headers)
                    return
                if path.endswith("/health"):
                    data = {"version": "2.10.0"}
                elif path.endswith("/config"):
                    data = {
                        "modes": ["array_compare"]
                        if state["advertised"]
                        else ["detection"],
                        "max_total_bases": 100000,
                        "max_record_bases": 1000,
                        "max_records": 20,
                        "max_request_bytes": 100000,
                    }
                elif path.endswith("/jobs") and request.method == "POST":
                    payload = request.post_data_json
                    assert (
                        payload["mode"] == "array_compare"
                        and json.loads(payload["sequence"]) == packet
                    )
                    state["submissions"].append(payload)
                    data = {
                        "job_id": fixture["job_id"],
                        "access_token": TOKEN,
                        "job": state["job"],
                    }
                elif "/jobs/" in path:
                    assert (
                        request.headers.get("authorization") == "Bearer " + TOKEN
                        and TOKEN not in request.url
                    )
                    if path.endswith("/artifacts/sharing"):
                        route.fulfill(
                            status=200,
                            headers=headers,
                            content_type="text/tab-separated-values",
                            body=TABLE,
                        )
                        return
                    data = state["job"]
                else:
                    raise AssertionError(path)
                route.fulfill(
                    status=202 if request.method == "POST" else 200,
                    headers=headers,
                    content_type="application/json",
                    body=json.dumps(data),
                )

            page.route("https://analysis.example.org/**", api)
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            check_offline_example(page, state, link='Compare observed arrays', input_label='Observed array JSON', sequence=json.dumps(packet), submit='Compute')
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            state["advertised"] = True
            page.reload()
            page.wait_for_load_state("networkidle")
            page.get_by_role("link", name="Compare observed arrays", exact=True).click()
            field = page.get_by_label("Observed array JSON", exact=True)
            field.fill(">genome\nACGT")
            expect(
                page.get_by_role("button", name="Compute", exact=True)
            ).to_be_disabled()
            bad = copy.deepcopy(packet)
            bad["arrays"][0]["spacers"][0]["sequence"] = None
            field.fill(json.dumps(bad))
            expect(
                page.get_by_role("button", name="Compute", exact=True)
            ).to_be_disabled()
            field.fill(json.dumps(packet))
            page.get_by_role("button", name="Compute", exact=True).click()
            page.get_by_role("tab", name="Array comparison", exact=True).click()
            panel = page.get_by_role("tabpanel", name="Array comparison", exact=True)
            expect(
                panel.get_by_role("group", name=re.compile("sharing graph"))
            ).to_be_visible()
            expect(panel.get_by_text("0.333", exact=True)).to_be_visible()
            expect(panel.get_by_text("×", exact=True)).to_have_count(2)
            expect(panel.get_by_text("Unknown", exact=True).first).to_be_visible()
            panel.get_by_role("button", name=re.compile("Highlight control-A")).click()
            expect(panel.locator(".observed-array.selected")).to_have_count(1)
            panel.locator(".observed-array").first.locator("summary").click()
            assert page.evaluate(
                "document.documentElement.scrollWidth <= window.innerWidth + 1"
            )
            page.evaluate("window.scrollTo(0, 0)")
            page.screenshot(path=output / "desktop.png", full_page=True)
            page.reload()
            page.wait_for_load_state("networkidle")
            page.get_by_role("tab", name="Files & methods", exact=True).click()
            with page.expect_download() as download:
                page.get_by_role("button", name=re.compile("array-sharing.tsv")).click()
            download.value.save_as(output / "array-sharing.tsv")
            assert (output / "array-sharing.tsv").read_bytes() == TABLE
            page.set_viewport_size({"width": 390, "height": 844})
            page.get_by_role("tab", name="Array comparison", exact=True).click()
            panel.locator(".observed-array").first.locator("summary").click()
            assert page.evaluate(
                "document.documentElement.scrollWidth <= window.innerWidth + 1"
            )
            assert panel.locator(".comparison-table").evaluate(
                "element => element.scrollWidth <= element.clientWidth + 1"
            )
            page.evaluate("window.scrollTo(0, 0)")
            page.screenshot(path=output / "mobile.png", full_page=True)
            # A truncated preview must never imply globally absent sharing.
            state["job"] = copy.deepcopy(fixture)
            state["job"]["summary"]["array_compare"]["edges"] = []
            state["job"]["summary"]["array_compare"]["edges_truncated"] = True
            page.reload()
            page.wait_for_load_state("networkidle")
            page.get_by_role("tab", name="Array comparison", exact=True).click()
            expect(
                panel.get_by_text(
                    "No sharing links in this preview; download the full comparison.",
                    exact=True,
                )
            ).to_be_visible()
            state["job"] = copy.deepcopy(fixture)
            state["job"]["summary"]["array_compare"] = None
            page.reload()
            page.wait_for_load_state("networkidle")
            expect(
                page.get_by_role("tabpanel", name="Overview", exact=True).get_by_text(
                    "A completed observed-array comparison was not reported."
                )
            ).to_be_visible()
            assert not errors, errors
            assert not page.evaluate("window.cspViolations")
            assert len(state["submissions"]) == 1
            (output / "acceptance.json").write_text(
                json.dumps(
                    {
                        "status": "passed",
                        "api": "intercepted real native CCTK public evidence",
                        "checks": [
                            "advertisement",
                            "array JSON admission",
                            "incomplete rejection",
                            "source IDs and unknown orientation",
                            "duplicate and deletion slots",
                            "sharing graph/table",
                            "interactive node",
                            "reload recovery",
                            "authenticated download",
                            "390px overflow",
                            "390px pairwise metrics visible without scrolling",
                            "truncated preview semantics",
                            "missing evidence",
                            "no JS/CSP errors",
                        ],
                    },
                    indent=2,
                )
                + "\n"
            )
            print("CCTK browser acceptance passed")
        finally:
            browser.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:4198")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check(args.origin.rstrip("/"), args.output)
