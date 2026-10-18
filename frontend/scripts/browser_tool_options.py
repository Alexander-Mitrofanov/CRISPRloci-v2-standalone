#!/usr/bin/env python3
"""Verify compact CLI controls in a local browser with intercepted API requests."""

import argparse
import json
from pathlib import Path
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]
CATALOG = json.loads((ROOT / "src/contracts/tool-options-v1.json").read_text())
MODES = sorted({mode for tool in CATALOG["tools"] for mode in tool["modes"]})
FIXED_MODES = {
    "locus_compare", "prophage_context", "insertion_context", "defense_context", "rna_family",
}


def check(origin, output):
    output.mkdir(parents=True, exist_ok=True)
    submissions = []
    errors = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            page = browser.new_page(viewport={"width": 1440, "height": 1050})
            page.on("pageerror", lambda error: errors.append(str(error)))

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
                    data = {"version": "2.6.0"}
                elif path.endswith("/config"):
                    data = {
                        "modes": MODES,
                        "max_total_bases": 20000,
                        "max_record_bases": 20000,
                        "max_records": 1000,
                        "max_request_bytes": 100000,
                    }
                elif path.endswith("/jobs") and request.method == "POST":
                    submissions.append(request.post_data_json)
                    route.fulfill(
                        status=503,
                        headers=headers,
                        content_type="application/json",
                        body=json.dumps({"detail": "Acceptance fixture captured the request."}),
                    )
                    return
                else:
                    raise AssertionError(f"Unexpected request: {path}")
                route.fulfill(
                    status=200,
                    headers=headers,
                    content_type="application/json",
                    body=json.dumps(data),
                )

            page.route("https://analysis.example.org/**", api)

            def visit(mode):
                page.goto(f"{origin}/?method={mode}")
                page.wait_for_load_state("networkidle")
                if mode in FIXED_MODES:
                    expect(page.locator(".advanced")).to_have_count(0)
                    return False
                expect(page.locator(".advanced")).not_to_have_attribute("open", "")
                expect(page.locator(".tool-options-tool[open]")).to_have_count(0)
                return True

            def open_tool(name):
                root = page.locator(".advanced > summary")
                root.click()
                page.locator(".tool-options-tool > summary").filter(has_text=name).click()

            for mode in MODES:
                if visit(mode):
                    assert page.locator(".tool-options-tool").count() > 0
            visit("detection")
            page.locator("#fasta-input").fill(">a\nACGTACGT\n")
            # Native disclosures must also work with the keyboard.
            page.locator(".advanced > summary").focus()
            page.keyboard.press("Enter")
            tool_summary = page.locator(".tool-options-tool > summary").filter(
                has_text="CRISPRidentify"
            )
            tool_summary.focus()
            page.keyboard.press("Enter")
            repeat_count = page.get_by_label("Minimum repeat count --min_repeats")
            expect(repeat_count).to_be_visible()
            repeat_count.fill("1")
            expect(repeat_count).to_have_attribute("aria-invalid", "true")
            expect(page.get_by_role("button", name="Compute", exact=True)).to_be_disabled()
            repeat_count.fill("5")
            expect(page.get_by_text("1 tool flag changed", exact=True)).to_be_visible()
            page.get_by_label("Minimum repeat length --min_len_rep").fill("60")
            expect(page.get_by_text("Must not exceed --max_len_rep (55).")).to_be_visible()
            expect(page.get_by_role("button", name="Compute", exact=True)).to_be_disabled()
            page.get_by_role("button", name="Reset --min_len_rep to default").click()
            page.get_by_label("Fast search --fast_run", exact=True).select_option("false")
            expect(page.locator("#tool-crispridentify-fast_run_seed_profile")).to_have_count(0)
            page.get_by_role("button", name="Reset --fast_run to default").click()
            expect(page.locator("#tool-crispridentify-fast_run_seed_profile")).to_be_visible()
            page.locator(".advanced").screenshot(path=output / "tool-options-desktop.png")
            with page.expect_response(lambda response: response.url.endswith("/jobs")):
                page.get_by_role("button", name="Compute", exact=True).click()
            assert submissions[-1]["tool_options"] == {
                "crispridentify": {"--min_repeats": 5}
            }

            visit("repeat_context")
            open_tool("CRISPRrepeat")
            page.get_by_label("Local window length --window").fill("50")
            expect(page.get_by_text("Must not exceed --window (50).")).to_be_visible()
            page.get_by_label("Maximum pairing span --span").fill("40")
            expect(page.get_by_text("Must not exceed --window (50).")).to_have_count(0)

            visit("repeats")
            open_tool("CRISPRrepeat")
            expect(page.locator("#tool-crisprrepeat-window")).to_have_count(0)
            page.get_by_label("Folding temperature (°C) --temperature").fill("25")
            page.set_viewport_size({"width": 390, "height": 844})
            page.locator(".advanced").scroll_into_view_if_needed()
            assert page.evaluate("document.documentElement.scrollWidth <= innerWidth + 1")
            page.screenshot(path=output / "tool-options-mobile.png", full_page=True)
            page.locator(".advanced > summary").click()
            expect(page.get_by_text("1 tool flag changed", exact=True)).to_be_visible()
            page.screenshot(path=output / "tool-options-collapsed-mobile.png", full_page=True)

            visit("repeat_map")
            open_tool("CRISPRmap")
            page.get_by_label("Maximum edit distance --max-distance").fill("0")
            page.get_by_label("Displayed reference hits --top-k").fill("5")
            page.get_by_role("button", name="Reset CRISPRmap flags", exact=True).click()
            expect(page.get_by_label("Maximum edit distance --max-distance")).to_have_value("3")
            expect(page.get_by_label("Displayed reference hits --top-k")).to_have_value("20")
            assert not errors, errors
            (output / "acceptance.json").write_text(
                json.dumps(
                    {
                        "status": "passed",
                        "api": "intercepted local fixtures; no scientific jobs submitted",
                        "checks": [
                            "collapsed configurable menus and no menu for fixed workflows",
                            "keyboard disclosure operation",
                            "exact CLI flags and changed defaults",
                            "range and cross-field validation",
                            "conditional and mode-specific controls",
                            "submitted CLI override payload",
                            "per-field and per-tool reset",
                            "390px layout and collapsed changes count",
                            "no JavaScript exceptions",
                        ],
                    },
                    indent=2,
                ) + "\n"
            )
            print("Tool options browser acceptance passed:", output)
        finally:
            browser.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:5178")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check(args.origin.rstrip("/"), args.output)
