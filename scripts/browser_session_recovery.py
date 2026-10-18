#!/usr/bin/env python3
"""Check session navigation and outage recovery against an isolated local test API.

Build the frontend with https://analysis.example.org as its API origin. Start a
preview and a disposable API/fixture worker before running this script. Only the
test HTTPS origin is forwarded to the supplied local API URL.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright


def main(origin: str, api_url: str, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    checks = []
    errors = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            context = browser.new_context(
                viewport={"width": 1440, "height": 1000}, accept_downloads=True,
            )
            context.on("page", lambda page: page.on("pageerror", lambda error: errors.append(str(error))))
            context.add_init_script(
                "window.cspViolations=[]; document.addEventListener('securitypolicyviolation', "
                "e => window.cspViolations.push(e.violatedDirective));"
            )

            def proxy(route):
                parts = urlsplit(route.request.url)
                target = api_url.rstrip("/") + parts.path + ("?" + parts.query if parts.query else "")
                route.fulfill(response=route.fetch(url=target))

            context.route("https://analysis.example.org/**", proxy)
            page = context.new_page()
            page.goto(origin + "/?method=orientation")
            page.wait_for_load_state("networkidle")
            expect(page.get_by_text("Service ready", exact=True)).to_be_visible()
            page.get_by_role("button", name="Load example", exact=True).click()
            page.get_by_role("button", name="View precomputed result", exact=True).click()
            expect(page.get_by_role("heading", name="Example results", exact=True)).to_be_visible()

            page.get_by_role("link", name="Methods", exact=True).click()
            page.get_by_role("link", name="Orientation-aware evolution", exact=True).click()
            history_length = page.evaluate("history.length")
            page.get_by_role("button", name="View precomputed result", exact=True).click()
            expect(page.get_by_role("heading", name="Example results", exact=True)).to_be_visible()
            assert page.evaluate("history.length") == history_length
            checks.append("reopening an existing example reveals its results without adding history")

            page.get_by_role("link", name="Methods", exact=True).click()
            page.get_by_role("link", name="Detect arrays", exact=True).click()
            page.locator("#fasta-input").fill(">a\nACGT\n>b\nACGA\n")
            history_length = page.evaluate("history.length")
            with page.expect_response(
                lambda response: response.request.method == "POST"
                and response.url.endswith("/api/v1/jobs")
            ) as submitted:
                page.get_by_role("button", name="Compute", exact=True).click()
            assert submitted.value.status == 202
            expect(page.locator("#job-heading")).to_be_visible()
            expect(page.locator("#fasta-input")).to_have_count(0)
            assert page.evaluate("history.length") == history_length
            assert "#job=" in page.url
            recovery = page.url
            expect(page.get_by_role("heading", name="Analysis complete", exact=True)).to_be_visible(timeout=30000)
            assert not page.evaluate("window.cspViolations")
            checks.append("submission after example opens job status and retains recovery without adding history")
            page.screenshot(path=str(output / "job-after-example.png"), full_page=True)

            page.get_by_role("tab", name="Files & methods", exact=True).click()
            with page.expect_download() as download:
                page.get_by_role("button", name="Download complete result bundle", exact=False).click()
            download.value.save_as(output / "fixture-result.zip")
            page.reload()
            expect(page.get_by_role("heading", name="Analysis complete", exact=True)).to_be_visible()
            assert page.url == recovery
            checks.append("authenticated download and reload recovery")
            page.get_by_role("link", name="Methods", exact=True).click()
            page.get_by_role("button", name="Return to current results", exact=True).click()
            expect(page.locator("#job-heading")).to_be_visible()
            checks.append("explicit return to current results")

            offline = context.new_page()
            requests = []
            healthy = False

            def health(route):
                requests.append(route.request.url)
                if healthy:
                    proxy(route)
                else:
                    route.fulfill(
                        status=503, json={"detail": "temporary outage"},
                        headers={"Access-Control-Allow-Origin": origin},
                    )

            offline.route("**/api/v1/health", health)
            offline.clock.install()
            offline.goto(origin + "/?method=detection")
            offline.wait_for_load_state("networkidle")
            expect(offline.get_by_text("Service unavailable", exact=True)).to_be_visible()
            draft = ">a\nACGT"
            offline.locator("#fasta-input").fill(draft)
            expect(offline.get_by_role("button", name="Compute", exact=True)).to_be_disabled()
            healthy = True
            offline.clock.fast_forward(10_001)
            expect(offline.get_by_text("Service ready", exact=True)).to_be_visible()
            expect(offline.get_by_role("button", name="Compute", exact=True)).to_be_enabled()
            expect(offline.locator("#fasta-input")).to_have_value(draft)
            assert len(requests) == 2
            offline.clock.fast_forward(60_000)
            assert len(requests) == 2
            checks.append("automatic recovery within ten seconds preserves draft and stops retrying after success")
            for width in (1440, 390):
                offline.set_viewport_size({"width": width, "height": 900})
                assert offline.evaluate("document.documentElement.scrollWidth <= innerWidth")
                offline.screenshot(path=str(output / f"recovered-{width}.png"), full_page=True)
            assert not errors, errors
            assert not page.evaluate("window.cspViolations")
            assert not offline.evaluate("window.cspViolations")
            checks.extend(["390px layout", "no JavaScript exceptions or production CSP violations"])
        finally:
            browser.close()
    (output / "acceptance.json").write_text(json.dumps({"status": "passed", "checks": checks}, indent=2) + "\n")
    print("Session and service recovery browser checks passed:", output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:14173")
    parser.add_argument("--api-url", default="http://127.0.0.1:18080")
    parser.add_argument("--output", type=Path, default=Path("output/playwright/session-recovery"))
    args = parser.parse_args()
    main(args.origin.rstrip("/"), args.api_url, args.output)
