#!/usr/bin/env python3
"""Check compiled UI and its real CSP, forwarding only its test API origin locally."""

from pathlib import Path
import json
from urllib.parse import urlsplit
from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]


def main():
    failures = []
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1100})
        page.on("pageerror", lambda error: failures.append(str(error)))
        page.add_init_script(
            "window.cspViolations=[]; document.addEventListener('securitypolicyviolation', e => window.cspViolations.push(e.violatedDirective));"
        )

        def proxy(route):
            path = urlsplit(route.request.url).path
            response = route.fetch(url="http://127.0.0.1:8000" + path)
            route.fulfill(response=response)

        page.route("https://analysis.example.org/**", proxy)
        page.goto("http://127.0.0.1:4173")
        page.wait_for_load_state("networkidle")
        expect(page.get_by_text("Service ready", exact=True)).to_be_visible()
        policy = page.locator(
            'meta[http-equiv="Content-Security-Policy"]'
        ).get_attribute("content")
        assert "style-src-attr 'none'" in policy and "script-src 'self'" in policy
        page.get_by_role("link", name="Orientation-aware evolution", exact=True).click()
        page.get_by_role("button", name="Load example").click()
        page.get_by_role("button", name="View precomputed result").click()
        page.get_by_role("tab", name="Orientation", exact=True).click()
        expect(page.locator("#orientation-heading")).to_be_visible(timeout=30000)
        page.get_by_role("tab", name="History", exact=True).click()
        expect(page.locator("#reconstruction-heading")).to_be_visible()
        assert not page.evaluate("window.cspViolations"), page.evaluate(
            "window.cspViolations"
        )
        # Exercise the new annotation component under the same compiled policy.
        page.goto("http://127.0.0.1:4173")
        page.wait_for_load_state("networkidle")
        page.get_by_role("link", name="Annotate a CRISPR locus", exact=True).click()
        page.locator("input[type=file]").set_input_files(
            ROOT / "tools/CRISPRtracrRNA-v3/examples/complete.fasta"
        )
        page.get_by_role("button", name="Compute", exact=True).click()
        expect(
            page.get_by_role("heading", name="Locus annotations", exact=True)
        ).to_be_visible(timeout=180000)
        page.get_by_role("tab", name="Cas systems", exact=True).click()
        assert (
            page.get_by_role("region", name="Cas gene table").get_by_role("row").count()
            > 1
        )
        assert not page.evaluate("window.cspViolations"), page.evaluate(
            "window.cspViolations"
        )
        assert not failures, failures
        page.get_by_role("tab", name="Overview", exact=True).click()
        page.locator("#annotation-overview-heading").evaluate(
            "e=>e.scrollIntoView({block:'start'})"
        )
        page.screenshot(
            path=str(ROOT / ".runtime/browser-evidence/production-locus.png")
        )
        browser.close()
    (ROOT / ".runtime/browser-evidence/production-csp.json").write_text(
        json.dumps(
            {
                "status": "passed",
                "api_transport": "Test-only Playwright forwarding from the configured HTTPS origin to the local API",
                "checks": [
                    "production CSP",
                    "cached evolutionary results",
                    "real locus submission and results",
                    "no CSP violations",
                    "no JavaScript exceptions",
                ],
            },
            indent=2,
        )
        + "\n"
    )
    print("Production browser/CSP checks passed")


if __name__ == "__main__":
    main()
