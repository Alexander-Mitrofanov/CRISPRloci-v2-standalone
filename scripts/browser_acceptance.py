#!/usr/bin/env python3
"""Browser acceptance against a locally running frontend and gated backend."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from urllib.parse import urlsplit
from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]


def main(origin, output, api_url=None):
    output.mkdir(parents=True, exist_ok=True)
    failures = []
    completed = {}
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        context = browser.new_context(
            viewport={"width": 1440, "height": 1100}, accept_downloads=True
        )
        if api_url:
            def proxy(route):
                parts = urlsplit(route.request.url)
                target = api_url.rstrip("/") + parts.path + ("?" + parts.query if parts.query else "")
                route.fulfill(response=route.fetch(url=target))

            context.route("https://analysis.example.org/**", proxy)
        page = context.new_page()
        page.on("pageerror", lambda error: failures.append(str(error)))

        def capture(response):
            if response.request.method == "GET" and re.search(r"/api/v1/jobs/[0-9a-f]{32}$", response.url):
                value = response.json()
                if value.get("status") == "completed":
                    completed[value["job_id"]] = value

        page.on("response", capture)
        page.goto(origin)
        page.wait_for_load_state("networkidle")
        expect(page.get_by_role("heading", name="Choose a method")).to_be_visible()
        expect(page.get_by_text("Service ready", exact=True)).to_be_visible(
            timeout=30000
        )
        page.screenshot(path=str(output / "home-desktop.png"), full_page=True)
        page.get_by_role("link", name="Annotate a CRISPR locus", exact=True).click()
        page.locator("input[type=file]").set_input_files(
            ROOT / "tools/CRISPRtracrRNA-v3/examples/complete.fasta"
        )
        button = page.get_by_role("button", name="Compute", exact=True)
        expect(button).to_be_enabled()
        with page.expect_response(
            lambda response: (
                response.request.method == "POST"
                and response.url.endswith("/api/v1/jobs")
            )
        ) as submitted:
            button.click()
        assert submitted.value.status == 202
        job_id = submitted.value.json()["job"]["job_id"]
        expect(
            page.get_by_role("heading", name="Locus annotations", exact=True)
        ).to_be_visible(timeout=180000)
        summary = completed[job_id]["summary"]
        cas_count = len(summary["casandra"]["genes"])
        rna_count = len(summary["tracrrna"]["candidates"])
        assert cas_count > 0 and rna_count > 0
        assert summary["tracrrna"]["truncated"] is False
        assert rna_count == summary["tracrrna"]["prediction_count"]
        expect(page.get_by_text("Cas protein calls", exact=True).first).to_be_visible()
        page.get_by_role("tab", name="Cas systems", exact=True).click()
        expect(
            page.get_by_role("region", name="Cas gene table").get_by_role("row")
        ).to_have_count(cas_count + 1)
        page.get_by_role("tab", name="tracrRNA", exact=True).click()
        expect(
            page.get_by_role("region", name="tracrRNA candidate table").get_by_role(
                "row"
            )
        ).to_have_count(rna_count + 1)
        page.get_by_role("tab", name="Files & methods", exact=True).click()
        page.locator("summary").filter(has_text="Methods & provenance").click()
        expect(
            page.get_by_role("link", name=re.compile("CasAndra revision"))
        ).to_have_attribute("href", re.compile("/commit/[0-9a-f]{40}$"))
        page.get_by_role("tab", name="Overview", exact=True).click()
        recovery_url = page.url
        assert "#job=" in recovery_url
        page.get_by_role(
            "heading", name="Locus annotations", exact=True
        ).scroll_into_view_if_needed()
        page.screenshot(path=str(output / "locus-desktop.png"))
        page.get_by_role("tab", name="Files & methods", exact=True).click()
        with page.expect_download() as archive:
            page.get_by_role(
                "button", name=re.compile("Download complete result bundle")
            ).click()
        assert archive.value.suggested_filename.endswith(".zip")
        archive.value.save_as(output / "browser-result.zip")
        page.reload()
        expect(
            page.get_by_role("heading", name="Locus annotations", exact=True)
        ).to_be_visible(timeout=30000)
        assert page.url == recovery_url
        page.set_viewport_size({"width": 390, "height": 844})
        page.get_by_role(
            "heading", name="Locus annotations", exact=True
        ).scroll_into_view_if_needed()
        assert page.evaluate(
            "document.documentElement.scrollWidth <= window.innerWidth"
        )
        page.screenshot(path=str(output / "locus-mobile.png"))
        fresh = context.new_page()
        fresh.set_viewport_size({"width": 390, "height": 844})
        fresh.goto(origin)
        fresh.wait_for_load_state("networkidle")
        assert fresh.evaluate(
            "document.documentElement.scrollWidth <= window.innerWidth"
        )
        fresh.screenshot(path=str(output / "home-mobile.png"), full_page=True)
        assert not failures, failures
        browser.close()
    (output / "acceptance.json").write_text(
        json.dumps(
            {
                "status": "passed",
                "job_id": job_id,
                "cas_rows": cas_count,
                "tracrrna_rows": rna_count,
                "checks": [
                    "real locus submission",
                    "Cas rows match validated API evidence",
                    "tracrRNA rows match validated API evidence",
                    "repository provenance links",
                    "authenticated browser download",
                    "reload recovery",
                    "390px overflow",
                    "no JavaScript exceptions",
                ],
            },
            indent=2,
        )
        + "\n"
    )
    print("Browser acceptance passed:", output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:5173")
    parser.add_argument("--api-url", help="Forward the test HTTPS API origin to this isolated local API")
    parser.add_argument(
        "--output", type=Path, default=ROOT / ".runtime/browser-evidence"
    )
    args = parser.parse_args()
    main(args.origin, args.output, args.api_url)
