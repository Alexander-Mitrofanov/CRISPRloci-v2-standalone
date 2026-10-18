#!/usr/bin/env python3
"""Compare stored example views with second independent real worker results."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlsplit
import zipfile

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]


def check(origin, work, output):
    output.mkdir(parents=True, exist_ok=True)
    assets = sorted((ROOT / "frontend/public/examples").glob("*/result.json"))
    errors, submissions, checks = [], [], []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            for result in assets:
                key = result.parent.name
                snapshot = json.loads(result.read_text())
                job = json.loads((work / f"{key}-run-2.json").read_text())
                assert job["mode"] == snapshot["job"]["mode"]
                bundle = work / "data/jobs" / job["job_id"] / "exports/attempt-001/results.zip"
                context = browser.new_context(accept_downloads=True)

                def api(route):
                    path = urlsplit(route.request.url).path
                    headers = {"Access-Control-Allow-Origin": origin,
                               "Access-Control-Allow-Headers": "authorization,content-type"}
                    if route.request.method == "OPTIONS":
                        route.fulfill(status=204, headers=headers)
                        return
                    if route.request.method == "POST":
                        submissions.append(key)
                        raise AssertionError("Precomputed examples must never submit a job")
                    if path.endswith("/health"):
                        value = {"version": "2.11.1"}
                    elif path.endswith("/config"):
                        value = {"modes": [], "max_records": 20, "max_total_bases": 2_000_000,
                                 "max_record_bases": 2_000_000, "max_request_bytes": 2_000_000}
                    elif path.endswith("/jobs/" + job["job_id"]):
                        assert route.request.headers["authorization"] == "Bearer " + "a" * 43
                        value = job
                    elif f"/jobs/{job['job_id']}/artifacts/" in path:
                        assert route.request.headers["authorization"] == "Bearer " + "a" * 43
                        artifact = next(item for item in job["artifacts"]
                                        if path.endswith("/" + item["artifact_id"]))
                        with zipfile.ZipFile(bundle) as archive:
                            payload = archive.read("results/" + artifact["name"])
                        route.fulfill(headers=headers, content_type=artifact["media_type"], body=payload)
                        return
                    else:
                        raise AssertionError("Unexpected API read: " + path)
                    route.fulfill(headers=headers, content_type="application/json", body=json.dumps(value))

                context.route("https://analysis.example.org/**", api)
                stored, live = context.new_page(), context.new_page()
                for page in (stored, live):
                    page.on("pageerror", lambda error: errors.append(str(error)))
                stored.goto(f"{origin}/?method={job['mode']}")
                stored.wait_for_load_state("networkidle")
                if key == "repeats-rna":
                    stored.locator('input[name="repeat-molecule"][value="RNA"]').check()
                if key == "tracrrna-v":
                    stored.locator("details.advanced > summary").click()
                    stored.locator('input[name="tracr-model"][value="V"]').check()
                stored.get_by_role("button", name="Load example", exact=True).click()
                view = stored.get_by_role("button", name="View precomputed result", exact=True)
                expect(view).to_be_enabled()
                view.click()
                expect(stored.locator(".results")).to_be_visible()
                assert "#job=" not in stored.url
                expect(stored.get_by_role("link", name="Download example input")).to_be_visible()
                live.goto(f"{origin}/?method={job['mode']}#job={job['job_id']}.{'a'*43}")
                expect(live.locator(".results")).to_be_visible()
                live.wait_for_load_state("networkidle")
                for page in (stored, live):
                    page.evaluate("document.fonts.ready")
                    # These controls require a private live job capability; compare
                    # the scientific presentation rather than authentication UI.
                    page.locator(".results .downloads, .results .sequence-handoff").evaluate_all(
                        "items => items.forEach(item => item.remove())"
                    )
                    # Example notices and live progress have different heights.
                    # Capture the unchanged result at a common layout origin so
                    # subpixel clipping does not shift rasterized text by a pixel.
                    page.locator(".results").evaluate(
                        "element => document.getElementById('main-content').replaceChildren(element)"
                    )
                tabs = stored.get_by_role("tab").all_text_contents()
                assert tabs == live.get_by_role("tab").all_text_contents()
                for width, height in ((1440, 1000), (390, 844)):
                    for page in (stored, live):
                        page.set_viewport_size({"width": width, "height": height})
                    for index, label in enumerate(tabs):
                        for page in (stored, live):
                            page.get_by_role("tab", name=label.strip(), exact=True).click()
                            page.locator('.result-tab-panel:not([hidden]) details').evaluate_all(
                                "items => items.forEach(item => item.open = true)"
                            )
                            # The actual duration remains in the product. Only this
                            # test normalizes nondeterministic timing cells.
                            page.locator('.reconstruction-table td:last-child').evaluate_all(
                                "items => items.forEach(item => item.textContent = '[run duration]')"
                            )
                            page.evaluate("document.activeElement?.blur()")
                        expected = stored.locator(".results").inner_text()
                        actual = live.locator(".results").inner_text()
                        assert expected == actual, f"Scientific text differs: {key}/{label}/{width}"
                        paths = [output / f"{key}-{index}-{width}-{kind}.png" for kind in ("stored", "independent")]
                        images = [page.locator(".results").screenshot(path=str(path), animations="disabled")
                                  for page, path in zip((stored, live), paths)]
                        assert images[0] == images[1], f"Scientific pixels differ: {key}/{label}/{width}; inspect {output}"
                        checks.append({"example": key, "tab": label.strip(), "width": width,
                                       "sha256": hashlib.sha256(images[0]).hexdigest()})
                    for page in (stored, live):
                        assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth"), f"Overflow: {key}/{width}"
                context.close()
                print(f"PASS {key}: independent-run scientific pixels, all tabs, desktop/mobile", flush=True)
            assert not errors, errors
            assert not submissions, submissions
        finally:
            browser.close()
    (output / "acceptance.json").write_text(json.dumps({
        "status": "passed", "examples": len(assets), "comparisons": checks,
        "exclusions": ["private authenticated download/continuation controls", "execution duration cells"],
        "checks": ["real second independent job", "pixel-identical scientific tabs", "390px overflow",
                   "no example job submissions", "no JavaScript exceptions"],
    }, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:5173")
    parser.add_argument("--work", type=Path, default=ROOT / ".runtime/example-generation")
    parser.add_argument("--output", type=Path, default=ROOT / ".runtime/browser-examples")
    args = parser.parse_args()
    check(args.origin.rstrip("/"), args.work, args.output)
