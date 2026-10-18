#!/usr/bin/env python3
"""Local JSON upload/download acceptance; intercepted API, no live jobs."""

import argparse
import json
from pathlib import Path
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]


def check(origin, output):
    output.mkdir(parents=True, exist_ok=True)
    packet = {
        "schema": {"name": "crisprloci.sequence-set", "version": "1.0.0"},
        "kind": "repeat",
        "molecule": "DNA",
        "complete": True,
        "records": [{"id": "repeat_1", "sequence": "ACGTCGATGCTAGCTAGTCGATCG"}],
    }
    payload = json.dumps(packet)
    compact_payload = (
        ROOT / "docs/examples/module-results/repeat_map.json"
    ).read_text()
    job = json.loads((ROOT / "frontend/tests/fixtures/repeat-map.json").read_text())
    job["artifacts"] = [
        {
            "artifact_id": "a" * 24,
            "name": "input-sequences.json",
            "media_type": "application/json",
            "size_bytes": len(payload),
        },
        {
            "artifact_id": "b" * 24,
            "name": "result.json",
            "media_type": "application/json",
            "size_bytes": len(compact_payload),
        },
    ]
    submissions, errors = [], []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            context = browser.new_context(
                viewport={"width": 1440, "height": 1000}, accept_downloads=True
            )
            page = context.new_page()
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
                        "modes": [
                            "repeat_map",
                            "repeats",
                            "viral_search",
                            "protospacer",
                            "array_compare",
                        ],
                        "max_total_bases": 20000,
                        "max_record_bases": 20000,
                        "max_records": 100,
                        "max_request_bytes": 100000,
                    }
                elif path.endswith("/jobs") and request.method == "POST":
                    submissions.append(request.post_data_json)
                    assert request.post_data_json["sequence"] == payload
                    data = {"job": job, "access_token": "a" * 43}
                elif "/jobs/" in path:
                    assert request.headers.get("authorization") == "Bearer " + "a" * 43
                    if "/artifacts/" in path:
                        route.fulfill(
                            status=200,
                            headers=headers,
                            content_type="application/json",
                            body=compact_payload
                            if path.endswith("b" * 24)
                            else payload,
                        )
                        return
                    data = job
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
            page.get_by_role("link", name="Compare observed arrays", exact=True).click()
            arrays = json.loads((ROOT / "frontend/tests/fixtures/array-compare-input.json").read_text())
            for width in (1440, 390):
                page.set_viewport_size({"width": width, "height": 1000})
                for kind in ("array", "spacer"):
                    malformed = json.loads(json.dumps(arrays))
                    if kind == "array":
                        malformed["arrays"][0] = None
                    else:
                        malformed["arrays"][0]["spacers"][0] = None
                    page.get_by_role("textbox").fill(json.dumps(malformed))
                    expect(page.locator("#fasta-errors")).to_be_visible()
                    expect(page.get_by_role("button", name="Compute", exact=True)).to_be_disabled()
                    page.get_by_role("textbox").fill(json.dumps(arrays))
                    expect(page.get_by_role("button", name="Compute", exact=True)).to_be_enabled()
            page.set_viewport_size({"width": 1440, "height": 1000})
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            page.get_by_role(
                "link", name="Map repeats to reference families", exact=True
            ).click()
            upload = page.locator("input[type=file]")
            upload.set_input_files(
                {
                    "name": "spacers.json",
                    "mimeType": "application/json",
                    "buffer": json.dumps({**packet, "kind": "spacer"}).encode(),
                }
            )
            expect(
                page.locator("#fasta-errors").get_by_text(
                    "This analysis requires repeat sequences, not spacer sequences."
                )
            ).to_be_visible()
            expect(
                page.get_by_role("button", name="Compute", exact=True)
            ).to_be_disabled()
            page.get_by_role("textbox").fill(json.dumps({**packet, "complete": False}))
            expect(
                page.locator("#fasta-errors").get_by_text(
                    "Incomplete sequence sets cannot be used as module inputs."
                )
            ).to_be_visible()
            upload.set_input_files(
                {
                    "name": "repeats.json",
                    "mimeType": "application/json",
                    "buffer": payload.encode(),
                }
            )
            expect(
                page.get_by_role("button", name="Compute", exact=True)
            ).to_be_enabled()
            page.screenshot(path=output / "upload-desktop.png", full_page=True)
            page.get_by_role("button", name="Compute", exact=True).click()
            page.get_by_role("tab", name="Files & methods", exact=True).click()
            expect(
                page.get_by_text("Continue an analysis:", exact=False)
            ).to_be_visible()
            expect(page.locator(".artifact-grid button").first).to_contain_text(
                "result.json"
            )
            with page.expect_download() as compact_download:
                page.get_by_role(
                    "button", name="result.json Compact scientific result · JSON"
                ).click()
            compact_saved = output / "result.json"
            compact_download.value.save_as(compact_saved)
            assert compact_saved.read_text() == compact_payload
            with page.expect_download() as download:
                page.get_by_role(
                    "button", name="input-sequences.json application/json"
                ).click()
            saved = output / "input-sequences.json"
            download.value.save_as(saved)
            assert saved.read_text() == payload
            page.set_viewport_size({"width": 390, "height": 844})
            assert page.evaluate(
                "document.documentElement.scrollWidth <= window.innerWidth + 1"
            )
            page.screenshot(path=output / "downloads-mobile.png", full_page=True)
            previous_hash = page.evaluate("location.hash")
            page.get_by_label("Sequence file and next analysis").select_option("a" * 24 + ":repeats")
            page.get_by_role("button", name="Use as analysis input").click()
            expect(page.get_by_role("textbox")).to_have_value(payload)
            assert page.evaluate("location.hash") == previous_hash
            assert len(submissions) == 1
            expect(page.get_by_role("button", name="Compute", exact=True)).to_be_enabled()
            page.screenshot(path=output / "continue-mobile.png", full_page=True)
            assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth + 1")
            page.get_by_role("button", name="Compute", exact=True).click()
            page.get_by_role("tab", name="Files & methods", exact=True).wait_for()
            assert [s["mode"] for s in submissions] == ["repeat_map", "repeats"]
            assert not errors, errors
            print(
                "JSON browser acceptance passed: malformed array recovery, upload, incompatible/incomplete rejection, submission, compact result, authenticated downloads, direct continuation, recovery preservation, mobile layout"
            )
        finally:
            browser.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check(args.origin.rstrip("/"), args.output)
