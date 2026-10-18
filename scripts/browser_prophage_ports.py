#!/usr/bin/env python3
"""Browser acceptance for validated-derived Phigaro typed ports; intercepted API only."""
import argparse
import copy
import json
from pathlib import Path
from urllib.parse import urlsplit
from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]
TOKEN = "a" * 43
JOB = "d" * 32


def check(origin, output):
    output.mkdir(parents=True, exist_ok=True)
    documents = json.loads((ROOT / "frontend/tests/fixtures/prophage-ports.json").read_text())
    payloads = {name: json.dumps(doc) for name, doc in documents.items()}
    for name in ("repeats.json", "spacers.json"):
        assert documents[name]["complete"] is True
        assert all(row["origin"]["array_id"] == "native_array_1" for row in documents[name]["records"])
    job = {"job_id": JOB, "mode": "prophage_context", "status": "completed", "options": {}, "summary": {},
           "artifacts": [{"artifact_id": name, "name": name, "size_bytes": len(payload), "media_type": "application/json"} for name, payload in payloads.items()]}
    errors, submissions, downloads = [], [], []
    recovery = origin + "/#job=" + JOB + "." + TOKEN
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 1440, "height": 1000}, accept_downloads=True)
            page.on("pageerror", lambda error: errors.append(str(error)))
            def api(route):
                req = route.request
                path = urlsplit(req.url).path
                headers = {"Access-Control-Allow-Origin": origin, "Access-Control-Allow-Headers": "authorization,content-type"}
                if req.method == "OPTIONS":
                    route.fulfill(status=204, headers=headers)
                    return
                if path.endswith("/health"):
                    data = {"version": "2.13.0"}
                elif path.endswith("/config"):
                    data = {"modes": ["prophage_context", "array_compare", "repeats", "repeat_map", "repeat_type", "viral_search", "spacer_association"], "max_total_bases": 100000, "max_record_bases": 10000, "max_records": 100, "max_request_bytes": 100000}
                elif path.endswith("/jobs") and req.method == "POST":
                    submissions.append(req.post_data_json)
                    raise AssertionError("Continuation must prepare a draft without enqueuing")
                elif "/jobs/" in path:
                    assert req.headers.get("authorization") == "Bearer " + TOKEN and TOKEN not in req.url
                    if "/artifacts/" in path:
                        name = path.rsplit("/", 1)[-1]
                        downloads.append(name)
                        route.fulfill(status=200, headers=headers, content_type="application/json", body=payloads[name])
                        return
                    data = job
                else:
                    raise AssertionError(path)
                route.fulfill(status=200, headers=headers, content_type="application/json", body=json.dumps(data))
            page.route("https://analysis.example.org/**", api)
            for name, mode in [("arrays.json", "array_compare"), ("repeats.json", "repeat_type"), ("spacers.json", "viral_search")]:
                page.goto(recovery)
                page.get_by_role("tab", name="Files & methods", exact=True).click()
                options = page.get_by_label("Sequence file and next analysis")
                assert options.locator(f'option[value="{name}:{mode}"]').count() == 1
                assert options.locator(f'option[value="{name}:rna_family"]').count() == 0
                with page.expect_download() as event:
                    page.get_by_role("button", name=name + " application/json", exact=True).click()
                saved = output / name
                event.value.save_as(saved)
                assert saved.read_text() == payloads[name]
                options.select_option(name + ":" + mode)
                page.get_by_role("button", name="Use as analysis input", exact=True).click()
                expect(page.get_by_role("textbox")).to_have_value(payloads[name])
                assert page.evaluate("location.hash") == "#job=" + JOB + "." + TOKEN
                assert not submissions
                # Source search accepts 16 bases; viral eligibility stays an independent 18–80 nt rule.
                if name == "spacers.json":
                    assert len(documents[name]["records"][0]["sequence"]) == 16
                    expect(page.get_by_text("18–80", exact=False).first).to_be_visible()
                page.set_viewport_size({"width": 390, "height": 844})
                assert page.evaluate("document.documentElement.scrollWidth <= innerWidth + 1")
                page.screenshot(path=output / (name + ".png"), full_page=True)
                page.set_viewport_size({"width": 1440, "height": 1000})
            incomplete = copy.deepcopy(documents["repeats.json"])
            incomplete["complete"] = False
            payloads["repeats.json"] = json.dumps(incomplete)
            page.goto(recovery)
            page.get_by_role("tab", name="Files & methods", exact=True).click()
            page.get_by_label("Sequence file and next analysis").select_option("repeats.json:repeat_type")
            page.get_by_role("button", name="Use as analysis input", exact=True).click()
            expect(page.get_by_role("alert").filter(has_text="Incomplete sequence sets")).to_be_visible()
            assert not submissions and not errors
            assert set(downloads) == {"arrays.json", "repeats.json", "spacers.json"}
            (output / "acceptance.json").write_text(json.dumps({"status": "passed", "evidence": "typed exports built from both-gate-validated prophage-spacer fixture", "checks": ["protected_downloads", "origin_and_complete_flags", "compatible_destinations", "ineligible_viral_query_remains_short", "draft_without_enqueue", "incomplete_rejection", "mobile"]}, indent=2) + "\n")
        finally:
            browser.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check(args.origin.rstrip("/"), args.output)
