#!/usr/bin/env python3
"""Browser acceptance with synthetic jobs and the real catalogue API via TestClient.

Run against a local frontend; this never submits jobs or starts scientific tools.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sqlite3
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from urllib.parse import urlsplit

from fastapi import FastAPI
from fastapi.testclient import TestClient
from playwright.sync_api import expect, sync_playwright

from catalog_snapshot import add_indexes, create_tables
from crispr_web.http.routes.catalog import build_catalog_router
from crispr_web.summaries.detection import _array_view

ROOT = Path(__file__).resolve().parents[1]
JOB_ID = "0123456789abcdef0123456789abcdef"


def fixture(directory):
    report = json.loads((ROOT / "backend/tests/fixtures/repeats/identify.json").read_text())
    array = _array_view(report["arrays"][0], report["source"]["id"])
    repeat = array["repeats"][0]["sequence"]
    spacer = array["spacers"][0]["sequence"]
    # A genuine non-exact neighbour plus a reverse-complement spacer.
    near = ("A" if repeat[0] != "A" else "C") + repeat[1:]
    reverse = spacer.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    path = directory / "catalog.sqlite3"
    with sqlite3.connect(path) as connection:
        create_tables(connection)
        summary = {
            "available": True, "snapshot_id": "browser-sequence-search",
            "counts": {"assemblies": 1, "arrays": 1, "repeats": 26, "spacers": 1,
                       "cas_systems": 0, "cas_genes": 0},
        }
        connection.execute("INSERT INTO export_metadata VALUES ('summary',?)", (json.dumps(summary),))
        connection.execute("INSERT INTO assemblies VALUES (1,'GCA_000000001.1','complete','valid',1,1,0,0)")
        connection.execute("INSERT INTO sequences(sequence_pk,assembly_pk,sequence_id) VALUES (1,1,'contig-1')")
        connection.execute(
            "INSERT INTO arrays(array_pk,sequence_pk,tool,call_id,repeat_consensus) "
            "VALUES (1,1,'CRISPRidentify','array-1',?)", (repeat,),
        )
        connection.executemany(
            "INSERT INTO repeats(repeat_pk,array_pk,ordinal,start,end,sequence) VALUES (?,1,?,100,129,?)",
            [(index, index, near) for index in range(1, 27)],
        )
        connection.execute(
            "INSERT INTO spacers(spacer_pk,array_pk,ordinal,kind,sequence) VALUES (1,1,1,'observed',?)", (reverse,),
        )
        add_indexes(connection)
    app = FastAPI()
    app.include_router(build_catalog_router(SimpleNamespace(catalog_database_path=path, catalog_search_index_path=None)))
    job = {
        "job_id": JOB_ID, "status": "completed", "mode": "detection", "progress": 1,
        "summary": {"detection": {"arrays": [array], "category_counts": {"Bona-fide": 1}}},
    }
    return app, job, near, reverse


def main(origin, output):
    output.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="crispr-sequence-browser-") as directory:
        app, job, near, reverse = fixture(Path(directory))
        with TestClient(app) as api, sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=True)
            errors, searches = [], []
            page = browser.new_page(viewport={"width": 1440, "height": 1000})
            page.emulate_media(reduced_motion="reduce")
            page.on("pageerror", lambda error: errors.append(str(error)))
            unavailable = False

            def route_api(route):
                nonlocal unavailable
                parsed = urlsplit(route.request.url)
                if parsed.path.startswith("/api/v1/catalog"):
                    if "match=similar" in parsed.query:
                        searches.append(parsed.query)
                        assert "authorization" not in route.request.headers
                        if unavailable:
                            unavailable = False
                            route.fulfill(status=503, json={"detail": "Reference database unavailable"})
                            return
                    response = api.get(parsed.path + ("?" + parsed.query if parsed.query else ""))
                    route.fulfill(status=response.status_code, json=response.json())
                elif parsed.path.startswith("/api/v1/jobs/"):
                    route.fulfill(json=job)
                else:
                    route.fulfill(json={"status": "ready", "api_version": "v1"})

            page.route("**/api/v1/**", route_api)
            page.goto(origin)
            page.wait_for_load_state("networkidle")
            entry = page.get_by_role("link", name="Search the database", exact=True)
            expect(entry).to_be_visible()
            page.screenshot(path=str(output / "home-desktop.png"), full_page=True)
            page.set_viewport_size({"width": 390, "height": 844})
            expect(entry).to_be_visible()
            assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
            page.screenshot(path=str(output / "home-mobile.png"), full_page=True)
            entry.focus()
            page.keyboard.press("Enter")
            expect(page.get_by_role("heading", name="CRISPR–Cas database", exact=True)).to_be_visible()
            expect(page.get_by_role("searchbox")).to_be_visible()
            page.get_by_role("searchbox").fill("GCA_000000001.1")
            page.get_by_role("button", name="Search", exact=True).click()
            expect(page.get_by_role("button", name="Browse genome GCA_000000001.1")).to_be_visible()
            page.go_back()
            expect(entry).to_be_visible()
            page.go_forward()
            expect(page.get_by_role("searchbox")).to_have_value("GCA_000000001.1")
            page.set_viewport_size({"width": 1440, "height": 1000})
            recovery = origin + "/#job=" + JOB_ID + "." + "a" * 43
            page.goto(recovery)
            page.wait_for_load_state("networkidle")
            expect(page.get_by_role("heading", name="Analysis results", exact=True)).to_be_visible()
            page.get_by_role("tab", name="Arrays", exact=True).click()
            page.locator("summary").filter(has_text="Repeats & spacers").click()
            repeat_buttons = page.get_by_role("button", name=re.compile("Search in DB: repeat"))
            spacer_buttons = page.get_by_role("button", name=re.compile("Search in DB: spacer"))
            assert repeat_buttons.count() == len(job["summary"]["detection"]["arrays"][0]["repeats"])
            assert spacer_buttons.count() == len(job["summary"]["detection"]["arrays"][0]["spacers"])
            assert not searches
            repeat_buttons.first.focus()
            page.keyboard.press("Enter")
            table = page.get_by_role("region", name="Similar sequence candidates")
            expect(table.get_by_text(near, exact=True).first).to_be_visible()
            expect(table.get_by_text("96.7% identity").first).to_be_visible()
            expect(table.get_by_text(re.compile("1 substitution" )).first).to_be_visible()
            page.get_by_role("button", name="Next candidates", exact=True).click()
            expect(page.get_by_text("Page 2 · 1 candidate")).to_be_visible()
            page.get_by_role("button", name="Previous candidates", exact=True).click()
            expect(page.get_by_text("Page 1 · 25 candidates")).to_be_visible()
            page.get_by_role("heading", name=re.compile("Similar candidates:")).evaluate(
                "element => element.scrollIntoView({block: 'start', behavior: 'instant'})"
            )
            page.screenshot(path=str(output / "repeat-candidates-desktop.png"))
            page.get_by_role("button", name="Close search", exact=True).click()
            expect(repeat_buttons.first).to_be_focused()
            spacer_buttons.first.click()
            expect(table.get_by_text(reverse, exact=True)).to_be_visible()
            expect(table.get_by_text(re.compile("Reverse complement"))).to_be_visible()
            assert page.url == recovery
            for width in (390, 1440):
                page.set_viewport_size({"width": width, "height": 900})
                page.get_by_role("heading", name=re.compile("Similar candidates:")).evaluate(
                    "element => element.scrollIntoView({block: 'start', behavior: 'instant'})"
                )
                assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
                page.screenshot(path=str(output / f"spacer-candidates-{width}.png"))
            page.get_by_role("button", name="Close search", exact=True).click()
            # The second spacer has no nearby reference sequence.
            spacer_buttons.nth(1).click()
            expect(page.get_by_text(re.compile("No candidates on this page within one substitution"))).to_be_visible()
            page.get_by_role("button", name="Close search", exact=True).click()
            unavailable = True
            repeat_buttons.first.click()
            expect(page.get_by_role("alert")).to_contain_text("Reference database unavailable")
            page.get_by_role("button", name="Retry search", exact=True).click()
            expect(table).to_be_visible()
            # Database array details provide the same action.
            page.get_by_role("link", name="Database", exact=True).click()
            page.get_by_role("searchbox").fill("GCA_000000001.1")
            page.get_by_role("button", name="Search", exact=True).click()
            page.get_by_role("tab", name="Arrays", exact=True).click()
            page.get_by_role("button", name="Open array array-1", exact=True).click()
            page.get_by_role("button", name="Search in DB: repeat 1", exact=True).click()
            expect(table).to_be_visible()
            assert not errors, errors
            browser.close()
    (output / "acceptance.json").write_text(json.dumps({
        "status": "passed", "search_requests": len(searches),
        "checks": ["main-page database link", "mobile home navigation", "database browser history",
                   "native detector summary", "every repeat and spacer has a lookup", "real catalogue API",
                   "one substitution", "reverse complement", "pagination", "empty results", "503 retry",
                   "keyboard focus", "recovery URL retained", "database array details", "390px layout"],
    }, indent=2) + "\n")
    print("Sequence database search browser acceptance passed:", output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:5179")
    parser.add_argument("--output", type=Path, default=ROOT / "output/playwright/sequence-search")
    args = parser.parse_args()
    main(args.origin, args.output)
