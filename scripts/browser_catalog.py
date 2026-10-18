#!/usr/bin/env python3
"""Check explicit database searches against a running local frontend and catalog."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright

LONG_SPACER = (
    "ACCAAGCTCATAGATTAAAGAATGAATGATTTATAGATGAAAAAGAAAGCTTATAGATTAGAGAATGAATTAATT"
    "ATAGGTAAAAAACAAGCTAATAGATTAAAAAATGAATGATTTATAGATTAAGAAGT"
)


def main(origin: str, output: Path, accession: str, skip_global: bool) -> None:
    output.mkdir(parents=True, exist_ok=True)
    errors: list[str] = []
    record_requests: list[str] = []
    checks: list[str] = []
    fixtures: dict[str, str] = {"accession": accession}
    global_results: dict[str, dict[str, int | bool]] = {}
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1440, "height": 1000}, accept_downloads=True)
        page.emulate_media(reduced_motion="reduce")
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.on("request", lambda request: record_requests.append(request.url)
                if urlsplit(request.url).path.startswith("/api/v1/catalog/") else None)
        page.goto(origin.rstrip("/") + "/")
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
        page.set_viewport_size({"width": 1440, "height": 1000})
        checks.append("main-page database entry, desktop/mobile layout and keyboard navigation")
        page.wait_for_load_state("networkidle")
        selector = page.get_by_role("combobox", name="Search by")
        query = page.get_by_role("searchbox")
        search = page.get_by_role("button", name="Search", exact=True)
        expect(selector).to_be_visible()
        expect(selector.locator("option")).to_have_text(["Accession number", "Repeat", "Spacer", "Cas gene"])
        expect(search).to_be_disabled()
        expect(page.get_by_role("table")).to_have_count(0)
        expect(page.get_by_role("tablist")).to_have_count(0)
        expect(page.get_by_role("complementary", name="Dataset coverage and provenance")).to_have_count(0)
        expect(page.get_by_role("region", name="Search results")).to_have_count(0)
        assert record_requests == [], record_requests
        page.screenshot(path=str(output / "search-desktop.png"), full_page=True)
        page.set_viewport_size({"width": 390, "height": 844})
        page.screenshot(path=str(output / "search-mobile.png"), full_page=True)
        assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
        page.set_viewport_size({"width": 1440, "height": 1000})
        query.fill("GCA_000")
        page.wait_for_load_state("networkidle")
        assert record_requests == []
        search.click()
        expect(page.get_by_role("alert")).to_contain_text("complete genome accession")
        assert record_requests == []
        selector.select_option("repeats")
        query.fill("ACGT")
        page.wait_for_load_state("networkidle")
        expect(page.get_by_role("region", name="Search results")).to_have_count(0)
        assert record_requests == []
        checks.append("no record requests or results before explicit valid search, including typing and changing type")

        def submit(kind: str, value: str) -> None:
            selector.select_option(kind)
            query.fill(value)
            search.click()

        submit("assemblies", accession.lower())
        expect(page.get_by_role("button", name=f"Browse genome {accession}")).to_be_visible()
        expect(page.get_by_role("tablist", name="Genome record type")).to_be_visible()
        expect(page.get_by_role("table").locator("tbody tr")).to_have_count(1)
        assert "assembly_accession=" + accession in record_requests[-1]
        page.screenshot(path=str(output / "genome-desktop.png"), full_page=True)
        count = len(record_requests)
        query.fill(accession + "0")
        page.wait_for_load_state("networkidle")
        expect(page.get_by_role("region", name="Search results")).to_have_count(0)
        assert len(record_requests) == count
        query.fill(accession)
        search.click()
        expect(page.get_by_role("button", name=f"Browse genome {accession}")).to_be_visible()
        checks.append("exact accession and clearing results on input edits")

        page.get_by_role("tab", name="Arrays", exact=True).click()
        first_array = page.get_by_role("button", name=re.compile(r"^Open array ")).first
        expect(first_array).to_be_visible()
        fixtures["array_call"] = first_array.inner_text()
        first_array.click()
        detail = page.get_by_role("region", name="Array details", exact=True)
        expect(detail).to_be_visible()
        expect(detail).to_contain_text("Unknown")
        repeat_table = detail.get_by_role("region", name="Repeat occurrences table")
        spacer_table = detail.get_by_role("region", name="Spacer occurrences table")
        expect(repeat_table).to_be_visible()
        expect(spacer_table).to_be_visible()
        fixtures["repeat"] = repeat_table.locator("code").first.inner_text()
        fixtures["spacer"] = next(value for value in spacer_table.locator("code").all_text_contents()
                                  if re.fullmatch(r"[ACGTRYSWKMBDHVN]+", value))
        page.screenshot(path=str(output / "array-desktop.png"), full_page=True)
        for kind in ("repeat", "spacer"):
            button = detail.get_by_role("button", name=re.compile(f"^Search in DB: {kind} ")).filter(visible=True).first
            expect(button).to_be_enabled()
            button.click()
            candidates = detail.get_by_role("region", name="Similar sequence candidates")
            expect(candidates).to_be_visible(timeout=30000)
            expect(candidates).to_contain_text("identity")
            assert any("match=similar" in request for request in record_requests)
            detail.get_by_role("button", name="Close search", exact=True).click()
            expect(button).to_be_focused()
        checks.append("repeat and spacer Search in DB return similar candidates through the public API")
        detail.get_by_role("button", name="Close array").click()
        with page.expect_download() as csv_download:
            page.get_by_role("button", name="Download this page as CSV", exact=True).click()
        csv_download.value.save_as(output / "arrays-page.csv")
        with page.expect_download() as fasta_download:
            page.get_by_role("button", name="Download consensus FASTA", exact=True).click()
        fasta_download.value.save_as(output / "consensus-page.fasta")
        assert (output / "consensus-page.fasta").read_text().startswith(">arrays|")
        checks.append("scoped arrays, repeat/spacer detail and bounded CSV/FASTA downloads")
        page.get_by_text("About this dataset", exact=True).click()
        coverage = page.get_by_role("complementary", name="Dataset coverage and provenance")
        expect(coverage).to_be_visible()
        coverage_text = coverage.inner_text()
        assert "3,331,847 genome records available here" in coverage_text
        page.get_by_role("tab", name="Cas genes", exact=True).click()
        table = page.get_by_role("region", name="Cas genes records").get_by_role("table")
        expect(table).to_be_visible()
        annotation_column = table.locator("thead th").all_text_contents().index("Annotation")
        annotations = table.locator(f"tbody td:nth-child({annotation_column + 1})").all_text_contents()
        fixtures["cas_annotation"] = next(value for value in annotations if value != "Not reported")
        evidence = table.locator("summary").filter(has_text="View evidence").first
        evidence.click()
        expect(page.get_by_text("Protein length (aa)", exact=True).first).to_be_visible()
        page.screenshot(path=str(output / "cas-genes-desktop.png"), full_page=True)
        checks.append("scoped Cas evidence and provenance available after search")
        count = len(record_requests)
        selector.select_option("spacers")
        page.wait_for_load_state("networkidle")
        expect(page.get_by_role("region", name="Search results")).to_have_count(0)
        assert len(record_requests) == count
        submit("assemblies", accession)
        expect(page.get_by_role("button", name=f"Browse genome {accession}")).to_be_visible()
        count = len(record_requests)
        page.get_by_role("button", name="New search", exact=True).click()
        page.wait_for_load_state("networkidle")
        expect(query).to_have_value("")
        expect(query).to_be_focused()
        expect(page.get_by_role("region", name="Search results")).to_have_count(0)
        assert len(record_requests) == count
        checks.append("changing type and New search clear results without record requests")

        if not skip_global:
            fixtures["long_spacer"] = LONG_SPACER
            fixtures["long_spacer_source"] = "GCA_000003215.1 / array 1000000008 / spacer 1000000019"
            assert len(LONG_SPACER) == 131
            for kind, fixture_key in (("repeats", "repeat"), ("spacers", "spacer"), ("spacers", "long_spacer")):
                submit(kind, fixtures[fixture_key].lower())
                expect(query).to_have_value(fixtures[fixture_key].lower())
                result_table = page.get_by_role("table")
                expect(result_table).to_be_visible(timeout=30000)
                expect(page.get_by_role("tablist")).to_have_count(0)
                sequences = result_table.locator("code").all_text_contents()
                assert sequences and all(value == fixtures[fixture_key] for value in sequences)
                page.screenshot(path=str(output / f"{fixture_key}-global-desktop.png"), full_page=True)
                if fixture_key == "long_spacer":
                    page.set_viewport_size({"width": 390, "height": 844})
                    assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
                    page.screenshot(path=str(output / "long-spacer-global-mobile.png"), full_page=True)
                    page.set_viewport_size({"width": 1440, "height": 1000})
                next_page = page.get_by_role("button", name="Next page", exact=True)
                paged = bool(next_page.count() and next_page.is_enabled())
                global_results[fixture_key] = {
                    "query_length": len(fixtures[fixture_key]), "first_page_rows": len(sequences), "next_page_checked": paged,
                }
                if paged:
                    next_page.click()
                    expect(page.get_by_role("navigation", name="Database pagination")).to_contain_text("Page 2")
                    assert all(value == fixtures[fixture_key] for value in result_table.locator("code").all_text_contents())
            submit("cas_genes", fixtures["cas_annotation"].lower())
            expect(page.get_by_role("table")).to_be_visible(timeout=30000)
            expect(page.get_by_role("tablist")).to_have_count(0)
            expect(page.get_by_role("table")).to_contain_text(fixtures["cas_annotation"])
            submit("cas_genes", "Cas9")
            expect(page.get_by_role("table")).to_be_visible(timeout=30000)
            expect(page.get_by_role("table")).to_contain_text(re.compile(r"cas[\W_]*9", re.IGNORECASE))
            page.screenshot(path=str(output / "cas9-global-desktop.png"), full_page=True)
            page.set_viewport_size({"width": 390, "height": 844})
            assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
            page.screenshot(path=str(output / "cas9-global-mobile.png"), full_page=True)
            page.set_viewport_size({"width": 1440, "height": 1000})
            checks.append("global exact repeat/spacer searches, available pagination, original Cas annotation and Cas9 alias")
            checks.append("real 131-nucleotide spacer submitted without truncation")
            checks.append("390px global long-sequence and Cas9 result layout")

        submit("assemblies", accession)
        expect(page.get_by_role("button", name=f"Browse genome {accession}")).to_be_visible()
        page.get_by_role("link", name="Methods", exact=True).click()
        expect(page.get_by_role("heading", name="Choose a method")).to_be_visible()
        page.go_back()
        expect(page.get_by_role("heading", name="CRISPR–Cas database")).to_be_visible()
        page.set_viewport_size({"width": 390, "height": 844})
        page.get_by_role("tab", name="Arrays", exact=True).click()
        expect(page.get_by_role("button", name=re.compile(r"^Open array ")).first).to_be_visible()
        page.screenshot(path=str(output / "arrays-mobile.png"), full_page=True)
        assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
        expect(page.get_by_role("link", name="Database", exact=True)).to_be_visible()
        page.get_by_role("button", name=re.compile(r"^Open array ")).first.click()
        expect(page.get_by_role("region", name="Array details", exact=True)).to_be_visible()
        assert page.evaluate("document.documentElement.scrollWidth <= window.innerWidth")
        page.screenshot(path=str(output / "array-mobile.png"), full_page=True)
        checks.extend(["Back navigation", "390px search/results layout", "no JavaScript exceptions"])
        assert not errors, errors
        browser.close()
    (output / "search-fixtures.json").write_text(json.dumps(fixtures, indent=2) + "\n")
    (output / "acceptance.json").write_text(json.dumps({
        "status": "passed", "global_search_verified": not skip_global, "coverage": coverage_text,
        "record_request_count": len(record_requests), "checks": checks,
        "global_results": global_results,
    }, indent=2) + "\n")
    print("Catalog search browser acceptance passed:", output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:5173")
    parser.add_argument("--output", type=Path, default=Path("output/playwright/catalog"))
    parser.add_argument("--accession", default="GCA_000008365.1")
    parser.add_argument("--skip-global", action="store_true", help="Verify accession/UI flows before global indexes are available")
    args = parser.parse_args()
    main(args.origin, args.output, args.accession, args.skip_global)
