#!/usr/bin/env python3
"""Browser acceptance for new context modes using intercepted API fixtures only."""

import argparse
import json
from pathlib import Path
from urllib.parse import urlsplit

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[1]
MODES = {
    "locus_compare": ("Compare gene neighborhoods", "Gene neighborhoods"),
    "prophage_context": ("Inspect prophage and spacer context", "Prophage context"),
    "insertion_context": ("Inspect insertion-sequence context", "Insertion sequences"),
    "defense_context": ("Survey the defense repertoire", "Defense repertoire"),
    "rna_family": ("Compare an RNA family", "RNA family"),
}


def fixtures():
    locus = json.loads((ROOT / "tests/fixtures/locus-compare.json").read_text())[
        "summary"
    ]["locus_compare"]
    return {
        "locus_compare": locus,
        "prophage_context": {
            "available": True,
            "counts": {
                "sources": 1,
                "assessed_sources": 1,
                "regions": 1,
                "genes": 2,
                "spacer_matches": 1,
            },
            "sources": [{"id": "one", "length": 20000, "status": "assessed"}],
            "source_inventory_complete": True,
            "source_arrays": [
                {
                    "id": "array1",
                    "native_id": "native_array_1",
                    "source_id": "one",
                    "start": 1,
                    "end": 30,
                }
            ],
            "source_spacers": [
                {
                    "id": "q1",
                    "array_id": "array1",
                    "source_id": "one",
                    "ordinal": 1,
                    "start": 8,
                    "end": 23,
                    "sequence": "ACGTTGCAAGTCGTCA",
                    "status": "eligible",
                    "boundary": None,
                }
            ],
            "regions": [
                {
                    "id": "region1",
                    "source_id": "one",
                    "start": 1,
                    "end": 200,
                    "gene_ids": ["g1", "g2"],
                    "supporting_gene_ids": ["g1"],
                    "peak_score": 0.8,
                    "first_gene_at_record_edge": True,
                    "last_gene_at_record_edge": False,
                }
            ],
            "spacer_matches": [
                {
                    "query_id": "q1",
                    "source_id": "one",
                    "array_id": "array1",
                    "spacer_index": 1,
                    "target_source_id": "one",
                    "start": 40,
                    "end": 70,
                    "strand": "+",
                    "region_ids": ["region1"],
                    "match_class": "sequence_match",
                    "identity": 1,
                }
            ],
        },
        "insertion_context": {
            "available": True,
            "counts": {"elements": 1},
            "elements": [
                {
                    "id": "is1",
                    "source_id": "one",
                    "start": 10,
                    "end": 100,
                    "family": "IS1",
                    "native_type": "p",
                    "contig_edge": False,
                    "transposase_orf": None,
                    "evalue_this_copy": None,
                    "terminal_inverted_repeats": None,
                }
            ],
        },
        "defense_context": {
            "available": True,
            "counts": {
                "sources": 2,
                "genes": 2,
                "systems": 1,
                "components": 1,
                "hits": 1,
            },
            "systems": [
                {
                    "id": "sys1",
                    "source_id": "one",
                    "type": "RM",
                    "subtype": "RM_I",
                    "model": "RM",
                    "gene_ids": ["g1"],
                    "wholeness": 0.5,
                    "score": 2,
                    "start": 20,
                    "end": 90,
                }
            ],
            "genes": [
                {"id": "g1", "source_id": "one", "start": 20, "end": 90, "strand": "+"}
            ],
            "components": [
                {
                    "system_id": "sys1",
                    "gene_id": "g1",
                    "profile": "example",
                    "status": "mandatory",
                    "score": 4,
                }
            ],
            "hits": [{"gene_id": "g1", "profile": "example", "score": 4}],
        },
        "rna_family": {
            "available": True,
            "counts": {
                "sources": 2,
                "unique_sequences": 2,
                "alignment_columns": 32,
                "mfe_pairs": 0,
            },
            "fold": {
                "consensus": "ACGU" * 8,
                "structure": "." * 32,
                "mfe_score_kcal_mol": 0,
                "ensemble_free_energy_kcal_mol": -0.1,
                "temperature_celsius": 37,
            },
            "alignment": [
                {"id": "one", "member_ids": ["one"], "sequence": "ACGU" * 8},
                {"id": "two", "member_ids": ["two"], "sequence": "GCGU" * 8},
            ],
            "pairs": [],
            "warnings": [],
        },
    }


def run(origin, output):
    output.mkdir(parents=True, exist_ok=True)
    evidence = fixtures()
    errors, submissions = [], []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            for mode, (title, tab) in MODES.items():
                context = browser.new_context(viewport={"width": 1440, "height": 1000})
                page = context.new_page()
                page.on("pageerror", lambda error: errors.append(str(error)))
                job = {
                    "job_id": "b" * 32,
                    "status": "completed",
                    "mode": mode,
                    "summary": {mode: evidence[mode]},
                    "artifacts": [],
                }

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
                        data = {"version": "2.12.0"}
                    elif path.endswith("/config"):
                        data = {
                            "modes": list(MODES),
                            "max_records": 100,
                            "max_record_bases": 2000000,
                            "max_total_bases": 2000000,
                            "max_request_bytes": 3000000,
                        }
                    elif path.endswith("/jobs") and request.method == "POST":
                        submissions.append(request.post_data_json)
                        assert request.post_data_json["mode"] == mode
                        data = {"job": job, "access_token": "a" * 43}
                    elif "/jobs/" in path:
                        assert (
                            request.headers.get("authorization") == "Bearer " + "a" * 43
                        )
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
                # Observe rendered link labels before using the declared workflow selector.
                names = page.get_by_role("link").all_text_contents()
                assert any(title in name for name in names), (mode, names)
                page.get_by_role("link", name=title, exact=True).click()
                expect(
                    page.get_by_role("button", name="Load example", exact=True)
                ).to_have_count(0)
                textbox = page.get_by_role("textbox")
                compute = page.get_by_role("button", name="Compute", exact=True)
                if mode == "locus_compare":
                    textbox.fill(">one\nACGT\n")
                    expect(compute).to_be_disabled()
                if mode == "rna_family":
                    page.get_by_role(
                        "radio",
                        name="RNA — family members in transcribed 5′→3′ order",
                        exact=True,
                    ).check()
                    textbox.fill(">one\n" + "ACGU" * 8 + "\n>two\n" + "ACGU" * 8 + "\n")
                    expect(compute).to_be_disabled()
                molecule = "U" if mode == "rna_family" else "T"
                textbox.fill(
                    f">one\n{('ACG' + molecule) * 8}\n>two\n{('GCG' + molecule) * 8}\n"
                )
                expect(compute).to_be_enabled()
                page.screenshot(path=output / f"{mode}-input.png", full_page=True)
                compute.click()
                page.get_by_role("tab", name=tab, exact=True).click()
                expect(
                    page.get_by_role("tabpanel", name=tab, exact=True)
                ).to_be_visible()
                if mode == "prophage_context":
                    page.locator("summary").filter(
                        has_text="Trace source spacer sequences"
                    ).click()
                    expect(
                        page.get_by_text(
                            "array1 / native_array_1 · spacer 1", exact=True
                        )
                    ).to_be_visible()
                page.screenshot(path=output / f"{mode}-desktop.png", full_page=True)
                page.set_viewport_size({"width": 390, "height": 844})
                assert page.evaluate(
                    "document.documentElement.scrollWidth <= innerWidth + 1"
                ), mode
                page.screenshot(path=output / f"{mode}-mobile.png", full_page=True)
                context.close()
            assert not errors, errors
            assert [item["mode"] for item in submissions] == list(MODES)
            assert submissions[-1]["molecule"] == "RNA"
            print(
                json.dumps(
                    {
                        "modes_checked": list(MODES),
                        "submissions": len(submissions),
                        "page_errors": errors,
                        "mobile_overflow": False,
                        "api": "intercepted fixtures; no live jobs submitted",
                    }
                )
            )
        finally:
            browser.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--origin", default="http://127.0.0.1:8491")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.origin, args.output)
