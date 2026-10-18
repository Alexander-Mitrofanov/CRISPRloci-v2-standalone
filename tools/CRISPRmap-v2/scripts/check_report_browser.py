#!/usr/bin/env python3
"""Smoke-test a generated local HTML report in Chromium at desktop/mobile sizes."""
from pathlib import Path
import argparse


def main():
    from playwright.sync_api import sync_playwright

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--screenshot-dir", type=Path)
    parser.add_argument("--chromium", type=Path, help="optional already-installed Chromium executable")
    args = parser.parse_args()
    with sync_playwright() as playwright:
        kwargs = {"executable_path": str(args.chromium)} if args.chromium else {}
        browser = playwright.chromium.launch(headless=True, **kwargs)
        errors, remote_requests = [], []
        for label, width, height in (("desktop", 1440, 960), ("mobile", 390, 844)):
            page = browser.new_page(viewport={"width": width, "height": height})
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.on("request", lambda request: remote_requests.append(request.url)
                    if request.url.startswith(("http:", "https:")) else None)
            page.goto(args.report.resolve().as_uri(), wait_until="networkidle")
            assert page.get_by_role("heading", name="CRISPRmap-v2", exact=True).count() == 1
            assert page.locator("tbody tr").count() > 0
            detail = page.locator("details").first
            detail.locator("summary").click()
            assert detail.evaluate("element => element.open")
            assert page.evaluate("document.documentElement.scrollWidth <= innerWidth"), "Page overflows viewport"
            if args.screenshot_dir:
                args.screenshot_dir.mkdir(parents=True, exist_ok=True)
                page.screenshot(path=str(args.screenshot_dir / f"{label}.png"), full_page=True)
            page.close()
        browser.close()
    assert not errors, errors
    assert not remote_requests, remote_requests
    print("Report browser check: desktop/mobile rendering, alignment details, no remote requests, no page errors")


if __name__ == "__main__":
    main()
