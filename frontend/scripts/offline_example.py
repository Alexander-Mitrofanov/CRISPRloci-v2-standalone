"""Browser acceptance for visible examples on services without live mode support."""
from playwright.sync_api import expect


def check_offline_example(page, state, *, link, input_label, sequence, submit="Compute"):
    expect(page.get_by_role("link", name=link, exact=True)).to_be_visible()
    page.get_by_role("link", name=link, exact=True).click()
    page.get_by_label(input_label, exact=True).fill(sequence)
    expect(page.get_by_role("button", name=submit, exact=True)).to_be_disabled()
    assert state["submissions"] == []
    page.get_by_role("button", name="Load example", exact=True).click()
    preview = page.get_by_role("button", name="View precomputed result", exact=True)
    expect(preview).to_be_enabled()
    preview.click()
    expect(page.get_by_role("tab", name="Overview", exact=True)).to_be_visible()
    assert state["submissions"] == []
