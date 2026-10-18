import { fireEvent, render, screen, waitFor, within } from "@testing-library/vue";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import CatalogView from "../src/components/catalog/CatalogView.vue";
import App from "../src/App.vue";
import { api } from "../src/api.js";
import {
  catalogArrayFixture,
  catalogPageFixture,
  catalogSummaryFixture,
} from "./support/catalogFixture.js";

beforeEach(() => {
  vi.spyOn(api, "catalogSummary").mockResolvedValue(catalogSummaryFixture);
  vi.spyOn(api, "catalogPage").mockImplementation(async (entity) => catalogPageFixture(entity));
  vi.spyOn(api, "catalogArray").mockResolvedValue(catalogArrayFixture);
  Object.defineProperty(Element.prototype, "scrollIntoView", {
    configurable: true,
    value: vi.fn(),
  });
});
afterEach(() => {
  vi.restoreAllMocks();
  window.history.replaceState(null, "", "/");
});

async function submit(type = "assemblies", value = "GCA_000001.1") {
  await fireEvent.update(await screen.findByRole("combobox", { name: "Search by" }), type);
  await fireEvent.update(screen.getByRole("searchbox"), value);
  await fireEvent.click(screen.getByRole("button", { name: "Search", exact: true }));
}

describe("search-only public database", () => {
  it("shows no records, counts or empty-result state until a valid search is submitted", async () => {
    render(CatalogView);
    await screen.findByRole("combobox", { name: "Search by" });
    expect(api.catalogPage).not.toHaveBeenCalled();
    expect(screen.queryByRole("table")).not.toBeInTheDocument();
    expect(screen.queryByRole("tablist")).not.toBeInTheDocument();
    expect(screen.queryByRole("region", { name: "Search results" })).not.toBeInTheDocument();
    expect(screen.queryByText("No matching records found.")).not.toBeInTheDocument();
    expect(screen.queryByText("About this dataset")).not.toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Search", exact: true })).toBeDisabled();
    await fireEvent.update(screen.getByRole("searchbox"), "GCA_000");
    expect(api.catalogPage).not.toHaveBeenCalled();
    await fireEvent.click(screen.getByRole("button", { name: "Search", exact: true }));
    expect(screen.getByRole("alert")).toHaveTextContent("complete genome accession");
    expect(api.catalogPage).not.toHaveBeenCalled();
    await fireEvent.update(screen.getByRole("combobox"), "repeats");
    expect(screen.queryByRole("alert")).not.toBeInTheDocument();
    expect(api.catalogPage).not.toHaveBeenCalled();
  });

  it("uses an exact accession and permits related genome records after submission", async () => {
    render(CatalogView);
    await submit("assemblies", "gca_000001.1");
    await screen.findByRole("button", { name: "Browse genome GCA_000001.1" });
    expect(api.catalogPage).toHaveBeenCalledWith(
      "assemblies",
      expect.objectContaining({ assembly_accession: "GCA_000001.1", q: "" }),
      expect.anything(),
    );
    await fireEvent.click(screen.getByRole("tab", { name: "Arrays", exact: true }));
    await fireEvent.click(await screen.findByRole("button", { name: "Open array array-1" }));
    const detail = await screen.findByRole("region", { name: "Array details" });
    expect(api.catalogPage).toHaveBeenCalledWith(
      "arrays",
      expect.objectContaining({ assembly_accession: "GCA_000001.1", q: "" }),
      expect.anything(),
    );
    expect(detail).toHaveTextContent("StrandUnknown");
    expect(detail).toHaveTextContent("Caller scoreNot reported");
    expect(detail).toHaveTextContent("Boundary 140 (0-based)");
    expect(detail).toHaveTextContent("No sequence (deletion)");
    await fireEvent.click(within(detail).getByRole("button", { name: "Close array" }));
    expect(screen.queryByRole("region", { name: "Array details" })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByText("About this dataset"));
    const coverage = screen.getByRole("complementary", { name: "Dataset coverage and provenance" });
    expect(coverage).toHaveTextContent("2 genome records available here");
    expect(coverage).toHaveTextContent("campaign catalog of 100");
  });

  it.each([
    ["repeats", "acgt acgt", "ACGTACGT"],
    ["spacers", "gattaca", "GATTACA"],
    ["spacers", "acgt".repeat(35), "ACGT".repeat(35)],
    ["cas_genes", "Cas9", "Cas9"],
  ])("searches %s across the catalogue only after submission", async (type, input, expected) => {
    render(CatalogView);
    await submit(type, input);
    await screen.findByRole("table");
    expect(api.catalogPage).toHaveBeenCalledWith(
      type,
      expect.objectContaining({ q: expected, assembly_accession: "", limit: 25 }),
      expect.anything(),
    );
    expect(screen.queryByRole("tablist")).not.toBeInTheDocument();
    await fireEvent.update(screen.getByRole("searchbox"), "");
    expect(screen.queryByRole("table")).not.toBeInTheDocument();
    expect(screen.queryByRole("region", { name: "Search results" })).not.toBeInTheDocument();
    expect(api.catalogPage).toHaveBeenCalledTimes(1);
  });

  it("preserves the submitted query across bounded pages and retries", async () => {
    api.catalogPage.mockImplementation(async (entity, filters) => ({
      ...catalogPageFixture(entity),
      next_cursor: !filters.cursor ? "next-page" : null,
    }));
    render(CatalogView);
    await submit("repeats", "ACGTACGT");
    await fireEvent.click(await screen.findByRole("button", { name: "Next page" }));
    await waitFor(() =>
      expect(api.catalogPage).toHaveBeenLastCalledWith(
        "repeats",
        expect.objectContaining({ q: "ACGTACGT", cursor: "next-page" }),
        expect.anything(),
      ),
    );
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "Previous page" })).toBeEnabled(),
    );
    api.catalogPage.mockRejectedValueOnce(new Error("Please retry this search."));
    await fireEvent.click(screen.getByRole("button", { name: "Previous page" }));
    expect(await screen.findByRole("alert")).toHaveTextContent("Please retry this search.");
    await fireEvent.click(screen.getByRole("button", { name: "Retry search" }));
    await screen.findByRole("table");
    expect(api.catalogPage).toHaveBeenLastCalledWith(
      "repeats",
      expect.objectContaining({ q: "ACGTACGT", cursor: null }),
      expect.anything(),
    );
  });

  it("cancels stale responses when editing, changing type or starting a new search", async () => {
    let resolveOld;
    api.catalogPage.mockImplementationOnce(
      () =>
        new Promise((resolve) => {
          resolveOld = resolve;
        }),
    );
    render(CatalogView);
    await submit();
    await fireEvent.update(screen.getByRole("searchbox"), "GCA_000002.1");
    resolveOld(catalogPageFixture("assemblies"));
    await waitFor(() =>
      expect(screen.queryByRole("region", { name: "Search results" })).not.toBeInTheDocument(),
    );
    await submit("spacers", "GATTACA");
    await screen.findByRole("table");
    await fireEvent.update(screen.getByRole("combobox"), "repeats");
    expect(screen.queryByRole("table")).not.toBeInTheDocument();
    expect(api.catalogPage).toHaveBeenCalledTimes(2);
    await submit();
    await screen.findByRole("table");
    await fireEvent.click(screen.getByRole("button", { name: "New search" }));
    expect(screen.queryByRole("region", { name: "Search results" })).not.toBeInTheDocument();
    expect(screen.getByRole("searchbox")).toHaveValue("");
    expect(screen.getByRole("searchbox")).toHaveFocus();
    expect(api.catalogPage).toHaveBeenCalledTimes(3);
  });

  it("rejects invalid sequences and distinguishes bounded empty pages from exhausted results", async () => {
    render(CatalogView);
    await submit("repeats", "not a sequence!");
    expect(screen.getByRole("alert")).toHaveTextContent("nucleotide sequence");
    expect(api.catalogPage).not.toHaveBeenCalled();
    await submit("repeats", "A".repeat(1001));
    expect(screen.getByRole("alert")).toHaveTextContent("1,000 characters or fewer");
    expect(screen.getByRole("searchbox")).toHaveValue("A".repeat(1001));
    expect(api.catalogPage).not.toHaveBeenCalled();
    api.catalogPage.mockResolvedValueOnce({ items: [], next_cursor: "later-shards", limit: 25 });
    await submit("repeats", "ACGT");
    await screen.findByText("No exact matches on this page. Continue to the next page.");
    api.catalogPage.mockResolvedValueOnce({ items: [], next_cursor: null, limit: 25 });
    await fireEvent.click(screen.getByRole("button", { name: "Next page" }));
    await screen.findByText("No matching records found.");
    expect(screen.getByRole("button", { name: "Next page" })).toBeDisabled();
  });

  it("handles unavailable data and retries without browsing the catalogue", async () => {
    api.catalogSummary.mockRejectedValueOnce(new Error("The database is offline."));
    const { unmount } = render(CatalogView);
    expect(await screen.findByRole("alert")).toHaveTextContent("The database is offline.");
    await fireEvent.click(screen.getByRole("button", { name: "Retry database" }));
    await screen.findByRole("searchbox");
    expect(api.catalogPage).not.toHaveBeenCalled();
    unmount();
    api.catalogSummary.mockResolvedValue({ available: false });
    render(CatalogView);
    await screen.findByRole("heading", { name: "Database is not available yet" });
    expect(screen.queryByRole("searchbox")).not.toBeInTheDocument();
  });

  it("supports direct database navigation and returning to methods", async () => {
    vi.spyOn(api, "health").mockResolvedValue({ version: "1.0.0" });
    vi.spyOn(api, "config").mockResolvedValue({ api_version: "v1" });
    window.history.replaceState(null, "", "/?view=database#unrelated-fragment");
    render(App);
    await screen.findByRole("searchbox");
    expect(api.catalogPage).not.toHaveBeenCalled();
    expect(window.location.hash).toBe("#unrelated-fragment");
    await fireEvent.click(screen.getByRole("link", { name: "Methods", exact: true }));
    expect(screen.getByRole("heading", { name: "Choose a method" })).toBeInTheDocument();
    await fireEvent.click(screen.getByRole("link", { name: "Database", exact: true }));
    expect(window.location.search).toBe("?view=database");
    expect(screen.getByRole("heading", { name: "CRISPR–Cas database" })).toHaveFocus();
  });
});
