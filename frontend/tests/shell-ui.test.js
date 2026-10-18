import { fireEvent, render, screen, waitFor, within } from "@testing-library/vue";
import { afterEach, describe, expect, it, vi } from "vitest";

import App from "../src/App.vue";
import { api } from "../src/api.js";
import { ANALYSIS_MODES } from "../src/science.js";
import { exampleAssetsFor } from "../src/example.js";

const configured = api.configured;
const renderApp = () => {
  api.configured = true;
  vi.spyOn(api, "health").mockResolvedValue({ version: "1.0.0" });
  vi.spyOn(api, "config").mockResolvedValue({
    api_version: "1.0.0",
    modes: ANALYSIS_MODES.map((mode) => mode.id),
  });
  Object.defineProperty(Element.prototype, "scrollIntoView", {
    configurable: true,
    value: vi.fn(),
  });
  return render(App);
};

afterEach(() => {
  api.configured = configured;
  vi.restoreAllMocks();
  window.history.replaceState(null, "", "/");
});

describe("method-first navigation", () => {
  it("offers the database and available methods without an upload form", async () => {
    const { container } = renderApp();
    await screen.findByRole("link", { name: "Spacer searches", exact: true });
    const picker = screen.getByRole("region", { name: "Choose a method" });
    await waitFor(() =>
      expect(within(picker).getAllByRole("link")).toHaveLength(ANALYSIS_MODES.length),
    );
    expect(within(picker).getByRole("link", { name: "Search the database" })).toHaveAttribute(
      "href",
      "?view=database",
    );
    for (const method of ANALYSIS_MODES.filter(
      (item) => !["protospacer", "viral_search"].includes(item.id),
    )) {
      expect(within(picker).getByRole("link", { name: method.title })).toHaveAttribute(
        "href",
        `?method=${method.id}`,
      );
    }
    expect(screen.queryByRole("textbox")).not.toBeInTheDocument();
    expect(container.querySelector('input[type="file"]')).toBeNull();
    expect(screen.queryByRole("img", { name: /Genomic locus/i })).not.toBeInTheDocument();
    expect(container.querySelector("footer")).toBeNull();
  });

  it("opens the database from the main page and restores it through browser history", async () => {
    vi.spyOn(api, "catalogSummary").mockResolvedValue({ available: true });
    const records = vi
      .spyOn(api, "catalogPage")
      .mockResolvedValue({ items: [], next_cursor: null });
    renderApp();
    await fireEvent.click(screen.getByRole("link", { name: "Search the database" }));
    expect(await screen.findByRole("heading", { name: "CRISPR–Cas database" })).toHaveFocus();
    expect(window.location.search).toBe("?view=database");
    expect(await screen.findByRole("combobox", { name: "Search by" })).toBeVisible();
    expect(records).not.toHaveBeenCalled();
    window.history.back();
    await waitFor(() =>
      expect(screen.getByRole("link", { name: "Search the database" })).toBeVisible(),
    );
    window.history.forward();
    await waitFor(() => expect(screen.getByRole("searchbox")).toBeVisible());
  });

  it.each(ANALYSIS_MODES)("opens the corresponding upload page for $id", async (method) => {
    renderApp();
    const spacer = ["protospacer", "viral_search"].includes(method.id);
    await fireEvent.click(
      await screen.findByRole("link", { name: spacer ? "Spacer searches" : method.title }),
    );
    if (method.id === "viral_search")
      await fireEvent.click(
        screen.getByRole("button", { name: "Find viruses for my spacers", exact: true }),
      );
    expect(
      screen.getByRole("heading", { level: 1, name: spacer ? "Spacer searches" : method.title }),
    ).toHaveFocus();
    expect(
      screen.getByRole("button", {
        name: method.id === "array_compare" ? "Upload arrays.json" : "Upload FASTA or JSON",
        exact: true,
      }),
    ).toBeInTheDocument();
    expect(
      screen.getByLabelText(
        method.id === "array_compare" ? "Upload array JSON file" : "Upload FASTA file",
      ),
    ).toBeInTheDocument();
    expect(screen.queryByRole("region", { name: "Choose a method" })).not.toBeInTheDocument();
    expect(window.location.search).toBe(`?method=${method.id}`);
    expect(Boolean(screen.queryByRole("button", { name: "Load example", exact: true }))).toBe(
      Boolean(exampleAssetsFor(method.id)),
    );
  });

  it("restores a direct method link and falls back safely for an unknown method", () => {
    window.history.replaceState(null, "", "/?method=cas");
    const { unmount } = renderApp();
    expect(screen.getByRole("heading", { level: 1, name: "Find Cas systems" })).toBeInTheDocument();
    unmount();
    window.history.replaceState(null, "", "/?method=unknown");
    renderApp();
    expect(screen.getByRole("heading", { name: "Choose a method" })).toBeInTheDocument();
    expect(screen.queryByRole("textbox")).toBeNull();
  });

  it("preserves method-specific drafts and supports browser Back and Forward", async () => {
    renderApp();
    await fireEvent.click(screen.getByRole("link", { name: /^Find Cas systems/ }));
    await fireEvent.update(screen.getByRole("textbox"), ">cas\nACGT\n");
    window.history.back();
    await waitFor(() =>
      expect(screen.getByRole("heading", { name: "Choose a method" })).toBeInTheDocument(),
    );
    window.history.forward();
    await waitFor(() => expect(screen.getByRole("textbox")).toHaveValue(">cas\nACGT\n"));
    await fireEvent.click(screen.getByRole("button", { name: "All methods" }));
    await fireEvent.click(screen.getByRole("link", { name: /^Detect arrays/ }));
    expect(screen.getByRole("textbox")).toHaveValue("");
    await fireEvent.click(screen.getByRole("button", { name: "All methods" }));
    await fireEvent.click(screen.getByRole("link", { name: /^Find Cas systems/ }));
    expect(screen.getByRole("textbox")).toHaveValue(">cas\nACGT\n");
  });
});
