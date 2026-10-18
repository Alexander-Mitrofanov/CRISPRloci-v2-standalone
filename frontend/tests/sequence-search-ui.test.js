import { fireEvent, render, screen, waitFor, within } from "@testing-library/vue";
import { afterEach, describe, expect, it, vi } from "vitest";
import { api } from "../src/api.js";
import ResultsView from "../src/components/results/ResultsView.vue";
import SequenceSearchPanel from "../src/components/catalog/SequenceSearchPanel.vue";
import { sequenceComparison } from "../src/features/catalog/index.js";

const candidate = (sequence = "ACGA", id = 1) => ({
  id,
  array_id: 42,
  accession: "GCA_000000001.1",
  sequence_id: "contig-1",
  ordinal: id,
  start: 100,
  end: 103,
  sequence,
});
const page = (items = [candidate()], next_cursor = null) => ({ items, next_cursor, limit: 25 });
afterEach(() => vi.restoreAllMocks());

describe("array sequence database search", () => {
  it("offers a lookup for every observed repeat/spacer and reports real comparison values", async () => {
    const lookup = vi.spyOn(api, "catalogPage").mockResolvedValue(page());
    render(ResultsView, {
      props: {
        job: {
          job_id: "test",
          mode: "detection",
          status: "completed",
          summary: {
            detection: {
              arrays: [
                {
                  source_id: "record-1",
                  array_id: "array-1",
                  repeats: [{ sequence: "ACGT" }, { sequence: "ACGA" }],
                  spacers: [{ sequence: "TTT" }, { sequence: null }],
                },
              ],
            },
          },
        },
      },
    });
    await fireEvent.click(screen.getByRole("tab", { name: "Arrays", exact: true }));
    const disclosure = screen.getByText("Repeats & spacers — record-1 / array-1");
    await fireEvent.click(disclosure);
    // jsdom does not implement the native details click toggle.
    disclosure.closest("details").open = true;
    const buttons = screen.getAllByRole("button", { name: /^Search in DB:/u });
    expect(buttons).toHaveLength(4);
    expect(buttons[3]).toBeDisabled();
    expect(lookup).not.toHaveBeenCalled();
    const originalUrl = window.location.href;
    await fireEvent.click(buttons[0]);
    expect(await screen.findByText("75.0% identity")).toBeVisible();
    expect(screen.getByText("1 substitution · Both")).toBeVisible();
    expect(screen.getByText("GCA_000000001.1")).toBeVisible();
    expect(lookup).toHaveBeenCalledWith(
      "repeats",
      expect.objectContaining({ q: "ACGT", match: "similar" }),
      expect.anything(),
    );
    await fireEvent.click(screen.getByRole("button", { name: "Close search" }));
    expect(buttons[0]).toHaveFocus();
    lookup.mockResolvedValue(page([candidate("AAA")]));
    await fireEvent.click(buttons[2]);
    expect(await screen.findByText("100.0% identity")).toBeVisible();
    expect(screen.getByText("0 substitutions · Reverse complement")).toBeVisible();
    expect(lookup).toHaveBeenLastCalledWith(
      "spacers",
      expect.objectContaining({ q: "TTT" }),
      expect.anything(),
    );
    expect(window.location.href).toBe(originalUrl);
  });

  it("keeps empty continuation pages distinct from completion and retries failed pages", async () => {
    const lookup = vi
      .spyOn(api, "catalogPage")
      .mockResolvedValueOnce(page([], "next"))
      .mockRejectedValueOnce(new Error("Reference database unavailable"))
      .mockResolvedValueOnce(page([candidate()]));
    render(SequenceSearchPanel, {
      props: { kind: "repeats", sequence: "ACGT", label: "Repeat 1" },
    });
    expect(await screen.findByText(/Continue to check the next database records/u)).toBeVisible();
    await fireEvent.click(screen.getByRole("button", { name: "Next candidates" }));
    expect(await screen.findByRole("alert")).toHaveTextContent("Reference database unavailable");
    await fireEvent.click(screen.getByRole("button", { name: "Retry search" }));
    expect(await screen.findByText("75.0% identity")).toBeVisible();
    expect(screen.getByText("Page 2 · 1 candidate")).toBeVisible();
    expect(lookup).toHaveBeenLastCalledWith(
      "repeats",
      expect.objectContaining({ cursor: "next" }),
      expect.anything(),
    );
    expect(screen.getByRole("button", { name: "Next candidates" })).toBeDisabled();
  });

  it("ignores stale replies and aborts when the chosen sequence changes or closes", async () => {
    let resolveOld;
    const lookup = vi
      .spyOn(api, "catalogPage")
      .mockImplementationOnce(
        () =>
          new Promise((resolve) => {
            resolveOld = resolve;
          }),
      )
      .mockResolvedValue(page([candidate("AAA")]));
    const { rerender, unmount } = render(SequenceSearchPanel, {
      props: { kind: "repeats", sequence: "ACGT", label: "Repeat 1" },
    });
    const firstSignal = lookup.mock.calls[0][2].signal;
    await rerender({ kind: "spacers", sequence: "TTT", label: "Spacer 1" });
    expect(await screen.findByText("100.0% identity")).toBeVisible();
    expect(firstSignal.aborted).toBe(true);
    resolveOld(page());
    await waitFor(() => expect(screen.queryByText("75.0% identity")).toBeNull());
    const candidates = screen.getByRole("region", { name: "Similar sequence candidates" });
    expect(within(candidates).getByText("AAA")).toBeVisible();
    const currentSignal = lookup.mock.calls[1][2].signal;
    unmount();
    expect(currentSignal.aborted).toBe(true);
  });

  it("compares ambiguity symbols literally, retains missingness and rejects gaps", () => {
    expect(sequenceComparison("ACGN", "ACGT")?.substitutions).toBe(1);
    expect(sequenceComparison("ACNN", "ACGT")).toBeNull();
    expect(sequenceComparison("ACGT", "ACG-T")).toBeNull();
    expect(sequenceComparison("ACGT", null)).toBeNull();
    expect(sequenceComparison("ACGT", "ACGT")?.strand).toBe("Both");
  });
});
