import { afterEach, describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/vue";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import { api } from "../src/api.js";
import { inspectSequenceInput } from "../src/features/submission/sequenceInput.js";
import { normalizePublicResult } from "../src/features/results/index.js";
import packet from "./fixtures/array-compare-input.json";
import fixture from "./fixtures/array-compare.json";
const service = { state: "online", modes: ["array_compare"] };
const limits = {
  maxRecords: 20,
  maxBases: 1000000,
  maxRecordBases: 1000,
  maxRequestBytes: 1000000,
};
afterEach(() => vi.restoreAllMocks());
describe("Observed array comparison", () => {
  it("submits actual array JSON with explicit equality and respects occurrence capacity", async () => {
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "c".repeat(32), access_token: "a".repeat(43) });
    const view = render(AnalysisForm, { props: { initialMode: "array_compare", service, limits } });
    const input = screen.getByLabelText("Observed array JSON");
    const button = screen.getByRole("button", { name: "Compute", exact: true });
    await fireEvent.update(input, ">genome\nACGT");
    expect(button).toBeDisabled();
    await fireEvent.update(input, JSON.stringify(packet));
    expect(button).toBeEnabled();
    await fireEvent.click(screen.getByText("Analysis options").closest("summary"));
    await fireEvent.click(document.querySelector(".tool-options-tool > summary"));
    await fireEvent.update(screen.getByLabelText(/DNA equality/), "reverse-complement");
    await fireEvent.click(button);
    expect(submit).toHaveBeenCalledWith(
      expect.objectContaining({
        mode: "array_compare",
        sequence: JSON.stringify(packet),
        tool_options: { cctk: { "--identity": "reverse-complement" } },
      }),
    );
    await view.rerender({ limits: { ...limits, maxRecords: 5 } });
    expect(button).toBeDisabled();
  });
  it("rejects missing/ambiguous spacer evidence without imposing boundary or orientation certainty", () => {
    const inspect = (doc) => inspectSequenceInput(JSON.stringify(doc), { mode: "array_compare" });
    expect(inspect(packet).valid).toBe(true);
    const doc = structuredClone(packet);
    doc.arrays[0].spacers[0].sequence = "NNNN";
    expect(inspect(doc).valid).toBe(false);
    doc.arrays[0].spacers[0].sequence = null;
    expect(inspect(doc).valid).toBe(false);
  });
  it.each(["array", "spacer"])("rejects a null %s and keeps the form usable", async (kind) => {
    const doc = structuredClone(packet);
    if (kind === "array") doc.arrays[0] = null;
    else doc.arrays[0].spacers[0] = null;
    render(AnalysisForm, { props: { initialMode: "array_compare", service, limits } });
    const input = screen.getByLabelText("Observed array JSON");
    const button = screen.getByRole("button", { name: "Compute", exact: true });
    await fireEvent.update(input, JSON.stringify(doc));
    expect(button).toBeDisabled();
    await fireEvent.update(input, JSON.stringify(packet));
    expect(button).toBeEnabled();
  });
  it("shows distinct gaps, duplicates, unknown strand, original IDs and a sharing graph", async () => {
    render(ResultsView, { props: { job: structuredClone(fixture) } });
    expect(screen.queryByRole("tab", { name: "Arrays", exact: true })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Array comparison", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Array comparison", exact: true }));
    expect(panel.getAllByText("Unknown", { exact: true }).length).toBeGreaterThan(0);
    expect(panel.getByRole("group", { name: /sharing graph/ })).toBeInTheDocument();
    expect(panel.getByText(/no observed DNA evidence/)).toBeInTheDocument();
    expect(panel.getAllByText("×").length).toBe(2);
    expect(panel.getByText("0.333")).toBeInTheDocument();
    const node = panel.getAllByRole("button", { name: /Highlight/ })[0];
    await fireEvent.click(node);
    expect(document.querySelector(".observed-array.selected")).not.toBeNull();
  });
  it("distinguishes an empty preview from globally absent sharing", async () => {
    const doc = structuredClone(fixture);
    doc.summary.array_compare.edges = [];
    doc.summary.array_compare.edges_truncated = true;
    render(ResultsView, { props: { job: doc } });
    await fireEvent.click(screen.getByRole("tab", { name: "Array comparison", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Array comparison", exact: true }));
    expect(panel.getByText(/No sharing links in this preview/)).toBeInTheDocument();
    expect(panel.queryByText(/No links meet the selected threshold/)).not.toBeInTheDocument();
  });
  it("never converts missing orientation or malformed metrics into claims", () => {
    const doc = structuredClone(fixture.summary);
    doc.array_compare.edges[0].jaccard = true;
    const result = normalizePublicResult(doc).array_compare;
    expect(result.edges).toEqual([]);
    expect(result.arrays[0].strand).toBeNull();
  });
});
