import { afterEach, describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/vue";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import { api } from "../src/api.js";
import { ANALYSIS_MODES, analysisModeAvailable, stagesForMode } from "../src/science.js";
import { normalizePublicResult } from "../src/features/results/index.js";
import { inspectSequenceInput } from "../src/features/submission/sequenceInput.js";
import locusFixture from "./fixtures/locus-compare.json";

const modes = [
  "locus_compare",
  "prophage_context",
  "insertion_context",
  "defense_context",
  "rna_family",
];
const service = { state: "online", modes };
const limits = {
  maxRecords: 100,
  maxBases: 2_000_000,
  maxRecordBases: 2_000_000,
  maxRequestBytes: 2_000_000,
};
afterEach(() => vi.restoreAllMocks());

describe("connected context workflows", () => {
  it.each(modes)("advertises %s only when installed and submits its exact mode", async (mode) => {
    expect(ANALYSIS_MODES.find((item) => item.id === mode).requiresAdvertisement).toBe(true);
    expect(analysisModeAvailable(mode, { state: "online", modes: [] })).toBe(false);
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "c".repeat(32), access_token: "a".repeat(43) });
    render(AnalysisForm, { props: { initialMode: mode, service, limits } });
    expect(
      screen.queryByRole("button", { name: "Load example", exact: true }),
    ).not.toBeInTheDocument();
    expect(screen.queryByText("Analysis options")).not.toBeInTheDocument();
    await fireEvent.update(
      screen.getByRole("textbox"),
      `>one\n${"ACGT".repeat(8)}\n>two\n${"GCGT".repeat(8)}\n`,
    );
    const button = screen.getByRole("button", { name: "Compute", exact: true });
    expect(button).toBeEnabled();
    await fireEvent.click(button);
    expect(submit).toHaveBeenCalledWith(expect.objectContaining({ mode }));
    expect(submit.mock.calls[0][0]).not.toHaveProperty("tool_options");
    expect(stagesForMode(mode).at(-1).id).toBe("package_results");
  });

  it("enforces locus number and per-record size before submission", async () => {
    render(AnalysisForm, { props: { initialMode: "locus_compare", service, limits } });
    const input = screen.getByRole("textbox");
    const button = screen.getByRole("button", { name: "Compute", exact: true });
    await fireEvent.update(input, ">one\nACGT\n");
    expect(button).toBeDisabled();
    await fireEvent.update(
      input,
      Array.from({ length: 9 }, (_, i) => `>record${i}\nACGT\n`).join(""),
    );
    expect(button).toBeDisabled();
    await fireEvent.update(input, `>one\n${"A".repeat(50001)}\n>two\nACGT\n`);
    expect(button).toBeDisabled();
    expect(screen.getByText("Limit: 50,000 bases per record")).toBeInTheDocument();
  });

  it("preserves RNA molecules and blocks duplicates, short and ambiguous family members", async () => {
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "c".repeat(32), access_token: "a".repeat(43) });
    render(AnalysisForm, { props: { initialMode: "rna_family", service, limits } });
    await fireEvent.click(screen.getByRole("radio", { name: /RNA — family members/ }));
    const input = screen.getByRole("textbox");
    const button = screen.getByRole("button", { name: "Compute", exact: true });
    for (const sequence of [
      ">a\nACGU\n>b\nGCGU\n",
      `>a\n${"ACGU".repeat(8)}\n>b\n${"ACGU".repeat(8)}\n`,
      `>a\n${"ACGU".repeat(8)}\n>b\n${"NNNN".repeat(8)}\n`,
    ]) {
      await fireEvent.update(input, sequence);
      expect(button).toBeDisabled();
    }
    await fireEvent.update(input, `>a\n${"ACGU".repeat(8)}\n>b\n${"GCGU".repeat(8)}\n`);
    expect(button).toBeEnabled();
    await fireEvent.click(button);
    expect(submit).toHaveBeenCalledWith(
      expect.objectContaining({ mode: "rna_family", molecule: "RNA" }),
    );
    const packet = {
      schema: { name: "crisprloci.sequence-set", version: "1.0.0" },
      kind: "rna",
      molecule: "RNA",
      complete: true,
      records: [{ id: "one", sequence: "ACGU".repeat(8) }],
    };
    expect(
      inspectSequenceInput(JSON.stringify(packet), { mode: "rna_family", molecule: "RNA" }).valid,
    ).toBe(true);
    expect(
      inspectSequenceInput(JSON.stringify(packet), { mode: "repeats", molecule: "RNA" }).valid,
    ).toBe(false);
  });

  it("renders verified gene maps and alignment evidence without a fabricated array or history tab", async () => {
    render(ResultsView, { props: { job: structuredClone(locusFixture) } });
    expect(screen.queryByRole("tab", { name: "Arrays", exact: true })).not.toBeInTheDocument();
    expect(screen.queryByRole("tab", { name: "History", exact: true })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Gene neighborhoods" }));
    const panel = within(screen.getByRole("tabpanel", { name: "Gene neighborhoods" }));
    expect(panel.getByRole("img", { name: "Gene map for one" })).toBeInTheDocument();
    expect(panel.getAllByText("100.0%")).toHaveLength(2);
    expect(panel.getByText(/not establish orthology/)).toBeInTheDocument();
  });

  it("keeps unknown metrics missing and rejects unbound protein links", () => {
    const summary = structuredClone(locusFixture.summary);
    summary.locus_compare.counts.genes = true;
    summary.locus_compare.links[0].target_gene = "absent";
    const result = normalizePublicResult(summary).locus_compare;
    expect(result.counts.genes).toBeNull();
    expect(result.links).toEqual([]);
    const missing = normalizePublicResult({});
    for (const mode of modes) expect(missing[mode].available).toBe(false);
  });

  it("distinguishes unavailable and truncated prophage evidence from a negative result", async () => {
    const job = {
      job_id: "c".repeat(32),
      status: "completed",
      mode: "prophage_context",
      summary: {},
    };
    const view = render(ResultsView, { props: { job } });
    const panel = within(screen.getByRole("tabpanel", { name: "Overview", exact: true }));
    expect(
      panel.getByText("Validated prophage-context evidence is unavailable."),
    ).toBeInTheDocument();
    await view.rerender({
      job: {
        ...job,
        summary: {
          prophage_context: {
            available: true,
            regions: [],
            spacer_matches: [],
            regions_truncated: true,
            spacer_matches_truncated: true,
            counts: {},
          },
        },
      },
    });
    expect(panel.getByText("No candidate regions in this preview.")).toBeInTheDocument();
    expect(panel.getByText("No spacer matches in this preview.")).toBeInTheDocument();
    expect(panel.queryByText(/This does not establish absence/)).not.toBeInTheDocument();
  });

  it.each([
    ["prophage_context", "No candidate regions in this preview."],
    ["insertion_context", "No insertion sequences in this preview."],
    ["defense_context", "No system calls in this preview."],
  ])("does not treat missing %s counts as zero calls", (mode, message) => {
    render(ResultsView, {
      props: {
        job: {
          job_id: "c".repeat(32),
          status: "completed",
          mode,
          summary: { [mode]: { available: true } },
        },
      },
    });
    const panel = within(screen.getByRole("tabpanel", { name: "Overview", exact: true }));
    expect(panel.getByText(message)).toBeInTheDocument();
  });
});
