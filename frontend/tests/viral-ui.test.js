import { afterEach, describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/vue";
import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import fixture from "./fixtures/viral.json";
const service = { state: "online", modes: ["protospacer", "viral_search"] };
const limits = { maxRecords: 20, maxBases: 10000, maxRecordBases: 10000, maxRequestBytes: 20000 };
const fasta = ">spacer\nACGTTGCAACGATTCGAGTC\n";
afterEach(() => vi.restoreAllMocks());
describe("bidirectional spacer search", () => {
  it("exposes two clearly labelled directions and emits navigation", async () => {
    const view = render(AnalysisForm, { props: { initialMode: "protospacer", service, limits } });
    expect(screen.getByRole("button", { name: "Find spacers in my sequence" })).toHaveAttribute(
      "aria-pressed",
      "true",
    );
    await fireEvent.click(screen.getByRole("button", { name: "Find viruses for my spacers" }));
    expect(view.emitted()["change-mode"]).toEqual([["viral_search"]]);
  });
  it.each([0, 1, 2])("submits the chosen substitution policy %i", async (mismatches) => {
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "c".repeat(32), access_token: "a".repeat(43) });
    render(AnalysisForm, { props: { initialMode: "viral_search", service, limits } });
    await fireEvent.update(screen.getByLabelText("Spacer sequences"), fasta);
    await fireEvent.update(screen.getByLabelText(/Maximum substitutions/), mismatches);
    await fireEvent.click(screen.getByRole("button", { name: "Search references", exact: true }));
    expect(submit).toHaveBeenCalledOnce();
    expect(submit.mock.calls[0][0]).toMatchObject({
      mode: "viral_search",
      sequence: fasta,
      viral_max_mismatches: mismatches,
    });
    expect(submit.mock.calls[0][0]).not.toHaveProperty("reference");
  });
  it("disables a withdrawn reference and rejects oversized spacer records", async () => {
    const submit = vi.spyOn(api, "submit");
    const view = render(AnalysisForm, { props: { initialMode: "viral_search", service, limits } });
    await fireEvent.update(screen.getByLabelText("Spacer sequences"), ">q\n" + "A".repeat(81));
    expect(screen.getByRole("button", { name: "Search references", exact: true })).toBeDisabled();
    await fireEvent.update(screen.getByLabelText("Spacer sequences"), fasta);
    expect(screen.getByRole("button", { name: "Search references", exact: true })).toBeEnabled();
    await view.rerender({ service: { state: "online", modes: ["protospacer"] } });
    expect(screen.getByRole("button", { name: "Find viruses for my spacers" })).toBeEnabled();
    expect(screen.getByRole("button", { name: "Search references", exact: true })).toBeDisabled();
    expect(screen.getByText(/precomputed example available/)).toBeInTheDocument();
    await fireEvent.submit(document.getElementById("analysis-form"));
    expect(submit).not.toHaveBeenCalled();
  });
});
describe("viral sequence evidence", () => {
  it("shows distinct support, matches, filtering and reference identity", async () => {
    render(ResultsView, { props: { job: structuredClone(fixture) } });
    expect(screen.queryByRole("tab", { name: "Arrays", exact: true })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Viral matches", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Viral matches", exact: true }));
    expect(panel.getByText(/do not add independent sequence support/)).toBeInTheDocument();
    expect(
      panel.getByText(fixture.summary.viral_search.reference.reference_sha256),
    ).toBeInTheDocument();
    expect(within(panel.getByRole("table")).getAllByRole("row")).toHaveLength(4);
    await fireEvent.update(panel.getByLabelText(/Filter displayed matches/), "no-such-accession");
    expect(panel.getByText("No displayed matches contain that text.")).toBeInTheDocument();
  });
  it.each([
    [null, /A completed viral search was not reported/],
    [{ available: true, execution_completed: false }, /A completed viral search was not reported/],
    [
      { available: true, execution_completed: true, outcome: "no_eligible_queries" },
      /No eligible spacers to search/,
    ],
    [
      { available: true, execution_completed: true, outcome: "no_raw_hsps" },
      /No full-length matches passed/,
    ],
  ])("keeps missing, incomplete, skipped and negative evidence distinct", (evidence, message) => {
    const job = structuredClone(fixture);
    job.summary.viral_search = evidence;
    render(ResultsView, { props: { job } });
    expect(
      within(screen.getByRole("tabpanel", { name: "Overview", exact: true })).getByText(message),
    ).toBeInTheDocument();
  });
});
