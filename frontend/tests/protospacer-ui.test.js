import { fireEvent, render, screen, within } from "@testing-library/vue";
import { afterEach, describe, expect, it, vi } from "vitest";
import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ModeSelector from "../src/components/submission/ModeSelector.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import { normalizePublicResult } from "../src/features/results/index.js";
import { stagesForMode } from "../src/science.js";
import fixture from "./fixtures/protospacer.json";

const service = { state: "online", modes: ["protospacer"] };
const limits = { maxRecords: 10, maxBases: 1000, maxRecordBases: 1000, maxRequestBytes: 10000 };
const fasta = ">target_A\nAAACGTTGCAACGATTCGAGTCTTGACTCGAATCGTTGCAACGTCC";
const job = () => structuredClone(fixture);
afterEach(() => vi.restoreAllMocks());

describe("advertised protospacer analysis", () => {
  it("keeps example access available when the live reference is absent", async () => {
    const view = render(ModeSelector);
    expect(
      screen.queryByRole("link", { name: "Spacer searches", exact: true }),
    ).toBeInTheDocument();
    await view.rerender({ service });
    expect(screen.getByRole("link", { name: "Spacer searches", exact: true })).toHaveAttribute(
      "href",
      "?method=protospacer",
    );
    await view.rerender({ service: { ...service, state: "offline" } });
    expect(
      screen.queryByRole("link", { name: "Spacer searches", exact: true }),
    ).toBeInTheDocument();
  });

  it("blocks direct-link submission without the configured capability", async () => {
    const submit = vi.spyOn(api, "submit");
    render(AnalysisForm, {
      props: { initialMode: "protospacer", service: { state: "online" }, limits },
    });
    await fireEvent.update(screen.getByLabelText("Contigs or small genomes"), fasta);
    expect(screen.getByRole("button", { name: "Search references", exact: true })).toBeDisabled();
    expect(screen.getByText("Selected analysis unavailable on this service")).toBeInTheDocument();
    await fireEvent.submit(document.getElementById("analysis-form"));
    expect(submit).not.toHaveBeenCalled();
  });

  it("submits target DNA without reference paths or irrelevant controls", async () => {
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "a".repeat(32), access_token: "a".repeat(43) });
    render(AnalysisForm, { props: { initialMode: "protospacer", service, limits } });
    expect(
      screen.getByText("Analysis options", { exact: true }).closest("details"),
    ).not.toHaveAttribute("open");
    expect(
      screen.getByText("CRISPRspacer", { selector: ".tool-options-tool strong" }),
    ).toBeInTheDocument();
    expect(document.querySelector(".tool-option-field")).not.toBeInTheDocument();
    expect(screen.queryByRole("radio")).not.toBeInTheDocument();
    await fireEvent.update(screen.getByLabelText("Contigs or small genomes"), fasta);
    await fireEvent.click(screen.getByRole("button", { name: "Search references", exact: true }));
    expect(submit).toHaveBeenCalledOnce();
    expect(submit.mock.calls[0][0]).toEqual({
      sequence: fasta,
      filename: "pasted-input.fasta",
      mode: "protospacer",
      category_policy: "bona_fide_possible",
      spacer_distance: 1,
      bias_corrections: true,
    });
  });

  it("stops submission if the advertised capability is withdrawn", async () => {
    const submit = vi.spyOn(api, "submit");
    const view = render(AnalysisForm, { props: { initialMode: "protospacer", service, limits } });
    await fireEvent.update(screen.getByLabelText("Contigs or small genomes"), fasta);
    expect(screen.getByRole("button", { name: "Search references", exact: true })).toBeEnabled();
    await view.rerender({ service: { state: "online", modes: [] } });
    await fireEvent.submit(document.getElementById("analysis-form"));
    expect(submit).not.toHaveBeenCalled();
  });

  it("uses only the reference-search progress stage", () => {
    expect(stagesForMode("protospacer").map((row) => row.id)).toEqual([
      "queued",
      "validate_input",
      "crisprspacer",
      "package_results",
    ]);
    expect(stagesForMode("orientation").some((row) => row.id === "crisprspacer")).toBe(false);
  });
});

describe("protospacer result evidence", () => {
  it("shows coordinates, relative strands, reference hashes and limitations without array inference", async () => {
    render(ResultsView, { props: { job: job() } });
    expect(screen.queryByRole("tab", { name: "Arrays", exact: true })).not.toBeInTheDocument();
    expect(screen.queryByRole("tab", { name: "History", exact: true })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Protospacer matches", exact: true }));
    const panel = within(
      screen.getByRole("tabpanel", { name: "Protospacer matches", exact: true }),
    );
    expect(panel.getByText("3–22")).toBeInTheDocument();
    expect(panel.getByText("25–44")).toBeInTheDocument();
    expect(panel.getByText("Forward (+)")).toBeInTheDocument();
    expect(panel.getByText("Reverse (−)")).toBeInTheDocument();
    expect(panel.getAllByText("ACGTTGCAACGATTCGAGTC")).toHaveLength(2);
    expect(panel.getByText("d".repeat(64))).toBeInTheDocument();
    expect(panel.getByText("2017-fixture")).toBeInTheDocument();
    expect(panel.getByText(/CRISPR-array overlap is unknown/)).toBeInTheDocument();
    expect(panel.getByText(/1-based and inclusive/)).toBeInTheDocument();
  });

  it.each([
    [undefined, /Protospacer evidence was not reported/],
    [{ available: true, counts: { matches: 0 } }, /Search completion was not confirmed/],
    [
      {
        available: true,
        search_complete: true,
        outcome: "no_assessable_windows",
        counts: { matches: 0 },
      },
      /No target windows could be assessed/,
    ],
    [
      { available: true, search_complete: true, outcome: "no_matches", counts: { matches: 0 } },
      /No matches to eligible reference spacers/,
    ],
  ])(
    "keeps missing, incomplete, unassessable and completed no-hit results distinct",
    (evidence, message) => {
      const value = job();
      value.summary.protospacer = evidence;
      render(ResultsView, { props: { job: value } });
      const panel = within(screen.getByRole("tabpanel", { name: "Overview", exact: true }));
      expect(panel.getByText(message)).toBeInTheDocument();
      if (evidence?.outcome !== "no_matches")
        expect(panel.queryByText(/No matches to eligible/)).not.toBeInTheDocument();
    },
  );

  it("caps the preview independently of complete search counts", async () => {
    const value = job();
    value.summary.protospacer.counts.matches = 120;
    value.summary.protospacer.matches = Array.from({ length: 110 }, () => ({
      ...value.summary.protospacer.matches[0],
    }));
    render(ResultsView, { props: { job: value } });
    await fireEvent.click(screen.getByRole("tab", { name: "Protospacer matches", exact: true }));
    const panel = within(
      screen.getByRole("tabpanel", { name: "Protospacer matches", exact: true }),
    );
    expect(panel.getByText(/preview shows 100 of 120/)).toBeInTheDocument();
    expect(within(panel.getByRole("table")).getAllByRole("row")).toHaveLength(101);
  });

  it("preserves missing counts, strips unsafe source URLs and renders identifiers as text", async () => {
    const value = job();
    value.summary.protospacer.reference.source_url = "javascript:alert(1)";
    value.summary.protospacer.matches[0].spacer_id = "<img src=x onerror=alert(1)>";
    const normalized = normalizePublicResult({
      protospacer: { counts: { matches: null, targets: NaN } },
    }).protospacer;
    expect(normalized.counts.matches).toBeNull();
    expect(normalized.counts.targets).toBeNull();
    expect(normalized.search_complete).toBeNull();
    render(ResultsView, { props: { job: value } });
    await fireEvent.click(screen.getByRole("tab", { name: "Protospacer matches", exact: true }));
    const panel = within(
      screen.getByRole("tabpanel", { name: "Protospacer matches", exact: true }),
    );
    expect(panel.getByText("<img src=x onerror=alert(1)>")).toBeInTheDocument();
    expect(panel.queryByRole("img")).not.toBeInTheDocument();
    expect(
      panel.queryByRole("link", { name: "Research spacer collection" }),
    ).not.toBeInTheDocument();
  });
});
