import { fireEvent, render, screen, within } from "@testing-library/vue";
import { afterEach, describe, expect, it, vi } from "vitest";
import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import { inspectFasta } from "../src/fasta.js";
import { normalizePublicResult } from "../src/features/results/index.js";
import { stagesForMode } from "../src/science.js";
import { buildSubmission } from "../src/submission.js";

const limits = {
  maxRecords: 2000,
  maxBases: 2_000_000,
  maxRecordBases: 2_000_000,
  maxRequestBytes: 2_000_000,
};
const form = (initialMode = "repeats") =>
  render(AnalysisForm, {
    props: { initialMode, service: { state: "online" }, limits },
  });
const instance = (orientation = "forward") => ({
  id: orientation,
  source_id: "contig_A",
  array_id: "array_1",
  repeat_index: 1,
  start: 100,
  end: 109,
  orientation,
  source_sequence: "GGGTTTCCC",
  rna_sequence: "GGGUUUCCC",
  terminal: true,
  isolated: {
    status: "completed",
    structure: "(((...)))",
    mfe_kcal_mol: -1.2,
    ensemble_kcal_mol: -1.39,
    motif_pairs: [
      [0, 8],
      [1, 7],
      [2, 6],
    ],
    pairs: [
      [0, 8, 0.75],
      [1, 7, 0.76],
      [2, 6, 0.74],
    ],
    unpaired: Array(9).fill(0.1),
  },
  comparison: {
    status: "completed",
    motif_origin: "isolated_mfe_hypothesis",
    motif_pairs: [
      [0, 8],
      [1, 7],
      [2, 6],
    ],
    isolated_pair_support: 0.75,
    context_pair_support: 0.6,
    support_change: -0.15,
  },
  context: {
    status: "completed",
    effective_window: 150,
    effective_span: 80,
    // Local pair averages may exceed one when summed for a base. Preserve native accessibility.
    intrarepeat_pairs: [
      [0, 8, 0.7],
      [0, 7, 0.8],
    ],
    unpaired_probabilities: Array(9).fill(0.25),
    competitors: [{ i: 4, j: 14, probability: 0.3 }],
    competitor_count: 25,
    competitors_truncated: true,
  },
  references: { status: "unsupported", reason: "reference_not_configured", hits: [] },
});
const job = (mode = "repeat_context") => ({
  job_id: "repeat-test",
  mode,
  status: "completed",
  options: { molecule: "DNA" },
  summary: {
    repeats: {
      status: "completed",
      molecule: "DNA",
      counts: { sources: 1, arrays: 1, instances: 2, contexts: 2 },
      model: {
        engine: "ViennaRNA",
        version: "2.7.0",
        parameters: "RNA_Turner2004",
        temperature: 37,
      },
      classification: {
        strand: { status: "unsupported", reason: "no_calibrated_strand_model" },
        subtype: { status: "unsupported", reason: "no_validated_subtype_model" },
        family: { status: "unsupported", reason: "no_validated_family_model" },
      },
      instances: [instance(), instance("reverse")],
      truncated: true,
    },
  },
});

afterEach(() => vi.restoreAllMocks());

describe("repeat input contract", () => {
  it("separates DNA/RNA alphabets and rejects mixed T/U before normalization", () => {
    expect(inspectFasta(">rna\nacgun", { molecule: "RNA" }).valid).toBe(true);
    expect(inspectFasta(">rna\nacgun").valid).toBe(false);
    expect(inspectFasta(">mixed\nACGTU", { molecule: "RNA" }).valid).toBe(false);
    expect(inspectFasta(">mixed\nACGTU").valid).toBe(false);
    expect(inspectFasta(">rna\nAſGU", { molecule: "RNA" }).valid).toBe(false);
  });

  it("submits declared RNA for repeats and preserves DNA-only genomic modes", async () => {
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "a".repeat(32), access_token: "a".repeat(43) });
    form();
    await fireEvent.click(screen.getByRole("radio", { name: /RNA repeats/ }));
    await fireEvent.update(screen.getByLabelText("Repeat sequences"), ">repeat\nGGGUUUCCC");
    await fireEvent.click(screen.getByRole("button", { name: "Compute", exact: true }));
    expect(submit).toHaveBeenCalledWith(
      expect.objectContaining({ mode: "repeats", molecule: "RNA" }),
    );
    expect(
      buildSubmission({ mode: "repeat_context", options: { molecule: "RNA" } }),
    ).not.toHaveProperty("molecule");
  });

  it("blocks supplied repeats longer than 200 nt", async () => {
    form();
    await fireEvent.update(
      screen.getByLabelText("Repeat sequences"),
      `>repeat\n${"A".repeat(201)}`,
    );
    expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeDisabled();
    expect(screen.getByText("Limit: 200 bases per record")).toBeInTheDocument();
  });

  it("caps supplied repeat records at 1000 while retaining stricter service limits", async () => {
    form();
    await fireEvent.update(
      screen.getByLabelText("Repeat sequences"),
      Array.from({ length: 1001 }, (_, index) => `>r${index}\nA`).join("\n"),
    );
    expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeDisabled();
    expect(screen.getByText("Limit: 1000 records")).toBeInTheDocument();
  });

  it("runs the repeat stage only for repeat workflows", () => {
    expect(stagesForMode("repeats").map((stage) => stage.id)).toEqual([
      "queued",
      "validate_input",
      "crisprrepeat",
      "package_results",
    ]);
    expect(stagesForMode("repeat_context").some((stage) => stage.id === "detect_arrays")).toBe(
      true,
    );
    expect(stagesForMode("orientation").some((stage) => stage.id === "crisprrepeat")).toBe(false);
  });
});

describe("repeat evidence results", () => {
  it("shows individual orientation alternatives, exact matrix values and native accessibility", async () => {
    render(ResultsView, { props: { job: job() } });
    await fireEvent.click(screen.getByRole("tab", { name: "Repeat evidence", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Repeat evidence", exact: true }));
    expect(panel.getByText("Forward hypothesis")).toBeInTheDocument();
    expect(panel.getByText("Reverse-complement hypothesis")).toBeInTheDocument();
    expect(panel.getByText(/preview is limited to 2/)).toBeInTheDocument();
    expect(panel.getAllByText("[100, 109)")).toHaveLength(2);
    expect(panel.getByText(/No calibrated strand model/)).toBeInTheDocument();
    expect(
      panel.getByRole("img", { name: /Isolated and array-context pair-support matrix/ }),
    ).toBeInTheDocument();
    const selected = within(panel.getByRole("region", { name: "Selected repeat detail" }));
    expect(selected.getByText("0.75")).toBeInTheDocument();
    expect(selected.getByText("0.7")).toBeInTheDocument();
    await fireEvent.click(selected.getByText("Single-base accessibility and sequence positions"));
    const unpaired = within(selected.getByRole("region", { name: "Unpaired probabilities table" }));
    expect(unpaired.getAllByText("0.2500")).toHaveLength(9);
    await fireEvent.click(panel.getByRole("radio", { name: /Reverse-complement hypothesis/ }));
    expect(panel.getByRole("radio", { name: /Reverse-complement hypothesis/ })).toBeChecked();
    expect(screen.queryByRole("tab", { name: "Orientation", exact: true })).not.toBeInTheDocument();
  });

  it("keeps repeat-only empty motifs and missing context valid", async () => {
    const value = job("repeats");
    value.summary.repeats.molecule = "RNA";
    value.summary.repeats.instances = [
      {
        ...instance(),
        array_id: null,
        repeat_index: null,
        start: null,
        end: null,
        context: null,
        isolated: {
          status: "completed",
          structure: ".........",
          mfe_kcal_mol: 0,
          ensemble_kcal_mol: 0,
          pairs: [],
          motif_pairs: [],
          unpaired: Array(9).fill(1),
        },
        comparison: {
          status: "unsupported",
          reason: "no_reference_pairs",
          motif_pairs: [],
          isolated_pair_support: null,
          context_pair_support: null,
          support_change: null,
        },
      },
    ];
    render(ResultsView, { props: { job: value } });
    expect(screen.queryByRole("tab", { name: "Arrays", exact: true })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Repeat evidence", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Repeat evidence", exact: true }));
    expect(panel.getByText(/No reference pairs.*Folding still completed/)).toBeInTheDocument();
    expect(panel.getByText("As supplied (RNA)")).toBeInTheDocument();
    expect(panel.getByText("No genomic interval")).toBeInTheDocument();
    expect(panel.getByText(/lower triangle is unavailable/)).toBeInTheDocument();
    expect(
      within(panel.getByRole("region", { name: "Selected repeat detail" })).getByText(
        "Support change",
      ).parentElement,
    ).toHaveTextContent("Not available");
    expect(panel.getByRole("region", { name: "Selected repeat detail" })).toHaveTextContent(
      "Not available",
    );
  });

  it("does not convert missing counts and motif support into zero", () => {
    const normalized = normalizePublicResult({
      repeats: { status: "failed", instances: [{ isolated: {}, comparison: {} }] },
    });
    expect(normalized.repeats.counts.instances).toBeNull();
    expect(normalized.repeats.instances[0].comparison.motif_pair_count).toBeNull();
    expect(normalized.repeats.instances[0].comparison.support_change).toBeNull();
    expect(normalized.repeats.instances[0].isolated.pairs_available).toBe(false);
    const value = job();
    value.summary.repeats = { status: "failed" };
    render(ResultsView, { props: { job: value } });
    const panel = within(screen.getByRole("tabpanel", { name: "Overview", exact: true }));
    expect(panel.getByText("Analysis failed")).toBeInTheDocument();
    expect(panel.getAllByText("Not available").length).toBeGreaterThan(0);
    expect(panel.queryByText(/No repeat instances/)).not.toBeInTheDocument();
  });
});
