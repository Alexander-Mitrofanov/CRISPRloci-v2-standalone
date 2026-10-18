import { describe, expect, it } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/vue";
import { normalizePublicResult } from "../src/features/results/index.js";
import { normalizeRepeatMap } from "../src/features/results/model/repeatMap.js";
import ResultsView from "../src/components/results/ResultsView.vue";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import DefenseContext from "../src/components/results/annotations/DefenseContext.vue";
import ViralResults from "../src/components/results/protospacer/ViralResults.vue";
import fixture from "./fixtures/repeat-map.json";

describe("council scientific contract regressions", () => {
  it.each([
    [0, 1],
    [1, 50],
    [3, 100],
  ])("retains effective map distance %s and %s reported ties", async (distance, topK) => {
    const job = structuredClone(fixture);
    const doc = job.summary.repeat_map;
    doc.policy = { max_distance: distance, top_k: topK };
    doc.results = [
      {
        ...doc.results[0],
        best_hit_count: 120,
        reported_hit_count: topK,
        hits_truncated: true,
        hits: Array.from({ length: topK }, (_, i) => ({
          reference_id: `ref_${i}`,
          orientation: "forward",
          query_aligned: "ACGT",
          reference_aligned: "ACGT",
        })),
      },
    ];
    const normalized = normalizeRepeatMap(doc);
    expect(normalized).toMatchObject({ maxDistance: distance, topK });
    expect(normalized.results[0]).toMatchObject({
      hitCount: 120,
      reportedHitCount: topK,
      truncated: true,
      previewTruncated: false,
    });
    expect(normalized.results[0].hits).toHaveLength(topK);
    render(ResultsView, { props: { job } });
    await fireEvent.click(screen.getByRole("tab", { name: "Repeat matches", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Repeat matches", exact: true }));
    expect(panel.getByText(new RegExp(`within ${distance} global`))).toBeInTheDocument();
    expect(panel.getByText(new RegExp(`^ref_${topK - 1} ·`))).toBeInTheDocument();
  });
  it("distinguishes browser preview loss from native tie truncation", () => {
    const doc = structuredClone(fixture.summary.repeat_map);
    doc.results[0].reported_hit_count = 101;
    doc.results[0].hits_truncated = false;
    doc.results[0].hits = Array.from({ length: 101 }, () => ({ reference_id: "r" }));
    const result = normalizeRepeatMap(doc).results[0];
    expect(result.hits).toHaveLength(100);
    expect(result.previewTruncated).toBe(true);
    expect(result.truncated).toBe(false);
  });
  it("preserves numeric and existing string defense strands and missingness", () => {
    const doc = normalizePublicResult({
      defense_context: {
        components_truncated: true,
        genes: [1, -1, "+", "-", null, 0].map((strand) => ({ strand })),
      },
    }).defense_context;
    expect(doc.genes.map((g) => g.strand)).toEqual(["+", "-", "+", "-", null, null]);
    render(DefenseContext, { props: { ...doc, available: true } });
    expect(screen.getByText(/Component preview is limited/)).toBeInTheDocument();
  });
  it("retains bounded RNA pair-type objects without inventing absent counts", () => {
    const doc = normalizePublicResult({
      rna_family: {
        pairs: [
          { pair_types: { GC: 2, AU: 0, "G-": 1, "--": 1 } },
          {},
          { pair_types: ["GC"] },
          { pair_types: { XX: 1 } },
          { pair_types: { GU: -1 } },
          { pair_types: { GU: 51 } },
        ],
      },
    }).rna_family;
    expect(doc.pairs.map((p) => p.pair_types)).toEqual([
      { GC: 2, AU: 0, "G-": 1, "--": 1 },
      null,
      null,
      null,
      null,
      null,
    ]);
  });
  it("states joint training and eligibility before prophage submission", () => {
    render(AnalysisForm, {
      props: {
        initialMode: "prophage_context",
        service: { state: "online", modes: ["prophage_context"] },
        limits: {},
      },
    });
    expect(screen.getByText(/Prodigal trains jointly/i)).toBeInTheDocument();
    expect(screen.getByText(/unrelated.*unchanged/i)).toBeInTheDocument();
    expect(screen.getByText(/20,000 nt.*eligibility/i)).toBeInTheDocument();
  });
  it("keeps viral overview concise while retaining interpretation and detailed evidence", async () => {
    const props = {
      available: true,
      completed: true,
      compact: true,
      outcome: "matches",
      reference: { name: "Ref", source_release: "r1", record_count: 10 },
      counts: { eligible_query_occurrences: 2, query_occurrences: 3 },
      substitutions: 3,
      truncated: true,
      matches: [{ query_id: "q", subject_accession: "ref", subject_title: "title" }],
      candidates: [{ accession: "ref" }],
      queries: [{ id: "skip", eligible: false, reasons: ["ambiguous_query"] }],
    };
    const view = render(ViralResults, { props });
    expect(screen.getByText(/BLAST retrieval is heuristic/)).toBeInTheDocument();
    expect(screen.getByText(/Ref · release r1/)).toBeInTheDocument();
    expect(screen.queryByRole("table")).not.toBeInTheDocument();
    expect(screen.queryByText("Evidence by viral reference")).not.toBeInTheDocument();
    await view.rerender({ compact: false });
    expect(screen.getByRole("table")).toBeInTheDocument();
    expect(screen.getByText("Evidence by viral reference")).toBeInTheDocument();
  });
});
