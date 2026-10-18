import { afterEach, describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/vue";
import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import { normalizePublicResult } from "../src/features/results/index.js";
import fixture from "./fixtures/spacer-association.json";
const service = { state: "online", modes: ["spacer_association"] };
const limits = {
  maxRecords: 1000,
  maxBases: 1000000,
  maxRecordBases: 10000,
  maxRequestBytes: 1000000,
};
afterEach(() => vi.restoreAllMocks());
describe("SpacePHARER", () => {
  it("requires an explicit grouping and a provisioned panel", async () => {
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "c".repeat(32), access_token: "a".repeat(43) });
    const view = render(AnalysisForm, {
      props: { initialMode: "spacer_association", service, limits },
    });
    const input = screen.getByLabelText("Spacer sequences");
    const button = screen.getByRole("button", { name: "Compute", exact: true });
    await fireEvent.update(input, ">r\nACGT\n");
    expect(button).toBeDisabled();
    await fireEvent.update(screen.getByLabelText("Spacer grouping"), "per_source");
    expect(button).toBeDisabled();
    await fireEvent.update(screen.getByLabelText("Spacer grouping"), "single_set");
    expect(button).toBeEnabled();
    await fireEvent.click(button);
    expect(submit.mock.calls[0][0]).toMatchObject({
      mode: "spacer_association",
      association_grouping: "single_set",
    });
    await fireEvent.update(input, ">long\n" + "A".repeat(201));
    expect(button).toBeDisabled();
    await fireEvent.update(input, ">r\nACGT\n");
    await view.rerender({ service: { state: "online", modes: ["detection"] } });
    expect(button).toBeDisabled();
  });
  it("accepts complete origin annotations only for per-source grouping", async () => {
    render(AnalysisForm, { props: { initialMode: "spacer_association", service, limits } });
    await fireEvent.update(screen.getByLabelText("Spacer grouping"), "per_source");
    const packet = {
      schema: { name: "crisprloci.sequence-set", version: "1.0.0" },
      kind: "spacer",
      molecule: "DNA",
      complete: true,
      records: [
        {
          id: "s1",
          sequence: "ACGT",
          origin: {
            source_id: "contigA",
            source_sha256: "a".repeat(64),
            array_id: "a",
            ordinal: 1,
            start: 1,
            end: 4,
          },
        },
      ],
    };
    await fireEvent.update(screen.getByLabelText("Spacer sequences"), JSON.stringify(packet));
    expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeEnabled();
    delete packet.records[0].origin;
    await fireEvent.update(screen.getByLabelText("Spacer sequences"), JSON.stringify(packet));
    expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeDisabled();
  });
  it("shows actual native candidates and retains missing FDR", async () => {
    render(ResultsView, { props: { job: structuredClone(fixture) } });
    expect(screen.queryByRole("tab", { name: "Arrays", exact: true })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Associations", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Associations", exact: true }));
    expect(panel.getByText("FDR calibration unavailable.")).toBeInTheDocument();
    expect(panel.getByText("source1 → AY369265.2")).toBeInTheDocument();
    expect(panel.getAllByText(/Native FDR: Not reported/)).toHaveLength(2);
    expect(panel.getByText(/Phages absent from it cannot be detected/)).toBeInTheDocument();
    const normalized = normalizePublicResult(fixture.summary).spacer_association;
    expect(normalized.associations[0].fdr).toBeNull();
    const doc = structuredClone(fixture.summary);
    doc.spacer_association.associations[0].native_fdr = 0;
    expect(normalizePublicResult(doc).spacer_association.associations[0].fdr).toBeNull();
  });
  it("shows a broad database's reported counts and source without a genome list", async () => {
    const job = structuredClone(fixture);
    job.summary.spacer_association.reference = {
      name: "dated-phage-reference",
      engine_version: "2.1.0",
      control_method: "reversed_translated_orfs",
      manifest_sha256: "b".repeat(64),
      genome_count: 32510,
      total_bases: 2147483648,
      scope: "broad_phage_reference",
      source: {
        collection: "NCBI GenBank and RefSeq",
        release: "Frozen snapshot",
        retrieved_at: "2026-10-10T12:00:00Z",
        selection_policy: "Complete eligible DNA records; exclusions recorded in the snapshot.",
      },
    };
    render(ResultsView, { props: { job } });
    await fireEvent.click(screen.getByRole("tab", { name: "Associations", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Associations", exact: true }));
    expect(
      panel.getByText(/32,510 genomes · 2,147,483,648 bases · Broad phage reference/),
    ).toBeInTheDocument();
    expect(panel.getByText("Collection: NCBI GenBank and RefSeq")).toBeInTheDocument();
    expect(panel.getByText("Release: Frozen snapshot")).toBeInTheDocument();
    expect(panel.getByText("Retrieved: 2026-10-10T12:00:00Z")).toBeInTheDocument();
    expect(panel.getByText(/Selection policy: Complete eligible DNA records/)).toBeInTheDocument();
    expect(panel.queryByText(/at most 100 genomes/)).not.toBeInTheDocument();
    const normalized = normalizePublicResult(job.summary).spacer_association.reference;
    expect(normalized.genomes).toEqual([]);
    expect(normalized.genomeCount).toBe(32510);
  });
  it("leaves missing database counts, scope and provenance unreported", async () => {
    const job = structuredClone(fixture);
    job.summary.spacer_association.reference = { name: "reference-with-missing-metadata" };
    render(ResultsView, { props: { job } });
    await fireEvent.click(screen.getByRole("tab", { name: "Associations", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Associations", exact: true }));
    expect(
      panel.getByText(/Not reported genomes · Not reported bases · Scope not reported/),
    ).toBeInTheDocument();
    for (const label of ["Collection", "Release", "Retrieved", "Selection policy"]) {
      expect(panel.getByText(`${label}: Not reported`)).toBeInTheDocument();
    }
    expect(normalizePublicResult(job.summary).spacer_association.reference).toMatchObject({
      genomeCount: null,
      totalBases: null,
      scope: null,
      source: { collection: null, release: null, retrievedAt: null, selectionPolicy: null },
    });
  });
  it("preserves historical reference inventory and derives only its reported totals", async () => {
    const job = structuredClone(fixture);
    job.summary.spacer_association.reference = {
      name: "historical-panel",
      engine_version: "2.0.2",
      control_method: "reversed_translated_orfs",
      manifest_sha256: "c".repeat(64),
      genomes: [
        { id: "legacyA", length: 1234, sha256: "a".repeat(64) },
        { id: "legacyB", length: 4321, sha256: "b".repeat(64) },
      ],
    };
    render(ResultsView, { props: { job } });
    await fireEvent.click(screen.getByRole("tab", { name: "Associations", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Associations", exact: true }));
    expect(panel.getByText(/2 genomes · 5,555 bases/)).toBeInTheDocument();
    expect(panel.getByText(/legacyA · 1,234 bases · SHA-256/)).toBeInTheDocument();
    expect(panel.getByText("Collection: Not reported")).toBeInTheDocument();
    const normalized = normalizePublicResult(job.summary).spacer_association.reference;
    expect(normalized.genomes).toHaveLength(2);
    expect(normalized.source.collection).toBeNull();
    delete job.summary.spacer_association.reference.genomes[1].length;
    expect(normalizePublicResult(job.summary).spacer_association.reference.totalBases).toBeNull();
  });
  it("retains explicit zero counts and rejects invalid count or provenance values", () => {
    const doc = structuredClone(fixture.summary);
    doc.spacer_association.reference = {
      genome_count: 0,
      total_bases: 0,
      scope: "targeted_panel",
      source: { collection: "", release: "   ", retrieved_at: false, selection_policy: null },
    };
    expect(normalizePublicResult(doc).spacer_association.reference).toMatchObject({
      genomeCount: 0,
      totalBases: 0,
      scope: "targeted_panel",
      source: { collection: null, release: null, retrievedAt: null, selectionPolicy: null },
    });
    doc.spacer_association.reference.genome_count = -1;
    doc.spacer_association.reference.total_bases = Number.MAX_SAFE_INTEGER + 1;
    expect(normalizePublicResult(doc).spacer_association.reference).toMatchObject({
      genomeCount: null,
      totalBases: null,
    });
    doc.spacer_association.reference = { genomes: [] };
    expect(normalizePublicResult(doc).spacer_association.reference).toMatchObject({
      genomeCount: null,
      totalBases: null,
    });
  });
});
