import { fireEvent, render, screen, within } from "@testing-library/vue";
import { describe, expect, it } from "vitest";
import ResultsView from "../src/components/results/ResultsView.vue";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import CasContextRelations from "../src/components/results/annotations/CasContextRelations.vue";
import { normalizePublicResult } from "../src/features/results/index.js";

describe("context evidence scope", () => {
  it.each(["prophage_context", "insertion_context", "defense_context"])(
    "applies the fixed 20-record budget for %s even when the service permits more",
    async (mode) => {
      render(AnalysisForm, {
        props: {
          initialMode: mode,
          service: { state: "online", modes: [mode] },
          limits: {
            maxRecords: 100,
            maxRecordBases: 5_000_000,
            maxBases: 5_000_000,
            maxRequestBytes: 6_000_000,
          },
        },
      });
      await fireEvent.update(
        screen.getByRole("textbox"),
        Array.from({ length: 21 }, (_, i) => `>source${i}\nACGT\n`).join(""),
      );
      expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeDisabled();
      expect(screen.getByText("Limit: 20 records")).toBeInTheDocument();
    },
  );
  it("retains both-strand palindrome evidence", () => {
    const result = normalizePublicResult({
      prophage_context: { spacer_matches: [{ strand: "both" }] },
    });
    expect(result.prophage_context.matches[0].strand).toBe("both");
  });
  it("distinguishes skipped short records from assessed records with no prophage calls", () => {
    render(ResultsView, {
      props: {
        job: {
          job_id: "c".repeat(32),
          status: "completed",
          mode: "prophage_context",
          summary: {
            prophage_context: {
              available: true,
              counts: { sources: 1, assessed_sources: 0, regions: 0, genes: 0, spacer_matches: 0 },
              sources: [{ id: "short", length: 500, status: "ineligible_short_sequence" }],
              regions: [],
              spacer_matches: [],
            },
          },
        },
      },
    });
    const panel = within(screen.getByRole("tabpanel", { name: "Overview", exact: true }));
    expect(panel.getByText(/No source records were assessed for prophages/)).toBeInTheDocument();
    expect(panel.queryByText(/No prophage candidates were reported/)).not.toBeInTheDocument();
    expect(panel.getByText(/0 of 1 source records assessed/)).toBeInTheDocument();
  });

  it("preserves truncated nearest-cassette ties without implying missing cassettes", () => {
    const summary = normalizePublicResult({
      defense_context: {
        cas_context: [
          {
            interval_id: "system",
            source_id: "source",
            nearest_cassette_count: 33,
            nearest_cassettes_truncated: true,
            nearest_cassettes: [
              {
                cassette_id: "cassette",
                start: 1,
                end: 10,
                overlap_bases: 0,
                intervening_bases: 100,
              },
            ],
          },
        ],
      },
    });
    render(CasContextRelations, { props: { relations: summary.defense_context.casContext } });
    expect(screen.getByText(/Nearest-cassette ties are limited/)).toHaveTextContent("33 total");
    expect(screen.queryByText(/No Cas cassette was reported/)).not.toBeInTheDocument();
  });
});
