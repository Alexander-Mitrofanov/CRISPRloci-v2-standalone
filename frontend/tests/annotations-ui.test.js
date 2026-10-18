import { fireEvent, render, screen, within } from "@testing-library/vue";
import { describe, expect, it } from "vitest";
import ResultsView from "../src/components/results/ResultsView.vue";
import { buildSubmission } from "../src/submission.js";

describe("CRISPRloci annotation workflows", () => {
  it.each(["loci", "tracrrna"])("submits the selected Type V-K models for %s", (mode) => {
    const payload = buildSubmission({
      sequence: ">a\nACGT\n",
      mode,
      options: { tracrModelType: "V" },
    });
    expect(payload.tracr_model_type).toBe("V");
    expect(payload.mode).toBe(mode);
  });

  it("presents annotation evidence, genomic context, null values and provenance", async () => {
    const job = {
      mode: "loci",
      status: "completed",
      summary: {
        detection: {
          arrays: [
            { source_id: "contig", array_id: "array-1", start: 10, end: 30, category: "Bona-fide" },
          ],
        },
        casandra: {
          status: "completed",
          cas_gene_count: 1,
          cassette_count: 0,
          cassettes: [],
          genes: [
            {
              id: "cas-1",
              source_id: "contig",
              start: 40,
              end: 100,
              strand: "-",
              profile: "Cas9",
              score_margin: null,
            },
          ],
        },
        tracrrna: {
          status: "completed",
          mode: "complete",
          model_type: "V",
          prediction_count: 1,
          sources: [{ id: "contig", length: 1000 }],
          candidates: [
            {
              id: "rna-1",
              source_id: "contig",
              start: 120,
              end: 160,
              strand: "+",
              score: 0,
              evalue: null,
            },
          ],
        },
        provenance: {
          repositories: [
            {
              display_name: "CasAndra",
              repository: "https://github.com/Alexander-Mitrofanov/CasAndra.git",
              commit: "a".repeat(40),
            },
          ],
        },
      },
    };
    render(ResultsView, { props: { job } });
    expect(screen.getByRole("heading", { name: "Locus annotations" })).toBeInTheDocument();
    expect(screen.getByRole("img", { name: "Genomic features on contig" })).toBeInTheDocument();
    expect(screen.getByText(/Type V-K/)).toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "tracrRNA", exact: true }));
    const table = screen.getByRole("region", { name: "tracrRNA candidate table" });
    expect(within(table).getByText("0")).toBeInTheDocument();
    expect(within(table).getByText("Not available")).toBeInTheDocument();
    expect(screen.queryByRole("heading", { name: /Reconstruction/ })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Files & methods", exact: true }));
    expect(screen.getByRole("link", { name: /CasAndra/ })).toHaveAttribute(
      "href",
      `https://github.com/Alexander-Mitrofanov/CasAndra/commit/${"a".repeat(40)}`,
    );
  });
});
