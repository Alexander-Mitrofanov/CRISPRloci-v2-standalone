import { render, screen } from "@testing-library/vue";
import { describe, expect, it } from "vitest";
import InsertionContextResults from "../src/components/results/InsertionContextResults.vue";

describe("insertion context evidence", () => {
  it("keeps absent ORF and missing copy E-value explicit", () => {
    render(InsertionContextResults, {
      props: {
        elements: [
          {
            id: "is1",
            source_id: "source1",
            start: 1,
            end: 12,
            family: "IS1",
            native_type: "p",
            contig_edge: true,
            transposase_orf: null,
            evalue_this_copy: null,
            terminal_inverted_repeats: null,
          },
        ],
        counts: { elements: 1 },
      },
    });
    expect(screen.getByText("Not identified")).toBeTruthy();
    expect(screen.getByText("Partial prediction")).toBeTruthy();
    expect(screen.getByText(/Target-site duplications were not assessed/)).toBeTruthy();
  });
  it("does not interpret zero calls as absence", () => {
    render(InsertionContextResults, { props: { counts: { elements: 0 } } });
    expect(screen.getByText(/does not establish their absence/)).toBeTruthy();
  });
});
