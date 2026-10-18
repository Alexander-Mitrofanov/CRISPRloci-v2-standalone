import { describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen } from "@testing-library/vue";
import { inspectSequenceInput } from "../src/features/submission/sequenceInput.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import { api } from "../src/api.js";

const packet = {
  schema: { name: "crisprloci.sequence-set", version: "1.0.0" },
  kind: "repeat",
  molecule: "DNA",
  complete: true,
  records: [{ id: "repeat_1", sequence: "ACGTCGAT" }],
};

describe("sequence JSON handoffs", () => {
  it("previews the same repeat port for folding and mapping", () => {
    for (const mode of ["repeats", "repeat_map"]) {
      const result = inspectSequenceInput(JSON.stringify(packet), { mode });
      expect(result.valid).toBe(true);
      expect(result.baseCount).toBe(8);
      expect(result.recordCount).toBe(1);
    }
  });
  it.each([
    { complete: false },
    { kind: "target" },
    { schema: { version: "2.0.0" } },
    { records: [] },
    { molecule: "RNA" },
    { records: [{ id: "bad\n>header", sequence: "ACGT" }] },
    { records: [packet.records[0], packet.records[0]] },
  ])("rejects incompatible and incomplete ports: %j", (change) => {
    expect(
      inspectSequenceInput(JSON.stringify({ ...packet, ...change }), { mode: "repeats" }).valid,
    ).toBe(false);
  });
  it("keeps spacer queries distinct from target DNA", () => {
    const text = JSON.stringify({ ...packet, kind: "spacer" });
    expect(inspectSequenceInput(text, { mode: "viral_search" }).valid).toBe(true);
    expect(inspectSequenceInput(text, { mode: "protospacer" }).valid).toBe(false);
  });
  it("submits the original JSON so the server validates it and retains origins", async () => {
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "a".repeat(32), access_token: "a".repeat(43) });
    try {
      render(AnalysisForm, {
        props: {
          initialMode: "repeats",
          service: { state: "online", modes: ["repeats"] },
          limits: { maxRecords: 10, maxBases: 1000, maxRecordBases: 200, maxRequestBytes: 10000 },
        },
      });
      const text = JSON.stringify(packet);
      await fireEvent.update(screen.getByRole("textbox"), text);
      expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeEnabled();
      await fireEvent.click(screen.getByRole("button", { name: "Compute", exact: true }));
      expect(submit).toHaveBeenCalledWith(
        expect.objectContaining({ mode: "repeats", sequence: text }),
      );
    } finally {
      submit.mockRestore();
    }
  });
});
