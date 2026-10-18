import { describe, expect, it } from "vitest";

import { fastaLines, inspectFasta } from "../src/fasta.js";

const LINE_BOUNDARIES = [
  ["CRLF", "\r\n"],
  ["LF", "\n"],
  ["CR", "\r"],
  ["vertical tab", "\v"],
  ["form feed", "\f"],
  ["file separator", "\u001C"],
  ["group separator", "\u001D"],
  ["record separator", "\u001E"],
  ["next line", "\u0085"],
  ["line separator", "\u2028"],
  ["paragraph separator", "\u2029"],
];

describe("FASTA inspection", () => {
  it.each(["repeat_", "repeat-", `${"x".repeat(63)}_suffix`])(
    "preserves canonical identifiers on repeated normalization: %s",
    (id) => {
      const first = inspectFasta(`>${id}\nACGT\n`).records[0].normalizedIdentifier;
      expect(inspectFasta(`>${first}\nACGT\n`).records[0].normalizedIdentifier).toBe(first);
    },
  );
  it("accepts multiple records with unique identifiers and IUPAC DNA", () => {
    const result = inspectFasta(">isolate_A note\nACGTRYSW\n>isolate_B\nNNACGT\n");
    expect(result.valid).toBe(true);
    expect(result.recordCount).toBe(2);
    expect(result.baseCount).toBe(14);
  });

  it.each(LINE_BOUNDARIES)("matches the backend for the %s boundary", (_name, boundary) => {
    expect([...fastaLines(`first${boundary}second`)]).toEqual(["first", "second"]);
    const result = inspectFasta(`>isolate_A${boundary}ac${boundary}>isolate_B${boundary}gt`);
    expect(result.valid).toBe(true);
    expect(result.records.map((record) => record.sequence)).toEqual(["AC", "GT"]);
  });

  it("accepts one FASTA record with every supported boundary mixed", () => {
    const lines = [">mixed", "A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B"];
    const text = lines.map((line, index) => `${line}${LINE_BOUNDARIES[index]?.[1] || ""}`).join("");
    expect([...fastaLines(text)]).toEqual(lines);
    const result = inspectFasta(text);
    expect(result.valid).toBe(true);
    expect(result.records[0].sequence).toBe("ACGTRYSWKMB");
  });

  it("rejects duplicated first-token identifiers", () => {
    const result = inspectFasta(">same first\nACGT\n>same second\nACGT\n");
    expect(result.valid).toBe(false);
    expect(result.errors.join(" ")).toMatch(/duplicated/i);
  });

  it("rejects protein symbols", () => {
    const result = inspectFasta(">x\nMPEPTIDE\n");
    expect(result.valid).toBe(false);
    expect(result.errors.join(" ")).toMatch(/unsupported DNA symbol/i);
  });

  it("accepts raw DNA using the backend-compatible fallback identifier", () => {
    const result = inspectFasta("ACGTRYSWKMBDHVN");
    expect(result.valid).toBe(true);
    expect(result.records[0].normalizedIdentifier).toBe("web_input");
  });

  it("rejects Unicode confusables before ASCII case normalization", () => {
    const result = inspectFasta(">confusable\nA\u017fGT\n");
    expect(result.valid).toBe(false);
    expect(result.errors.join(" ")).toMatch(/unsupported DNA symbol/i);
    expect(result.records[0].sequence).toBe("A\u017fGT");
  });

  it("rejects headers that are too long, contain controls, or collide after sanitization", () => {
    const long = inspectFasta(`>${"x".repeat(201)}\nACGT`);
    expect(long.errors.join(" ")).toMatch(/exceeds 200 characters/i);

    const control = inspectFasta(">bad\u0001header\nACGT");
    expect(control.errors.join(" ")).toMatch(/control characters/i);

    const collision = inspectFasta(">sample/a\nACGT\n>sample?a\nACGT");
    expect(collision.errors.join(" ")).toMatch(/collide after safe filename normalization/i);
  });
});
