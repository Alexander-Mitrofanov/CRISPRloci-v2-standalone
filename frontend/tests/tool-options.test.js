import { describe, expect, it } from "vitest";
import { buildSubmission } from "../src/submission.js";
import {
  selectedToolOptions,
  toolOptionError,
  toolsForMode,
  validToolOptions,
} from "../src/toolOptions.js";

describe("shipped tool flag payloads", () => {
  it("keeps only changed, known flags used by the selected workflow", () => {
    const result = selectedToolOptions("detection", {
      crispridentify: {
        "--min_repeats": 5,
        "--fast_run": true,
        "--cpu": 999,
        "--unknown": "arbitrary input",
      },
      crisprmap: { "--max-distance": 1 },
      arbitrary: { "--exec": "shell" },
    });
    expect(result).toEqual({ crispridentify: { "--min_repeats": 5 } });
  });

  it("omits flags that do not affect the selected CLI mode", () => {
    expect(
      selectedToolOptions("repeats", {
        crisprrepeat: { "--temperature": 25, "--window": 100, "--span": 50 },
      }),
    ).toEqual({ crisprrepeat: { "--temperature": 25 } });
    expect(
      toolsForMode("repeat_context")
        .find((tool) => tool.id === "crisprrepeat")
        .flags.map((field) => field.flag),
    ).toEqual(expect.arrayContaining(["--window", "--span"]));
  });

  it("removes dependent flags from the payload when their feature is disabled", () => {
    expect(
      selectedToolOptions("detection", {
        crispridentify: { "--fast_run": false, "--fast_run_seed_profile": "sensitive" },
      }),
    ).toEqual({ crispridentify: { "--fast_run": false } });
  });

  it.each(["", 1.5, 0, 100, Number.NaN, Number.POSITIVE_INFINITY])(
    "rejects an invalid repeat count %s before submission",
    (value) => {
      expect(validToolOptions("detection", { crispridentify: { "--min_repeats": value } })).toBe(
        false,
      );
    },
  );

  it("preserves numeric choice types and rejects arbitrary choices", () => {
    const numericChoice = { type: "enum", default: 11, choices: [1, 4, 11] };
    expect(toolOptionError(numericChoice, 4)).toBe("");
    expect(toolOptionError(numericChoice, "4")).not.toBe("");
    expect(validToolOptions("cas", { casandra: { "--gene-mode": "made-up" } })).toBe(false);
  });

  it("checks related limits against both explicit values and effective defaults", () => {
    expect(validToolOptions("detection", { crispridentify: { "--min_len_rep": 60 } })).toBe(false);
    expect(validToolOptions("detection", { crispridentify: { "--max_len_spacer": 10 } })).toBe(
      false,
    );
    expect(validToolOptions("repeat_context", { crisprrepeat: { "--window": 50 } })).toBe(false);
    expect(
      validToolOptions("repeat_context", { crisprrepeat: { "--window": 50, "--span": 40 } }),
    ).toBe(true);
  });

  it("puts actual CLI flag overrides into the public job request", () => {
    const request = buildSubmission({
      sequence: ">a\nACGT\n",
      filename: "input.fasta",
      mode: "repeat_map",
      options: {
        toolOptions: { crisprmap: { "--max-distance": 0, "--top-k": 4 } },
      },
    });
    expect(request.tool_options).toEqual({ crisprmap: { "--max-distance": 0, "--top-k": 4 } });
  });
});
