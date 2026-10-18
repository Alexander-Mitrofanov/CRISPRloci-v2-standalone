import { describe, expect, it } from "vitest";

import { buildSubmission } from "../src/submission.js";

describe("public job submission", () => {
  it("emits only the backend's bounded, allowlisted scientific fields", () => {
    const payload = buildSubmission({
      sequence: ">a\nACGT\n>b\nACGT\n",
      filename: "isolates.fasta",
      mode: "orientation",
      options: {
        categoryPolicy: "bona_fide_only",
        spacerEditDistance: 2,
        biasCorrection: false,
        ignoredArbitraryArgs: ["--unsafe"],
      },
    });

    expect(payload).toEqual({
      sequence: ">a\nACGT\n>b\nACGT\n",
      filename: "isolates.fasta",
      mode: "orientation",
      category_policy: "bona_fide_only",
      spacer_distance: 2,
      bias_corrections: false,
    });
    expect(payload).not.toHaveProperty("args");
    expect(payload).not.toHaveProperty("tree_policy");
  });
});
