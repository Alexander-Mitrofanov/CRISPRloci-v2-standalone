import { describe, expect, it } from "vitest";

import { normalizePublicResult } from "../src/features/results/index.js";

describe("public result model", () => {
  it("normalizes legacy aliases once for all presentation slices", () => {
    const result = normalizePublicResult({
      categories: { "Bona-fide": 1 },
      detected_arrays: [
        {
          Category: "Bona-fide",
          "Confidence score": 0.75,
          End: 20,
          Name: "array-1",
          "Number of spacers": 3,
          Start: 5,
          Strand: "+",
          record_id: "record-1",
        },
      ],
      orientation_evidence: {
        groups: [
          {
            delta_lnL: 8,
            forward_ln_likelihood: -10,
            group_id: "group-1",
            reverse_ln_likelihood: -18,
            threshold: 5,
          },
        ],
        selected_reconstructions: [
          {
            gains: 7,
            group: "group-1",
            leaf_count: 4,
            losses: 2,
            log_likelihood: -10,
            runtime_seconds: 1.5,
            unique_spacers: 9,
          },
        ],
      },
      preflight: {
        eligible_groups: 1,
        excluded_arrays: 2,
        retained_arrays: 4,
        unknown_strand_arrays: 1,
      },
      spacerplacer: {
        results: [{ group: "group-1", insertions: 7, losses: 2 }],
      },
    });

    expect(result.detection).toMatchObject({
      category_counts: { "Bona-fide": 1 },
      arrays: [
        {
          array_id: "array-1",
          category: "Bona-fide",
          end: 20,
          model_score: 0.75,
          source_id: "record-1",
          spacer_count: 3,
          start: 5,
          strand: "+",
        },
      ],
    });
    expect(result.adapter).toMatchObject({
      emitted_array_count: 4,
      emitted_group_count: 1,
      skipped_array_count: 2,
      unknown_strand_excluded_count: 1,
    });
    expect(result.orientation.comparisons[0]).toMatchObject({
      confidence_threshold: 5,
      forward_ln_likelihood_bdm: -10,
      forward_minus_reverse_ln_likelihood_bdm: 8,
      group: "group-1",
      reverse_ln_likelihood_bdm: -18,
    });
    expect(result.orientation.selected_reconstructions[0]).toMatchObject({
      acquisitions: 7,
      bdm_log_likelihood: -10,
      deletions: 2,
      leaf_count: 4,
      name: "group-1",
      runtime_seconds: 1.5,
      unique_spacers: 9,
    });
    expect(result.reconstruction.results[0]).toMatchObject({
      acquisitions: 7,
      deletions: 2,
      name: "group-1",
    });
  });

  it("preserves missing and null scientific values instead of inventing zeroes", () => {
    const result = normalizePublicResult({
      detection: { arrays: [{ array_id: "array-1", model_score: null }] },
      orientation: {
        selected_reconstructions: [{ name: "group-1", losses: null }],
      },
    });

    expect(result.detection.arrays[0].model_score).toBeNull();
    expect(result.orientation.selected_reconstructions[0].deletions).toBeNull();
    expect(result.orientation.selected_reconstructions[0].acquisitions).toBeUndefined();
    expect(result.reconstruction).toBeNull();
  });

  it("normalizes job-level provenance and warnings before presentation", () => {
    const result = normalizePublicResult(
      {
        provenance: { versions: { SpacerPlacer: "1.2.3" } },
        warnings: [{ code: "summary-warning" }],
      },
      {
        options: { category_policy: "bona-fide", warnings: [{ code: "option-warning" }] },
        request: { options: { category_policy: "fallback" } },
        warnings: [{ code: "job-warning" }],
      },
    );

    expect(result.provenance).toMatchObject({
      parameters: { category_policy: "bona-fide", warnings: [{ code: "option-warning" }] },
      tool_versions: { SpacerPlacer: "1.2.3" },
    });
    expect(result.warnings.map((warning) => warning.code)).toEqual([
      "summary-warning",
      "option-warning",
      "job-warning",
    ]);
  });

  it("is stable when canonical normalized data is normalized again", () => {
    const once = normalizePublicResult({
      adapter: { emitted_array_count: 2, groups: [] },
      detection: { arrays: [], category_counts: {} },
      orientation: { comparisons: [] },
    });

    expect(normalizePublicResult(once)).toEqual(once);
  });
});
