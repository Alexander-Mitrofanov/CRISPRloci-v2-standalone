import { describe, expect, it } from "vitest";

import publicApiContract from "../src/contracts/public-api-v1.json";
import {
  ANALYSIS_MODES,
  CATEGORY_POLICIES,
  DEFAULT_CATEGORY_POLICY,
  TERMINAL_STATUSES,
} from "../src/science.js";

describe("public API contract", () => {
  it("covers every frontend analysis mode and terminal status", () => {
    expect(new Set(ANALYSIS_MODES.map((mode) => mode.id))).toEqual(
      new Set(publicApiContract.enums.analysis_modes),
    );
    expect([...TERMINAL_STATUSES]).toEqual(publicApiContract.enums.terminal_statuses);
    expect(CATEGORY_POLICIES).toEqual(publicApiContract.enums.category_policies);
    expect(CATEGORY_POLICIES).toEqual(["bona_fide_possible", "bona_fide_only"]);
    expect(DEFAULT_CATEGORY_POLICY).toBe(publicApiContract.enums.category_policies[0]);
  });
});
