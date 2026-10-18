import { fireEvent, render, screen } from "@testing-library/vue";
import { describe, expect, it } from "vitest";

import ResultsView from "../src/components/results/ResultsView.vue";
import { cloneResultJob, resultCredential } from "./support/resultFixture.js";

describe("scientific reconstruction presentation", () => {
  it("keeps absent reconstruction metrics absent and accepts documented aliases", async () => {
    const missingJob = cloneResultJob();
    const missing = missingJob.summary.orientation.selected_reconstructions[0];
    for (const key of [
      "nb of reconstructed insertions",
      "nb of reconstructed deletions",
      "nb of unique spacers",
      "nb of spacers in alignment",
      "nb of unique spacer arrays",
      "nb of leafs (after combining non-uniques)",
      "test_statistic (-2*ln_lh_ratio)",
      "chi2_quantile",
      "ln_lh_idm",
      "ln_lh_bdm",
      "Deletion model preferred by LRT",
      "nb of reconstructed duplications",
      "nb of reconstructed rearrangements",
      "nb of reconstructed reacquisitions",
      "nb of reconstructed independent gains",
    ])
      delete missing[key];
    const { unmount } = render(ResultsView, {
      props: { job: missingJob, credential: resultCredential },
    });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    expect(document.querySelector(".spacerplacer-verdict strong")).toHaveTextContent(
      "— acquisitions · — deletions",
    );
    expect(document.querySelector(".event-graphic .model-unavailable")).toHaveTextContent(
      "totals were not reported",
    );
    expect(document.querySelector(".model-selection .model-unavailable")).toHaveTextContent(
      "Likelihood-ratio statistic not reported",
    );
    expect(document.querySelector(".event-ribbon")).not.toBeInTheDocument();
    unmount();

    const aliasJob = cloneResultJob();
    const alias = aliasJob.summary.orientation.selected_reconstructions[0];
    for (const key of [
      "nb of reconstructed insertions",
      "nb of reconstructed deletions",
      "nb of unique spacers",
      "nb of spacers in alignment",
      "nb of unique spacer arrays",
      "nb of leafs (after combining non-uniques)",
      "test_statistic (-2*ln_lh_ratio)",
      "chi2_quantile",
      "ln_lh_idm",
      "ln_lh_bdm",
      "Deletion model preferred by LRT",
    ])
      delete alias[key];
    Object.assign(alias, {
      gains: 7,
      losses: 2,
      unique_spacers: 9,
      aligned_spacers: 10,
      unique_arrays: 3,
      leaf_count: 4,
      likelihood_ratio_statistic: 6,
      model_selection_cutoff: 3,
      idm_log_likelihood: -10,
      log_likelihood: -7,
      preferred_model: "BDM",
    });
    render(ResultsView, { props: { job: aliasJob, credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    expect(document.querySelector(".spacerplacer-verdict strong")).toHaveTextContent(
      "7 acquisitions · 2 deletions",
    );
    expect(screen.getByRole("img", { name: "7 acquisitions and 2 deletions" })).toBeInTheDocument();
    expect(
      screen.getByRole("img", {
        name: /likelihood-ratio statistic 6, cutoff 3.*Preferred model BDM/i,
      }),
    ).toBeInTheDocument();
    expect(document.querySelector(".inventory-counts")).toHaveTextContent("10 aligned positions");
  });
});
