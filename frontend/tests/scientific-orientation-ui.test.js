import { fireEvent, render, screen } from "@testing-library/vue";
import { describe, expect, it } from "vitest";

import ResultsView from "../src/components/results/ResultsView.vue";
import { cloneResultJob, resultCredential } from "./support/resultFixture.js";

describe("scientific orientation and history presentation", () => {
  it("shows orientation boundary semantics and both likelihood hypotheses", async () => {
    render(ResultsView, { props: { job: cloneResultJob(), credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    await fireEvent.click(screen.getByRole("tab", { name: "Orientation", exact: true }));
    expect(screen.getByText("1 decisive · 0 unresolved")).toBeInTheDocument();
    expect(screen.getByText(/threshold is an evidence rule/i)).toHaveTextContent(
      /not a p-value or probability/i,
    );
    expect(screen.getAllByText("Input spacer order").length).toBeGreaterThan(0);
    expect(screen.getAllByText("Reversed spacer order").length).toBeGreaterThan(0);
    const orientationPlot = screen.getByRole("img", { name: /Delta log likelihood 15\.50/i });
    expect(orientationPlot.tagName.toLowerCase()).toBe("svg");
    expect(
      Number(orientationPlot.querySelector(".orientation-marker").getAttribute("x1")),
    ).toBeGreaterThan(50);
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    const eventRibbon = screen.getByRole("img", { name: "42 acquisitions and 4 deletions" });
    expect(
      Number(eventRibbon.querySelector(".event-ribbon-gains").getAttribute("width")),
    ).toBeCloseTo((42 / 46) * 100, 5);
    expect(
      Number(eventRibbon.querySelector(".event-ribbon-losses").getAttribute("width")),
    ).toBeCloseTo((4 / 46) * 100, 5);
    const modelGauge = screen.getByRole("img", {
      name: /likelihood-ratio statistic 0, cutoff 3\.841.*Preferred model IDM/i,
    });
    expect(modelGauge.querySelector("svg")).toBeInTheDocument();
    expect(modelGauge.querySelector(".model-value")).toHaveAttribute("x1", "0");
    expect(Number(modelGauge.querySelector(".model-cutoff").getAttribute("x1"))).toBeGreaterThan(
      60,
    );
    expect(screen.getByText("LRT statistic").parentElement).toHaveTextContent("0");
    expect(screen.getByText("IDM lnL").parentElement).toHaveTextContent(
      "−23.745".replace("−", "-"),
    );
    expect(document.querySelectorAll("[style]")).toHaveLength(0);
  });

  it("switches structured input and reverse histories and exposes node events", async () => {
    render(ResultsView, { props: { job: cloneResultJob(), credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    const input = screen.getByRole("button", { name: /Input spacer order.*lnL/i });
    const reverse = screen.getByRole("button", { name: /Reversed spacer order.*lnL/i });
    expect(input).toHaveAttribute("aria-pressed", "true");
    await fireEvent.click(reverse);
    expect(reverse).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByText("45", { selector: ".history-summary strong" })).toBeInTheDocument();
    const node = screen.getAllByRole("button", { name: /Inspect node example_record_04/i })[0];
    await fireEvent.click(node);
    expect(document.querySelector(".ancestor-state strong")).toHaveTextContent("example_record_04");
    expect(
      screen.getByText("Comparison history", { selector: ".history-summary strong" }),
    ).toBeInTheDocument();
  });

  it("preserves topology, canvas, exact-branch, and loss-count inspection", async () => {
    render(ResultsView, { props: { job: cloneResultJob(), credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    const topology = screen.getByRole("button", { name: "Readable topology" });
    const branchScale = screen.getByRole("button", { name: "Shared branch scale" });
    expect(topology).toHaveAttribute("aria-pressed", "true");
    const branchNode = screen.getAllByRole("button", {
      name: /Inspect node example_record_03/i,
    })[0];
    const topologyX = Number(branchNode.querySelector("circle").getAttribute("cx"));
    await fireEvent.click(branchScale);
    expect(branchScale).toHaveAttribute("aria-pressed", "true");
    expect(Number(branchNode.querySelector("circle").getAttribute("cx"))).not.toBe(topologyX);

    await fireEvent.click(screen.getByRole("button", { name: "Fit overview" }));
    expect(document.querySelector(".history-canvas")).toHaveClass("is-fit");
    await fireEvent.click(branchNode);
    expect(document.querySelector(".ancestor-state small")).toHaveTextContent("−3 losses");

    await fireEvent.click(screen.getByText("Exact node and branch data"));
    expect(screen.getByRole("columnheader", { name: "Parent" })).toBeInTheDocument();
    expect(screen.getByRole("columnheader", { name: "Branch length" })).toBeInTheDocument();
    expect(
      screen.getByRole("row", { name: /example_record_03.*Inner3.*Observed leaf.*0\.00049/i }),
    ).toBeInTheDocument();
  });

  it("exposes every special-event category in exact node data", async () => {
    const job = cloneResultJob();
    const root = job.summary.orientation.reconstructions
      .find((entry) => entry.hypothesis === "input")
      .nodes.find((node) => node.name === "Inner4");
    Object.assign(root, {
      contradictions: [7],
      duplications: [8],
      rearrangements: [9],
      reacquisitions: [10],
      independent_gains: [11],
      other_duplication_events: [12],
    });
    render(ResultsView, { props: { job, credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    await fireEvent.click(screen.getByText("Exact node and branch data"));
    const rootRow = screen.getByRole("row", { name: /Inner4.*Inferred root/i });
    expect(rootRow).toHaveTextContent("contradictions: 7");
    expect(rootRow).toHaveTextContent("duplications: 8");
    expect(rootRow).toHaveTextContent("rearrangements: 9");
    expect(rootRow).toHaveTextContent("reacquisitions: 10");
    expect(rootRow).toHaveTextContent("independent gains: 11");
    expect(rootRow).toHaveTextContent("other duplications: 12");
  });

  it("caps oversized histories and falls back for malformed structured trees", async () => {
    const oversizedJob = cloneResultJob();
    const input = oversizedJob.summary.orientation.reconstructions.find(
      (entry) => entry.hypothesis === "input",
    );
    const leaves = Array.from({ length: 41 }, (_, index) => `leaf_${index + 1}`);
    input.newick = `(${leaves.map((name) => `${name}:0.1`).join(",")})Root:0;`;
    input.nodes = [
      { name: "Root", spacers: [], gains: [], loss_blocks: [] },
      ...leaves.map((name) => ({ name, spacers: [], gains: [], loss_blocks: [] })),
    ];
    const { unmount } = render(ResultsView, {
      props: { job: oversizedJob, credential: resultCredential },
    });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    expect(
      screen.getByText("Structured history available in the result bundle"),
    ).toBeInTheDocument();
    expect(
      screen.getByText(/inline browser is capped at 40 leaves and 120 nodes/i),
    ).toBeInTheDocument();
    unmount();

    const malformedJob = cloneResultJob();
    malformedJob.summary.orientation.reconstructions.find(
      (entry) => entry.hypothesis === "input",
    ).newick = "(not-valid";
    render(ResultsView, { props: { job: malformedJob, credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    expect(
      screen.getByRole("img", {
        name: /Supported history model tree.*repeat_cfad14d7184c_group_001/i,
      }),
    ).toBeInTheDocument();
  });

  it("never substitutes an available comparison history for a missing supported history", async () => {
    const job = cloneResultJob();
    const comparison = job.summary.orientation.comparisons[0];
    comparison.forward_ln_likelihood_bdm = -39.24813645021139;
    comparison.reverse_ln_likelihood_bdm = -23.745096131593357;
    comparison.forward_minus_reverse_ln_likelihood_bdm = -15.503040318618034;
    comparison.prediction = "Reverse";
    comparison.recommended_reverse = true;
    job.summary.orientation.reconstructions = job.summary.orientation.reconstructions.filter(
      (entry) => entry.hypothesis === "input",
    );

    render(ResultsView, { props: { job, credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));

    expect(
      screen.getByText(
        "Reversed spacer order is supported, but its structured reconstruction is unavailable.",
      ),
    ).toBeInTheDocument();
    expect(
      screen.getByText(/available hypothesis is shown for inspection only/i),
    ).toHaveTextContent(/not substituted for the missing reported history/i);
    expect(
      screen.getByText("Inspection only", { selector: ".history-summary strong" }),
    ).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /Reversed spacer order/i })).toBeDisabled();
    expect(
      screen.queryByText("Supported history", { selector: ".history-summary strong" }),
    ).not.toBeInTheDocument();
  });
});
