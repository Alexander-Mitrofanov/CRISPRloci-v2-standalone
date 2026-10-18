import { fireEvent, render, screen, within } from "@testing-library/vue";
import { describe, expect, it } from "vitest";

import ResultsView from "../src/components/results/ResultsView.vue";
import { cloneResultJob, resultCredential } from "./support/resultFixture.js";

describe("scientific result overview", () => {
  it("supports keyboard navigation without changing the recovery URL", async () => {
    const job = cloneResultJob();
    const { rerender } = render(ResultsView, { props: { job, credential: resultCredential } });
    const originalUrl = window.location.href;
    const overview = screen.getByRole("tab", { name: "Overview", exact: true });
    overview.focus();
    await fireEvent.keyDown(overview, { key: "ArrowLeft" });
    const files = screen.getByRole("tab", { name: "Files & methods", exact: true });
    expect(files).toHaveFocus();
    expect(files).toHaveAttribute("aria-selected", "true");
    await fireEvent.keyDown(files, { key: "Home" });
    expect(overview).toHaveFocus();
    await fireEvent.keyDown(overview, { key: "ArrowRight" });
    const arrays = screen.getByRole("tab", { name: "Arrays", exact: true });
    expect(arrays).toHaveFocus();
    await rerender({ job: { ...job } });
    expect(arrays).toHaveAttribute("aria-selected", "true");
    await fireEvent.keyDown(arrays, { key: "End" });
    expect(files).toHaveFocus();
    await fireEvent.keyDown(files, { key: "ArrowRight" });
    expect(overview).toHaveFocus();
    expect(window.location.href).toBe(originalUrl);
  });

  it.each([
    ["leader", ["Overview", "Leader context", "Arrays", "Files & methods"]],
    ["detection", ["Overview", "Arrays", "Files & methods"]],
    ["cas", ["Overview", "Cas systems", "Files & methods"]],
    ["tracrrna", ["Overview", "tracrRNA", "Files & methods"]],
    [
      "loci",
      ["Overview", "Cas systems", "tracrRNA", "Leader context", "Arrays", "Files & methods"],
    ],
    ["reconstruction", ["Overview", "Arrays", "History", "Files & methods"]],
  ])("shows only relevant result tabs for %s", (mode, labels) => {
    render(ResultsView, { props: { job: { ...cloneResultJob(), mode } } });
    expect(screen.getAllByRole("tab").map((tab) => tab.textContent.trim())).toEqual(labels);
  });

  it("preserves history selections when switching result tabs", async () => {
    render(ResultsView, { props: { job: cloneResultJob() } });
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    await fireEvent.click(screen.getByRole("button", { name: /Reversed spacer order.*lnL/i }));
    await fireEvent.click(screen.getByRole("tab", { name: "Arrays", exact: true }));
    expect(screen.queryByRole("button", { name: /Reversed spacer order.*lnL/i })).toBeNull();
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    expect(screen.getByRole("button", { name: /Reversed spacer order.*lnL/i })).toHaveAttribute(
      "aria-pressed",
      "true",
    );
  });

  it("shows one result panel at a time with accessible tab relationships", async () => {
    render(ResultsView, { props: { job: cloneResultJob(), credential: resultCredential } });
    const navigation = screen.getByRole("tablist", { name: "Result sections" });
    for (const name of ["Overview", "Arrays", "Orientation", "History", "Files & methods"]) {
      const tab = within(navigation).getByRole("tab", { name, exact: true });
      await fireEvent.click(tab);
      expect(tab).toHaveAttribute("aria-selected", "true");
      const panels = screen.getAllByRole("tabpanel");
      expect(panels).toHaveLength(1);
      expect(panels[0]).toHaveAttribute("id", tab.getAttribute("aria-controls"));
      expect(panels[0]).toHaveAttribute("aria-labelledby", tab.id);
    }
  });

  it("keeps model scores optional and never formats them as probabilities", async () => {
    render(ResultsView, { props: { job: cloneResultJob(), credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "Arrays", exact: true }));
    expect(
      screen.queryByRole("columnheader", { name: "Raw CRISPRidentify Model score" }),
    ).toBeNull();
    await fireEvent.click(screen.getByRole("checkbox", { name: "Show model scores" }));
    expect(screen.getByText("Raw CRISPRidentify Model score")).toBeInTheDocument();
    expect(screen.getByText(/not a calibrated probability/i)).toBeInTheDocument();
    expect(document.body.textContent).not.toMatch(/81\.67%/);
  });

  it("connects detector members, canonical repeat, and evolutionary outputs", async () => {
    render(ResultsView, { props: { job: cloneResultJob(), credential: resultCredential } });
    expect(screen.getAllByText("example_record_01").length).toBeGreaterThan(0);
    await fireEvent.click(screen.getByRole("tab", { name: "Files & methods", exact: true }));
    await fireEvent.click(screen.getByText("Filtering & evidence", { selector: "summary" }));
    expect(screen.getByRole("img", { name: /Canonical repeat CGGTTCAT/i })).toHaveAttribute(
      "tabindex",
      "0",
    );
    expect(screen.getAllByText("42 acquisitions · 4 deletions").length).toBeGreaterThan(0);
  });

  it("surfaces preflight counts, no-deletion caveats, warnings, and tree policy", () => {
    const job = cloneResultJob();
    job.summary.warnings = [
      {
        code: "rho_bias_fit_unavailable_used_uncorrected",
        title: "Bias fallback",
        message: "uncorrected estimates retained",
      },
    ];
    job.summary.orientation.selected_reconstructions[0]["nb of reconstructed deletions"] = 0;
    render(ResultsView, { props: { job, credential: resultCredential } });
    expect(screen.getByText("Unknown strand excluded")).toBeInTheDocument();
    expect(screen.getByText(/No deletion events were reconstructed/i)).toBeInTheDocument();
    expect(screen.getAllByText(/estimated separately/i).length).toBeGreaterThan(0);
    expect(screen.getByText("Bias fallback")).toBeInTheDocument();
  });

  it("treats no eligible groups as informative success", () => {
    render(ResultsView, {
      props: {
        job: { ...cloneResultJob(), status: "completed_no_eligible_groups" },
        credential: resultCredential,
      },
    });
    expect(screen.getByText("No eligible evolutionary groups")).toBeInTheDocument();
    expect(screen.getByText(/Detection completed/i)).toBeInTheDocument();
    expect(screen.queryByRole("heading", { name: /Spacer orientation/i })).not.toBeInTheDocument();
  });
});
