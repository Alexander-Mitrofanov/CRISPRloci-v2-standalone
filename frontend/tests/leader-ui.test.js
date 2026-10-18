import { fireEvent, render, screen, within } from "@testing-library/vue";
import { describe, expect, it } from "vitest";
import ResultsView from "../src/components/results/ResultsView.vue";
import AdvancedOptions from "../src/components/submission/AdvancedOptions.vue";
import { buildSubmission } from "../src/submission.js";
import { normalizePublicResult } from "../src/features/results/index.js";
import { stagesForMode } from "../src/science.js";

const job = (mode = "leader") => ({
  job_id: "leader-example",
  mode,
  status: "completed",
  summary: {
    crisprleader: {
      status: "completed",
      array_count: 1,
      context_count: 2,
      flank_length: 37,
      truncated: true,
      contexts: [
        {
          id: "array:left",
          source_id: "contig",
          array_id: "array",
          side: "left",
          strand: "+",
          segments: [],
          available_length: 0,
          requested_length: 37,
          sequence: "",
          observed_repeat: "ACGT",
          truncated_reason: "unknown_topology_boundary",
        },
        {
          id: "array:right",
          source_id: "contig",
          array_id: "array",
          side: "right",
          strand: "-",
          segments: [{ start: 100, end: 137 }],
          available_length: 37,
          requested_length: 37,
          sequence: "ACGT",
          observed_repeat: "ACGT",
          truncated_reason: null,
        },
      ],
    },
  },
});

describe("CRISPRleader context", () => {
  it.each(["leader", "loci"])("submits and displays the window setting for %s", async (mode) => {
    const payload = buildSubmission({
      sequence: ">a\nACGT\n",
      mode,
      options: { leaderFlankLength: 37 },
    });
    expect(payload.leader_flank_length).toBe(37);
    const { emitted } = render(AdvancedOptions, {
      props: { mode, modelValue: { leaderFlankLength: 500 } },
    });
    await fireEvent.update(screen.getByLabelText(/Leader context window \(nt per side\)/), "37");
    expect(emitted()["update:modelValue"][0][0].leaderFlankLength).toBe(37);
    expect(screen.getByText(/Leader prediction is unavailable/)).toBeInTheDocument();
    expect(stagesForMode(mode).some((stage) => stage.id === "extract_leader_context")).toBe(true);
  });

  it.each(["leader", "loci"])("shows context-only evidence in a dedicated %s tab", async (mode) => {
    render(ResultsView, { props: { job: job(mode) } });
    await fireEvent.click(screen.getByRole("tab", { name: "Leader context", exact: true }));
    const panel = screen.getByRole("tabpanel", { name: "Leader context", exact: true });
    expect(within(panel).getByText(/Leader prediction is unavailable/)).toBeInTheDocument();
    const table = within(panel).getByRole("region", { name: "Leader context table" });
    expect(within(table).getByText("[100, 137)")).toBeInTheDocument();
    expect(within(table).getByText("0 / 37")).toBeInTheDocument();
    expect(within(panel).getByText(/preview is limited/)).toBeInTheDocument();
    expect(within(panel).getByText(/Empty windows remain/)).toBeInTheDocument();
    expect(screen.queryByRole("tab", { name: "Orientation", exact: true })).not.toBeInTheDocument();
  });

  it("does not turn missing counts into zero or claim leader absence", async () => {
    const value = job();
    value.summary.crisprleader = { status: "completed" };
    render(ResultsView, { props: { job: value } });
    const panel = screen.getByRole("tabpanel", { name: "Overview", exact: true });
    expect(within(panel).getAllByText(/Not available/).length).toBeGreaterThan(0);
    expect(screen.queryByText(/No accepted arrays were available/)).not.toBeInTheDocument();
    const normalized = normalizePublicResult({
      crisprleader: { contexts: [null, { available_length: null }] },
    });
    expect(normalized.crisprleader.contexts[0].available_length).toBeNull();
  });
});

it("keeps scientific warnings visible and retains routine publication notes with the files", async () => {
  const value = job();
  value.summary.warnings = [
    "No supported transcription direction",
    "Worker log omitted because it was empty: attempt-001.stderr.log",
    "Bundle file omitted because it exceeds the size limit: important.json",
  ];
  render(ResultsView, { props: { job: value } });
  const warnings = screen.getByRole("list", { name: "Analysis warnings" });
  expect(within(warnings).getAllByRole("listitem")).toHaveLength(2);
  expect(within(warnings).getByText(/No supported transcription direction/)).toBeInTheDocument();
  expect(within(warnings).getByText(/exceeds the size limit/)).toBeInTheDocument();
  await fireEvent.click(screen.getByRole("tab", { name: "Files & methods", exact: true }));
  await fireEvent.click(screen.getByText("Publication notes (1)", { exact: true }));
  expect(screen.getByText(/Worker log omitted/)).toBeVisible();
});
