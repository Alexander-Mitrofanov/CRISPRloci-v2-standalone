import { fireEvent, render, screen, waitFor } from "@testing-library/vue";
import { afterEach, describe, expect, it, vi } from "vitest";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import snapshot from "../public/example-result.json";
import { EXAMPLE_FASTA_PATH, exampleAssetsFor } from "../src/example.js";
import { ANALYSIS_MODES } from "../src/science.js";

const service = { state: "online", expiresHours: 72, modes: ANALYSIS_MODES.map((mode) => mode.id) };
const limits = {
  maxRecords: 20,
  maxBases: 2_000_000,
  maxRecordBases: 2_000_000,
  maxRequestBytes: 2_000_000,
  maxHeaderCharacters: 200,
};
const renderForm = (initialMode) =>
  render(AnalysisForm, { props: { service, limits, initialMode } });
const mockSubmit = () =>
  vi.spyOn(api, "submit").mockResolvedValue({
    job_id: "0123456789abcdef0123456789abcdef",
    access_token: "a".repeat(43),
    status: "queued",
  });
const openTool = async (name) => {
  await fireEvent.click(screen.getByText("Analysis options").closest("summary"));
  await fireEvent.click(
    screen.getByText(name, { selector: ".tool-options-tool strong" }).closest("summary"),
  );
};

afterEach(() => {
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

describe("compact per-tool CLI settings", () => {
  it.each(ANALYSIS_MODES)("offers a collapsed menu for the $id workflow", ({ id }) => {
    renderForm(id);
    if (
      [
        "locus_compare",
        "prophage_context",
        "insertion_context",
        "defense_context",
        "rna_family",
      ].includes(id)
    ) {
      expect(screen.queryByText("Analysis options")).not.toBeInTheDocument();
      return;
    }
    const root = screen.getByText("Analysis options").closest("details");
    expect(root).not.toHaveAttribute("open");
    expect(root.querySelectorAll(".tool-options-tool").length).toBeGreaterThan(0);
    for (const section of root.querySelectorAll(".tool-options-tool")) {
      expect(section).not.toHaveAttribute("open");
    }
  });

  it("shows exact flag names, effects and defaults, and submits a changed value", async () => {
    const submit = mockSubmit();
    renderForm("detection");
    await openTool("CRISPRidentify");
    const input = screen.getByLabelText("Minimum repeat count --min_repeats");
    expect(input).toHaveValue(3);
    expect(input.closest(".tool-option-field")).toHaveTextContent("Default: 3");
    expect(input).toHaveAccessibleDescription(/Minimum repeats required in a candidate array/);
    await fireEvent.update(input, "5");
    expect(screen.getByText("1 tool flag changed")).toBeInTheDocument();
    await fireEvent.update(screen.getByRole("textbox"), ">a\nACGT\n");
    await fireEvent.click(screen.getByRole("button", { name: "Compute", exact: true }));
    expect(submit).toHaveBeenCalledWith(
      expect.objectContaining({ tool_options: { crispridentify: { "--min_repeats": 5 } } }),
    );
  });

  it("resets a flag and removes the override from submission", async () => {
    const submit = mockSubmit();
    renderForm("repeat_map");
    await openTool("CRISPRmap");
    const input = screen.getByLabelText("Maximum edit distance --max-distance");
    await fireEvent.update(input, "0");
    await fireEvent.click(screen.getByRole("button", { name: "Reset --max-distance to default" }));
    expect(input).toHaveValue(3);
    expect(screen.queryByText("1 tool flag changed")).not.toBeInTheDocument();
    await fireEvent.update(screen.getByRole("textbox"), ">a\nACGT\n");
    await fireEvent.click(screen.getByRole("button", { name: "Compute", exact: true }));
    expect(submit.mock.calls[0][0]).not.toHaveProperty("tool_options");
  });

  it("blocks invalid numerical values with a specific accessible error", async () => {
    const submit = mockSubmit();
    renderForm("detection");
    await openTool("CRISPRidentify");
    await fireEvent.update(screen.getByRole("textbox"), ">a\nACGT\n");
    const input = screen.getByLabelText("Minimum repeat count --min_repeats");
    await fireEvent.update(input, "1");
    expect(input).toHaveAttribute("aria-invalid", "true");
    expect(input).toHaveAccessibleDescription(/Use 2 or more/);
    expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeDisabled();
    expect(submit).not.toHaveBeenCalled();
  });

  it("shows dependent seed settings only when the fast search uses them", async () => {
    renderForm("detection");
    await openTool("CRISPRidentify");
    expect(
      screen.getByLabelText("Fast-search seed profile --fast_run_seed_profile"),
    ).toBeInTheDocument();
    await fireEvent.update(screen.getByLabelText("Fast search --fast_run"), "false");
    expect(
      screen.queryByLabelText("Fast-search seed profile --fast_run_seed_profile"),
    ).not.toBeInTheDocument();
  });

  it("submits the selected CLI mode value", async () => {
    const submit = mockSubmit();
    renderForm("cas");
    await openTool("CasAndra");
    await fireEvent.update(screen.getByLabelText("Gene calling mode --gene-mode"), "meta");
    await fireEvent.update(screen.getByRole("textbox"), ">a\nACGT\n");
    await fireEvent.click(screen.getByRole("button", { name: "Compute", exact: true }));
    expect(submit).toHaveBeenCalledWith(
      expect.objectContaining({ tool_options: { casandra: { "--gene-mode": "meta" } } }),
    );
  });

  it("runs a fresh analysis when a verified example's tool flags change", async () => {
    const example = readFileSync(resolve(process.cwd(), "public", EXAMPLE_FASTA_PATH), "utf8");
    vi.stubGlobal(
      "fetch",
      vi.fn(async (url) =>
        String(url).endsWith(exampleAssetsFor("orientation").input)
          ? { ok: true, text: async () => example }
          : { ok: true, text: async () => JSON.stringify(snapshot) },
      ),
    );
    const submit = mockSubmit();
    renderForm("orientation");
    await fireEvent.click(screen.getByRole("button", { name: "Load example" }));
    await waitFor(() => expect(screen.getByRole("textbox")).toHaveValue(example));
    expect(screen.getByRole("button", { name: "View precomputed result" })).toBeEnabled();
    await openTool("evOr");
    await fireEvent.update(
      screen.getByLabelText("Orientation decision boundary --orientation_decision_boundary"),
      "8",
    );
    expect(
      screen.queryByRole("button", { name: "View precomputed result" }),
    ).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("button", { name: "Compute", exact: true }));
    expect(submit).toHaveBeenCalledWith(
      expect.objectContaining({
        tool_options: { crispr_evor: { "--orientation_decision_boundary": 8 } },
      }),
    );
  });
});
