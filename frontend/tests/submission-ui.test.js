import { fireEvent, render, screen } from "@testing-library/vue";
import { afterEach, describe, expect, it, vi } from "vitest";

import arrayPacket from "./fixtures/array-compare-input.json";
import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import { ANALYSIS_MODES } from "../src/science.js";

const service = { state: "online", expiresHours: 72, modes: ANALYSIS_MODES.map((mode) => mode.id) };
const limits = {
  maxRecords: 20,
  maxBases: 2_000_000,
  maxRecordBases: 2_000_000,
  maxRequestBytes: 2_000_000,
  maxHeaderCharacters: 200,
};
const renderForm = (props = {}) => render(AnalysisForm, { props: { service, limits, ...props } });

afterEach(() => vi.restoreAllMocks());

describe("Vue submission policy", () => {
  it.each(ANALYSIS_MODES)("submits the method selected for the $id page", async ({ id }) => {
    const submitSpy = vi.spyOn(api, "submit").mockResolvedValue({
      job_id: "0123456789abcdef0123456789abcdef",
      access_token: "a".repeat(43),
      status: "queued",
    });
    renderForm({ initialMode: id });
    await fireEvent.update(
      screen.getByRole("textbox"),
      id === "array_compare"
        ? JSON.stringify(arrayPacket)
        : id === "rna_family"
          ? `>a\n${"ACGT".repeat(8)}\n>b\n${"GCGT".repeat(8)}\n`
          : ">a\nACGT\n>b\nACGT\n",
    );
    if (id === "spacer_association") {
      const compute = screen.getByRole("button", { name: "Compute", exact: true });
      expect(compute).toBeDisabled();
      await fireEvent.update(screen.getByLabelText("Spacer grouping"), "per_record");
    }
    await fireEvent.click(
      screen.getByRole("button", {
        name: ["protospacer", "viral_search"].includes(id) ? "Search references" : "Compute",
        exact: true,
      }),
    );
    expect(submitSpy).toHaveBeenCalledWith(expect.objectContaining({ mode: id }));
  });

  it("exposes a concise form heading hierarchy", () => {
    renderForm();

    expect(
      screen.getByRole("heading", { level: 1, name: "Annotate a CRISPR locus" }),
    ).toBeInTheDocument();
    expect(screen.getByRole("heading", { level: 3, name: "Sequences" })).toBeInTheDocument();
    expect(screen.getByRole("heading", { level: 3, name: "Options" })).toBeInTheDocument();
    expect(
      screen.queryByRole("radio", { name: /Orientation-aware evolution/ }),
    ).not.toBeInTheDocument();
  });

  it("presents advanced analysis policy as a native disclosure button", async () => {
    renderForm();
    const disclosure = screen.getByText("Analysis options").closest("summary");
    const details = disclosure.closest("details");

    expect(disclosure.tagName).toBe("SUMMARY");
    expect(details).not.toHaveAttribute("open");
    expect(disclosure).toHaveTextContent(/Show options/i);
    await fireEvent.click(disclosure);
    expect(details).toHaveAttribute("open");
    expect(disclosure).toHaveTextContent(/Hide options/i);
  });

  it("puts operator visibility, token scope, retention, and sensitive-data limits before submit", () => {
    renderForm();
    const notice = screen.getByRole("note", { name: /sequence privacy and retention/i });
    expect(notice).toHaveTextContent(/service operator/i);
    expect(notice).toHaveTextContent(/not end-to-end encryption/i);
    expect(notice).toHaveTextContent(/results are available as soon as the analysis completes/i);
    expect(notice).not.toHaveTextContent(/72 hours/i);
    expect(notice).toHaveTextContent(/each job shows its expiry time/i);
    expect(notice).toHaveTextContent(
      /personal, clinical, controlled, or unpublished sensitive sequences/i,
    );
    expect(notice).toHaveTextContent(/institutionally approved private route/i);
  });

  it("explains the cohort gate, bounded service, and correction fallback", async () => {
    renderForm({ initialMode: "orientation" });
    expect(
      screen.getByText(/public-service cohort policy requires 2 records/i),
    ).toBeInTheDocument();
    expect(screen.getByText(/at least two comparable detected arrays/i)).toBeInTheDocument();
    expect(
      screen.getByText(/larger genomes require an institutional batch route/i),
    ).toBeInTheDocument();
    const explanation = screen.getByText(/refine α\/ρ deletion-parameter estimates/i);
    expect(explanation).toHaveTextContent(/orientation ΔlnL is unchanged/i);
    expect(explanation).toHaveTextContent(/records a warning instead of failing the job/i);
    expect(explanation).toHaveTextContent(/strict all-corrections-or-error/i);
  });

  it("blocks a second submission while a credential is active", () => {
    renderForm({ hasActiveJob: true });
    expect(screen.getByText(/another job is open/i).closest("div")).toHaveAttribute(
      "role",
      "status",
    );
    expect(screen.getByRole("button", { name: /current job still open/i })).toBeDisabled();
  });

  it("synchronously rejects rapid duplicate submissions", async () => {
    const submitSpy = vi.spyOn(api, "submit").mockReturnValue(new Promise(() => {}));
    renderForm();
    await fireEvent.update(
      screen.getByLabelText(/contigs or small genomes/i),
      ">a\nACGT\n>b\nACGT\n",
    );
    const form = screen.getByRole("button", { name: /compute/i }).closest("form");
    form.dispatchEvent(new Event("submit", { bubbles: true, cancelable: true }));
    form.dispatchEvent(new Event("submit", { bubbles: true, cancelable: true }));
    form.dispatchEvent(new Event("submit", { bubbles: true, cancelable: true }));
    expect(submitSpy).toHaveBeenCalledTimes(1);
  });

  it("rejects an oversized file before reading it", async () => {
    const text = vi.fn();
    renderForm({ limits: { ...limits, maxRequestBytes: 512 } });
    await fireEvent.change(screen.getByLabelText("Upload FASTA file"), {
      target: { files: [{ name: "large.fasta", size: 1_000, text }] },
    });
    expect(text).not.toHaveBeenCalled();
    expect(screen.getByRole("alert")).toHaveTextContent(/selected file is/i);
  });

  it("enforces the JSON request limit for pasted input", async () => {
    renderForm({ limits: { ...limits, maxRequestBytes: 120 } });
    await fireEvent.update(
      screen.getByLabelText(/contigs or small genomes/i),
      `>a\n${"A".repeat(80)}`,
    );
    expect(screen.getByText(/request fits the upload limit/i).closest("li")).toHaveClass(
      "not-ready",
    );
  });
});
