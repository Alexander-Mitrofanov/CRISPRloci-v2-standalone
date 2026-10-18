import { afterEach, describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/vue";
import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import { normalizePublicResult } from "../src/features/results/index.js";
import fixture from "./fixtures/repeat-type.json";
const service = { state: "online", modes: ["repeat_type"] };
const limits = {
  maxRecords: 2000,
  maxBases: 1000000,
  maxRecordBases: 10000,
  maxRequestBytes: 1000000,
};
afterEach(() => vi.restoreAllMocks());
describe("RepeatTyper", () => {
  it("submits bounded DNA repeat inputs only when provisioned", async () => {
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "c".repeat(32), access_token: "a".repeat(43) });
    const view = render(AnalysisForm, { props: { initialMode: "repeat_type", service, limits } });
    const input = screen.getByLabelText("Repeat sequences");
    const button = screen.getByRole("button", { name: "Compute", exact: true });
    await fireEvent.update(input, ">r\n" + "A".repeat(201));
    expect(button).toBeDisabled();
    await fireEvent.update(
      input,
      Array.from({ length: 1001 }, (_, i) => `>r${i}\nACGT\n`).join(""),
    );
    expect(button).toBeDisabled();
    await fireEvent.update(input, ">r\nACGT\n");
    expect(button).toBeEnabled();
    await fireEvent.click(button);
    expect(submit.mock.calls[0][0]).toMatchObject({
      mode: "repeat_type",
      sequence: ">r\nACGT\n",
    });
    await view.rerender({ service: { state: "online", modes: ["detection"] } });
    expect(button).toBeDisabled();
  });
  it("retains model winner, unknowns, unsupported bases and filters", async () => {
    render(ResultsView, { props: { job: structuredClone(fixture) } });
    expect(screen.queryByRole("tab", { name: "Arrays", exact: true })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Repeat subtypes", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Repeat subtypes", exact: true }));
    expect(panel.getByText(/not calibrated confidence/)).toBeInTheDocument();
    expect(panel.getByText(/winning label falls below/)).toBeInTheDocument();
    expect(panel.getByText(/no prediction was made/)).toBeInTheDocument();
    await fireEvent.update(panel.getByLabelText("Filter displayed repeats"), "ambiguous");
    expect(panel.getByText("ambiguous · unsupported")).toBeInTheDocument();
    expect(panel.queryByText("duplicate · predicted")).not.toBeInTheDocument();
  });
  it("never converts missing or invalid scientific scores into zeros", () => {
    const doc = structuredClone(fixture);
    doc.summary.repeat_type.results[0].native_score = true;
    const normalized = normalizePublicResult(doc.summary, doc).repeat_type;
    expect(normalized.results[0].nativeScore).toBeNull();
    expect(normalized.results[3].nativeScore).toBeNull();
    render(ResultsView, { props: { job: { ...fixture, summary: { repeat_type: null } } } });
    expect(
      within(screen.getByRole("tabpanel", { name: "Overview", exact: true })).getByText(
        "A completed repeat classification was not reported.",
      ),
    ).toBeInTheDocument();
  });
});
