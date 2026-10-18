import { afterEach, describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/vue";
import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import fixture from "./fixtures/repeat-map.json";
const service = { state: "online", modes: ["repeat_map"] };
const limits = {
  maxRecords: 1000,
  maxBases: 100000,
  maxRecordBases: 10000,
  maxRequestBytes: 200000,
};
afterEach(() => vi.restoreAllMocks());
describe("repeat reference mapping", () => {
  it("submits DNA repeats and enforces bounds and reference availability", async () => {
    const submit = vi
      .spyOn(api, "submit")
      .mockResolvedValue({ job_id: "c".repeat(32), access_token: "a".repeat(43) });
    const view = render(AnalysisForm, { props: { initialMode: "repeat_map", service, limits } });
    const input = screen.getByLabelText("Repeat sequences");
    const button = screen.getByRole("button", { name: "Compute", exact: true });
    await fireEvent.update(input, ">r\n" + "A".repeat(201));
    expect(button).toBeDisabled();
    await fireEvent.update(input, Array.from({ length: 101 }, (_, i) => `>r${i}\nACGT\n`).join(""));
    expect(button).toBeDisabled();
    await fireEvent.update(input, ">r\nACGT\n");
    expect(button).toBeEnabled();
    await fireEvent.click(button);
    expect(submit.mock.calls[0][0]).toMatchObject({ mode: "repeat_map", sequence: ">r\nACGT\n" });
    expect(submit.mock.calls[0][0]).not.toHaveProperty("reference");
    await view.rerender({ service: { state: "online", modes: ["detection"] } });
    expect(button).toBeDisabled();
  });
  it("preserves conflicting, missing, candidate and no-hit evidence", async () => {
    render(ResultsView, { props: { job: structuredClone(fixture) } });
    expect(screen.queryByRole("tab", { name: "Arrays", exact: true })).not.toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Repeat matches", exact: true }));
    const panel = within(screen.getByRole("tabpanel", { name: "Repeat matches", exact: true }));
    expect(panel.getAllByText(/Annotation evidence includes every best tie/)).toHaveLength(3);
    expect(panel.getAllByText(/Family: Unassigned/).length).toBeGreaterThan(0);
    expect(panel.getByText(/No hit does not establish biological novelty/)).toBeInTheDocument();
    await fireEvent.update(panel.getByLabelText("Filter displayed repeats"), "negative");
    expect(panel.getByText("negative · no_hit")).toBeInTheDocument();
    expect(panel.queryByText("exact · exact")).not.toBeInTheDocument();
  });
  it("does not turn missing output into zero matches", () => {
    render(ResultsView, { props: { job: { ...fixture, summary: { repeat_map: null } } } });
    expect(
      within(screen.getByRole("tabpanel", { name: "Overview", exact: true })).getByText(
        "A completed repeat mapping was not reported.",
      ),
    ).toBeInTheDocument();
  });
});
