import { webcrypto } from "node:crypto";
import { fireEvent, render, screen, waitFor, within } from "@testing-library/vue";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import ResultsView from "../src/components/results/ResultsView.vue";
import { api } from "../src/api.js";
import { createGenomeBrowser } from "../src/features/genome/renderer.js";
import { genomeFixture, genomeArtifacts } from "./support/genomeFixture.js";
import { resultCredential } from "./support/resultFixture.js";

vi.mock("../src/features/genome/renderer.js", () => ({ createGenomeBrowser: vi.fn() }));
let viewer;
let read;
let fixture;
function job() {
  return {
    job_id: resultCredential.jobId,
    mode: "loci",
    status: "completed",
    artifacts: genomeArtifacts,
    summary: { tracrrna: { status: "completed", sources: fixture.result.data.sources } },
  };
}
beforeEach(() => {
  vi.stubGlobal("crypto", webcrypto);
  fixture = genomeFixture();
  viewer = {
    search: vi.fn().mockResolvedValue(true),
    visible: vi.fn(),
    dispose: vi.fn(),
    zoomIn: vi.fn(),
    zoomOut: vi.fn(),
  };
  vi.mocked(createGenomeBrowser).mockResolvedValue(viewer);
  read = vi
    .spyOn(api, "readSequenceArtifact")
    .mockImplementation(async (_job, id) =>
      JSON.stringify(id === "result-artifact" ? fixture.result : fixture.input),
    );
});
afterEach(() => {
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
  vi.clearAllMocks();
});

describe("genome result integration", () => {
  it("loads only on demand, reuses one viewer across tabs, and connects feature selection to navigation", async () => {
    const view = render(ResultsView, { props: { job: job(), credential: resultCredential } });
    expect(read).not.toHaveBeenCalled();
    expect(screen.getByText("Genomic context", { exact: true })).toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Genome & locus" }));
    await waitFor(() => expect(screen.getByRole("button", { name: "Zoom in" })).toBeEnabled());
    expect(screen.queryByText("Genomic context", { exact: true })).not.toBeInTheDocument();
    expect(read).toHaveBeenCalledTimes(2);
    expect(read).toHaveBeenCalledWith(
      resultCredential.jobId,
      "result-artifact",
      resultCredential.accessToken,
      expect.objectContaining({ signal: expect.any(AbortSignal), maxBytes: 20_000_000 }),
    );
    const invocation = vi.mocked(createGenomeBrowser).mock.calls[0];
    expect(JSON.stringify(invocation.slice(1, 3))).not.toContain(resultCredential.accessToken);
    await fireEvent.update(screen.getByLabelText("Feature", { exact: true }), "feature_0");
    expect(viewer.search).toHaveBeenLastCalledWith("record_1:342-660");
    const detail = screen.getByRole("region", { name: "Selected genomic feature" });
    expect(detail).toHaveTextContent("422–580 (159 nt)");
    expect(detail).toHaveTextContent("Not established");
    await fireEvent.click(screen.getByRole("button", { name: "Open evidence tables" }));
    expect(screen.getByRole("tab", { name: "Arrays" })).toHaveAttribute("aria-selected", "true");
    await fireEvent.click(screen.getByRole("tab", { name: "Genome & locus" }));
    expect(createGenomeBrowser).toHaveBeenCalledTimes(1);
    expect(read).toHaveBeenCalledTimes(2);
    expect(
      within(screen.getByRole("region", { name: "Selected genomic feature" })).getByRole("heading"),
    ).toHaveTextContent("array-1");
    view.unmount();
    expect(viewer.dispose).toHaveBeenCalledOnce();
  });

  it("does not offer a browser for spacer-only modes or another job's capability", () => {
    render(ResultsView, {
      props: { job: { ...job(), mode: "viral_search" }, credential: resultCredential },
    });
    expect(screen.queryByRole("tab", { name: "Genome & locus" })).not.toBeInTheDocument();
  });

  it("reports incomplete results instead of rendering empty tracks", async () => {
    fixture.result.complete = false;
    render(ResultsView, { props: { job: job(), credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "Genome & locus" }));
    expect(await screen.findByRole("alert")).toHaveTextContent(
      "Complete genome evidence is unavailable",
    );
    expect(createGenomeBrowser).not.toHaveBeenCalled();
    await fireEvent.click(screen.getByRole("tab", { name: "Overview", exact: true }));
    expect(screen.getByText("Genomic context", { exact: true })).toBeVisible();
  });

  it("keeps the existing map if the renderer cannot initialize", async () => {
    vi.mocked(createGenomeBrowser).mockRejectedValueOnce(new Error("Initialization failed"));
    render(ResultsView, { props: { job: job(), credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "Genome & locus" }));
    expect(await screen.findByRole("alert")).toHaveTextContent("could not start");
    await fireEvent.click(screen.getByRole("tab", { name: "Overview", exact: true }));
    expect(screen.getByText("Genomic context", { exact: true })).toBeVisible();
  });

  it("labels source DNA as interval reference when occurrence evidence is unverified", async () => {
    Object.assign(fixture.result.data.detection.arrays[0].spacers[0], {
      sequence: null,
      status: "not_source_verified",
    });
    render(ResultsView, { props: { job: job(), credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "Genome & locus" }));
    await waitFor(() => expect(screen.getByRole("button", { name: "Zoom in" })).toBeEnabled());
    await fireEvent.update(screen.getByLabelText("Feature", { exact: true }), "feature_2");
    const detail = screen.getByRole("region", { name: "Selected genomic feature" });
    expect(detail).toHaveTextContent("Occurrence sequence is not source-verified");
    expect(detail).toHaveTextContent("Reference DNA at this interval (source-forward)");
    expect(detail).not.toHaveTextContent("Observed occurrence");
  });

  it("aborts artifact loading when the job changes and ignores late results", async () => {
    let resolve;
    read.mockImplementationOnce(
      () =>
        new Promise((finish) => {
          resolve = finish;
        }),
    );
    const view = render(ResultsView, { props: { job: job(), credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "Genome & locus" }));
    const signal = read.mock.calls[0][3].signal;
    await view.rerender({ job: { ...job(), job_id: "other" } });
    expect(signal.aborted).toBe(true);
    resolve(JSON.stringify(fixture.result));
    await Promise.resolve();
    expect(read).toHaveBeenCalledTimes(1);
    expect(createGenomeBrowser).not.toHaveBeenCalled();
  });

  it("disposes a viewer that finishes initializing after its job was closed", async () => {
    let resolve;
    vi.mocked(createGenomeBrowser).mockImplementationOnce(
      () =>
        new Promise((finish) => {
          resolve = finish;
        }),
    );
    const view = render(ResultsView, { props: { job: job(), credential: resultCredential } });
    await fireEvent.click(screen.getByRole("tab", { name: "Genome & locus" }));
    await waitFor(() => expect(createGenomeBrowser).toHaveBeenCalledOnce());
    view.unmount();
    resolve(viewer);
    await waitFor(() => expect(viewer.dispose).toHaveBeenCalledOnce());
  });
});
