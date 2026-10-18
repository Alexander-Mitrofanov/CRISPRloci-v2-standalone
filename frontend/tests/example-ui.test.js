import { readFileSync } from "node:fs";
import { resolve } from "node:path";

import { fireEvent, render, screen, waitFor } from "@testing-library/vue";
import { afterEach, describe, expect, it, vi } from "vitest";

import snapshot from "../public/example-result.json";
import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import { EXAMPLE_FASTA_PATH, exampleAssetsFor } from "../src/example.js";
import { inspectFasta } from "../src/fasta.js";

const exampleFasta = readFileSync(resolve(process.cwd(), "public", EXAMPLE_FASTA_PATH), "utf8");
const service = { state: "online", expiresHours: 72 };
const limits = {
  maxRecords: 20,
  maxBases: 2_000_000,
  maxRecordBases: 2_000_000,
  maxHeaderCharacters: 200,
  maxRequestBytes: 2_000_000,
};

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe("verified local example", () => {
  it("loads bound assets and reveals the result only after explicit viewing", async () => {
    const fetchMock = vi.fn(async (url) =>
      String(url).endsWith(exampleAssetsFor("orientation").input)
        ? { ok: true, text: async () => exampleFasta }
        : String(url).endsWith(exampleAssetsFor("orientation").result)
          ? { ok: true, text: async () => JSON.stringify(snapshot) }
          : Promise.reject(new Error("Unexpected URL")),
    );
    vi.stubGlobal("fetch", fetchMock);
    const submitSpy = vi.spyOn(api, "submit");
    const onExampleLoaded = vi.fn();
    const onSubmitted = vi.fn();
    render(AnalysisForm, {
      props: { service, limits, initialMode: "orientation", onExampleLoaded, onSubmitted },
    });
    await fireEvent.click(screen.getByRole("button", { name: "Increase spacer edit distance" }));
    await fireEvent.click(screen.getByRole("button", { name: "Load example" }));
    const input = screen.getByLabelText(/contigs or small genomes/i);
    await waitFor(() => expect(input).toHaveValue(exampleFasta));
    expect(onExampleLoaded).toHaveBeenLastCalledWith(null);
    expect(inspectFasta(input.value)).toMatchObject({
      valid: true,
      recordCount: 5,
      baseCount: 8_380,
    });
    expect(document.getElementById("edit-distance")).toHaveTextContent("1");
    await fireEvent.click(screen.getByRole("button", { name: "View precomputed result" }));
    await waitFor(() => expect(onExampleLoaded).toHaveBeenLastCalledWith(snapshot));
    expect(submitSpy).not.toHaveBeenCalled();
    expect(onSubmitted).not.toHaveBeenCalled();
  });

  it("uses the live submission path after the example is edited", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn(async (url) =>
        String(url).endsWith(exampleAssetsFor("orientation").input)
          ? { ok: true, text: async () => exampleFasta }
          : { ok: true, text: async () => JSON.stringify(snapshot) },
      ),
    );
    const submitSpy = vi.spyOn(api, "submit").mockResolvedValue({
      job_id: "0123456789abcdef0123456789abcdef",
      access_token: "a".repeat(43),
      status: "queued",
    });
    const onSubmitted = vi.fn();
    render(AnalysisForm, { props: { service, limits, initialMode: "orientation", onSubmitted } });
    await fireEvent.click(screen.getByRole("button", { name: "Load example" }));
    const input = screen.getByLabelText(/contigs or small genomes/i);
    await waitFor(() => expect(input).toHaveValue(exampleFasta));
    await fireEvent.update(input, `${exampleFasta.trimEnd()}A\n`);
    await fireEvent.click(screen.getByRole("button", { name: "Compute" }));
    await waitFor(() => expect(onSubmitted).toHaveBeenCalledOnce());
    expect(submitSpy).toHaveBeenCalledOnce();
  });

  it("fails closed when the stored input does not match its result", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn(async (url) =>
        String(url).endsWith(exampleAssetsFor("orientation").input)
          ? { ok: true, text: async () => `${exampleFasta}\n` }
          : { ok: true, text: async () => JSON.stringify(snapshot) },
      ),
    );
    render(AnalysisForm, { props: { service, limits, initialMode: "orientation" } });
    await fireEvent.click(screen.getByRole("button", { name: "Load example" }));
    expect(await screen.findByText(/does not match its precomputed result/i)).toBeInTheDocument();
    expect(screen.getByLabelText(/contigs or small genomes/i)).toHaveValue("");
  });

  it("renders the cached result through the ordinary Vue result components", async () => {
    render(ResultsView, { props: { job: snapshot.job, exampleSnapshot: snapshot } });
    expect(screen.getByText("Completed")).toBeInTheDocument();
    expect(screen.getAllByText(/example_record_/).length).toBeGreaterThanOrEqual(5);
    await fireEvent.click(screen.getByRole("tab", { name: "Files & methods", exact: true }));
    await fireEvent.click(screen.getByText("Filtering & evidence", { selector: "summary" }));
    expect(screen.getByRole("heading", { name: /What reached the model/i })).toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "History", exact: true }));
    expect(screen.getByRole("heading", { name: "Spacer history" })).toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Files & methods", exact: true }));
    expect(
      screen.getByRole("heading", { name: /How detections became evolutionary evidence/i }),
    ).toBeInTheDocument();
    await fireEvent.click(screen.getByRole("tab", { name: "Overview", exact: true }));
    expect(screen.getByRole("heading", { name: "Run summary" })).toBeInTheDocument();
    expect(screen.getAllByText("Input order supported").length).toBeGreaterThan(0);
    expect(
      screen.queryByRole("button", { name: /Download complete result bundle/i }),
    ).not.toBeInTheDocument();
  });
});
