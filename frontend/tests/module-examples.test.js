import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { fireEvent, render, screen, waitFor } from "@testing-library/vue";
import { afterEach, describe, expect, it, vi } from "vitest";

import { api } from "../src/api.js";
import AnalysisForm from "../src/components/submission/AnalysisForm.vue";
import ResultsView from "../src/components/results/ResultsView.vue";
import {
  EXAMPLE_ASSETS,
  exampleAssetsFor,
  exampleSubmissionMatches,
  validateExampleInput,
  validateExampleSnapshot,
} from "../src/example.js";
import { ANALYSIS_MODES } from "../src/science.js";

const fixtures = Object.values(EXAMPLE_ASSETS).map((asset) => ({
  ...asset,
  resultText: readFileSync(resolve("public", asset.result), "utf8"),
  sequence: readFileSync(resolve("public", asset.input), "utf8"),
  snapshot: JSON.parse(readFileSync(resolve("public", asset.result), "utf8")),
}));
const limits = {
  maxRecords: 20,
  maxBases: 2_000_000,
  maxRecordBases: 2_000_000,
  maxRequestBytes: 2_000_000,
};

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe("source-verified shipped examples", () => {
  it("offers the existing fifteen methods and both alternate input families only", () => {
    const newModes = new Set([
      "locus_compare",
      "prophage_context",
      "insertion_context",
      "defense_context",
      "rna_family",
    ]);
    expect(Object.keys(EXAMPLE_ASSETS).sort()).toEqual(
      [
        ...ANALYSIS_MODES.map((mode) => mode.id).filter((id) => !newModes.has(id)),
        "repeats-rna",
        "tracrrna-v",
      ].sort(),
    );
    expect(Object.keys(EXAMPLE_ASSETS)).toHaveLength(17);
    for (const mode of newModes) expect(exampleAssetsFor(mode)).toBeNull();
    expect(exampleAssetsFor("repeats", { molecule: "RNA" }).key).toBe("repeats-rna");
    expect(exampleAssetsFor("tracrrna", { tracrModelType: "V" }).key).toBe("tracrrna-v");
    expect(exampleAssetsFor("unknown")).toBeNull();
  });

  it.each(fixtures)(
    "binds $key to its exact independently reproduced input",
    async ({ sequence, snapshot }) => {
      expect(validateExampleSnapshot(snapshot)).toBe(snapshot);
      await expect(validateExampleInput(snapshot, sequence)).resolves.toMatchObject({ snapshot });
      await expect(validateExampleInput(snapshot, `${sequence}\n`)).rejects.toThrow(
        /does not match/,
      );
      const submission = {
        sequence,
        filename: snapshot.example.input.filename,
        mode: snapshot.job.mode,
        ...snapshot.example.submission_options,
      };
      expect(exampleSubmissionMatches(snapshot, submission)).toBe(true);
      for (const [key, value] of Object.entries(snapshot.example.submission_options)) {
        expect(
          exampleSubmissionMatches(snapshot, {
            ...submission,
            [key]:
              typeof value === "boolean"
                ? !value
                : typeof value === "number"
                  ? value + 1
                  : "changed",
          }),
        ).toBe(false);
      }
      expect(
        exampleSubmissionMatches(snapshot, { ...submission, tool_options: { changed: {} } }),
      ).toBe(false);
    },
  );

  it.each(fixtures.filter((fixture) => !fixture.key.includes("-")))(
    "loads $key without an online API or a submitted job",
    async ({ input, result, resultText, sequence, snapshot }) => {
      vi.stubGlobal(
        "fetch",
        vi.fn(async (url) => ({
          ok: true,
          text: async () =>
            String(url).endsWith(input)
              ? sequence
              : String(url).endsWith(result)
                ? resultText
                : "unexpected",
        })),
      );
      const submit = vi.spyOn(api, "submit");
      const onExampleLoaded = vi.fn();
      render(AnalysisForm, {
        props: {
          initialMode: snapshot.job.mode,
          service: { state: "offline", modes: [] },
          limits,
          onExampleLoaded,
        },
      });
      await fireEvent.click(screen.getByRole("button", { name: "Load example", exact: true }));
      await waitFor(() => expect(screen.getByRole("textbox")).toHaveValue(sequence));
      const view = screen.getByRole("button", { name: "View precomputed result", exact: true });
      expect(view).toBeEnabled();
      expect(screen.getByText("Ready to view example")).toBeInTheDocument();
      expect(screen.queryByText("Review input requirements")).not.toBeInTheDocument();
      await fireEvent.click(view);
      await waitFor(() => expect(onExampleLoaded).toHaveBeenLastCalledWith(snapshot));
      expect(submit).not.toHaveBeenCalled();
      await fireEvent.update(screen.getByRole("textbox"), `${sequence}\n`);
      expect(screen.queryByRole("button", { name: "View precomputed result" })).toBeNull();
    },
  );

  it.each(fixtures)(
    "renders $key with identical scientific content for stored and ordinary jobs",
    async ({ snapshot }) => {
      const stored = render(ResultsView, {
        props: { job: snapshot.job, exampleSnapshot: snapshot },
      });
      const content = stored.container.textContent;
      const tabs = screen.getAllByRole("tab").map((tab) => tab.textContent);
      stored.unmount();
      const live = render(ResultsView, { props: { job: snapshot.job } });
      expect(live.container.textContent).toBe(content);
      expect(screen.getAllByRole("tab").map((tab) => tab.textContent)).toEqual(tabs);
    },
  );

  it("rejects a capability or local path in a public result", () => {
    const snapshot = structuredClone(fixtures[0].snapshot);
    snapshot.job.access_token = "private";
    expect(() => validateExampleSnapshot(snapshot)).toThrow();
    delete snapshot.job.access_token;
    snapshot.example.source_note = "/home/operator/private-result";
    expect(() => validateExampleSnapshot(snapshot)).toThrow();
  });

  it("does not open a stale result when input changes during its hash verification", async () => {
    const { input, result, resultText, sequence, snapshot } = fixtures.find(
      (row) => row.key === "orientation",
    );
    vi.stubGlobal(
      "fetch",
      vi.fn(async (url) => ({
        ok: true,
        text: async () =>
          String(url).endsWith(input)
            ? sequence
            : String(url).endsWith(result)
              ? resultText
              : "unexpected",
      })),
    );
    const onExampleLoaded = vi.fn();
    render(AnalysisForm, {
      props: { initialMode: "orientation", service: { state: "offline" }, limits, onExampleLoaded },
    });
    await fireEvent.click(screen.getByRole("button", { name: "Load example", exact: true }));
    await waitFor(() => expect(screen.getByRole("textbox")).toHaveValue(sequence));
    const digest = globalThis.crypto.subtle.digest.bind(globalThis.crypto.subtle);
    let finish;
    const held = new Promise((resolve) => {
      finish = resolve;
    });
    vi.spyOn(globalThis.crypto.subtle, "digest").mockImplementationOnce(async (...args) => {
      await held;
      return digest(...args);
    });
    await fireEvent.click(
      screen.getByRole("button", { name: "View precomputed result", exact: true }),
    );
    await fireEvent.update(screen.getByRole("textbox"), `${sequence}\n`);
    finish();
    expect(
      await screen.findByText(
        "The example input or options changed. Review them before continuing.",
      ),
    ).toBeInTheDocument();
    expect(onExampleLoaded).toHaveBeenLastCalledWith(null);
    expect(onExampleLoaded).not.toHaveBeenCalledWith(snapshot);
  });
});
