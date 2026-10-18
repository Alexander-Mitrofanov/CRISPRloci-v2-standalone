import { readFileSync } from "node:fs";
import { resolve } from "node:path";

import { fireEvent, render, screen, waitFor } from "@testing-library/vue";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

import App from "../src/App.vue";
import { api } from "../src/api.js";
import snapshot from "../public/example-result.json";
import { EXAMPLE_FASTA_PATH, exampleAssetsFor } from "../src/example.js";

const exampleFasta = readFileSync(resolve("public", EXAMPLE_FASTA_PATH), "utf8");
const job = { job_id: "0123456789abcdef0123456789abcdef", status: "queued", mode: "detection" };
const token = "a".repeat(43);
const configured = api.configured;

beforeEach(() => {
  api.configured = true;
  vi.spyOn(api, "health").mockResolvedValue({ version: "2.1.0" });
  vi.spyOn(api, "config").mockResolvedValue({ api_version: "v1" });
  vi.spyOn(api, "submit").mockResolvedValue({ job, access_token: token });
  vi.spyOn(api, "getJob").mockResolvedValue(job);
  vi.stubGlobal(
    "fetch",
    vi.fn(async (url) =>
      String(url).endsWith(exampleAssetsFor("orientation").input)
        ? { ok: true, text: async () => exampleFasta }
        : { ok: true, text: async () => JSON.stringify(snapshot) },
    ),
  );
  Object.defineProperty(Element.prototype, "scrollIntoView", {
    configurable: true,
    value: vi.fn(),
  });
});

afterEach(() => {
  api.configured = configured;
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
  window.history.replaceState(null, "", "/");
});

async function openExample() {
  render(App);
  await fireEvent.click(screen.getByRole("link", { name: "Orientation-aware evolution" }));
  await fireEvent.click(screen.getByRole("button", { name: "Load example" }));
  await fireEvent.click(await screen.findByRole("button", { name: "View precomputed result" }));
  await screen.findByRole("heading", { name: "Analysis results" });
}

describe("navigation from an existing example session", () => {
  it("opens a newly submitted job and preserves its recovery fragment", async () => {
    await openExample();
    await fireEvent.click(screen.getByRole("link", { name: "Methods" }));
    await fireEvent.click(screen.getByRole("link", { name: "Detect arrays" }));
    await fireEvent.update(screen.getByRole("textbox"), ">a\nACGT\n");
    const historyLength = window.history.length;
    await fireEvent.click(screen.getByRole("button", { name: "Compute", exact: true }));

    expect(await screen.findByRole("heading", { name: "Queued" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "Cancel job" })).toBeInTheDocument();
    expect(screen.queryByRole("heading", { name: "Analysis results" })).not.toBeInTheDocument();
    expect(screen.queryByRole("textbox", { name: /contigs/i })).not.toBeInTheDocument();
    expect(window.location.hash).toBe(`#job=${job.job_id}.${token}`);
    expect(window.history.length).toBe(historyLength);
    await waitFor(() => expect(api.submit).toHaveBeenCalledOnce());
  });

  it("can view the same prepared example again without reloading it", async () => {
    await openExample();
    await fireEvent.click(screen.getByRole("link", { name: "Methods" }));
    await fireEvent.click(screen.getByRole("link", { name: "Orientation-aware evolution" }));
    const historyLength = window.history.length;
    await fireEvent.click(screen.getByRole("button", { name: "View precomputed result" }));

    expect(await screen.findByRole("heading", { name: "Analysis results" })).toBeInTheDocument();
    expect(api.submit).not.toHaveBeenCalled();
    expect(window.history.length).toBe(historyLength);
  });
});
