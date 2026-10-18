import { effectScope, reactive } from "vue";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, waitFor } from "@testing-library/vue";
import { useSequenceHandoff } from "../src/features/results/index.js";
import ports from "../src/contracts/module-ports-v1.json";
import { stagesForMode } from "../src/science.js";
import { validAssociationGrouping } from "../src/features/submission/association.js";
import App from "../src/App.vue";
import { api } from "../src/api.js";

const packet = {
  schema: { name: "crisprloci.sequence-set", version: "1.0.0" },
  kind: "repeat",
  molecule: "DNA",
  complete: true,
  records: [{ id: "repeat_", sequence: "ACGTACGT" }],
};
const credential = { jobId: "a".repeat(32), accessToken: "b".repeat(43) };
const scopes = [];
const configured = api.configured;
function setup(mode = "detection", text = JSON.stringify(packet)) {
  const props = reactive({
    job: {
      job_id: credential.jobId,
      status: "completed",
      mode,
      artifacts: ["input-sequences.json", "repeats.json", "spacers.json", "arrays.json"].map(
        (name, i) => ({ name, artifact_id: String(i), size_bytes: 200 }),
      ),
    },
    credential,
    service: { state: "online", modes: Object.keys(ports) },
    limits: {},
  });
  const emit = vi.fn();
  const client = { readSequenceArtifact: vi.fn().mockResolvedValue(text) };
  const scope = effectScope();
  scopes.push(scope);
  return { props, emit, client, result: scope.run(() => useSequenceHandoff(props, emit, client)) };
}
afterEach(() => {
  scopes.splice(0).forEach((scope) => scope.stop());
  api.configured = configured;
  vi.restoreAllMocks();
  window.history.replaceState(null, "", "/");
});

describe("complete output handoffs", () => {
  it.each(Object.entries(ports))("offers only compatible input consumers for %s", (mode, role) => {
    const { result } = setup(mode);
    const expected = role === "array" ? "spacer" : role;
    expect(
      result.choices.value
        .filter((c) => c.artifact.name === "input-sequences.json")
        .map((c) => c.mode)
        .sort(),
    ).toEqual(
      Object.keys(ports)
        .filter((m) => ports[m] === expected)
        .sort(),
    );
  });
  it("loads exact protected JSON into a draft without submitting a job", async () => {
    const { result, client, emit } = setup();
    result.selected.value = "1:repeat_type";
    await result.prepare();
    expect(client.readSequenceArtifact).toHaveBeenCalledWith(
      credential.jobId,
      "1",
      credential.accessToken,
      expect.objectContaining({ signal: expect.any(AbortSignal) }),
    );
    expect(emit).toHaveBeenCalledWith(
      "continue-analysis",
      expect.objectContaining({
        mode: "repeat_type",
        sequence: JSON.stringify(packet),
        sourceJobId: credential.jobId,
      }),
    );
  });
  it.each([{ complete: false }, { records: [] }, { kind: "target" }])(
    "blocks unusable exports: %j",
    async (change) => {
      const { result, emit } = setup("detection", JSON.stringify({ ...packet, ...change }));
      result.selected.value = "1:repeats";
      await result.prepare();
      expect(emit).not.toHaveBeenCalled();
      expect(result.error.value).not.toBe("");
    },
  );
  it("keeps RNA inputs as RNA and rejects DNA-only consumers", async () => {
    const { result, emit } = setup(
      "repeats",
      JSON.stringify({ ...packet, molecule: "RNA", records: [{ id: "r", sequence: "ACGU" }] }),
    );
    result.selected.value = "0:repeat_map";
    await result.prepare();
    expect(emit).not.toHaveBeenCalled();
    result.selected.value = "0:repeats";
    await result.prepare();
    expect(emit).toHaveBeenCalledWith(
      "continue-analysis",
      expect.objectContaining({ molecule: "RNA" }),
    );
  });
  it("discards late downloads after switching jobs", async () => {
    const { props, client, result, emit } = setup();
    let finish;
    client.readSequenceArtifact.mockImplementation(
      () =>
        new Promise((resolve) => {
          finish = resolve;
        }),
    );
    result.selected.value = "1:repeats";
    const work = result.prepare();
    props.job.job_id = "c".repeat(32);
    finish(JSON.stringify(packet));
    await work;
    expect(emit).not.toHaveBeenCalled();
    expect(client.readSequenceArtifact.mock.calls[0][3].signal.aborted).toBe(true);
  });
  it("does not expose connections from failed jobs or unavailable reference modes", () => {
    const { props, result } = setup();
    props.service.modes = ["repeats"];
    expect(result.choices.value.some((c) => c.mode === "repeat_type")).toBe(false);
    props.job.status = "failed";
    expect(result.choices.value).toEqual([]);
  });
  it("accepts source grouping for a UTF-8 BOM port", () => {
    expect(
      validAssociationGrouping(
        "\ufeff" +
          JSON.stringify({
            records: [{ origin: { source_id: "s", source_sha256: "a".repeat(64) } }],
          }),
        "per_source",
      ),
    ).toBe(true);
  });
  it.each(["detection", "orientation", "reconstruction"])(
    "keeps CCTK out of the %s pipeline",
    (mode) => {
      expect(stagesForMode(mode).some((stage) => stage.id === "cctk")).toBe(false);
      expect(stagesForMode("array_compare").map((stage) => stage.id)).toEqual([
        "queued",
        "validate_input",
        "cctk",
        "package_results",
      ]);
    },
  );

  it("continues a completed job and preserves the old recovery link until explicit submission", async () => {
    Object.defineProperty(Element.prototype, "scrollIntoView", {
      configurable: true,
      value: vi.fn(),
    });
    api.configured = true;
    const original = {
      job_id: credential.jobId,
      status: "completed",
      mode: "detection",
      summary: {},
      artifacts: [{ name: "repeats.json", artifact_id: "repeat-port", size_bytes: 200 }],
    };
    const next = { job_id: "c".repeat(32), status: "queued", mode: "repeats" };
    vi.spyOn(api, "health").mockResolvedValue({ version: "2.11.0" });
    vi.spyOn(api, "config").mockResolvedValue({ modes: Object.keys(ports) });
    vi.spyOn(api, "getJob").mockImplementation((id) =>
      Promise.resolve(id === original.job_id ? original : next),
    );
    vi.spyOn(api, "readSequenceArtifact").mockResolvedValue(JSON.stringify(packet));
    vi.spyOn(api, "submit").mockResolvedValue({ job: next, access_token: "d".repeat(43) });
    window.history.replaceState(null, "", `/#job=${credential.jobId}.${credential.accessToken}`);
    render(App);
    await fireEvent.click(await screen.findByRole("tab", { name: "Files & methods" }));
    await fireEvent.update(
      screen.getByLabelText("Sequence file and next analysis"),
      "repeat-port:repeats",
    );
    await fireEvent.click(screen.getByRole("button", { name: "Use as analysis input" }));
    expect(await screen.findByRole("textbox")).toHaveValue(JSON.stringify(packet));
    expect(window.location.hash).toContain(credential.jobId);
    expect(api.submit).not.toHaveBeenCalled();
    await waitFor(() =>
      expect(screen.getByRole("button", { name: "Compute", exact: true })).toBeEnabled(),
    );
    await fireEvent.click(screen.getByRole("button", { name: "Compute", exact: true }));
    await waitFor(() =>
      expect(api.submit).toHaveBeenCalledWith(
        expect.objectContaining({ mode: "repeats", sequence: JSON.stringify(packet) }),
      ),
    );
    expect(window.location.hash).toBe(`#job=${next.job_id}.${"d".repeat(43)}`);
  });
});
