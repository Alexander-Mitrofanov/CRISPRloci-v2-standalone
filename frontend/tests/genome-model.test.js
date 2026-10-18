import { webcrypto } from "node:crypto";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { buildGenomeModel } from "../src/features/genome/index.js";
import { genomeFixture } from "./support/genomeFixture.js";

beforeEach(() => vi.stubGlobal("crypto", webcrypto));
afterEach(() => vi.unstubAllGlobals());
const build = ({ result, input }) => buildGenomeModel(result, input, result.module);

describe("genome evidence contract", () => {
  it("rejects track density beyond the renderer's row limit", async () => {
    const fixture = genomeFixture();
    fixture.result.data.cas.genes = Array.from({ length: 1001 }, (_, index) => ({
      id: `gene-${index}`,
      source_id: "alpha",
      start: 1,
      end: 300,
      strand: "+",
    }));
    await expect(build(fixture)).rejects.toThrow(/1,000 features overlap/);
  });
  it("aligns inclusive arrays with half-open flanks, preserving occurrence identity and unknowns", async () => {
    const model = await build(genomeFixture());
    const arrays = model.tracks.find((t) => t.id === "arrays").features;
    expect(arrays.map((f) => [f.chr, f.start, f.end, f.strand])).toEqual([
      ["record_1", 421, 580, null],
      ["record_2", 421, 580, null],
    ]);
    expect(arrays[0].end - arrays[0].start).toBe(159);
    expect(arrays[0].details).toContainEqual(["Raw model score", "0"]);
    expect(arrays[1].details).toContainEqual(["Raw model score", "Not reported"]);
    const flanks = model.tracks.find((t) => t.id === "flanks").features;
    expect(flanks.map((f) => [f.start, f.end, f.strand])).toEqual([
      [384, 421, null],
      [580, 617, null],
    ]);
    expect(flanks[1].details).toContainEqual(["Sequence extraction strand", "-"]);
    expect(model.tracks.find((t) => t.id === "cas").features[0].strand).toBe("-");
    expect(new Set(model.features.map((f) => f.key)).size).toBe(model.features.length);
  });

  it("keeps single-base features, deletion boundaries and clipped empty windows distinct", async () => {
    const fixture = genomeFixture();
    fixture.result.data.cas.genes[0].start = 1;
    fixture.result.data.cas.genes[0].end = 1;
    fixture.result.data.detection.arrays[0].spacers.push(
      ...[0, 580, 2000].map((position, index) => ({
        ordinal: index + 3,
        sequence: null,
        status: "deletion",
        boundary_position: position,
      })),
    );
    fixture.result.data.flanks.contexts[0].segments = [{ start: 0, end: 0 }];
    const model = await build(fixture);
    expect(model.features.find((f) => f.trackId === "cas")).toMatchObject({ start: 0, end: 1 });
    expect(model.features.filter((f) => f.point).map((f) => [f.start, f.end])).toEqual([
      [0, 0],
      [580, 580],
      [2000, 2000],
    ]);
    expect(model.notes.join(" ")).toContain("1 empty or unlocated");
  });

  it("retains unverified occurrence evidence without asserting an observed sequence", async () => {
    const fixture = genomeFixture();
    Object.assign(fixture.result.data.detection.arrays[0].spacers[0], {
      sequence: null,
      status: "not_source_verified",
    });
    const model = await build(fixture);
    const spacer = model.features.find((f) => f.trackId === "spacers");
    expect(spacer).toMatchObject({ start: 451, end: 485 });
    expect(spacer.details).toContainEqual([
      "Evidence",
      "Occurrence sequence is not source-verified",
    ]);
    expect(spacer.details).toContainEqual(["Sequence status", "not source verified"]);
  });

  it.each([
    [
      "wrong role",
      (f) => {
        f.input.kind = "spacer";
      },
    ],
    [
      "RNA",
      (f) => {
        f.input.molecule = "RNA";
      },
    ],
    [
      "incomplete input",
      (f) => {
        f.input.complete = false;
      },
    ],
    [
      "incomplete result",
      (f) => {
        f.result.complete = false;
      },
    ],
    [
      "wrong module",
      (f) => {
        f.result.module = "viral_search";
      },
    ],
    [
      "unknown version",
      (f) => {
        f.result.schema.version = "2.0.0";
      },
    ],
    [
      "different source",
      (f) => {
        f.result.data.sources[0].sha256 = "f".repeat(64);
      },
    ],
    [
      "bad length",
      (f) => {
        f.result.data.sources[0].length--;
      },
    ],
    [
      "wrong coordinates",
      (f) => {
        f.result.data.detection.coordinate_system = "unknown";
      },
    ],
    [
      "outside source",
      (f) => {
        f.result.data.cas.genes[0].end = 2001;
      },
    ],
    [
      "unknown source",
      (f) => {
        f.result.data.cas.genes[0].source_id = "other";
      },
    ],
  ])("rejects %s instead of showing a misleading empty map", async (_name, change) => {
    const fixture = genomeFixture();
    change(fixture);
    await expect(build(fixture)).rejects.toThrow();
  });

  it("requires unambiguous record matching, even when two inputs contain identical DNA", async () => {
    const fixture = genomeFixture();
    fixture.input.records[1].sequence = fixture.input.records[0].sequence;
    fixture.result.data.sources[1].sha256 = fixture.result.data.sources[0].sha256;
    await expect(build(fixture)).resolves.toMatchObject({
      sources: [{ id: "alpha" }, { id: "beta" }],
    });
    fixture.result.data.sources[0].id = "renamed";
    await expect(build(fixture)).rejects.toThrow(/unambiguously/);
  });

  it("supports empty assessed Cas results using the complete input inventory", async () => {
    const fixture = genomeFixture();
    fixture.result.module = "cas";
    fixture.result.data = { ...fixture.result.data.cas, genes: [] };
    const model = await build(fixture);
    expect(model.sources).toHaveLength(2);
    expect(model.features).toEqual([]);
    expect(model.tracks.map((t) => t.id)).toEqual(["cassettes", "cas"]);
  });

  it("joins defense profile evidence to genomic genes without treating protein positions as genomic coordinates", async () => {
    const fixture = genomeFixture();
    fixture.result.module = "defense_context";
    fixture.result.data = {
      sources: fixture.result.data.sources,
      policy: { coordinates: "one_based_inclusive" },
      genes: [{ id: "gene", source_id: "alpha", start: 101, end: 400, strand: -1 }],
      systems: [{ id: "system", source_id: "alpha", start: 101, end: 400, type: "RM" }],
      components: [{ gene_id: "gene", system_id: "system", profile: "RM", score: 0 }],
      hits: [{ gene_id: "gene", profile: "candidate", protein_start: 4, protein_end: 20 }],
    };
    const model = await build(fixture);
    expect(model.features.find((f) => f.trackId === "hits")).toMatchObject({
      start: 100,
      end: 400,
      strand: "-",
    });
  });

  it("preserves palindromic match strand and search completeness", async () => {
    const fixture = genomeFixture();
    fixture.result.module = "protospacer";
    fixture.result.data = {
      coordinate_system: "1-based-inclusive-target-forward",
      targets: fixture.result.data.sources,
      search_complete: false,
      policy: { pam: "not_evaluated" },
      matches: [
        {
          spacer_id: "s",
          target_id: "alpha",
          target_start: 1,
          target_end: 16,
          relative_strand: "both",
        },
      ],
    };
    const { features } = await build(fixture);
    expect(features[0]).toMatchObject({ start: 0, end: 16, strand: null });
    expect(features[0].details).toContainEqual(["Relative strand", "both"]);
    expect(features[0].details).toContainEqual(["Search complete", "false"]);
  });
});
