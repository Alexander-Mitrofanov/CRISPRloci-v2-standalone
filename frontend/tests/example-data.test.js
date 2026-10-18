import { createHash } from "node:crypto";
import { existsSync, readFileSync } from "node:fs";
import { resolve } from "node:path";

import { describe, expect, it } from "vitest";

import snapshot from "../public/example-result.json";
import {
  EXAMPLE_FASTA_PATH,
  EXAMPLE_RESULT_PATH,
  EXAMPLE_SCHEMA_VERSION,
  validateExampleInput,
  validateExampleSnapshot,
} from "../src/example.js";
import { inspectFasta } from "../src/fasta.js";

const publicRoot = resolve(process.cwd(), "public");
const fasta = readFileSync(resolve(publicRoot, EXAMPLE_FASTA_PATH), "utf8");

function clone(value) {
  return JSON.parse(JSON.stringify(value));
}

describe("masked stored-input/precomputed-output contract", () => {
  it("binds the public FASTA to the precomputed snapshot by raw and normalized SHA-256", async () => {
    expect(validateExampleSnapshot(snapshot)).toBe(snapshot);
    await expect(validateExampleInput(snapshot, fasta)).resolves.toMatchObject({ snapshot });
    const inspection = inspectFasta(fasta);
    expect(inspection).toMatchObject({ valid: true, recordCount: 5, baseCount: 8_380 });
    expect(createHash("sha256").update(fasta).digest("hex")).toBe(
      snapshot.example.input.file_sha256,
    );
    expect(snapshot.example.input.normalized_sha256).toBe(
      "7657eff497d8ef6b4e5e97dda849f0f60a4488001b52e8201381bba91c6556e5",
    );
    expect(EXAMPLE_SCHEMA_VERSION).toBe("1.3.0");
    expect(snapshot.schema.version).toBe(EXAMPLE_SCHEMA_VERSION);
    expect(EXAMPLE_RESULT_PATH).toBe("example-result.json");
    expect(existsSync(resolve(publicRoot, EXAMPLE_RESULT_PATH))).toBe(true);
    expect(existsSync(resolve(publicRoot, "obsolete-example-result.json"))).toBe(false);
  });

  it("keeps every teaching claim equal to the completed job summary", () => {
    const { findings } = snapshot.example;
    const { detection, adapter, orientation } = snapshot.job.summary;
    const mainReconstruction = orientation.selected_reconstructions.find(
      (row) => row.name === findings.reconstruction.group,
    );
    const mainComparison = orientation.comparisons.find(
      (row) => row.group === findings.orientation.group,
    );
    expect(findings.detection.arrays).toBe(detection.array_count);
    expect(findings.detection.bona_fide).toBe(detection.category_counts["Bona-fide"]);
    expect(findings.detection.possible).toBe(detection.category_counts.Possible ?? 0);
    expect(findings.preflight.modeled_arrays).toBe(adapter.emitted_array_count);
    expect(findings.preflight.eligible_groups).toBe(adapter.emitted_group_count);
    expect(findings.preflight.excluded_arrays).toBe(adapter.skipped_array_count);
    expect(findings.reconstruction.insertions).toBe(
      mainReconstruction["nb of reconstructed insertions"],
    );
    expect(findings.reconstruction.deletions).toBe(
      mainReconstruction["nb of reconstructed deletions"],
    );
    expect(findings.reconstruction.duplications).toBe(
      mainReconstruction["nb of reconstructed duplications"],
    );
    expect(findings.orientation.delta_ln_likelihood).toBe(
      mainComparison.forward_minus_reverse_ln_likelihood_bdm,
    );
    expect(Math.abs(findings.orientation.delta_ln_likelihood)).toBeGreaterThan(
      findings.orientation.confidence_threshold,
    );
    expect(findings.orientation.decision).toBe("Input order supported");
    const mainGroup = adapter.groups.find((row) => row.name === findings.orientation.group);
    expect(mainGroup.arrays).toHaveLength(5);
    expect(mainGroup.arrays.map((row) => row.source_id)).toEqual(
      snapshot.example.records.map((row) => row.record_id),
    );
    expect(mainGroup.repeat_key).toMatch(/^[ACGT]+$/);
  });

  it("keeps both ancestral histories canonical and makes the supported history visibly more parsimonious", () => {
    const { orientation, adapter } = snapshot.job.summary;
    const group = snapshot.example.findings.orientation.group;
    const tree = orientation.trees.find((row) => row.group === group);
    const input = orientation.reconstructions.find(
      (row) => row.group === group && row.hypothesis === "input",
    );
    const reverse = orientation.reconstructions.find(
      (row) => row.group === group && row.hypothesis === "reverse",
    );
    const members = adapter.groups.find((row) => row.name === group).arrays;
    expect(orientation.trees_truncated).toBe(false);
    expect(orientation.reconstructions_truncated).toBe(false);
    expect(tree.selected_newick).toBe(tree.forward_newick);
    expect(input.newick).toBe(tree.forward_newick);
    expect(reverse.newick).toBe(tree.reverse_newick);
    expect(new Set(tree.forward_newick.match(/example_record_\d{2}(?=:)/g))).toEqual(
      new Set(members.map((row) => row.source_id)),
    );
    expect(input.nodes.map((node) => node.name)).toEqual(
      expect.arrayContaining([
        "Inner1",
        "Inner2",
        "Inner3",
        "Inner4",
        ...members.map((row) => row.source_id),
      ]),
    );
    expect(input.nodes).toHaveLength(9);
    expect(reverse.nodes).toHaveLength(9);
    expect(new Set(input.spacer_order)).toEqual(
      new Set(Array.from({ length: 42 }, (_, index) => index + 1)),
    );
    expect(new Set(reverse.spacer_order)).toEqual(new Set(input.spacer_order));
    for (const member of members) {
      const inputLeaf = input.nodes.find((node) => node.name === member.source_id);
      const reverseLeaf = reverse.nodes.find((node) => node.name === member.source_id);
      expect(new Set(reverseLeaf.spacers)).toEqual(new Set(inputLeaf.spacers));
      expect(inputLeaf.spacers).toHaveLength(member.spacer_count);
    }
    const gainTotal = (row) => row.nodes.reduce((sum, node) => sum + node.gains.length, 0);
    const lossTotal = (row) =>
      row.nodes.reduce(
        (sum, node) =>
          sum + node.loss_blocks.reduce((blockSum, block) => blockSum + block.length, 0),
        0,
      );
    expect([input.acquisition_count, input.deletion_count]).toEqual([42, 4]);
    expect([reverse.acquisition_count, reverse.deletion_count]).toEqual([42, 45]);
    expect([gainTotal(input), lossTotal(input)]).toEqual([42, 4]);
    expect([gainTotal(reverse), lossTotal(reverse)]).toEqual([42, 45]);
    expect(input.nodes.find((node) => node.name === "Inner4").gains).toHaveLength(1);
    expect(reverse.nodes.find((node) => node.name === "Inner4").gains).toHaveLength(42);
  });

  it("publishes only ordered masked record identifiers", () => {
    const inspection = inspectFasta(fasta);
    const expectedIds = [
      "example_record_01",
      "example_record_02",
      "example_record_03",
      "example_record_04",
      "example_record_05",
    ];
    expect(inspection.records.map((record) => record.identifier)).toEqual(expectedIds);
    expect(snapshot.example.records.map((record) => record.record_id)).toEqual(expectedIds);
    expect(snapshot.job.summary.detection.arrays.map((record) => record.source_id)).toEqual(
      expectedIds,
    );
    expect(
      snapshot.example.records.every((record) =>
        Object.keys(record).every((key) =>
          [
            "record_id",
            "sequence_length",
            "source_array_orientation",
            "expected_spacer_count",
          ].includes(key),
        ),
      ),
    ).toBe(true);
    expect(snapshot.example.records.reduce((sum, record) => sum + record.sequence_length, 0)).toBe(
      inspection.baseCount,
    );
    expect(JSON.stringify(snapshot)).not.toMatch(/(?:CP|FR|LN|LR|AP)\d{6}|ncbi\.nlm\.nih\.gov/i);
  });

  it("fails closed for byte changes, normalized changes, credentials, and inconsistent claims", async () => {
    await expect(validateExampleInput(snapshot, `${fasta}\n`)).rejects.toThrow(/does not match/i);
    const wrongHash = clone(snapshot);
    wrongHash.example.input.normalized_sha256 = "0".repeat(64);
    await expect(validateExampleInput(wrongHash, fasta)).rejects.toThrow(/does not match/i);
    const credentialLeak = clone(snapshot);
    credentialLeak.access_token = "not-public";
    expect(() => validateExampleSnapshot(credentialLeak)).toThrow(/incomplete or incompatible/i);
    const identityLeak = clone(snapshot);
    identityLeak.example.records[0].accession = "source-record";
    expect(() => validateExampleSnapshot(identityLeak)).toThrow(/incomplete or incompatible/i);
    const changedClaim = clone(snapshot);
    changedClaim.example.findings.reconstruction.insertions += 1;
    expect(() => validateExampleSnapshot(changedClaim)).toThrow(/incomplete or incompatible/i);
    const weakenedDecision = clone(snapshot);
    weakenedDecision.job.summary.orientation.comparisons[0].decisive = false;
    expect(() => validateExampleSnapshot(weakenedDecision)).toThrow(/incomplete or incompatible/i);
    const wrongSelectedTree = clone(snapshot);
    wrongSelectedTree.job.summary.orientation.trees[0].selected_newick =
      wrongSelectedTree.job.summary.orientation.trees[0].reverse_newick;
    expect(() => validateExampleSnapshot(wrongSelectedTree)).toThrow(/incomplete or incompatible/i);
    const changedEventTotal = clone(snapshot);
    changedEventTotal.job.summary.orientation.reconstructions[0].deletion_count += 1;
    expect(() => validateExampleSnapshot(changedEventTotal)).toThrow(/incomplete or incompatible/i);
    const invalidCanonicalSpacer = clone(snapshot);
    invalidCanonicalSpacer.job.summary.orientation.reconstructions[1].nodes[0].spacers[0] = 0;
    expect(() => validateExampleSnapshot(invalidCanonicalSpacer)).toThrow(
      /incomplete or incompatible/i,
    );
  });
});
