import { inspectFasta } from "../../fasta.js";
import { inspectSequenceInput } from "../../sequenceInput.js";
import { fail } from "./contract.js";
import { validateExampleSnapshot } from "./snapshot.js";

function normalizedFasta(inspection) {
  const lines = [];
  for (const record of inspection.records) {
    lines.push(`>${record.normalizedIdentifier}`);
    for (let offset = 0; offset < record.sequence.length; offset += 80) {
      lines.push(record.sequence.slice(offset, offset + 80));
    }
  }
  return `${lines.join("\n")}\n`;
}

async function sha256Hex(value) {
  if (!globalThis.crypto?.subtle) {
    fail("This browser cannot verify the stored example input.");
  }
  const digest = await globalThis.crypto.subtle.digest("SHA-256", new TextEncoder().encode(value));
  return [...new Uint8Array(digest)].map((byte) => byte.toString(16).padStart(2, "0")).join("");
}

export async function validateExampleInput(
  snapshotValue,
  sequence,
  { maxHeaderCharacters = 200 } = {},
) {
  const snapshot = validateExampleSnapshot(snapshotValue);
  if (snapshot.schema.name === "crisprloci.web-example") {
    const inspection = inspectSequenceInput(sequence, {
      mode: snapshot.job.mode,
      molecule: snapshot.job.options.molecule || "DNA",
      maxHeaderCharacters,
    });
    if (
      !inspection.valid ||
      inspection.recordCount !== snapshot.example.input.record_count ||
      inspection.baseCount !== snapshot.example.input.base_count ||
      (await sha256Hex(sequence)) !== snapshot.example.input.file_sha256 ||
      (snapshot.job.mode !== "array_compare" &&
        (await sha256Hex(normalizedFasta(inspection))) !== snapshot.example.input.normalized_sha256)
    )
      fail("The stored example input does not match its precomputed result.");
    return { snapshot, inspection };
  }
  const inspection = inspectFasta(sequence, { maxHeaderCharacters });
  const expectedIds = snapshot.example.records.map((record) => record.record_id);
  const observedIds = inspection.records.map((record) => record.identifier);
  if (
    !inspection.valid ||
    inspection.recordCount !== snapshot.example.input.record_count ||
    inspection.baseCount !== snapshot.example.input.base_count ||
    observedIds.length !== expectedIds.length ||
    observedIds.some((identifier, index) => identifier !== expectedIds[index])
  )
    fail("The stored example input does not match its precomputed result.");
  const [fileHash, normalizedHash] = await Promise.all([
    sha256Hex(sequence),
    sha256Hex(normalizedFasta(inspection)),
  ]);
  if (
    fileHash !== snapshot.example.input.file_sha256 ||
    normalizedHash !== snapshot.example.input.normalized_sha256
  )
    fail("The stored example input does not match its precomputed result.");
  return { snapshot, inspection };
}
