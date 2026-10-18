import { inspectArrayInput } from "./arrayInput.js";
import ports from "../../contracts/module-ports-v1.json" with { type: "json" };
import { inspectFasta } from "../../fasta.js";

// The server independently validates the entire contract. Preview only the
// compatible sequence port, keeping the original JSON in the submitted request.
export function inspectSequenceInput(text, { mode, molecule = "DNA", ...limits }) {
  if (mode === "array_compare") return inspectArrayInput(text, limits);
  if (!text.trimStart().startsWith("{")) return inspectFasta(text, { molecule, ...limits });
  const invalid = (message) => ({
    valid: false,
    records: [],
    recordCount: 0,
    baseCount: 0,
    errors: [message],
  });
  if (new TextEncoder().encode(text).byteLength > 20_000_000)
    return invalid("Sequence JSON exceeds the input byte limit.");
  let packet;
  try {
    packet = JSON.parse(text.replace(/^\uFEFF/, ""));
  } catch {
    return invalid("Invalid JSON. Upload a CRISPRloci sequence-set v1 file.");
  }
  if (packet?.schema?.name !== "crisprloci.sequence-set" || packet.schema.version !== "1.0.0")
    return invalid(
      "Upload a CRISPRloci sequence-set v1 file (repeats.json, spacers.json or input-sequences.json).",
    );
  if (packet.complete !== true)
    return invalid("Incomplete sequence sets cannot be used as module inputs.");
  if (packet.kind !== ports[mode])
    return invalid(
      `This analysis requires ${ports[mode]} sequences, not ${packet.kind} sequences.`,
    );
  if (packet.molecule !== molecule)
    return invalid(`Select ${packet.molecule} input to use this sequence set.`);
  if (!Array.isArray(packet.records) || !packet.records.length || packet.records.length > 100_000)
    return invalid("This sequence set has no sequences or exceeds the record limit.");
  const alphabet = molecule === "RNA" ? /^[ACGURYSWKMBDHVN]+$/ : /^[ACGTRYSWKMBDHVN]+$/;
  if (
    packet.records.some(
      (row) =>
        typeof row?.id !== "string" ||
        !/^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$/.test(row.id) ||
        typeof row.sequence !== "string" ||
        !alphabet.test(row.sequence),
    )
  )
    return invalid("Sequence JSON contains an invalid identifier or sequence.");
  return inspectFasta(packet.records.map((row) => `>${row.id}\n${row.sequence}\n`).join(""), {
    molecule,
    ...limits,
  });
}
