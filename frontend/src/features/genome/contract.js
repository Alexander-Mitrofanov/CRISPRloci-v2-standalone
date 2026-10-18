export const GENOME_MODES = new Set([
  "detection",
  "loci",
  "cas",
  "tracrrna",
  "leader",
  "protospacer",
  "prophage_context",
  "insertion_context",
  "defense_context",
]);
export const MAX_GENOME_BASES = 10_000_000;
export const MAX_GENOME_FEATURES = 50_000;

export function genomeArtifacts(artifacts = []) {
  const find = (name) =>
    (Array.isArray(artifacts) ? artifacts : []).find(
      (row) => row?.name === name || row?.filename === name,
    );
  return { result: find("result.json"), input: find("input-sequences.json") };
}

export function genomeAvailable(mode, artifacts) {
  const { result, input } = genomeArtifacts(artifacts);
  return GENOME_MODES.has(mode) && Boolean(result && input);
}

export function requireValue(
  condition,
  message = "The genome evidence has an unsupported format.",
) {
  if (!condition) throw new Error(message);
}

export function rows(value) {
  requireValue(Array.isArray(value));
  return value;
}

export function textId(value) {
  requireValue(typeof value === "string" && value.length > 0 && value.length <= 1024);
  return value;
}

export function coordinates(data, expected = "1-based-inclusive-source-forward") {
  requireValue(data?.coordinate_system === expected, "Unsupported genomic coordinate system.");
}

async function sha256(value) {
  const bytes = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(value));
  return Array.from(new Uint8Array(bytes), (v) => v.toString(16).padStart(2, "0")).join("");
}

export async function sourceRecords(input, declared) {
  requireValue(
    input?.schema?.name === "crisprloci.sequence-set" &&
      input.schema.version === "1.0.0" &&
      input.complete === true &&
      input.kind === "target" &&
      input.molecule === "DNA",
    "Complete genomic DNA is required. Repeat, spacer and RNA inputs cannot serve as a genome.",
  );
  const records = rows(input.records);
  requireValue(records.length > 0 && records.length <= 1000);
  let total = 0;
  const ids = new Set();
  const sources = [];
  for (const record of records) {
    const id = textId(record.id);
    requireValue(!ids.has(id));
    ids.add(id);
    requireValue(
      typeof record.sequence === "string" && /^[ACGTRYSWKMBDHVN]+$/.test(record.sequence),
    );
    total += record.sequence.length;
    requireValue(
      total <= MAX_GENOME_BASES,
      "This input exceeds the genome viewer's 10 million base limit. Complete downloads remain available.",
    );
    sources.push({
      id,
      inputId: id,
      chr: `record_${sources.length + 1}`,
      sequence: record.sequence,
      length: record.sequence.length,
      sha256: await sha256(record.sequence),
      topology: "unknown",
      assessment: null,
    });
  }
  if (declared === undefined) return sources; // Standalone Cas uses original input record IDs.
  const assigned = new Set();
  const declaredIds = new Set();
  for (const source of rows(declared)) {
    textId(source.id);
    requireValue(!declaredIds.has(source.id));
    declaredIds.add(source.id);
    requireValue(Number.isSafeInteger(source.length) && /^[a-f0-9]{64}$/.test(source.sha256));
    const candidates = sources.filter(
      (record) =>
        record.length === source.length && record.sha256 === source.sha256 && !assigned.has(record),
    );
    const exact = sources.find((candidate) => candidate.inputId === source.id);
    const record = exact
      ? candidates.includes(exact)
        ? exact
        : null
      : candidates.length === 1
        ? candidates[0]
        : null;
    requireValue(
      record,
      "The result cannot be matched unambiguously to its input DNA. Genome browsing is unavailable.",
    );
    assigned.add(record);
    Object.assign(record, {
      id: source.id,
      topology: source.topology ?? "unknown",
      assessment: source.status ?? null,
    });
  }
  requireValue(new Set(sources.map((s) => s.id)).size === sources.length);
  return sources;
}
