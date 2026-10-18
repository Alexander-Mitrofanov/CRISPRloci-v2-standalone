// Only modes with committed, source-verified example assets are offered.
// New live workflows must not create links to nonexistent example files.
const exampleModes = [
  "array_compare",
  "spacer_association",
  "repeat_type",
  "repeat_map",
  "viral_search",
  "loci",
  "leader",
  "cas",
  "tracrrna",
  "detection",
  "repeat_context",
  "repeats",
  "protospacer",
  "reconstruction",
  "orientation",
];

const asset = (key, json = false) =>
  Object.freeze({
    key,
    input: `examples/${key}/input.${json ? "json" : "fasta"}`,
    result: `examples/${key}/result.json`,
  });

export const EXAMPLE_ASSETS = Object.freeze(
  Object.fromEntries([
    ...exampleModes.map((mode) => [
      mode,
      asset(mode, ["array_compare", "spacer_association"].includes(mode)),
    ]),
    ["repeats-rna", asset("repeats-rna")],
    ["tracrrna-v", asset("tracrrna-v")],
  ]),
);

export function exampleAssetsFor(mode, options = {}) {
  const key =
    mode === "repeats" && options.molecule === "RNA"
      ? "repeats-rna"
      : mode === "tracrrna" && (options.tracrModelType || options.tracr_model_type) === "V"
        ? "tracrrna-v"
        : mode;
  return EXAMPLE_ASSETS[key] || null;
}

export function canonicalJson(value) {
  if (Array.isArray(value)) return `[${value.map(canonicalJson).join(",")}]`;
  if (value && typeof value === "object")
    return `{${Object.keys(value)
      .sort()
      .map((key) => `${JSON.stringify(key)}:${canonicalJson(value[key])}`)
      .join(",")}}`;
  return JSON.stringify(value);
}

export function exampleSubmissionMatches(snapshot, submission) {
  if (snapshot?.schema?.name !== "crisprloci.web-example") return false;
  const { sequence: _sequence, filename: _filename, mode, ...options } = submission;
  return (
    mode === snapshot.job.mode &&
    canonicalJson(options) === canonicalJson(snapshot.example.submission_options)
  );
}
