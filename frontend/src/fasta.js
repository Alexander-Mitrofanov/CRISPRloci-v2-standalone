export const IUPAC_DNA = new Set("ACGTRYSWKMBDHVN".split(""));
const IUPAC_DNA_INPUT = new Set("ACGTRYSWKMBDHVNacgtryswkmbdhvn".split(""));
const IUPAC_RNA_INPUT = new Set("ACGURYSWKMBDHVNacguryswkmbdhvn".split(""));
const LINE_BOUNDARIES = new Set([
  "\n",
  "\r",
  "\v",
  "\f",
  "\u001C",
  "\u001D",
  "\u001E",
  "\u0085",
  "\u2028",
  "\u2029",
]);
const PYTHON_WHITESPACE_CLASS = [
  "\\u0009-\\u000D",
  "\\u001C-\\u0020",
  "\\u0085",
  "\\u00A0",
  "\\u1680",
  "\\u2000-\\u200A",
  "\\u2028",
  "\\u2029",
  "\\u202F",
  "\\u205F",
  "\\u3000",
].join("");
const PYTHON_WHITESPACE = new RegExp(`[${PYTHON_WHITESPACE_CLASS}]+`, "g");
const PYTHON_EDGE_WHITESPACE = new RegExp(
  `^[${PYTHON_WHITESPACE_CLASS}]+|[${PYTHON_WHITESPACE_CLASS}]+$`,
  "g",
);

function stripPythonWhitespace(value) {
  return value.replace(PYTHON_EDGE_WHITESPACE, "");
}

export function* fastaLines(value) {
  let start = 0;
  let offset = 0;
  while (offset < value.length) {
    const character = value[offset];
    if (!LINE_BOUNDARIES.has(character)) {
      offset += 1;
      continue;
    }
    yield value.slice(start, offset);
    if (character === "\r" && value[offset + 1] === "\n") offset += 1;
    offset += 1;
    start = offset;
  }
  if (start < value.length) yield value.slice(start);
}

function safeIdentifier(value, fallback, maxLength = 64) {
  if (value.length <= maxLength && /^[A-Za-z0-9][A-Za-z0-9_-]*$/.test(value)) return value;
  const ascii = String(value)
    .normalize("NFKD")
    .replace(/[^\x00-\x7F]/g, "");
  return (
    ascii
      .replace(/[^A-Za-z0-9_-]+/g, "_")
      .replace(/^[_-]+|[_-]+$/g, "")
      .slice(0, maxLength) || fallback
  );
}

export function inspectFasta(text, { maxHeaderCharacters = 200, molecule = "DNA" } = {}) {
  let source = stripPythonWhitespace(String(text || "").replace(/^\uFEFF+/, ""));
  const errors = [];
  const records = [];
  const identifiers = new Set();
  const safeIdentifiers = new Set();
  let current = null;
  const alphabet = molecule === "RNA" ? IUPAC_RNA_INPUT : IUPAC_DNA_INPUT;

  if (!source) {
    return {
      valid: false,
      records,
      recordCount: 0,
      baseCount: 0,
      errors: ["Add at least one FASTA record."],
    };
  }

  if (!source.startsWith(">")) source = `>web_input\n${source}`;

  let lineNumber = 0;
  for (const rawLine of fastaLines(source)) {
    lineNumber += 1;
    const line = stripPythonWhitespace(rawLine);
    if (!line) continue;
    if (line.startsWith(">")) {
      const header = line.slice(1).trim();
      const identifier = header.split(PYTHON_WHITESPACE)[0];
      if (!identifier) errors.push(`Line ${lineNumber}: FASTA header is empty.`);
      if (header.length > maxHeaderCharacters)
        errors.push(`Line ${lineNumber}: FASTA header exceeds ${maxHeaderCharacters} characters.`);
      if (
        [...header].some(
          (character) => character.charCodeAt(0) < 32 || character.charCodeAt(0) === 127,
        )
      ) {
        errors.push(`Line ${lineNumber}: FASTA header contains control characters.`);
      }
      if (identifier && identifiers.has(identifier)) {
        errors.push(`Record identifier “${identifier}” is duplicated.`);
      }
      if (identifier) identifiers.add(identifier);
      const normalizedIdentifier = safeIdentifier(identifier, `sequence_${records.length + 1}`);
      if (safeIdentifiers.has(normalizedIdentifier)) {
        errors.push(
          `Record identifiers collide after safe filename normalization: “${normalizedIdentifier}”.`,
        );
      }
      safeIdentifiers.add(normalizedIdentifier);
      current = { header, identifier, normalizedIdentifier, sequence: "" };
      records.push(current);
      continue;
    }
    if (!current) {
      errors.push(`Line ${lineNumber}: sequence appears before the first FASTA header.`);
      continue;
    }
    const rawSequence = line.replace(PYTHON_WHITESPACE, "");
    const invalid = [...new Set(rawSequence)].filter((symbol) => !alphabet.has(symbol));
    if (invalid.length) {
      errors.push(
        `Line ${lineNumber}: unsupported ${molecule} symbol${invalid.length > 1 ? "s" : ""} ${invalid.join(", ")}.`,
      );
    }
    // Normalize only literal ASCII lowercase after validation. In particular,
    // never let Unicode case conversion turn a confusable into an IUPAC base.
    const sequence = rawSequence.replace(/[a-z]/g, (symbol) => symbol.toUpperCase());
    current.sequence += sequence;
  }

  for (const record of records) {
    if (!record.sequence)
      errors.push(`Record “${record.identifier || "unnamed"}” has no sequence.`);
  }

  const baseCount = records.reduce((total, record) => total + record.sequence.length, 0);
  return {
    valid: errors.length === 0 && records.length > 0,
    records,
    recordCount: records.length,
    baseCount,
    errors: [...new Set(errors)],
  };
}

export function readableBases(value) {
  const number = Number(value);
  if (!Number.isFinite(number)) return "—";
  if (number < 1_000) return `${number.toLocaleString()} bp`;
  if (number < 1_000_000) return `${(number / 1_000).toFixed(number < 10_000 ? 1 : 0)} kbp`;
  return `${(number / 1_000_000).toFixed(2)} Mbp`;
}
