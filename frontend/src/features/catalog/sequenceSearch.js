const dna = /^[ACGTRYSWKMBDHVN]+$/iu;
const complement = Object.fromEntries(
  [..."ACGTRYSWKMBDHVN"].map((base, index) => [base, "TGCAYRSWMKVHDBN"[index]]),
);

export function sequenceSearchIssue(sequence, status) {
  if (status === "deletion") return "No sequence (deletion).";
  if (status === "too_long" || sequence?.length > 1000)
    return "Database search supports sequences up to 1,000 bases.";
  if (!sequence) return "Sequence not reported.";
  if (!dna.test(sequence)) return "A DNA sequence is required for database search.";
  return "";
}

export function sequenceComparison(query, candidate) {
  if (
    sequenceSearchIssue(query) ||
    sequenceSearchIssue(candidate) ||
    candidate.length !== query.length
  )
    return null;
  const forward = query.toUpperCase();
  const reverse = [...forward]
    .reverse()
    .map((base) => complement[base])
    .join("");
  const sequence = candidate.toUpperCase();
  const substitutions = (oriented) =>
    [...oriented].reduce((count, base, index) => count + Number(base !== sequence[index]), 0);
  const same = substitutions(forward);
  const opposite = substitutions(reverse);
  const differences = Math.min(same, opposite);
  if (differences > 1) return null;
  return {
    substitutions: differences,
    identity: ((100 * (query.length - differences)) / query.length).toFixed(1),
    strand: same === opposite ? "Both" : same < opposite ? "Same" : "Reverse complement",
  };
}
