import { asList } from "./values.js";

// Keep scientific limitations prominent. Routine privacy/empty-log publication
// messages remain available in Files & methods without obscuring the results.
const routinePublication = (value) =>
  typeof value === "string" &&
  (/^Adapter manifest published with \d+ private field\(s\) removed\.$/.test(value) ||
    /^(?:Bundle file omitted because it is empty:|Worker log omitted because it was empty:|Worker log omitted because it contained a private path:)/.test(
      value,
    ) ||
    /^Public artifact omitted because it was invalid or contained a private path: logs\//.test(
      value,
    ));

export function normalizeWarnings(summary) {
  const warnings = asList(summary.warnings);
  return {
    warnings: warnings.filter((value) => !routinePublication(value)),
    publication_notes: [
      ...new Set([
        ...asList(summary.publication_notes).filter(routinePublication),
        ...warnings.filter(routinePublication),
      ]),
    ],
  };
}
