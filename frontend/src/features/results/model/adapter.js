import { normalizeDetectionArray } from "./detection.js";
import { asList, firstDefined, firstRecord } from "./values.js";

function normalizeGroup(value, index) {
  const group = firstRecord(value) || {};
  return {
    ...group,
    array_count: firstDefined(group, "array_count", "count"),
    arrays: asList(group.arrays).map(normalizeDetectionArray),
    name: String(firstDefined(group, "name", "group", "group_id", "id") || `group_${index + 1}`),
    repeat_key: firstDefined(group, "repeat_key", "canonical_repeat"),
  };
}

export function normalizeAdapter(summary) {
  const source = firstRecord(summary?.adapter, summary?.preflight, summary?.group_preflight);
  if (!source) return {};
  const emittedGroupCount = firstDefined(
    source,
    "emitted_group_count",
    "eligible_groups",
    "groups_count",
  );
  return {
    ...source,
    emitted_array_count: firstDefined(source, "emitted_array_count", "retained_arrays", "retained"),
    emitted_group_count:
      emittedGroupCount !== undefined
        ? emittedGroupCount
        : Array.isArray(source.groups)
          ? undefined
          : source.groups,
    groups: asList(source.groups).map(normalizeGroup),
    skipped_array_count: firstDefined(source, "skipped_array_count", "excluded_arrays", "excluded"),
    skipped_by_reason: firstDefined(
      source,
      "skipped_by_reason",
      "skip_reasons",
      "exclusion_reasons",
    ),
    unknown_strand_excluded_count:
      firstDefined(
        source,
        "unknown_strand_excluded_count",
        "unknown_strand",
        "unknown_strand_arrays",
      ) ?? 0,
  };
}
