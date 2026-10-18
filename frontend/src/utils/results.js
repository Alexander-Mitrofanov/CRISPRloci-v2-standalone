import { asArray, finiteMetric } from "./formatting.js";

const MAX_ADAPTER_MANIFEST_BYTES = 1_000_000;

export { MAX_ADAPTER_MANIFEST_BYTES };

export function groupIdentity(item, index = 0) {
  return String(item?.group || item?.name || `group_${index + 1}`);
}

export function comparisonDecisionFor(item, orientation = {}) {
  const rawThreshold = Number(item?.confidence_threshold ?? orientation?.confidence_threshold ?? 5);
  const threshold = Number.isFinite(rawThreshold) && rawThreshold >= 0 ? rawThreshold : 5;
  const reported = finiteMetric(item?.forward_minus_reverse_ln_likelihood_bdm);
  const forward = finiteMetric(item?.forward_ln_likelihood_bdm);
  const reverse = finiteMetric(item?.reverse_ln_likelihood_bdm);
  const delta = reported ?? (forward != null && reverse != null ? forward - reverse : null);
  if (item?.decisive === false || delta == null || Math.abs(delta) <= threshold) {
    return { label: "Unresolved", threshold, delta };
  }
  return {
    label: delta > threshold ? "Input order supported" : "Reverse input order supported",
    threshold,
    delta,
  };
}

function safeText(value, maxLength = 240) {
  if (value == null) return null;
  const text = String(value).trim();
  return text ? text.slice(0, maxLength) : null;
}

function safeCount(value) {
  const number = Number(value);
  return Number.isInteger(number) && number >= 0 && number <= 100_000 ? number : null;
}

export function sanitizeAdapterMembership(document) {
  if (!document || typeof document !== "object" || Array.isArray(document)) return [];
  return asArray(document.groups)
    .slice(0, 100)
    .map((group) => {
      if (!group || typeof group !== "object" || Array.isArray(group)) return null;
      const name = safeText(group.name);
      if (!name) return null;
      const rawArrays = asArray(group.arrays);
      const arrays = rawArrays
        .slice(0, 200)
        .map((array) => {
          if (!array || typeof array !== "object" || Array.isArray(array)) return null;
          return {
            source_id: safeText(array.source_id) || "unknown",
            array_id: safeText(array.array_id) || "unknown",
            category: safeText(array.category, 80) || "unknown",
            spacer_count: safeCount(array.spacer_count),
            strand: safeText(array.strand, 32) || "unknown",
            input_sequence_orientation: safeText(array.input_sequence_orientation, 32) || "unknown",
            ccdb_strand: safeText(array.ccdb_strand, 32) || "unknown",
          };
        })
        .filter(Boolean);
      return {
        name,
        array_count: safeCount(group.array_count),
        repeat_key: safeText(group.repeat_key, 512),
        arrays,
        arrays_truncated: rawArrays.length > arrays.length,
      };
    })
    .filter(Boolean);
}

export function adapterHasMembership(summary) {
  const groups = asArray(summary?.adapter?.groups);
  return (
    groups.length > 0 &&
    groups.every((group) => Number(group?.array_count) === 0 || asArray(group?.arrays).length > 0)
  );
}

export function mergeAdapterMembership(summary, membershipGroups) {
  if (!membershipGroups?.length || adapterHasMembership(summary)) return summary;
  const adapter = summary?.adapter;
  if (!adapter || typeof adapter !== "object") return summary;
  const membershipByName = new Map(membershipGroups.map((group) => [String(group.name), group]));
  const currentGroups = asArray(adapter.groups);
  const groups = currentGroups.length
    ? currentGroups.map((group) => {
        const membership = membershipByName.get(String(group?.name));
        return membership
          ? {
              ...group,
              repeat_key: group?.repeat_key || membership.repeat_key,
              arrays: membership.arrays,
              arrays_truncated: membership.arrays_truncated,
            }
          : group;
      })
    : membershipGroups;
  return { ...summary, adapter: { ...adapter, groups } };
}
