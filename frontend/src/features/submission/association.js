/** Preview only; the API independently validates origins and grouping. */
export function validAssociationGrouping(text, grouping) {
  if (["single_set", "per_record"].includes(grouping)) return true;
  if (grouping !== "per_source") return false;
  try {
    const packet = JSON.parse(text.replace(/^\uFEFF/, ""));
    const rows = packet.records;
    return (
      Array.isArray(rows) &&
      rows.length > 0 &&
      rows.every(
        (r) =>
          typeof r.origin?.source_id === "string" && /^[a-f0-9]{64}$/.test(r.origin?.source_sha256),
      ) &&
      new Set(rows.map((r) => JSON.stringify([r.origin.source_id, r.origin.source_sha256]))).size <=
        100
    );
  } catch {
    return false;
  }
}
