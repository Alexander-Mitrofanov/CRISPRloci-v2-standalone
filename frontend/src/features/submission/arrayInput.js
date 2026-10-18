import { inspectFasta } from "../../fasta.js";

/** Browser preview only; strict field/coordinate validation belongs to the API. */
export function inspectArrayInput(text, limits = {}) {
  const invalid = (message) => ({
    valid: false,
    records: [],
    recordCount: 0,
    baseCount: 0,
    errors: [message],
  });
  if (new TextEncoder().encode(text).byteLength > 20_000_000)
    return invalid("Array JSON exceeds the input byte limit.");
  let packet;
  try {
    packet = JSON.parse(text.replace(/^\uFEFF/, ""));
  } catch {
    return invalid("Upload arrays.json from a detection analysis.");
  }
  if (
    packet?.schema?.name !== "crisprloci.arrays" ||
    packet.schema.version !== "1.0.0" ||
    packet.coordinate_system !== "1-based-inclusive-source-forward" ||
    !Array.isArray(packet.arrays) ||
    !packet.arrays.length ||
    packet.arrays.length > 200
  )
    return invalid("Supply a CRISPRloci arrays v1 file with observed arrays.");
  const rows = [];
  const scopes = new Set();
  let slots = 0;
  for (const [index, array] of packet.arrays.entries()) {
    if (!array || typeof array !== "object" || Array.isArray(array))
      return invalid("Each array must contain observed array evidence.");
    const scope = JSON.stringify([array.source_id, array.id]);
    if (
      scopes.has(scope) ||
      ![null, "+", "-"].includes(array.strand) ||
      !Array.isArray(array.spacers) ||
      !array.spacers.length ||
      array.spacers.length > 500
    )
      return invalid("Array identities, orientation or spacer slots are invalid.");
    scopes.add(scope);
    slots += array.spacers.length;
    for (const [n, spacer] of array.spacers.entries()) {
      if (!spacer || typeof spacer !== "object" || Array.isArray(spacer))
        return invalid("Each spacer slot must contain observed DNA or an explicit deletion.");
      if (spacer.ordinal !== n + 1)
        return invalid("Spacer slots must preserve their supplied order.");
      if (spacer.status === "deletion" && spacer.sequence === null) continue;
      if (
        spacer.status != null ||
        typeof spacer.sequence !== "string" ||
        !/^[ACGT]{1,1000}$/.test(spacer.sequence) ||
        spacer.end - spacer.start + 1 !== spacer.sequence.length
      )
        return invalid(
          "Every observed spacer needs unambiguous A/C/G/T and its source interval; missing sequences cannot be skipped.",
        );
      rows.push(`>occ${index}_${n}\n${spacer.sequence}\n`);
    }
  }
  if (slots > 10000 || !rows.length)
    return invalid("At least one observed spacer is required; at most 10,000 slots are supported.");
  return {
    ...inspectFasta(rows.join(""), { ...limits, molecule: "DNA" }),
    arrayCount: packet.arrays.length,
  };
}
