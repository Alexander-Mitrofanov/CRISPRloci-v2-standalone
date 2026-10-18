import { asList, asRecord } from "./values.js";
const text = (v) => (typeof v === "string" ? v : null);
const number = (v) => (typeof v === "number" && Number.isFinite(v) && v >= 0 ? v : null);
const count = (v) => (Number.isSafeInteger(v) && v >= 0 ? v : null);
const sourceText = (v) => (typeof v === "string" && v.trim() ? v : null);
export function normalizeAssociation(value) {
  const doc = asRecord(value) || {};
  const ref = asRecord(doc.reference) || {};
  const source = asRecord(ref.source) || {};
  const legacyGenomes = Array.isArray(ref.genomes) && ref.genomes.length ? ref.genomes : null;
  const legacyBases =
    legacyGenomes && legacyGenomes.every((g) => asRecord(g) && count(g.length) !== null)
      ? count(legacyGenomes.reduce((total, g) => total + g.length, 0))
      : null;
  return {
    available: doc.available === true && doc.execution_completed === true,
    calibrationStatus: ["estimated", "unavailable", "not_run"].includes(doc.calibration_status)
      ? doc.calibration_status
      : null,
    calibrationCounts: {
      target: count(doc.calibration?.target_score_count),
      control: count(doc.calibration?.control_score_count),
    },
    fdr: number(doc.requested_fdr),
    reference: {
      name: text(ref.name),
      digest: text(ref.manifest_sha256),
      version: text(ref.engine_version),
      control: text(ref.control_method),
      genomeCount: count(ref.genome_count) ?? legacyGenomes?.length ?? null,
      totalBases: count(ref.total_bases) ?? legacyBases,
      scope: ["broad_phage_reference", "targeted_panel"].includes(ref.scope) ? ref.scope : null,
      source: {
        collection: sourceText(source.collection),
        release: sourceText(source.release),
        retrievedAt: sourceText(source.retrieved_at),
        selectionPolicy: sourceText(source.selection_policy),
      },
      genomes: asList(ref.genomes)
        .filter(asRecord)
        .slice(0, 100)
        .map((g) => ({ id: text(g.id), length: count(g.length), sha256: text(g.sha256) })),
    },
    groups: asList(doc.groups)
      .filter(asRecord)
      .slice(0, 100)
      .map((g) => ({
        id: text(g.id),
        members: asList(g.members)
          .filter((v) => typeof v === "string")
          .slice(0, 1000),
        unsupported: asList(g.unsupported_occurrence_ids)
          .filter((v) => typeof v === "string")
          .slice(0, 1000),
      })),
    counts: Object.fromEntries(
      ["records", "groups", "associations", "candidates", "support_hits", "unsupported"].map(
        (key) => [key, count(doc.counts?.[key])],
      ),
    ),
    truncated: doc.associations_truncated === true,
    associations: asList(doc.associations)
      .filter(asRecord)
      .slice(0, 100)
      .map((row) => ({
        groupId: text(row.group_id),
        targetId: text(row.target_id),
        score: number(row.combined_score),
        fdr:
          doc.calibration_status === "estimated" && number(row.native_fdr) <= 1
            ? number(row.native_fdr)
            : null,
        status: text(row.status),
        hitCount: count(row.native_hit_count),
        truncated: row.hits_truncated === true,
        hits: asList(row.hits)
          .filter(asRecord)
          .slice(0, 50)
          .map((hit) => ({
            occurrenceIds: asList(hit.occurrence_ids)
              .filter((v) => typeof v === "string")
              .slice(0, 1000),
            pBestHit: number(hit.p_best_hit),
            queryStart: count(hit.query_start),
            queryEnd: count(hit.query_end),
            targetStart: count(hit.target_start),
            targetEnd: count(hit.target_end),
            queryAlignment: text(hit.query_alignment),
            targetAlignment: text(hit.target_alignment),
          })),
      })),
  };
}
