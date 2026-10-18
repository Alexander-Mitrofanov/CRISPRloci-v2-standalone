import { normalizeArrayCompare } from "./arrayCompare.js";
import { normalizeLocusCompare } from "./locusCompare.js";
import { normalizeProphageContext } from "./prophageContext.js";
import { normalizeInsertionContext } from "./insertionContext.js";
import { normalizeDefenseContext } from "./defenseContext.js";
import { normalizeRnaFamily } from "./rnaFamily.js";
import { normalizeAssociation } from "./association.js";
import { normalizeAdapter } from "./adapter.js";
import { normalizeLeader } from "./leader.js";
import { normalizeRepeats } from "./repeats.js";
import { normalizeRepeatType } from "./repeatType.js";
import { normalizeRepeatMap } from "./repeatMap.js";
import { normalizeViral } from "./viral.js";
import { normalizeProtospacer } from "./protospacer.js";
import { normalizeDetection } from "./detection.js";
import { normalizeOrientation } from "./orientation.js";
import { normalizeProvenance } from "./provenance.js";
import { normalizeReconstruction } from "./reconstruction.js";
import { normalizeWarnings } from "./warnings.js";
import { firstRecord } from "./values.js";

export function normalizePublicResult(value, context = {}) {
  const summary = firstRecord(value) || {};
  const provenance = normalizeProvenance(summary, context);
  return {
    ...summary,
    ...provenance,
    ...normalizeWarnings({ ...summary, ...provenance }),
    adapter: normalizeAdapter(summary),
    detection: normalizeDetection(summary),
    crisprleader: normalizeLeader(summary.crisprleader),
    repeats: normalizeRepeats(summary.repeats),
    array_compare: normalizeArrayCompare(summary.array_compare),
    locus_compare: normalizeLocusCompare(summary.locus_compare),
    prophage_context: normalizeProphageContext(summary.prophage_context),
    insertion_context: normalizeInsertionContext(summary.insertion_context),
    defense_context: normalizeDefenseContext(summary.defense_context),
    rna_family: normalizeRnaFamily(summary.rna_family),
    spacer_association: normalizeAssociation(summary.spacer_association),
    repeat_type: normalizeRepeatType(summary.repeat_type),
    repeat_map: normalizeRepeatMap(summary.repeat_map),
    viral_search: normalizeViral(summary.viral_search),
    protospacer: normalizeProtospacer(summary.protospacer),
    orientation: normalizeOrientation(summary),
    reconstruction: normalizeReconstruction(summary),
  };
}
