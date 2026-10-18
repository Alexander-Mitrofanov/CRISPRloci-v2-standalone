import { MAX_GENOME_FEATURES, requireValue, textId } from "./contract.js";

export function featureCollector(sources) {
  const bySource = new Map(sources.map((source) => [source.id, source]));
  const tracks = [];
  const features = [];
  let unplaced = 0;

  function track(id, label, color) {
    const value = { id, label, color, features: [] };
    tracks.push(value);
    return value;
  }

  function add(
    target,
    row,
    {
      start = row.start,
      end = row.end,
      halfOpen = false,
      label = row.id,
      strand = row.strand,
      details = [],
      parentId = null,
      point = false,
    } = {},
  ) {
    const source = bySource.get(row.source_id);
    requireValue(source, "A feature refers to an unknown input record.");
    if (start == null && end == null) {
      unplaced++;
      return;
    }
    requireValue(Number.isSafeInteger(start) && Number.isSafeInteger(end));
    requireValue(halfOpen || end >= start);
    const from = halfOpen ? start : start - 1;
    requireValue(
      from >= 0 && end >= from && end <= source.length,
      "A feature interval falls outside its source sequence.",
    );
    if (from === end && !point) {
      unplaced++;
      return;
    }
    requireValue(
      features.length < MAX_GENOME_FEATURES,
      "This result exceeds the genome viewer's 50,000 feature limit. Complete downloads remain available.",
    );
    const direction =
      strand === 1 || strand === "+" ? "+" : strand === -1 || strand === "-" ? "-" : null;
    const feature = {
      key: `feature_${features.length}`,
      sourceId: source.id,
      chr: source.chr,
      start: from,
      end,
      point,
      strand: direction,
      label: textId(label),
      kind: target.label,
      trackId: target.id,
      parentId,
      details: details
        .filter((entry) => entry[1] !== undefined)
        .map(([key, value]) => [key, value === null ? "Not reported" : String(value)]),
    };
    target.features.push(feature);
    features.push(feature);
    return feature;
  }

  return {
    tracks,
    features,
    track,
    add,
    get unplaced() {
      return unplaced;
    },
  };
}
