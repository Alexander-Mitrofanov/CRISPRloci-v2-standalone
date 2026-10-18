import { GENOME_MODES, requireValue, sourceRecords } from "./contract.js";
import { featureCollector } from "./features.js";
import { arrayTracks, casTracks, tracrTracks, flankTracks } from "./annotations.js";
import { insertionTracks, defenseTracks, prophageTracks, protospacerTracks } from "./context.js";

function checkTrackDensity(tracks) {
  // IGV's static feature packer permits 1,000 overlapping rows. Fail explicitly
  // before that limit would silently leave features without drawable rows.
  for (const track of tracks) {
    const events = track.features
      .flatMap((feature) =>
        feature.point
          ? []
          : [
              { source: feature.sourceId, position: feature.start, change: 1 },
              { source: feature.sourceId, position: feature.end, change: -1 },
            ],
      )
      .sort(
        (a, b) =>
          a.source.localeCompare(b.source) || a.position - b.position || a.change - b.change,
      );
    let overlapping = 0;
    for (const event of events) {
      overlapping += event.change;
      requireValue(
        overlapping <= 1000,
        "More than 1,000 features overlap in a track. Use the complete downloads to inspect this dense result.",
      );
    }
  }
}

export async function buildGenomeModel(result, input, mode) {
  requireValue(
    GENOME_MODES.has(mode) &&
      result?.module === mode &&
      result.schema?.name === "crisprloci.module-result" &&
      result.schema.version === "1.0.0",
  );
  requireValue(
    result.complete === true && result.data && result.status === "completed",
    "Complete genome evidence is unavailable for this job. Consult its detailed downloads.",
  );
  const data = result.data;
  const declared = mode === "protospacer" ? data.targets : data.sources;
  requireValue(mode === "cas" || Array.isArray(declared));
  const sources = await sourceRecords(input, declared);
  const c = featureCollector(sources);
  if (mode === "loci") {
    arrayTracks(data.detection, c);
    casTracks(data.cas, c);
    tracrTracks(data.tracrrna, c);
    flankTracks(data.flanks, c);
  } else if (mode === "detection") arrayTracks(data, c);
  else if (mode === "cas") casTracks(data, c);
  else if (mode === "tracrrna") tracrTracks(data, c);
  else if (mode === "leader") flankTracks(data, c, true);
  else if (mode === "insertion_context") insertionTracks(data, c);
  else if (mode === "defense_context") defenseTracks(data, c);
  else if (mode === "prophage_context") prophageTracks(data, c);
  else if (mode === "protospacer") protospacerTracks(data, c);
  checkTrackDensity(c.tracks);
  const notes = [
    "Coordinates are 1-based inclusive on the submitted forward sequence. Spatial proximity does not establish functional association.",
    "Tracks show reported evidence for this analysis only. Empty tracks do not establish biological absence.",
    "Raw scores and classifications are model evidence, not calibrated probabilities or proof of activity.",
  ];
  if (typeof data.interpretation === "string") notes.push(data.interpretation);
  const searchComplete =
    mode === "protospacer" ? data.search_complete : data.search_policy?.search_complete;
  if (["protospacer", "prophage_context"].includes(mode) && searchComplete !== true) {
    notes.push(
      "Spacer search completeness is not established. The displayed matches do not cover an exhaustive search.",
    );
  }
  if (mode === "prophage_context" && data.source_inventory_complete !== true) {
    notes.push(
      "The source array/spacer inventory is incomplete; consult the reported search scope in the evidence tables.",
    );
  }
  if (["loci", "leader"].includes(mode))
    notes.push(
      "Both adjacent context windows are shown. Leader prediction was not performed; extraction strand is not transcription direction.",
    );
  if (mode === "prophage_context")
    notes.push(
      "Short records may be unassessed for prophages. Sequence matches do not establish targeting activity; PAM was not assessed.",
    );
  if (c.unplaced)
    notes.push(
      `${c.unplaced} empty or unlocated features cannot be drawn. Their evidence remains in the result tables and downloads.`,
    );
  return { sources, tracks: c.tracks, features: c.features, notes };
}
