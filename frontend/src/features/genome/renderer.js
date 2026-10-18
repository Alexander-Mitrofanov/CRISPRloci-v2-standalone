// This module receives only validated in-memory science, never API credentials
// or artifact URLs. Keep the large browser dependency behind this lazy import.
const escapeLabel = (value) =>
  String(value).replace(
    /[&<>"']/g,
    (character) =>
      ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[character],
  );

export async function createGenomeBrowser(element, sources, tracks, onSelect) {
  const { default: igv } = await import("igv/dist/igv.esm.js");
  const fasta = new File(
    sources.flatMap((source) => [`>${source.chr}\n`, source.sequence, "\n"]),
    "genome.fasta",
    { type: "text/plain" },
  );
  const browser = await igv.createBrowser(element, {
    loadDefaultGenomes: false,
    queryParametersSupported: false,
    search: false,
    showNavigation: false,
    showSVGButton: false,
    showIdeogram: false,
    showCircularView: false,
    showCircularViewButton: false,
    showChromosomeWidget: false,
    minimumBases: 1,
    reference: {
      id: "submitted-sequences",
      name: "Submitted DNA",
      fastaURL: fasta,
      wholeGenomeView: false,
    },
    locus: `${sources[0].chr}:1-${sources[0].length}`,
    tracks: tracks
      .filter((track) => track.features.length)
      .map((track) => ({
        id: track.id,
        type: "annotation",
        format: "bed",
        name: track.label,
        color: track.color,
        displayMode: "EXPANDED",
        height: 65,
        maxHeight: 180,
        removable: false,
        searchable: false,
        features: track.features.map((feature) => ({
          chr: feature.chr,
          start: feature.start,
          end: feature.end,
          strand: feature.strand ?? ".",
          name: escapeLabel(feature.label),
          featureKey: feature.key,
        })),
      })),
    listeners: {
      trackclick(_track, _data, _position, features) {
        const keys = (features || []).map((feature) => feature.featureKey).filter(Boolean);
        if (keys.length) onSelect(keys);
        return false; // Vue renders evidence as text; never use tool-generated popup HTML.
      },
    },
  });
  return {
    search: (locus) => browser.search(locus),
    zoomIn: () => browser.zoomIn(),
    zoomOut: () => browser.zoomOut(),
    visible: () => browser.visibilityChange(),
    loci: () => browser.currentLoci(),
    svg: () => browser.toSVG(),
    dispose: () => igv.removeBrowser(browser),
  };
}
