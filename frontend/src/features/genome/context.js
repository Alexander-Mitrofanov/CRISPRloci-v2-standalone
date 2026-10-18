import { coordinates, requireValue, rows } from "./contract.js";

export function insertionTracks(data, c) {
  coordinates(data);
  const elements = c.track("insertions", "Insertion candidates", "#98553d");
  const orfs = c.track("transposases", "Transposase ORFs", "#3158a5");
  const terminals = c.track("terminals", "Terminal inverted repeats", "#a76716");
  for (const row of rows(data.elements)) {
    c.add(elements, row, {
      strand: null,
      details: [
        ["Family", row.family],
        ["Native type", row.native_type],
        ["Contig edge", row.contig_edge],
        ["Target-site duplication", row.target_site_duplication_status],
      ],
    });
    if (row.transposase_orf)
      c.add(
        orfs,
        { ...row.transposase_orf, source_id: row.source_id },
        {
          label: `${row.id} · transposase`,
          parentId: row.id,
        },
      );
    for (const side of ["left", "right"]) {
      const terminal = row.terminal_inverted_repeats?.[side];
      if (terminal)
        c.add(
          terminals,
          { ...terminal, source_id: row.source_id },
          {
            label: `${row.id} · ${side} terminal repeat`,
            strand: null,
            parentId: row.id,
          },
        );
    }
  }
}

export function defenseTracks(data, c) {
  requireValue(data.policy?.coordinates === "one_based_inclusive");
  const genes = c.track("genes", "Predicted genes", "#65758b");
  const systems = c.track("defense", "Defense systems", "#0f766e");
  const components = c.track("components", "Defense components", "#3158a5");
  const hits = c.track("hits", "Profile hits", "#a76716");
  const byGene = new Map(rows(data.genes).map((gene) => [gene.id, gene]));
  for (const row of byGene.values()) c.add(genes, row);
  for (const row of rows(data.systems))
    c.add(systems, row, {
      label: `${row.subtype || row.type} · ${row.id}`,
      strand: null,
      details: [
        ["Model", row.model],
        ["Wholeness", row.wholeness],
        ["Raw score", row.score],
      ],
    });
  for (const [values, track] of [
    [data.components, components],
    [data.hits, hits],
  ]) {
    for (const row of rows(values)) {
      const gene = byGene.get(row.gene_id);
      requireValue(gene, "Defense evidence refers to an unknown gene.");
      c.add(track, gene, {
        label: `${row.profile} · ${row.gene_id}`,
        parentId: row.system_id,
        details: [
          ["Gene", row.gene_id],
          ["System", row.system_id ?? null],
          ["Profile", row.profile],
          ["Raw score", row.score],
          ["E-value", row.i_evalue],
          ["Interpretation", "Profile evidence; functional activity not assessed"],
        ],
      });
    }
  }
}

export function prophageTracks(data, c) {
  coordinates(data);
  const regions = c.track("prophage", "Prophage candidates", "#85629f");
  const genes = c.track("genes", "Predicted genes", "#65758b");
  for (const row of rows(data.regions))
    c.add(regions, row, {
      strand: null,
      details: [
        ["Peak model score", row.peak_score],
        ["Boundary", "Computational prediction"],
      ],
    });
  for (const row of rows(data.genes)) c.add(genes, row);
  const arrays = c.track("arrays", "Reported arrays", "#0f766e");
  for (const row of rows(data.source_arrays))
    c.add(arrays, row, { strand: row.orientation, details: [["Category", row.category]] });
  const matches = c.track("matches", "External spacer matches", "#a76716");
  for (const [index, row] of rows(data.spacer_matches).entries()) {
    c.add(
      matches,
      { ...row, source_id: row.target_source_id },
      {
        label: `${row.query_id} · match ${index + 1}`,
        strand: row.strand,
        details: [
          ["Query occurrence", row.query_id],
          ["Relative strand", row.strand],
          ["Array", row.array_id],
          ["Search complete", data.search_policy?.search_complete],
          ["Interpretation", "Exact sequence match; PAM and targeting activity not assessed"],
        ],
      },
    );
  }
}

export function protospacerTracks(data, c) {
  coordinates(data, "1-based-inclusive-target-forward");
  const matches = c.track("matches", "Observed-spacer matches", "#a76716");
  for (const [index, row] of rows(data.matches).entries())
    c.add(
      matches,
      { source_id: row.target_id, start: row.target_start, end: row.target_end },
      {
        label: `${row.spacer_id} · match ${index + 1}`,
        strand: row.relative_strand,
        details: [
          ["Spacer", row.spacer_id],
          ["Relative strand", row.relative_strand],
          ["Search complete", data.search_complete],
          ["PAM", data.policy?.pam],
          ["Array overlap", data.policy?.array_overlap],
        ],
      },
    );
}
