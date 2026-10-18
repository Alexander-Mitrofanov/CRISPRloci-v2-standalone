export const catalogSummaryFixture = {
  available: true,
  snapshot_id: "test-snapshot",
  title: "CRISPR-Cas database campaign",
  source_url: "https://example.org/campaign",
  generated_at: "2026-09-22T00:00:00Z",
  counts: { assemblies: 2, arrays: 1, repeats: 2, spacers: 2, cas_genes: 1, cas_systems: 1 },
  coverage: {
    catalog_assemblies: 100,
    successful_assemblies: 98,
    failed_assemblies: 1,
    unavailable_assemblies: 1,
    imported_assemblies: 2,
  },
  limitations: ["Strand inference was disabled.", "Only accepted arrays were retained."],
  provenance: { caller: "CRISPRidentify 2.3.2" },
};

export const catalogArrayFixture = {
  id: 101,
  accession: "GCA_000001.1",
  sequence_id: "NC_000001.1",
  call_id: "array-1",
  start: 100,
  end: 140,
  strand: null,
  repeat_consensus: "ACGTACGT",
  category: "Bona-fide",
  tool: "CRISPRidentify",
  tool_version: "2.3.2",
  certainty_score: null,
  evidence_level: null,
};

export const catalogRows = {
  assemblies: [
    {
      id: 1,
      accession: "GCA_000001.1",
      organism_name: "Example organism",
      status: "complete",
      array_count: 1,
      cas_gene_count: 1,
    },
  ],
  arrays: [catalogArrayFixture],
  repeats: [
    {
      id: 201,
      array_id: 101,
      accession: "GCA_000001.1",
      ordinal: 1,
      start: 100,
      end: 107,
      sequence: "ACGTACGT",
    },
    {
      id: 202,
      array_id: 101,
      accession: "GCA_000001.1",
      ordinal: 2,
      start: 133,
      end: 140,
      sequence: "ACGTACGT",
    },
  ],
  spacers: [
    {
      id: 301,
      array_id: 101,
      accession: "GCA_000001.1",
      ordinal: 1,
      kind: "observed",
      start: 108,
      end: 132,
      sequence: "GATTACAGATTACAGATTACAGATTA",
    },
    {
      id: 302,
      array_id: 101,
      accession: "GCA_000001.1",
      ordinal: 2,
      kind: "deletion",
      start: null,
      end: null,
      boundary_position: 140,
      sequence: null,
    },
  ],
  cas_genes: [
    {
      id: 401,
      accession: "GCA_000001.1",
      call_id: "gene-1",
      annotation: "cas1",
      start: 300,
      end: 450,
      strand: "+",
      tool: "CasAndra",
    },
  ],
  cas_systems: [
    {
      id: 501,
      accession: "GCA_000001.1",
      call_id: "system-1",
      start: 300,
      end: 600,
      cas_type: "I",
      subtype: "I-E",
      confidence: null,
      tool: "CasAndra",
    },
  ],
};

export const catalogPageFixture = (entity) => ({
  items: catalogRows[entity],
  next_cursor: null,
  limit: 25,
});
