import { asRecord } from "./values.js";
import {
  available,
  casRelations,
  count,
  counts,
  interval,
  metric,
  probability,
  rows,
  strand,
  text,
} from "./contextValues.js";

export function normalizeLocusCompare(value) {
  const doc = asRecord(value) || {};
  const loci = rows(doc.loci, 8).map((locus) => ({
    id: text(locus.id),
    sourceId: text(locus.source_id),
    ...interval(locus),
    atSourceStart: locus.at_source_start === true,
    atSourceEnd: locus.at_source_end === true,
    genes: rows(locus.genes, 120).map((gene) => ({
      id: text(gene.id),
      ...interval(gene),
      strand: strand(gene.strand),
      proteinLength: count(gene.protein_length),
      translationTable: count(gene.translation_table),
      partial5: typeof gene.partial_5prime === "boolean" ? gene.partial_5prime : null,
      partial3: typeof gene.partial_3prime === "boolean" ? gene.partial_3prime : null,
    })),
  }));
  const geneIds = new Map(
    loci.map((locus) => [locus.id, new Set(locus.genes.map((gene) => gene.id))]),
  );
  const links = rows(doc.links)
    .filter(
      (link) =>
        geneIds.get(link.query_locus)?.has(link.query_gene) &&
        geneIds.get(link.target_locus)?.has(link.target_gene) &&
        probability(link.identity) !== null &&
        probability(link.similarity) !== null,
    )
    .map((link) => ({
      queryLocus: text(link.query_locus),
      queryGene: text(link.query_gene),
      targetLocus: text(link.target_locus),
      targetGene: text(link.target_gene),
      identity: probability(link.identity),
      similarity: probability(link.similarity),
      score: metric(link.score),
    }));
  return {
    available: available(doc),
    loci,
    links,
    casContext: casRelations(doc.cas_context),
    counts: counts(doc.counts, ["loci", "genes", "links", "evaluated_protein_pairs"]),
    identityCutoff: probability(doc.options?.identity),
    truncated: doc.links_truncated === true,
  };
}
