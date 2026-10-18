<script setup>
defineProps({
  available: Boolean,
  calibrationStatus: { type: String, default: null },
  calibrationCounts: { type: Object, required: true },
  fdr: { type: Number, default: null },
  reference: { type: Object, required: true },
  groups: { type: Array, required: true },
  counts: { type: Object, required: true },
  associations: { type: Array, required: true },
  truncated: Boolean,
  compact: Boolean,
});
const show = (v) => v ?? "Not reported";
const showCount = (v) => (v == null ? "Not reported" : v.toLocaleString("en-US"));
const scope = (v) =>
  ({
    broad_phage_reference: "Broad phage reference",
    targeted_panel: "Targeted phage panel",
  })[v] || "Scope not reported";
const status = (v) =>
  ({
    association: "Accepted association",
    uncalibrated_candidate: "Uncalibrated candidate",
    above_cutoff_candidate: "Candidate above FDR cutoff",
  })[v] || "Status unavailable";
</script>
<template>
  <section
    class="result-section association-results"
    aria-label="Spacer group association evidence"
  >
    <h3>Spacer group association evidence</h3>
    <p v-if="!available">A completed association analysis was not reported.</p>
    <template v-else>
      <p>
        {{ show(counts.associations) }} accepted associations ·
        {{ show(counts.candidates) }} candidates · {{ show(counts.support_hits) }} supporting
        alignments · {{ show(counts.unsupported) }} unsupported occurrences
      </p>
      <p>
        Similarity supports a group–target association hypothesis. It does not establish active
        infection, cleavage, a PAM, or a complete host identity.
      </p>
      <p v-if="calibrationStatus !== 'estimated'">
        <strong>FDR calibration unavailable.</strong> Candidates are evidence for inspection; none
        are accepted predictions. Unsupported sequences are not negative evidence.
      </p>
      <p v-else>
        Native empirical FDR is estimated. Acceptance requires FDR at or below {{ show(fdr) }}.
      </p>
      <p>
        Calibration pools scores across all groups in this run; group membership and support remain
        separate. Run-wide pooled score populations: {{ show(calibrationCounts.target) }} target and
        {{ show(calibrationCounts.control) }} control scores. A combined score is not a probability.
      </p>
      <p>
        Reference database: <strong>{{ show(reference.name) }}</strong> ·
        {{ showCount(reference.genomeCount) }}
        {{ reference.genomeCount === 1 ? "genome" : "genomes" }} ·
        {{ showCount(reference.totalBases) }} bases · {{ scope(reference.scope) }}.
      </p>
      <p>
        Searches cover only this database's recorded scope. Phages absent from it cannot be
        detected; no match does not rule out a biological interaction.
      </p>
      <p v-if="compact">
        Open Associations to inspect group membership, support and reference provenance.
      </p>
      <template v-else>
        <p v-if="truncated">
          This preview is limited to 100 associations. Download complete JSON or TSV from Files
          &amp; methods.
        </p>
        <p v-if="!associations.length">
          No reportable association candidates. This does not establish absence of biological
          interactions.
        </p>
        <article
          v-for="row in associations"
          :key="`${row.groupId}/${row.targetId}`"
          class="repeat-instance"
        >
          <h4>{{ row.groupId }} → {{ row.targetId }}</h4>
          <p>
            <strong>{{ status(row.status) }}</strong> · Combined score: {{ show(row.score) }} ·
            Native FDR: {{ show(row.fdr) }} · {{ show(row.hitCount) }} alignments
          </p>
          <details>
            <summary>Supporting spacer occurrences and alignments</summary>
            <p>
              Coordinates are 1-based inclusive. Descending coordinates indicate reverse-complement
              orientation. Alignments include gaps; no flanking sequence is interpreted as a PAM.
            </p>
            <p v-if="row.truncated">
              First 50 alignments shown; the complete JSON retains all support.
            </p>
            <div
              v-for="(hit, index) in row.hits"
              :key="index"
              class="support-hit"
            >
              <p class="sequence-text">
                Occurrences: {{ hit.occurrenceIds.join(", ") }} · Best-hit P:
                {{ show(hit.pBestHit) }}
              </p>
              <p>
                Query {{ show(hit.queryStart) }}–{{ show(hit.queryEnd) }} · Target
                {{ show(hit.targetStart) }}–{{ show(hit.targetEnd) }}
              </p>
              <pre
                >{{ hit.queryAlignment }}
{{ hit.targetAlignment }}</pre>
            </div>
          </details>
        </article>
        <details>
          <summary>Explicit group membership</summary>
          <p>
            Source annotations supplied with an export are not independently verified host
            identities.
          </p>
          <p
            v-for="group in groups"
            :key="group.id"
            class="sequence-text"
          >
            <strong>{{ group.id }}</strong
            >: {{ group.members.join(", ") }}. Unsupported:
            {{ group.unsupported.join(", ") || "none" }}.
          </p>
        </details>
        <details>
          <summary>Reference provenance</summary>
          <p>SpacePHARER {{ show(reference.version) }} · Control: {{ show(reference.control) }}</p>
          <p class="sequence-text">Reference SHA-256: {{ show(reference.digest) }}</p>
          <p class="sequence-text">Collection: {{ show(reference.source.collection) }}</p>
          <p class="sequence-text">Release: {{ show(reference.source.release) }}</p>
          <p>Retrieved: {{ show(reference.source.retrievedAt) }}</p>
          <p class="sequence-text">
            Selection policy: {{ show(reference.source.selectionPolicy) }}
          </p>
          <p
            v-for="genome in reference.genomes"
            :key="genome.id"
            class="sequence-text"
          >
            {{ show(genome.id) }} · {{ showCount(genome.length) }} bases · SHA-256
            {{ show(genome.sha256) }}
          </p>
          <p>
            <a
              href="https://doi.org/10.1093/bioinformatics/btab222"
              target="_blank"
              rel="noopener noreferrer"
              >SpacePHARER methods paper</a
            >
          </p>
        </details>
      </template>
    </template>
  </section>
</template>
<style scoped>
.sequence-text {
  overflow-wrap: anywhere;
}
.support-hit {
  padding: 0.5rem 0;
}
.support-hit pre {
  max-width: 100%;
  overflow-x: auto;
  font-size: 0.85rem;
}
</style>
