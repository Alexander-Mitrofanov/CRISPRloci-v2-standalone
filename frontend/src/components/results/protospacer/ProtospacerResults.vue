<script setup>
import { useId } from "vue";

defineProps({
  available: Boolean,
  searchComplete: { type: Boolean, default: null },
  outcome: { type: String, default: null },
  reference: { type: Object, required: true },
  counts: { type: Object, required: true },
  matches: { type: Array, required: true },
  truncated: Boolean,
  compact: Boolean,
});
const heading = useId();
const count = (value) => (value == null ? "Not reported" : value.toLocaleString());
const interval = (row) =>
  row.target_start != null && row.target_end != null
    ? `${count(row.target_start)}–${count(row.target_end)}`
    : "Not reported";
const strand = (value) =>
  ({ "+": "Forward (+)", "-": "Reverse (−)", both: "Both (palindromic spacer)" })[value] ||
  "Not reported";
</script>

<template>
  <section
    class="protospacer-results"
    :aria-labelledby="heading"
  >
    <h3 :id="heading">Reference spacer matches</h3>
    <p class="annotation-lead">
      Full-length exact matches to reference spacers in the submitted linear DNA records, searched
      in both orientations. A match does not establish functional targeting, infection, immunity or
      a host association.
    </p>
    <p
      v-if="!available"
      class="annotation-empty"
    >
      Protospacer evidence was not reported. Missing evidence is not a zero-match result.
    </p>
    <p
      v-else-if="searchComplete !== true"
      class="annotation-empty"
    >
      Search completion was not confirmed. These results cannot establish a completed search with no
      matches.
    </p>
    <div class="annotation-counts">
      <div>
        <strong>{{ count(counts.matches) }}</strong
        ><span>Match occurrences</span>
      </div>
      <div>
        <strong>{{ count(counts.matched_spacers) }}</strong
        ><span>Matched reference spacers</span>
      </div>
      <div>
        <strong>{{ count(counts.targets) }}</strong
        ><span>Target records</span>
      </div>
      <div>
        <strong>{{ count(counts.target_bases) }}</strong
        ><span>Target bases</span>
      </div>
    </div>
    <p class="annotation-lead">
      CRISPR-array overlap is unknown; source-array copies may be included. PAM compatibility is not
      evaluated. Strand is relative to the submitted forward sequence, not a predicted transcription
      direction. Multiple occurrences are not independent biological confirmation.
    </p>
    <p
      v-if="available && searchComplete === true && outcome === 'no_assessable_windows'"
      class="annotation-empty"
    >
      No target windows could be assessed under this search policy. Short or ambiguous input is
      unassessed; this is not evidence that protospacers are absent.
    </p>
    <p
      v-else-if="
        available && searchComplete === true && outcome === 'no_matches' && counts.matches === 0
      "
      class="annotation-empty"
    >
      No matches to eligible reference spacers were reported in the submitted linear sequences under
      this exact-match policy. This does not establish absence of other protospacers.
    </p>
    <dl class="protospacer-reference">
      <div>
        <dt>Reference collection</dt>
        <dd>{{ reference.name || "Not reported" }}</dd>
      </div>
      <div>
        <dt>Source release</dt>
        <dd>{{ reference.source_release || "Not reported" }}</dd>
      </div>
      <div>
        <dt>Eligible reference spacers</dt>
        <dd>{{ count(reference.eligible_record_count) }} / {{ count(reference.record_count) }}</dd>
      </div>
    </dl>
    <template v-if="!compact">
      <p class="annotation-lead">
        Reference spacers with ambiguous bases or outside the supported length range are excluded.
        Ambiguous target windows are skipped; incomplete sequence remains unassessed.
      </p>
      <dl class="protospacer-reference">
        <div>
          <dt>Spacer length range</dt>
          <dd>{{ count(reference.min_length) }}–{{ count(reference.max_length) }} nt</dd>
        </div>
        <div>
          <dt>Excluded ambiguous spacers</dt>
          <dd>{{ count(reference.skipped_ambiguous) }}</dd>
        </div>
        <div>
          <dt>Excluded spacer lengths</dt>
          <dd>{{ count(reference.skipped_length) }}</dd>
        </div>
        <div>
          <dt>Ambiguous target seed windows</dt>
          <dd>{{ count(counts.ambiguous_seed_windows) }}</dd>
        </div>
        <div>
          <dt>Source file SHA-256</dt>
          <dd>
            <code>{{ reference.source_file_sha256 || "Not reported" }}</code>
          </dd>
        </div>
        <div>
          <dt>Database SHA-256</dt>
          <dd>
            <code>{{ reference.database_sha256 || "Not reported" }}</code>
          </dd>
        </div>
        <div v-if="reference.source_url">
          <dt>Reference source</dt>
          <dd>
            <a
              :href="reference.source_url"
              target="_blank"
              rel="noopener noreferrer"
              >Research spacer collection</a
            >
          </dd>
        </div>
      </dl>
      <p
        v-if="truncated"
        class="annotation-notice"
      >
        This preview shows {{ matches.length }} of {{ count(counts.matches) }} reported match
        occurrences. Download the complete match table from Files &amp; methods.
      </p>
      <p
        v-if="!matches.length && counts.matches !== 0 && available"
        class="annotation-empty"
      >
        No match preview was reported. Consult Files &amp; methods for the retained evidence.
      </p>
      <template v-if="matches.length">
        <p class="annotation-lead">
          Target coordinates are 1-based and inclusive on the submitted forward sequence. Spacer
          sequence is shown in the reference orientation.
        </p>
        <div
          class="annotation-scroll"
          tabindex="0"
          role="region"
          aria-label="Protospacer match table"
        >
          <table class="annotation-table protospacer-table">
            <thead>
              <tr>
                <th scope="col">Target record</th>
                <th scope="col">Coordinates</th>
                <th scope="col">Relative strand</th>
                <th scope="col">Reference spacer</th>
                <th scope="col">Length (nt)</th>
                <th scope="col">Spacer sequence</th>
              </tr>
            </thead>
            <tbody>
              <tr
                v-for="(row, index) in matches"
                :key="index"
              >
                <td>{{ row.target_id || "Not reported" }}</td>
                <td>{{ interval(row) }}</td>
                <td>{{ strand(row.relative_strand) }}</td>
                <td>{{ row.spacer_id || "Not reported" }}</td>
                <td>{{ count(row.query_length) }}</td>
                <td>
                  <code>{{ row.spacer_sequence || "Not reported" }}</code>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </template>
    </template>
  </section>
</template>

<style scoped>
.protospacer-results h3 {
  margin-top: 0;
}
.protospacer-reference {
  display: grid;
  gap: 0.9rem;
  margin-block: 1.5rem;
}
.protospacer-reference div {
  display: grid;
  grid-template-columns: minmax(10rem, 14rem) minmax(0, 1fr);
  gap: 1rem;
}
.protospacer-reference dt {
  font-weight: 600;
}
.protospacer-reference dd {
  margin: 0;
  overflow-wrap: anywhere;
}
.protospacer-reference code {
  font-size: 0.82rem;
}
.protospacer-table td:first-child,
.protospacer-table td:nth-child(4) {
  max-width: 20rem;
  overflow-wrap: anywhere;
  white-space: normal;
}
.protospacer-table code {
  display: block;
  min-width: 15rem;
  max-width: 28rem;
  white-space: normal;
  overflow-wrap: anywhere;
  line-height: 1.6;
}
@media (max-width: 650px) {
  .protospacer-reference div {
    grid-template-columns: 1fr;
    gap: 0.3rem;
  }
}
</style>
