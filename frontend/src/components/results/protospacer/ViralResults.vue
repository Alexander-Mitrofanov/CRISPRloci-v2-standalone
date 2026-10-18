<script setup>
import { computed, ref, useId } from "vue";
const props = defineProps({
  available: Boolean,
  compact: Boolean,
  completed: { type: Boolean, default: null },
  outcome: { type: String, default: null },
  reference: { type: Object, required: true },
  counts: { type: Object, required: true },
  matches: { type: Array, required: true },
  candidates: { type: Array, required: true },
  queries: { type: Array, required: true },
  substitutions: { type: Number, default: null },
  truncated: Boolean,
});
const heading = useId(),
  filterId = useId(),
  filter = ref("");
const count = (v) => (v == null ? "Not reported" : v.toLocaleString());
const rows = computed(() =>
  props.matches.filter((r) =>
    [r.query_id, r.subject_accession, r.subject_title]
      .join(" ")
      .toLowerCase()
      .includes(filter.value.toLowerCase()),
  ),
);
const skipped = computed(() => props.queries.filter((q) => q.eligible === false));
const accessionUrl = (a) => `https://www.ncbi.nlm.nih.gov/nuccore/${encodeURIComponent(a)}`;
</script>
<template>
  <section
    class="viral-results"
    :aria-labelledby="heading"
  >
    <h3 :id="heading">Candidate viral matches</h3>
    <p class="annotation-lead">
      Full-length, ungapped spacer alignments with up to {{ count(substitutions) }} substitutions.
      BLAST retrieval is heuristic. A sequence match does not establish infection, host range or
      functional targeting.
    </p>
    <p
      v-if="!available || completed !== true"
      class="annotation-empty"
    >
      A completed viral search was not reported. Missing evidence is not a zero-match result.
    </p>
    <template v-else>
      <div class="annotation-counts">
        <div>
          <strong>{{ count(counts.matched_subject_accessions) }}</strong
          ><span>Viral reference records</span>
        </div>
        <div>
          <strong>{{ count(counts.matched_query_occurrences) }}</strong
          ><span>Spacers with matches</span>
        </div>
        <div>
          <strong>{{ count(counts.accepted_hsps) }}</strong
          ><span>Accepted alignments</span>
        </div>
        <div>
          <strong
            >{{ count(counts.eligible_query_occurrences) }} /
            {{ count(counts.query_occurrences) }}</strong
          ><span>Eligible input spacers</span>
        </div>
      </div>
      <p
        v-if="['no_eligible_queries', 'no_raw_hsps', 'all_hsps_filtered'].includes(outcome)"
        class="annotation-empty"
      >
        {{
          outcome === "no_eligible_queries"
            ? "No eligible spacers to search. Use unambiguous DNA spacers of 18–80 nt."
            : "No full-length matches passed this search policy. This does not establish the absence of a viral relationship."
        }}
      </p>
      <p
        v-if="
          !['matches', 'no_eligible_queries', 'no_raw_hsps', 'all_hsps_filtered'].includes(outcome)
        "
        class="annotation-empty"
      >
        Search outcome was not reported. Missing evidence is not a zero-match result.
      </p>
      <p class="viral-reference">
        {{ reference.name || "Reference not reported" }} · release
        {{ reference.source_release || "unknown" }} · {{ count(reference.record_count) }} reference
        records
      </p>
      <details
        v-if="candidates.length && !compact"
        class="viral-candidates"
        open
      >
        <summary>Evidence by viral reference</summary>
        <p>
          Repeated copies of the same spacer do not add independent sequence support. Segments are
          separate records; records are not virus species counts.
        </p>
        <ul>
          <li
            v-for="item in candidates"
            :key="item.accession"
          >
            <div>
              <a
                :href="accessionUrl(item.accession)"
                target="_blank"
                rel="noopener noreferrer"
                >{{ item.accession }}</a
              ><span>{{ item.title }}</span>
            </div>
            <small
              >{{ count(item.unique_query_sequence_count) }} distinct spacer sequences ·
              {{ count(item.location_count) }} locations</small
            >
          </li>
        </ul>
      </details>
      <template v-if="matches.length && !compact"
        ><label :for="filterId">Filter displayed matches by spacer, accession or title</label
        ><input
          :id="filterId"
          v-model="filter"
          type="search"
          class="viral-filter"
        />
        <div
          class="annotation-table-scroll"
          tabindex="0"
          aria-label="Viral match table"
        >
          <table class="annotation-table">
            <thead>
              <tr>
                <th scope="col">Spacer</th>
                <th scope="col">Viral accession</th>
                <th scope="col">Coordinates</th>
                <th scope="col">Strand</th>
                <th scope="col">Substitutions</th>
              </tr>
            </thead>
            <tbody>
              <tr
                v-for="(hit, index) in rows"
                :key="index"
              >
                <td>
                  {{ hit.query_id }}<code>{{ hit.query_sequence }}</code>
                </td>
                <td>
                  <a
                    :href="accessionUrl(hit.subject_accession)"
                    target="_blank"
                    rel="noopener noreferrer"
                    >{{ hit.subject_accession }}</a
                  ><span class="viral-title">{{ hit.subject_title }}</span>
                </td>
                <td>{{ count(hit.subject_start) }}–{{ count(hit.subject_end) }}</td>
                <td>
                  {{
                    hit.relative_strand === "+"
                      ? "Forward (+)"
                      : hit.relative_strand === "-"
                        ? "Reverse (−)"
                        : "Not reported"
                  }}
                </td>
                <td>{{ count(hit.substitutions) }}</td>
              </tr>
            </tbody>
          </table>
        </div>
        <p
          v-if="!rows.length"
          role="status"
        >
          No displayed matches contain that text.
        </p>
      </template>
      <details v-if="skipped.length && !compact">
        <summary>Skipped spacers ({{ skipped.length }} shown)</summary>
        <ul>
          <li
            v-for="q in skipped"
            :key="q.id"
          >
            {{ q.id }}:
            {{
              q.reasons
                .map((r) =>
                  r === "ambiguous_query"
                    ? "ambiguous bases"
                    : "outside the supported length range",
                )
                .join(", ")
            }}
          </li>
        </ul>
      </details>
      <p
        v-if="truncated"
        class="annotation-notice"
      >
        This is a limited preview. Download viral-matches.tsv or viral-matches.json for the complete
        validated result.
      </p>
      <p v-if="compact">
        Open Viral matches for reference evidence, skipped-query reasons and alignments. Complete
        validated evidence is available in Files &amp; methods.
      </p>
      <details v-if="!compact">
        <summary>Reference identity and method</summary>
        <p>
          Coordinates are 1-based inclusive on the forward reference. Both strands are searched. PAM
          compatibility and array overlap are not evaluated. Terminal mismatches may be clipped by
          BLAST and fail the full-length filter.
        </p>
        <code class="viral-hash">{{ reference.reference_sha256 || "Digest not reported" }}</code>
      </details>
    </template>
  </section>
</template>
<style scoped>
.viral-results {
  margin-top: 1rem;
}
.viral-reference {
  color: var(--ink-soft);
  padding: 0.8rem 0;
  border-bottom: 1px solid var(--line);
}
.viral-candidates {
  margin: 1rem 0 1.5rem;
}
.viral-candidates ul {
  list-style: none;
  padding: 0;
  max-height: 28rem;
  overflow: auto;
}
.viral-candidates li {
  display: flex;
  justify-content: space-between;
  gap: 1rem;
  padding: 0.8rem 0;
  border-bottom: 1px solid var(--line);
}
.viral-candidates li div {
  display: flex;
  flex-direction: column;
  gap: 0.25rem;
}
.viral-candidates small {
  flex-shrink: 0;
  color: var(--ink-soft);
}
.viral-results summary {
  cursor: pointer;
  font-weight: 600;
  padding: 0.5rem 0;
}
.viral-filter {
  display: block;
  width: min(100%, 36rem);
  padding: 0.65rem;
  margin: 0.5rem 0 1rem;
  border: 1px solid var(--line-dark);
  border-radius: 0.35rem;
}
.viral-title {
  display: block;
  max-width: 40ch;
  font-size: 0.85rem;
  margin-top: 0.3rem;
}
.viral-results td code {
  display: block;
  margin-top: 0.35rem;
  font-size: 0.8rem;
  overflow-wrap: anywhere;
  max-width: 38ch;
}
.viral-hash {
  overflow-wrap: anywhere;
}
.annotation-table-scroll {
  overflow-x: auto;
}
@media (max-width: 680px) {
  .viral-candidates li {
    flex-direction: column;
  }
  .viral-candidates small {
    white-space: normal;
  }
}
</style>
