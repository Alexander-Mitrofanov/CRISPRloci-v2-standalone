<script setup>
import { useId } from "vue";
import { formatNumber } from "../../../utils/formatting.js";

defineProps({
  evidence: { type: Object, required: true },
  compact: { type: Boolean, default: false },
});
const headingId = useId();
const metric = (value) => (value == null ? "Not available" : formatNumber(value));
const location = (row) =>
  row.segments.length
    ? row.segments
        .map((segment) => `[${formatNumber(segment.start)}, ${formatNumber(segment.end)})`)
        .join(" + ")
    : row.available_length === 0
      ? "Empty at contig boundary"
      : "Not available";
</script>

<template>
  <section
    class="annotation-results"
    :aria-labelledby="headingId"
  >
    <h3 :id="headingId">Leader context</h3>
    <p class="annotation-notice">
      CRISPRleader v2 currently extracts context only. Leader prediction is unavailable. Neither
      side is selected; window length is an extraction setting, not a biological boundary.
    </p>
    <template v-if="evidence.status === 'completed'">
      <div class="annotation-counts">
        <div>
          <strong>{{ metric(evidence.array_count) }}</strong
          ><span>Accepted arrays</span>
        </div>
        <div>
          <strong>{{ metric(evidence.context_count) }}</strong
          ><span>Context windows</span>
        </div>
        <div>
          <strong>{{ metric(evidence.flank_length) }} nt</strong><span>Requested on each side</span>
        </div>
      </div>
      <template v-if="!compact">
        <p class="annotation-lead">
          Intervals are zero-based, half-open, on the submitted forward sequence. Left windows use
          the forward strand; right windows use its reverse complement so both sequences point
          toward the array. This does not infer transcription direction. Bona-fide and Possible
          arrays are included. Unknown topology clips windows at contig boundaries.
        </p>
        <div
          v-if="evidence.contexts.length"
          class="annotation-scroll"
          tabindex="0"
          role="region"
          aria-label="Leader context table"
        >
          <table class="annotation-table">
            <thead>
              <tr>
                <th scope="col">Contig / array</th>
                <th scope="col">Side</th>
                <th scope="col">Interval [start, end)</th>
                <th scope="col">Sequence strand</th>
                <th scope="col">Available / requested nt</th>
                <th scope="col">Sequence evidence</th>
              </tr>
            </thead>
            <tbody>
              <tr
                v-for="row in evidence.contexts"
                :key="`${row.source_id}:${row.id}`"
              >
                <td>
                  {{ row.source_id ?? "Not available" }} / {{ row.array_id ?? "Not available" }}
                </td>
                <td>{{ row.side ?? "Not available" }}</td>
                <td>{{ location(row) }}</td>
                <td>{{ row.strand ?? "Not available" }}</td>
                <td>
                  {{ metric(row.available_length) }} / {{ metric(row.requested_length)
                  }}<small v-if="row.truncated_reason">Clipped at contig boundary</small>
                </td>
                <td>
                  <details class="leader-sequence">
                    <summary>View sequences</summary>
                    <p>Context</p>
                    <code>{{
                      row.sequence === ""
                        ? "Empty at contig boundary"
                        : (row.sequence ?? "Not available")
                    }}</code>
                    <p>Observed terminal repeat</p>
                    <code>{{ row.observed_repeat ?? "Not available" }}</code>
                  </details>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
        <p
          v-else-if="evidence.context_count === 0"
          class="annotation-empty"
        >
          No accepted arrays were available for context extraction. Leader presence or absence was
          not evaluated.
        </p>
        <p
          v-if="evidence.truncated"
          class="annotation-notice"
        >
          This preview is limited. Download the complete context exports from Files &amp; methods.
        </p>
        <p class="annotation-lead">
          Downloads include JSON, TSV, FASTA, GFF3 and BED. JSON and BED use zero-based, half-open
          intervals; GFF3 uses one-based, inclusive intervals. Empty windows remain in JSON and TSV.
        </p>
      </template>
    </template>
    <p
      v-else
      class="annotation-empty"
    >
      No leader context result was reported for this job.
    </p>
  </section>
</template>
