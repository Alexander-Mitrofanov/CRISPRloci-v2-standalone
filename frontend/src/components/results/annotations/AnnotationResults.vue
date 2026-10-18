<script setup>
import { computed } from "vue";
import { formatNumber, finiteMetric } from "../../../utils/formatting.js";

const props = defineProps({
  cas: { type: Object, required: true },
  tracr: { type: Object, required: true },
  arrays: { type: Array, default: () => [] },
  view: { type: String, default: "all" },
  showMap: { type: Boolean, default: true },
});
const heading = computed(
  () =>
    ({ cas: "Cas systems", tracrrna: "tracrRNA candidates" })[props.view] || "Locus annotations",
);
const headingId = computed(() => `annotation-${props.view}-heading`);
const showOverview = computed(() => ["all", "overview"].includes(props.view));
const hasCas = computed(() => props.cas.status === "completed");
const hasTracr = computed(() => props.tracr.status === "completed");
const sources = computed(() => (props.tracr.sources || []).slice(0, 20));
const metric = (value) => (finiteMetric(value) == null ? "Not available" : formatNumber(value));
const location = (row) => `${formatNumber(row.start)}–${formatNumber(row.end)}`;
const featureX = (start, length) => 15 + ((Math.max(1, start) - 1) / length) * 770;
const featureWidth = (row, length) => Math.max(2, ((row.end - row.start + 1) / length) * 770);
const forSource = (rows, id) => (rows || []).filter((row) => row.source_id === id);
</script>

<template>
  <section
    class="annotation-results"
    :aria-labelledby="headingId"
  >
    <h3 :id="headingId">{{ heading }}</h3>
    <p
      v-if="
        (view === 'cas' && !hasCas) || (view === 'tracrrna' && !hasTracr) || (!hasCas && !hasTracr)
      "
      class="annotation-empty"
    >
      No annotation result was reported for this view.
    </p>
    <p
      v-if="!showOverview"
      class="annotation-lead"
    >
      Coordinates are 1-based, end-inclusive; strand is reported separately.
    </p>
    <div
      v-if="showOverview"
      class="annotation-counts"
    >
      <div v-if="hasCas">
        <strong>{{ metric(cas.cas_gene_count) }}</strong
        ><span>Cas protein calls</span>
      </div>
      <div v-if="hasCas">
        <strong>{{ metric(cas.cassette_count) }}</strong
        ><span>Cas cassettes</span>
      </div>
      <div v-if="hasTracr">
        <strong>{{ metric(tracr.prediction_count) }}</strong
        ><span>tracrRNA candidates</span>
      </div>
    </div>
    <details
      v-if="showMap && showOverview && sources.length"
      class="result-details"
    >
      <summary>Genomic context</summary>
      <p class="annotation-lead">
        Coordinates are 1-based and end-inclusive. Proximity does not establish a functional
        association.
      </p>
      <div class="locus-legend">
        <span class="feature-array">CRISPR array</span><span class="feature-cas">Cas gene</span
        ><span class="feature-tracr">tracrRNA candidate</span>
      </div>
      <figure
        v-for="source in sources"
        :key="source.id"
        class="locus-strip"
      >
        <strong>{{ source.id }}</strong>
        <div class="annotation-scroll">
          <svg
            viewBox="0 0 800 85"
            role="img"
            :aria-label="`Genomic features on ${source.id}`"
          >
            <path
              d="M15 45H785"
              stroke="#b7c5d5"
              stroke-width="2"
            />
            <rect
              v-for="row in forSource(arrays, source.id)"
              :key="`a-${row.array_id}`"
              :x="featureX(row.start, source.length)"
              y="29"
              :width="featureWidth(row, source.length)"
              height="30"
              fill="#0f766e"
            >
              <title>Array {{ row.array_id }}: {{ location(row) }}</title>
            </rect>
            <rect
              v-for="row in forSource(cas.genes, source.id)"
              :key="`g-${row.id}`"
              :x="featureX(row.start, source.length)"
              y="34"
              :width="featureWidth(row, source.length)"
              height="22"
              fill="#3158a5"
            >
              <title>{{ row.profile || "Cas gene" }}: {{ location(row) }} ({{ row.strand }})</title>
            </rect>
            <rect
              v-for="row in forSource(tracr.candidates, source.id)"
              :key="`t-${row.id}`"
              :x="featureX(row.start, source.length)"
              y="12"
              :width="featureWidth(row, source.length)"
              height="15"
              fill="#a76716"
            >
              <title>tracrRNA candidate {{ row.id }}: {{ location(row) }}</title>
            </rect>
            <text
              x="15"
              y="78"
              fill="#52647a"
              font-size="12"
            >
              1
            </text>
            <text
              x="785"
              y="78"
              text-anchor="end"
              fill="#52647a"
              font-size="12"
            >
              {{ formatNumber(source.length) }} nt
            </text>
          </svg>
        </div>
        <figcaption>
          Small features are enlarged. Exact intervals are in the result tables.
        </figcaption>
      </figure>
    </details>
    <template v-if="hasCas && ['all', 'cas'].includes(view)">
      <h4>Cas cassettes</h4>
      <p class="annotation-lead">
        Model: {{ cas.model || "Not reported" }}. Classification confidence is a model output, not
        proof of biological activity.
      </p>
      <div
        v-if="cas.cassettes?.length"
        class="annotation-scroll"
        tabindex="0"
        role="region"
        aria-label="Cas cassette table"
      >
        <table class="annotation-table">
          <thead>
            <tr>
              <th scope="col">Contig</th>
              <th scope="col">Interval</th>
              <th scope="col">Class</th>
              <th scope="col">Type / subtype</th>
              <th scope="col">Cas genes</th>
              <th scope="col">Method</th>
              <th scope="col">Confidence</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="row in cas.cassettes"
              :key="`${row.source_id}:${row.id}`"
            >
              <td>{{ row.source_id }}</td>
              <td>{{ location(row) }}</td>
              <td>{{ row.class ?? "Unassigned" }}</td>
              <td>{{ row.type ?? "Unassigned" }} / {{ row.subtype ?? "Unassigned" }}</td>
              <td>{{ metric(row.cas_gene_count) }}</td>
              <td>{{ row.method ?? "Not reported" }}</td>
              <td>{{ metric(row.confidence) }}</td>
            </tr>
          </tbody>
        </table>
      </div>
      <p
        v-else
        class="annotation-empty"
      >
        No Cas cassette passed the model's evidence criteria.
      </p>
      <h4>Cas protein calls</h4>
      <div
        v-if="cas.genes?.length"
        class="annotation-scroll"
        tabindex="0"
        role="region"
        aria-label="Cas gene table"
      >
        <table class="annotation-table">
          <thead>
            <tr>
              <th scope="col">Contig</th>
              <th scope="col">Protein</th>
              <th scope="col">Interval</th>
              <th scope="col">Strand</th>
              <th scope="col">Profile</th>
              <th scope="col">Type / subtype</th>
              <th scope="col">Score margin</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="row in cas.genes"
              :key="`${row.source_id}:${row.id}`"
            >
              <td>{{ row.source_id }}</td>
              <td>{{ row.id }}</td>
              <td>{{ location(row) }}</td>
              <td>{{ row.strand }}</td>
              <td>{{ row.profile ?? "Not reported" }}</td>
              <td>{{ row.type ?? "Unassigned" }} / {{ row.subtype ?? "Unassigned" }}</td>
              <td>{{ metric(row.score_margin) }}</td>
            </tr>
          </tbody>
        </table>
      </div>
      <p
        v-else
        class="annotation-empty"
      >
        No protein passed the Cas calling criteria.
      </p>
      <p
        v-if="cas.truncated"
        class="annotation-notice"
      >
        This preview is limited. Download the annotation exports for all validated results.
      </p>
    </template>
    <template v-if="hasTracr && ['all', 'tracrrna'].includes(view)">
      <h4>tracrRNA candidates</h4>
      <p class="annotation-lead">
        {{
          tracr.mode === "complete" ? "Complete evidence workflow" : "Covariance-model screening"
        }}
        · {{ tracr.model_type === "V" ? "Type V-K" : "Type II" }}. Predictions are hypotheses.
        {{
          tracr.mode === "complete"
            ? "Scores are rankings"
            : "Scores are covariance-model bit scores"
        }}, not probabilities. Intervals do not establish complete native transcript boundaries.
      </p>
      <div
        v-if="tracr.candidates?.length"
        class="annotation-scroll"
        tabindex="0"
        role="region"
        aria-label="tracrRNA candidate table"
      >
        <table class="annotation-table">
          <thead>
            <tr>
              <th scope="col">Contig</th>
              <th scope="col">Candidate</th>
              <th scope="col">Interval</th>
              <th scope="col">Strand</th>
              <th scope="col">Origin</th>
              <th scope="col">Evidence</th>
              <th scope="col">Score</th>
              <th scope="col">Model E-value</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="row in tracr.candidates"
              :key="row.id"
            >
              <td>{{ row.source_id }}</td>
              <td>{{ row.id }}</td>
              <td>{{ location(row) }}</td>
              <td>{{ row.strand }}</td>
              <td>{{ row.origin ?? "Not reported" }}</td>
              <td>{{ row.evidence_class ?? "Not reported" }}</td>
              <td>{{ metric(row.score) }}</td>
              <td>{{ row.evalue == null ? "Not available" : row.evalue.toExponential(2) }}</td>
            </tr>
          </tbody>
        </table>
      </div>
      <p
        v-else
        class="annotation-empty"
      >
        No tracrRNA candidate passed the selected workflow's criteria. This does not establish
        absence of tracrRNA.
      </p>
      <p
        v-if="tracr.truncated"
        class="annotation-notice"
      >
        The candidate preview is limited. Download the annotation exports for all candidates.
      </p>
    </template>
  </section>
</template>
