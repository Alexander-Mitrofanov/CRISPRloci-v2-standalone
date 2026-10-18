<script setup>
defineProps({
  elements: { type: Array, default: () => [] },
  counts: { type: Object, default: () => ({}) },
  truncated: { type: Boolean, default: false },
});
const metric = (value) =>
  value === null || value === undefined ? "Not reported" : Number(value).toExponential(2);
</script>

<template>
  <section class="result-panel">
    <h3>Insertion-sequence context</h3>
    <p>
      {{ counts.elements ?? "Not reported" }} predicted elements. These predictions support
      inspection of locus organization; an overlap does not establish gene disruption or loss of
      function.
    </p>
    <p v-if="!elements.length">
      {{
        truncated || counts.elements !== 0
          ? "No insertion sequences in this preview."
          : "No insertion sequences were reported with the pinned reference models. This does not establish their absence."
      }}
    </p>
    <div
      v-else
      class="table-scroll"
    >
      <table>
        <thead>
          <tr>
            <th>Source</th>
            <th>Interval</th>
            <th>Family</th>
            <th>Native call</th>
            <th>Transposase ORF</th>
            <th>Copy E-value</th>
            <th>Inverted repeats</th>
          </tr>
        </thead>
        <tbody>
          <tr
            v-for="element in elements"
            :key="element.id"
          >
            <td>{{ element.source_id }}</td>
            <td>
              {{ element.start }}–{{ element.end
              }}<span v-if="element.contig_edge"> · contig edge</span>
            </td>
            <td>{{ element.family }}</td>
            <td>
              {{
                element.native_type === "c"
                  ? "Complete prediction"
                  : element.native_type === "p"
                    ? "Partial prediction"
                    : "Not reported"
              }}
            </td>
            <td>
              {{
                element.transposase_orf
                  ? `${element.transposase_orf.start}–${element.transposase_orf.end} (${element.transposase_orf.strand})`
                  : "Not identified"
              }}
            </td>
            <td>{{ metric(element.evalue_this_copy) }}</td>
            <td>{{ element.terminal_inverted_repeats ? "Reported" : "Not reported" }}</td>
          </tr>
        </tbody>
      </table>
    </div>
    <p v-if="truncated">
      This preview is limited. Download the complete insertion-context report for all predictions.
    </p>
    <p>
      Coordinates are 1-based and inclusive. Strand describes the transposase ORF. Target-site
      duplications were not assessed; terminal inverted repeats are different evidence.
    </p>
  </section>
</template>
