<script setup>
import { computed } from "vue";
import { repeatStageText } from "../../../features/results/index.js";
import RepeatPairMatrix from "./RepeatPairMatrix.vue";

const props = defineProps({
  instance: { type: Object, required: true },
});
const metric = (value, digits = 4) => (value == null ? "Not available" : value.toFixed(digits));
const unpaired = computed(() =>
  Array.from({ length: props.instance.rna_sequence?.length || 0 }, (_, index) => ({
    index,
    base: props.instance.rna_sequence[index],
    isolated: props.instance.isolated.unpaired[index],
    context: props.instance.context.unpaired_probabilities[index],
  })),
);
</script>

<template>
  <section
    class="repeat-instance"
    aria-label="Selected repeat detail"
  >
    <h4>
      Selected repeat · {{ instance.source_id
      }}<template v-if="instance.repeat_index != null">
        / repeat {{ instance.repeat_index }}</template
      >
    </h4>
    <p
      v-if="instance.terminal"
      class="repeat-note"
    >
      Terminal repeat: local support may depend on the observed array boundary.
    </p>
    <dl class="repeat-sequences">
      <div>
        <dt>Source sequence</dt>
        <dd>
          <code>{{ instance.source_sequence ?? "Not available" }}</code>
        </dd>
      </div>
      <div>
        <dt>RNA hypothesis (5′→3′)</dt>
        <dd>
          <code>{{ instance.rna_sequence ?? "Not available" }}</code>
        </dd>
      </div>
      <div>
        <dt>Isolated MFE structure</dt>
        <dd>
          <code>{{ instance.isolated.structure ?? "Not available" }}</code>
        </dd>
      </div>
    </dl>
    <p
      v-if="instance.isolated.status !== 'completed'"
      class="annotation-empty"
    >
      {{ repeatStageText(instance.isolated) }}
    </p>
    <p
      v-else
      class="repeat-note"
    >
      MFE {{ metric(instance.isolated.mfe_kcal_mol, 2) }} kcal/mol · ensemble free energy
      {{ metric(instance.isolated.ensemble_kcal_mol, 2) }} kcal/mol. The isolated MFE motif is a
      structural hypothesis.
    </p>
    <dl
      class="repeat-comparison"
      aria-label="Motif support comparison"
    >
      <div>
        <dt>Isolated pair support</dt>
        <dd>{{ metric(instance.comparison.isolated_pair_support) }}</dd>
      </div>
      <div>
        <dt>Context pair support</dt>
        <dd>{{ metric(instance.comparison.context_pair_support) }}</dd>
      </div>
      <div>
        <dt>Support change</dt>
        <dd>{{ metric(instance.comparison.support_change) }}</dd>
      </div>
      <div>
        <dt>Motif pairs</dt>
        <dd>{{ instance.comparison.motif_pair_count ?? "Not available" }}</dd>
      </div>
    </dl>
    <p
      v-if="instance.comparison.status !== 'completed'"
      class="annotation-empty"
    >
      {{ repeatStageText(instance.comparison) }}
    </p>
    <p class="repeat-note">
      Mean support for the listed MFE pairs, not the probability that the entire motif forms. Values
      describe thermodynamic predictions, not cleavage or immune activity.
    </p>
    <p
      v-if="instance.context.status === 'completed'"
      class="repeat-note"
    >
      Local context: effective window {{ instance.context.effective_window ?? "not reported" }} nt;
      maximum span {{ instance.context.effective_span ?? "not reported" }} nt. Context uses pair
      averages, not a full-array MFE structure.
    </p>
    <RepeatPairMatrix
      v-if="instance.rna_sequence && instance.isolated.status === 'completed'"
      :sequence="instance.rna_sequence"
      :isolated-pairs="instance.isolated.pairs"
      :context-pairs="instance.context.intrarepeat_pairs"
      :motif-pairs="instance.isolated.motif_pairs"
      :isolated-available="
        instance.isolated.status === 'completed' && instance.isolated.pairs_available
      "
      :context-available="
        instance.context.status === 'completed' && instance.context.pairs_available
      "
    />
    <details
      v-if="instance.isolated.status === 'completed'"
      class="repeat-details"
    >
      <summary>Single-base accessibility and sequence positions</summary>
      <p>
        Zero-based positions within the selected RNA hypothesis. Context values are native unpaired
        probabilities, not one minus the row sum.
      </p>
      <div
        class="annotation-scroll"
        tabindex="0"
        role="region"
        aria-label="Unpaired probabilities table"
      >
        <table class="annotation-table">
          <thead>
            <tr>
              <th>Position</th>
              <th>Base</th>
              <th>Isolated unpaired</th>
              <th>Context unpaired</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="base in unpaired"
              :key="base.index"
            >
              <td>{{ base.index }}</td>
              <td>{{ base.base }}</td>
              <td>{{ metric(base.isolated) }}</td>
              <td>{{ metric(base.context) }}</td>
            </tr>
          </tbody>
        </table>
      </div>
    </details>
    <details
      v-if="instance.context.status === 'completed'"
      class="repeat-details"
    >
      <summary>
        Pairs to the surrounding array ({{
          instance.context.competitor_count ?? "count not reported"
        }})
      </summary>
      <p>
        Zero-based positions in the oriented array. Each row is an individual predicted interaction;
        values must not be summed.
      </p>
      <p
        v-if="instance.context.competitors_truncated"
        class="annotation-notice"
      >
        This interaction preview is limited. Complete pair exports are in Files &amp; methods.
      </p>
      <div
        v-if="instance.context.competitors.length"
        class="annotation-scroll"
        tabindex="0"
        role="region"
        aria-label="Array interactions table"
      >
        <table class="annotation-table">
          <thead>
            <tr>
              <th>Array position i</th>
              <th>Array position j</th>
              <th>Local pair support</th>
            </tr>
          </thead>
          <tbody>
            <tr
              v-for="pair in instance.context.competitors"
              :key="`${pair.i}:${pair.j}`"
            >
              <td>{{ pair.i }}</td>
              <td>{{ pair.j }}</td>
              <td>{{ metric(pair.probability) }}</td>
            </tr>
          </tbody>
        </table>
      </div>
      <p v-else-if="instance.context.competitor_count === 0">
        No outside-repeat pairs were reported under these settings.
      </p>
    </details>
    <details class="repeat-details">
      <summary>Reference similarity</summary>
      <p>{{ repeatStageText(instance.references) }}</p>
      <p v-if="instance.references.release">
        Reference release: {{ instance.references.release }}.
      </p>
      <p v-if="instance.references.total_hit_count != null">
        {{ instance.references.reported_hit_count ?? "Not reported" }} reported /
        {{ instance.references.total_hit_count }} matching records. Full alignments and origin
        evidence are in Files &amp; methods. Reference similarity does not establish subtype or
        query host.
      </p>
    </details>
  </section>
</template>
