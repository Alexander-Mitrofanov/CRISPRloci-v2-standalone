<script setup>
import { computed } from "vue";

import { categoryClass } from "../../../science.js";
import { asArray, formatNumber, signedNumber } from "../../../utils/formatting.js";
import { comparisonDecisionFor, groupIdentity } from "../../../utils/results.js";

const props = defineProps({
  summary: { type: Object, required: true },
  membershipStatus: { type: String, default: "inline" },
});
const groups = computed(() => asArray(props.summary?.adapter?.groups));
const rows = computed(() => {
  const orientation = props.summary.orientation || {};
  const comparisons = asArray(orientation.comparisons);
  const reconstruction = props.summary.reconstruction || {};
  const reconstructions = asArray(orientation.selected_reconstructions).length
    ? asArray(orientation.selected_reconstructions)
    : asArray(reconstruction.results);
  const comparisonByGroup = new Map(
    comparisons.map((item, index) => [groupIdentity(item, index), item]),
  );
  const reconstructionByGroup = new Map(
    reconstructions.map((item, index) => [groupIdentity(item, index), item]),
  );
  return groups.value.map((group, index) => {
    const name = groupIdentity(group, index);
    const comparison = comparisonByGroup.get(name);
    const selected = reconstructionByGroup.get(name);
    return {
      group,
      name,
      members: asArray(group?.arrays),
      comparison,
      selected,
      decision: comparison ? comparisonDecisionFor(comparison, orientation) : null,
    };
  });
});
</script>

<template>
  <section
    v-if="rows.length"
    class="result-section group-map-section"
    aria-labelledby="group-map-heading"
  >
    <div class="result-heading">
      <div>
        <p class="eyebrow">Connected evidence</p>
        <h3 id="group-map-heading">How detections became evolutionary evidence</h3>
      </div>
      <p>
        Follow calls through their canonical repeat into the evOr comparison and SpacerPlacer
        history.
      </p>
    </div>
    <div class="group-bridge-list">
      <article
        v-for="(row, index) in rows"
        :key="row.name"
        :class="['group-bridge', `group-tone-${index % 4}`]"
      >
        <div class="group-bridge-heading">
          <div>
            <small>Evolutionary group {{ index + 1 }}</small
            ><strong
              >{{ formatNumber(row.group?.array_count ?? row.members.length) }} connected
              arrays</strong
            >
          </div>
          <code :title="row.name">{{ row.name }}</code>
        </div>
        <div class="repeat-band">
          <span
            ><b>Grouping key · canonical repeat</b
            ><small v-if="row.group?.repeat_key"
              >{{ String(row.group.repeat_key).length }} nt</small
            ></span
          ><span
            class="repeat-sequence"
            role="img"
            tabindex="0"
            :aria-label="`Canonical repeat ${row.group?.repeat_key || 'not reported'}`"
            ><i
              v-for="(base, baseIndex) in String(row.group?.repeat_key || '')"
              :key="baseIndex"
              :class="['repeat-base', `repeat-base-${base.toLowerCase()}`]"
              >{{ base }}</i
            ><span
              v-if="!row.group?.repeat_key"
              class="repeat-unavailable"
              >Canonical repeat not reported</span
            ></span
          >
        </div>
        <div class="group-bridge-flow">
          <div class="group-members">
            <div class="flow-label"><span>1</span><strong>CRISPRidentify detections</strong></div>
            <div
              v-if="row.members.length"
              class="group-member-list"
            >
              <div
                v-for="(member, memberIndex) in row.members"
                :key="`${member?.source_id}:${member?.array_id}:${memberIndex}`"
                class="group-member"
              >
                <div>
                  <strong>{{ member?.source_id || "Unknown record" }}</strong
                  ><small>{{ member?.array_id || "Unknown array" }}</small>
                </div>
                <span :class="['category-pill', `category-${categoryClass(member?.category)}`]">{{
                  member?.category || "Unclassified"
                }}</span
                ><span class="member-strand">strand {{ member?.strand || "?" }}</span
                ><span
                  class="mini-spacer-array"
                  :aria-label="`${formatNumber(member?.spacer_count)} detected spacers`"
                  ><i
                    v-for="number in Math.min(22, Number(member?.spacer_count) || 0)"
                    :key="number"
                    class="mini-spacer-count"
                  /><b v-if="Number(member?.spacer_count) > 22"
                    >+{{ Number(member.spacer_count) - 22 }}</b
                  ></span
                >
              </div>
            </div>
            <div
              v-else
              class="membership-pending"
            >
              {{
                membershipStatus === "loading"
                  ? "Loading exact group members…"
                  : "Exact member mapping was not available."
              }}
            </div>
          </div>
          <div
            class="group-connector"
            aria-hidden="true"
          >
            <span>2</span><i /><strong>same repeat<br />shared spacers</strong><b>→</b>
          </div>
          <div class="group-outcomes">
            <div class="flow-label"><span>3</span><strong>Evolutionary results</strong></div>
            <div class="outcome-card outcome-evor">
              <small>CRISPR-evOr orientation</small
              ><template v-if="row.decision"
                ><strong>{{ row.decision.label }}</strong
                ><span
                  >Δ lnL {{ signedNumber(row.decision.delta, 2) }} · boundary ±{{
                    formatNumber(row.decision.threshold, 2)
                  }}</span
                ></template
              ><strong v-else>Not evaluated</strong>
            </div>
            <div class="outcome-card outcome-spacerplacer">
              <small>SpacerPlacer reported history</small
              ><template v-if="row.selected"
                ><strong
                  >{{ formatNumber(row.selected.acquisitions) }} acquisitions ·
                  {{ formatNumber(row.selected.deletions) }} deletions</strong
                ><span
                  >{{ formatNumber(row.selected.unique_spacers) }} unique spacers
                  reconstructed</span
                ></template
              ><strong v-else>Not reconstructed</strong>
            </div>
          </div>
        </div>
      </article>
    </div>
  </section>
</template>
