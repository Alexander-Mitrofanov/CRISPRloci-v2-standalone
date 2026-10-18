<script setup>
defineProps({ relations: { type: Array, default: () => [] } });
const number = (value) => (value == null ? "Not reported" : value.toLocaleString());
</script>

<template>
  <details
    v-if="relations.length"
    class="cas-context-relations"
  >
    <summary>Inspect coordinate relationships to Cas cassettes</summary>
    <p>
      Overlap and nearest-cassette distance describe submitted coordinates. Neither establishes gene
      disruption, functional interaction or shared ancestry.
    </p>
    <ul>
      <li
        v-for="relation in relations"
        :key="`${relation.sourceId}:${relation.intervalId}`"
      >
        <strong>{{ relation.sourceId }} · {{ relation.intervalId }}</strong>
        <p v-if="!relation.nearestCassettes.length">
          {{
            relation.nearestCassetteCount === 0
              ? "No Cas cassette was reported for this source."
              : "No Cas cassette relationships are shown in this preview."
          }}
        </p>
        <ul v-else>
          <li
            v-for="cassette in relation.nearestCassettes"
            :key="cassette.id"
          >
            {{ cassette.id }} · {{ number(cassette.start) }}–{{ number(cassette.end) }} ·
            {{ number(cassette.overlapBases) }} overlapping bases ·
            {{ number(cassette.interveningBases) }} intervening bases
          </li>
        </ul>
        <p v-if="relation.nearestCassettesTruncated">
          Nearest-cassette ties are limited in this preview ({{
            number(relation.nearestCassetteCount)
          }}
          total). Inspect the complete cassette annotations for all evidence.
        </p>
      </li>
    </ul>
  </details>
</template>

<style scoped>
li {
  overflow-wrap: anywhere;
  margin: 0.6rem 0;
}
p {
  max-width: 90ch;
}
</style>
