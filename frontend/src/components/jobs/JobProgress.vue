<script setup>
import { computed } from "vue";

import { TERMINAL_STATUSES, stagesForMode } from "../../science.js";
import { formatDate } from "../../utils/formatting.js";
import AppIcon from "../common/AppIcon.vue";
import RecoveryLink from "./RecoveryLink.vue";

const ACTIVE_STATUSES = new Set([
  "queued",
  "running",
  "validate_input",
  "detect_arrays",
  "adapt_arrays",
  "preflight_groups",
  "reconstruct_spacer_histories",
  "compare_orientations",
  "package_results",
]);
const props = defineProps({
  job: { type: Object, required: true },
  credential: { type: Object, required: true },
  cancelling: Boolean,
});
defineEmits(["cancel", "forget"]);

const statusCopy = (status) =>
  ({
    queued: "Queued",
    running: "Analysis running",
    validate_input: "Validating input",
    detect_arrays: "Detecting arrays",
    adapt_arrays: "Preparing arrays",
    preflight_groups: "Checking groups",
    reconstruct_spacer_histories: "Reconstructing histories",
    compare_orientations: "Comparing orientations",
    package_results: "Packaging results",
    completed: "Analysis complete",
    completed_no_eligible_groups: "Detection complete — no eligible evolutionary groups",
    failed: "Analysis failed",
    cancelled: "Analysis cancelled",
    expired: "Results expired",
  })[status] || "Waiting for status";
const stages = computed(() =>
  stagesForMode(props.job?.mode || props.job?.request?.mode || "orientation"),
);
const current = computed(
  () =>
    props.job?.stage ||
    (ACTIVE_STATUSES.has(props.job?.status) && props.job?.status !== "running"
      ? props.job.status
      : null),
);
const currentIndex = computed(() => stages.value.findIndex((stage) => stage.id === current.value));
const terminal = computed(() => TERMINAL_STATUSES.has(props.job?.status));
const successful = computed(() =>
  ["completed", "completed_no_eligible_groups"].includes(props.job?.status),
);
</script>

<template>
  <section
    :class="['job-panel', `job-${job?.status || 'queued'}`]"
    aria-labelledby="job-heading"
  >
    <div class="job-heading">
      <div>
        <h2
          id="job-heading"
          tabindex="-1"
        >
          {{ statusCopy(job?.status || "queued") }}
        </h2>
        <p class="job-id">
          Job <code>{{ credential.jobId }}</code> ·
          {{
            job?.expires_at || credential.expiresAt
              ? `expires ${formatDate(job?.expires_at || credential.expiresAt)}`
              : "retention starts when the run finishes"
          }}
        </p>
      </div>
      <span
        :class="['job-badge', successful ? 'success' : terminal ? 'terminal' : 'active']"
        role="status"
        ><i />{{ successful ? "Ready" : terminal ? statusCopy(job?.status) : "In progress" }}</span
      >
    </div>
    <RecoveryLink :credential="credential" />
    <ol
      v-if="!terminal"
      class="stage-list"
      aria-label="Analysis progress"
    >
      <li
        v-for="(stage, index) in stages"
        :key="stage.id"
        :class="
          currentIndex >= 0 && index < currentIndex
            ? 'complete'
            : index === currentIndex
              ? 'current'
              : 'pending'
        "
        :aria-current="index === currentIndex ? 'step' : undefined"
      >
        <span
          ><AppIcon
            v-if="currentIndex >= 0 && index < currentIndex"
            name="check"
            :size="15"
          /><template v-else>{{ index + 1 }}</template></span
        >
        <div>
          <strong>{{ stage.label }}</strong
          ><small>{{ stage.detail }}</small>
        </div>
      </li>
    </ol>
    <p
      v-if="!terminal && job?.status === 'running' && !job?.stage"
      class="stage-unavailable"
    >
      <AppIcon
        name="info"
        :size="15"
      />
      The workflow is running, but no reliable tool-level stage is exposed.
    </p>
    <div
      v-if="!terminal"
      class="queue-row"
    >
      <span>{{
        job?.queue_position
          ? `Queue position ${job.queue_position}`
          : "Keep this tab open, or copy the recovery link before closing it."
      }}</span
      ><button
        class="cancel-button"
        type="button"
        :disabled="cancelling"
        @click="$emit('cancel')"
      >
        <AppIcon
          name="stop"
          :size="16"
        />{{ cancelling ? "Cancelling…" : "Cancel job" }}
      </button>
    </div>
    <div
      v-if="job?.status === 'failed'"
      class="job-message error"
      role="alert"
    >
      <AppIcon name="warning" />
      <div>
        <strong>The workflow did not complete</strong>
        <p>{{ job.error?.message || job.error || "The service reported an analysis failure." }}</p>
      </div>
    </div>
    <div
      v-if="job?.status === 'cancelled'"
      class="job-message"
    >
      <AppIcon name="info" />
      <div>
        <strong>Job cancelled</strong>
        <p>Partial working files are not presented as completed results.</p>
      </div>
    </div>
    <div
      v-if="terminal"
      class="queue-row"
    >
      <span>Download the results and copy the recovery link before leaving this job.</span
      ><button
        class="leave-job-button"
        type="button"
        @click="$emit('forget')"
      >
        Leave this job and start another
      </button>
    </div>
  </section>
</template>
