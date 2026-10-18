<script setup>
defineProps({ name: { type: String, default: "info" }, size: { type: Number, default: 20 } });
</script>

<template>
  <svg
    class="icon"
    :width="size"
    :height="size"
    viewBox="0 0 24 24"
    aria-hidden="true"
    fill="none"
    stroke="currentColor"
    stroke-width="1.8"
    stroke-linecap="round"
    stroke-linejoin="round"
  >
    <template v-if="name === 'upload'">
      <path d="M12 16V4M7 9l5-5 5 5" />
      <path d="M5 15v4h14v-4" />
    </template>
    <path
      v-else-if="name === 'check'"
      d="m5 12 4 4L19 6"
    />
    <template v-else-if="name === 'arrow'"><path d="M5 12h14M14 7l5 5-5 5" /></template>
    <template v-else-if="name === 'download'">
      <path d="M12 4v12M7 11l5 5 5-5" />
      <path d="M5 20h14" />
    </template>
    <rect
      v-else-if="name === 'stop'"
      x="6"
      y="6"
      width="12"
      height="12"
      rx="2"
    />
    <template v-else-if="name === 'file'">
      <path d="M7 3h7l4 4v14H7z" />
      <path d="M14 3v5h5M10 13h5M10 17h5" />
    </template>
    <template v-else-if="name === 'shield'">
      <path d="M12 3 5 6v5c0 4.6 2.8 8.1 7 10 4.2-1.9 7-5.4 7-10V6z" />
      <path d="m9 12 2 2 4-5" />
    </template>
    <template v-else-if="name === 'warning'">
      <path d="M12 3 2.8 20h18.4z" />
      <path d="M12 9v5M12 17.5v.1" />
    </template>
    <template v-else-if="name === 'external'">
      <path d="M10 6H6v12h12v-4M13 5h6v6M19 5l-8 8" />
    </template>
    <template v-else-if="name === 'link'">
      <path d="M10 13a5 5 0 0 0 7.1.1l2-2a5 5 0 0 0-7.1-7.1l-1.1 1.1" />
      <path d="M14 11a5 5 0 0 0-7.1-.1l-2 2A5 5 0 0 0 12 20l1.1-1.1" />
    </template>
    <path
      v-else-if="name === 'plus'"
      d="M12 5v14M5 12h14"
    />
    <path
      v-else-if="name === 'minus'"
      d="M5 12h14"
    />
    <template v-else-if="name === 'tree'">
      <path d="M7 4v16M7 8h6M7 16h6" />
      <circle
        cx="16"
        cy="8"
        r="2"
      />
      <circle
        cx="16"
        cy="16"
        r="2"
      />
    </template>
    <template v-else>
      <circle
        cx="12"
        cy="12"
        r="9"
      />
      <path d="M12 11v6M12 7.5v.1" />
    </template>
  </svg>
</template>
