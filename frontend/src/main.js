import { createApp } from "vue";

import App from "./App.vue";
import { applyFrameBootPolicy } from "./frameGuard.js";
import "./styles/index.css";

if (applyFrameBootPolicy(window, document)) {
  createApp(App).mount("#root");
}
