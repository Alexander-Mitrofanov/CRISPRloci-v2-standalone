"use strict";
const $ = (id) => document.getElementById(id);
const token = window.location.hash.slice(1);
let options;
let selection;
let polling = false;
const modeNames = {
  sample: "Bundled example",
  download: "Download public reference",
  local: "Use local files",
  historical: "Download 2013 reference",
  import: "Import export bundle",
};

function error(message) {
  $("message").textContent = message;
  $("message").hidden = !message;
}
async function api(path, body) {
  const response = await fetch(`/api/${path}`, {
    method: body === undefined ? "GET" : "POST",
    headers: {
      Authorization: `Bearer ${token}`,
      "Content-Type": "application/json",
    },
    body: body === undefined ? undefined : JSON.stringify(body),
    cache: "no-store",
    credentials: "omit",
    referrerPolicy: "no-referrer",
  });
  const result = await response.json();
  if (!response.ok)
    throw new Error(
      result.error || "The installer could not complete this request.",
    );
  return result;
}
function screen(name, index) {
  for (const id of ["choose", "review-panel", "progress-panel"])
    $(id).hidden = id !== name;
  document.querySelectorAll("#navigation li").forEach((item, i) => {
    item.className = i === index ? "current" : "";
    if (i === index) item.setAttribute("aria-current", "step");
    else item.removeAttribute("aria-current");
  });
}
function bytes(value) {
  return value >= 1024 ** 3
    ? `${(value / 1024 ** 3).toFixed(1)} GiB`
    : `${(value / 1024 ** 2).toFixed(1)} MiB`;
}
function collect() {
  const value = {
    cpu: Number($("cpu").value),
    developer: $("developer").checked,
    references: {},
  };
  for (const key of Object.keys(options.references))
    value.references[key] = {
      mode: $(`${key}-mode`).value,
      path: $(`${key}-path`).value.trim(),
      sha256: $(`${key}-sha`).value.trim(),
    };
  return value;
}
function reference(key, spec, initial) {
  const row = document.createElement("div");
  row.className = "reference";
  const top = document.createElement("div");
  top.className = "reference-top";
  const label = document.createElement("label");
  label.htmlFor = `${key}-mode`;
  label.textContent = spec.title;
  const select = document.createElement("select");
  select.id = `${key}-mode`;
  for (const mode of spec.modes) {
    const option = document.createElement("option");
    option.value = mode;
    option.textContent = modeNames[mode];
    select.append(option);
  }
  select.value = initial.mode;
  top.append(label, select);
  row.append(top);
  const description = document.createElement("p");
  description.textContent = spec.description;
  row.append(description);
  const fields = document.createElement("div");
  fields.className = "path-fields";
  for (const [field, title] of [
    ["path", "Absolute local path"],
    ["sha", "Export bundle SHA-256"],
  ]) {
    const wrapper = document.createElement("div");
    wrapper.id = `${key}-${field}-field`;
    const name = document.createElement("label");
    name.htmlFor = `${key}-${field}`;
    name.textContent = title;
    const input = document.createElement("input");
    input.id = `${key}-${field}`;
    input.type = "text";
    input.value = field === "path" ? initial.path : initial.sha256;
    input.autocomplete = "off";
    input.spellcheck = false;
    wrapper.append(name, input);
    fields.append(wrapper);
  }
  function changed() {
    $(`${key}-path-field`).hidden = !["local", "import"].includes(select.value);
    $(`${key}-sha-field`).hidden = select.value !== "import";
  }
  select.addEventListener("change", changed);
  row.append(fields);
  $("references").append(row);
  changed();
}
function render(state) {
  if (state.status === "idle") return;
  screen("progress-panel", 2);
  const running = state.status === "running";
  const complete = state.status === "complete";
  $("status-title").textContent = complete
    ? "Your workspace is ready."
    : running
      ? "Preparing your workspace."
      : state.status === "cancelled"
        ? "Installation cancelled."
        : "Installation needs attention.";
  $("status-description").textContent = complete
    ? "The scientific release passed validation. Start the local interface to begin an analysis."
    : state.error ||
      "Completed downloads and verified environments are reused when you retry.";
  $("overall").max = state.total || 1;
  $("overall").value = state.completed || 0;
  $("overall-text").textContent =
    `${state.completed || 0} of ${state.total || 0} stages`;
  $("stage-title").textContent = complete
    ? "Installation complete"
    : state.stages.find((s) => s.id === state.current)?.title || "Preparing";
  const progress = state.progress;
  $("tool-progress").hidden = state.current !== "tools" || !state.tool_progress;
  if (state.tool_progress?.total) {
    $("tools").max = state.tool_progress.total;
    $("tools").value = state.tool_progress.completed;
    $("tool-text").textContent =
      `${state.tool_progress.completed} / ${state.tool_progress.total}`;
  }
  $("operation-title").textContent =
    progress?.label || (running ? "Working on the current stage…" : "Finished");
  if (progress?.total) {
    $("operation").max = progress.total;
    $("operation").value = progress.completed;
  } else if (complete) {
    $("operation").max = 1;
    $("operation").value = 1;
  } else $("operation").removeAttribute("value");
  $("operation-text").textContent = progress
    ? progress.unit === "bytes"
      ? `${bytes(progress.completed)}${progress.total ? ` / ${bytes(progress.total)}` : " downloaded"}`
      : `${progress.completed}${progress.total ? ` / ${progress.total}` : ""} ${progress.unit}`
    : "";
  $("stage-list").replaceChildren(
    ...state.stages.map((stage, index) => {
      const item = document.createElement("li");
      item.textContent = stage.title;
      item.className =
        index < state.completed
          ? "done"
          : stage.id === state.current && running
            ? "active"
            : "";
      return item;
    }),
  );
  $("log").textContent = state.logs.join("\n");
  $("cancel").hidden = !running;
  $("retry").hidden = running || complete;
  $("start").hidden = !complete;
  $("start-help").hidden = !complete;
}
async function poll() {
  if (polling) return;
  polling = true;
  try {
    const state = await api("state");
    render(state);
    if (state.status === "running") window.setTimeout(poll, 700);
  } catch (e) {
    error(e.message);
  } finally {
    polling = false;
  }
}
async function appPoll() {
  try {
    const state = await api("app");
    $("open").hidden = !state.ready;
    $("stop-app").hidden = !state.active;
    $("start").disabled = state.active;
    if (state.active && !state.ready) {
      $("status-description").textContent =
        "Starting the validated local services…";
      window.setTimeout(appPoll, 1000);
    } else if (state.ready)
      $("status-description").textContent =
        "The local application is running. Open CRISPRloci to begin.";
  } catch (e) {
    error(e.message);
  }
}
$("review").addEventListener("click", async () => {
  error("");
  try {
    const plan = await api("plan", collect());
    selection = plan.choices;
    $("plan").replaceChildren(
      ...plan.stages.map((stage) => {
        const item = document.createElement("li");
        item.textContent = stage.title;
        return item;
      }),
    );
    $("install").disabled = plan.system.errors.length > 0;
    if (plan.system.errors.length) error(plan.system.errors.join("\n"));
    screen("review-panel", 1);
  } catch (e) {
    error(e.message);
  }
});
$("back").addEventListener("click", () => {
  error("");
  screen("choose", 0);
});
$("retry").addEventListener("click", () => {
  error("");
  screen("choose", 0);
});
$("install").addEventListener("click", async () => {
  error("");
  $("install").disabled = true;
  try {
    render(await api("install", selection));
    poll();
  } catch (e) {
    error(e.message);
    $("install").disabled = false;
  }
});
$("cancel").addEventListener("click", async () => {
  try {
    await api("cancel", {});
    $("status-description").textContent =
      "Stopping the current task. Completed verified work is kept for a retry.";
  } catch (e) {
    error(e.message);
  }
});
$("start").addEventListener("click", async () => {
  try {
    await api("app/start", {});
    appPoll();
  } catch (e) {
    error(e.message);
  }
});
$("stop-app").addEventListener("click", async () => {
  try {
    await api("app/stop", {});
    $("status-description").textContent = "The local services have stopped.";
    appPoll();
  } catch (e) {
    error(e.message);
  }
});
(async () => {
  try {
    options = await api("options");
    $("platform").textContent =
      `${options.system.platform} ${options.system.architecture}`;
    $("storage").textContent = `${bytes(options.system.free_bytes)} free`;
    $("location").textContent = options.system.project;
    $("cpu").value = options.defaults.cpu;
    for (const [key, spec] of Object.entries(options.references))
      reference(key, spec, options.defaults.references[key]);
    if (options.system.errors.length) error(options.system.errors.join("\n"));
    await poll();
  } catch (e) {
    error(e.message);
  }
})();
