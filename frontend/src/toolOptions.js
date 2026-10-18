import catalog from "./contracts/tool-options-v1.json";

export const TOOL_OPTIONS = catalog.tools;

export function toolsForMode(mode, options = {}) {
  return TOOL_OPTIONS.filter((tool) => tool.modes.includes(mode)).map((tool) => ({
    ...tool,
    flags: tool.flags.filter(
      (field) =>
        (!field.modes || field.modes.includes(mode)) &&
        Object.entries(field.requires || {}).every(([flag, required]) => {
          const defaultValue = tool.flags.find((item) => item.flag === flag)?.default;
          return (options[tool.id]?.[flag] ?? defaultValue) === required;
        }),
    ),
  }));
}

export function selectedToolOptions(mode, options = {}) {
  const selected = {};
  for (const tool of toolsForMode(mode, options)) {
    const values = {};
    for (const field of tool.flags) {
      const value = options[tool.id]?.[field.flag];
      if (value !== undefined && value !== field.default) values[field.flag] = value;
    }
    if (Object.keys(values).length) selected[tool.id] = values;
  }
  return selected;
}

export function toolOptionCount(mode, options = {}) {
  return Object.values(selectedToolOptions(mode, options)).reduce(
    (count, values) => count + Object.keys(values).length,
    0,
  );
}

export function toolOptionError(field, value) {
  if (value === undefined || value === field.default) return "";
  if (field.type === "boolean") return typeof value === "boolean" ? "" : "Choose on or off.";
  if (field.type === "enum") {
    return field.choices.includes(value) ? "" : "Choose one of the listed values.";
  }
  if (typeof value !== "number" || !Number.isFinite(value)) return "Enter a number.";
  if (field.type === "integer" && !Number.isInteger(value)) return "Enter a whole number.";
  if (field.min !== undefined && value < field.min) return `Use ${field.min} or more.`;
  if (field.max !== undefined && value > field.max) return `Use ${field.max} or less.`;
  return "";
}

export function validToolOptions(mode, options = {}) {
  return Object.values(toolOptionErrors(mode, options)).every(
    (fields) => Object.keys(fields).length === 0,
  );
}

export function toolOptionErrors(mode, options = {}) {
  const errors = {};
  for (const tool of toolsForMode(mode, options)) {
    errors[tool.id] = {};
    for (const field of tool.flags) {
      const error = toolOptionError(field, options[tool.id]?.[field.flag]);
      if (error) errors[tool.id][field.flag] = error;
    }
    const constraints =
      tool.id === "crispridentify"
        ? [
            ["--min_len_rep", "--max_len_rep"],
            ["--min_len_spacer", "--max_len_spacer"],
          ]
        : tool.id === "crisprrepeat"
          ? [["--span", "--window"]]
          : [];
    for (const [lowerFlag, upperFlag] of constraints) {
      const lowerField = tool.flags.find((field) => field.flag === lowerFlag);
      const upperField = tool.flags.find((field) => field.flag === upperFlag);
      if (!lowerField || !upperField) continue;
      const lower = options[tool.id]?.[lowerFlag] ?? lowerField.default;
      const upper = options[tool.id]?.[upperFlag] ?? upperField.default;
      if (typeof lower === "number" && typeof upper === "number" && lower > upper) {
        errors[tool.id][lowerFlag] = `Must not exceed ${upperFlag} (${upper}).`;
      }
    }
  }
  return errors;
}

export function displayToolDefault(value) {
  if (value === null) return "automatic";
  if (typeof value === "boolean") return value ? "on" : "off";
  return String(value);
}
