import { asArray } from "./formatting.js";
import { layoutNewick, parseNewick } from "./newick.js";

export const MAX_HISTORY_SPACER_COLUMNS = 160;
export const MAX_HISTORY_LEAVES = 40;
export const MAX_HISTORY_NODES = 120;

export function spacerColorClass(spacer) {
  const number = Math.abs(Number(spacer));
  return `history-color-${Number.isFinite(number) ? number % 12 : 0}`;
}

export function publicNodeName(value) {
  const text = String(value || "unnamed node");
  return text.includes("__") ? text.split("__")[0] : text;
}

export function historyLossCount(node) {
  return asArray(node?.loss_blocks).reduce((total, block) => total + asArray(block).length, 0);
}

export function historyValueList(value) {
  const values = asArray(value)
    .flat(4)
    .filter((item) => item != null && String(item).trim());
  return values.length ? values.map(String).join(", ") : "—";
}

export function specialEventSummary(node) {
  return (
    [
      ["contradictions", node?.contradictions],
      ["duplications", node?.duplications],
      ["rearrangements", node?.rearrangements],
      ["reacquisitions", node?.reacquisitions],
      ["independent gains", node?.independent_gains],
      ["other duplications", node?.other_duplication_events],
    ]
      .filter(([, values]) => historyValueList(values) !== "—")
      .map(([label, values]) => `${label}: ${historyValueList(values)}`)
      .join("; ") || "—"
  );
}

export function entryRootGains(entry) {
  const root = parseNewick(entry?.newick);
  if (!root) return null;
  const rootData = asArray(entry?.nodes).find((node) => String(node?.name) === String(root.name));
  return rootData ? asArray(rootData.gains).length : null;
}

export function entryTreeHeight(entry) {
  return layoutNewick(entry?.newick, 100)?.maxDistance ?? null;
}
