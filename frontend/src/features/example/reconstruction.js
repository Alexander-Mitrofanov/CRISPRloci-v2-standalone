import { canonicalInteger, fail, finite, nonNegativeInteger, object } from "./contract.js";
import { parseNewick } from "./newick.js";

const RECONSTRUCTION_KEYS = new Set([
  "group",
  "hypothesis",
  "newick",
  "spacer_order",
  "nodes",
  "acquisition_count",
  "deletion_count",
]);
const RECONSTRUCTION_NODE_KEYS = new Set([
  "name",
  "spacers",
  "gains",
  "loss_blocks",
  "contradictions",
  "duplications",
  "rearrangements",
  "reacquisitions",
  "independent_gains",
  "other_duplication_events",
]);
const SPECIAL_EVENT_FIELDS = [
  "contradictions",
  "duplications",
  "rearrangements",
  "reacquisitions",
  "independent_gains",
  "other_duplication_events",
];

function validateSpacerIds(values, canonicalIds, { allowEmpty = true } = {}) {
  if (
    !Array.isArray(values) ||
    (!allowEmpty && values.length === 0) ||
    values.some((value) => !canonicalInteger(value) || !canonicalIds.has(value)) ||
    new Set(values).size !== values.length
  )
    fail();
}

function validateCanonicalReconstruction(reconstruction, treeNewick, group, members) {
  if (
    !object(reconstruction) ||
    Object.keys(reconstruction).some((key) => !RECONSTRUCTION_KEYS.has(key)) ||
    Object.keys(reconstruction).length !== RECONSTRUCTION_KEYS.size ||
    reconstruction.group !== group ||
    reconstruction.newick !== treeNewick ||
    !Array.isArray(reconstruction.spacer_order) ||
    !Array.isArray(reconstruction.nodes) ||
    !nonNegativeInteger(reconstruction.acquisition_count) ||
    !nonNegativeInteger(reconstruction.deletion_count)
  )
    fail();
  const canonicalIds = new Set(reconstruction.spacer_order);
  if (
    canonicalIds.size !== reconstruction.spacer_order.length ||
    reconstruction.spacer_order.some((value) => !canonicalInteger(value)) ||
    [...canonicalIds]
      .sort((left, right) => left - right)
      .some((value, index) => value !== index + 1)
  )
    fail();
  const tree = parseNewick(treeNewick);
  const treeNodeNames = new Set(tree.nodeNames);
  const nodesByName = new Map();
  let acquisitionCount = 0;
  let deletionCount = 0;
  const specialCounts = Object.fromEntries(SPECIAL_EVENT_FIELDS.map((field) => [field, 0]));
  for (const node of reconstruction.nodes) {
    if (
      !object(node) ||
      Object.keys(node).some((key) => !RECONSTRUCTION_NODE_KEYS.has(key)) ||
      Object.keys(node).length !== RECONSTRUCTION_NODE_KEYS.size ||
      typeof node.name !== "string" ||
      !treeNodeNames.has(node.name) ||
      nodesByName.has(node.name)
    )
      fail();
    validateSpacerIds(node.spacers, canonicalIds);
    validateSpacerIds(node.gains, canonicalIds);
    if (!Array.isArray(node.loss_blocks)) fail();
    const flattenedLosses = [];
    for (const block of node.loss_blocks) {
      validateSpacerIds(block, canonicalIds, { allowEmpty: false });
      flattenedLosses.push(...block);
    }
    if (new Set(flattenedLosses).size !== flattenedLosses.length) fail();
    for (const field of SPECIAL_EVENT_FIELDS) {
      validateSpacerIds(node[field], canonicalIds);
      specialCounts[field] += node[field].length;
    }
    acquisitionCount += node.gains.length;
    deletionCount += flattenedLosses.length;
    nodesByName.set(node.name, node);
  }
  if (
    nodesByName.size !== treeNodeNames.size ||
    reconstruction.nodes.length !== treeNodeNames.size ||
    acquisitionCount !== reconstruction.acquisition_count ||
    deletionCount !== reconstruction.deletion_count
  )
    fail();
  const memberBySource = new Map(members.map((member) => [String(member.source_id), member]));
  if (tree.leaves.length !== memberBySource.size) fail();
  const leafSpacers = new Map();
  for (const leaf of tree.leaves) {
    const member = memberBySource.get(leaf);
    const node = nodesByName.get(leaf);
    if (!member || !node || node.spacers.length !== Number(member.spacer_count)) fail();
    leafSpacers.set(
      leaf,
      [...node.spacers].sort((left, right) => left - right),
    );
  }
  return { canonicalIds, leafSpacers, specialCounts, tree };
}

export function validateEvolutionaryReconstructions(
  summary,
  mainGroup,
  comparison,
  selectedReconstruction,
) {
  const orientation = summary.orientation;
  const trees = orientation.trees;
  const reconstructions = orientation.reconstructions;
  if (
    !Array.isArray(trees) ||
    trees.length !== summary.adapter.groups.length ||
    orientation.trees_truncated !== false ||
    !Array.isArray(reconstructions) ||
    reconstructions.length !== summary.adapter.groups.length * 2 ||
    orientation.reconstructions_truncated !== false
  )
    fail();
  const treeRow = trees.find((item) => item?.group === mainGroup.name);
  const groupReconstructions = reconstructions.filter((item) => item?.group === mainGroup.name);
  const input = groupReconstructions.find((item) => item.hypothesis === "input");
  const reverse = groupReconstructions.find((item) => item.hypothesis === "reverse");
  const selectedHypothesis = comparison.prediction === "Reverse" ? "reverse" : "input";
  const expectedSelectedTree =
    selectedHypothesis === "reverse" ? treeRow?.reverse_newick : treeRow?.forward_newick;
  if (
    !object(treeRow) ||
    Object.keys(treeRow).some(
      (key) => !["group", "forward_newick", "reverse_newick", "selected_newick"].includes(key),
    ) ||
    Object.keys(treeRow).length !== 4 ||
    groupReconstructions.length !== 2 ||
    !input ||
    !reverse ||
    input.hypothesis !== "input" ||
    reverse.hypothesis !== "reverse" ||
    treeRow.selected_newick !== expectedSelectedTree
  )
    fail();
  const inputResult = validateCanonicalReconstruction(
    input,
    treeRow.forward_newick,
    mainGroup.name,
    mainGroup.arrays,
  );
  const reverseResult = validateCanonicalReconstruction(
    reverse,
    treeRow.reverse_newick,
    mainGroup.name,
    mainGroup.arrays,
  );
  if (
    inputResult.canonicalIds.size !== reverseResult.canonicalIds.size ||
    [...inputResult.canonicalIds].some((value) => !reverseResult.canonicalIds.has(value))
  )
    fail();
  for (const [leaf, spacers] of inputResult.leafSpacers) {
    if (spacers.join(",") !== reverseResult.leafSpacers.get(leaf)?.join(",")) fail();
  }
  const selected = selectedHypothesis === "reverse" ? reverse : input;
  const selectedResult = selectedHypothesis === "reverse" ? reverseResult : inputResult;
  const serializedTreeTolerance = (selectedResult.tree.nodeNames.length + 1) * 0.0000051;
  if (
    selected.newick !== treeRow.selected_newick ||
    selected.acquisition_count !==
      Number(selectedReconstruction["nb of reconstructed insertions"]) ||
    selected.deletion_count !== Number(selectedReconstruction["nb of reconstructed deletions"]) ||
    selectedResult.canonicalIds.size !== Number(selectedReconstruction["nb of unique spacers"]) ||
    selectedResult.tree.leaves.length !==
      Number(selectedReconstruction["nb of leafs (after combining non-uniques)"]) ||
    !finite(selectedReconstruction["tree length"]) ||
    !finite(selectedReconstruction["tree height"]) ||
    Math.abs(selectedResult.tree.totalLength - Number(selectedReconstruction["tree length"])) >
      serializedTreeTolerance ||
    Math.abs(selectedResult.tree.height - Number(selectedReconstruction["tree height"])) >
      serializedTreeTolerance
  )
    fail();
  const specialMetrics = {
    duplications: "nb of reconstructed duplications",
    rearrangements: "nb of reconstructed rearrangements",
    reacquisitions: "nb of reconstructed reacquisitions",
    independent_gains: "nb of reconstructed independent gains",
    other_duplication_events: "nb of reconstructed other dup. events",
  };
  for (const [eventField, metricField] of Object.entries(specialMetrics)) {
    if (selectedResult.specialCounts[eventField] !== Number(selectedReconstruction[metricField]))
      fail();
  }
}
