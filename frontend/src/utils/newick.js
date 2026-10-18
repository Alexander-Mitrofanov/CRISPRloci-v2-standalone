export function parseNewick(value) {
  const source = String(value || "").trim();
  if (!source || source.length > 10_000) return null;
  let position = 0;
  let nextId = 0;
  const skip = () => {
    while (position < source.length && /\s/.test(source[position])) position += 1;
  };
  const readLabel = (required) => {
    skip();
    const start = position;
    while (position < source.length && /[A-Za-z0-9_.+|\-]/.test(source[position])) position += 1;
    if (required && start === position) throw new Error("missing Newick label");
    return source.slice(start, position);
  };
  const readLength = () => {
    skip();
    if (source[position] !== ":") return 0;
    position += 1;
    skip();
    const start = position;
    while (position < source.length && /[0-9eE+\-.]/.test(source[position])) position += 1;
    const number = Number(source.slice(start, position));
    if (!Number.isFinite(number) || number < 0) throw new Error("invalid Newick length");
    return number;
  };
  const readNode = () => {
    skip();
    const node = { id: (nextId += 1), name: "", length: 0, children: [] };
    if (source[position] === "(") {
      position += 1;
      node.children.push(readNode());
      while (true) {
        skip();
        if (source[position] !== ",") break;
        position += 1;
        node.children.push(readNode());
      }
      skip();
      if (source[position] !== ")" || node.children.length < 2)
        throw new Error("invalid Newick branch");
      position += 1;
      node.name = readLabel(false);
    } else node.name = readLabel(true);
    node.length = readLength();
    return node;
  };
  try {
    const root = readNode();
    skip();
    if (source[position] !== ";") return null;
    position += 1;
    skip();
    return position === source.length ? root : null;
  } catch {
    return null;
  }
}

export function layoutNewick(
  value,
  height = 320,
  scaleMode = "branch",
  sharedDistance = null,
  treeWidth = 430,
) {
  const tree = parseNewick(value);
  if (!tree) return null;
  const nodes = [];
  const leaves = [];
  const visit = (node, depth = 0, distance = 0, parent = null) => {
    Object.assign(node, { depth, distance, parent });
    nodes.push(node);
    if (node.children.length)
      node.children.forEach((child) => visit(child, depth + 1, distance + child.length, node));
    else leaves.push(node);
  };
  visit(tree);
  const maxDepth = Math.max(1, ...nodes.map((node) => node.depth));
  const maxDistance = Math.max(0, ...nodes.map((node) => node.distance));
  const top = 34;
  const bottom = height - 34;
  leaves.forEach((leaf, index) => {
    leaf.y =
      leaves.length === 1
        ? height / 2
        : top + (index / Math.max(1, leaves.length - 1)) * (bottom - top);
  });
  const place = (node) => {
    if (!node.children.length) return node.y;
    const positions = node.children.map(place);
    node.y = positions.reduce((sum, item) => sum + item, 0) / positions.length;
    return node.y;
  };
  place(tree);
  const distanceExtent = Number(sharedDistance) > 0 ? Number(sharedDistance) : maxDistance;
  nodes.forEach((node) => {
    const measure =
      scaleMode === "branch" && distanceExtent > 0
        ? node.distance / distanceExtent
        : node.depth / maxDepth;
    node.x = 28 + measure * treeWidth;
  });
  return { tree, nodes, leaves, height, maxDistance };
}
