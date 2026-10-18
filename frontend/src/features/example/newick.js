import { fail } from "./contract.js";

export function parseNewick(value) {
  const source = typeof value === "string" ? value.trim() : "";
  if (!source || source.length > 10_000) fail();
  let position = 0;
  const skip = () => {
    while (position < source.length && /\s/.test(source[position])) position += 1;
  };
  const readLabel = () => {
    skip();
    const start = position;
    while (position < source.length && /[A-Za-z0-9_.+|\-]/.test(source[position])) position += 1;
    if (start === position) fail();
    return source.slice(start, position);
  };
  const readLength = () => {
    skip();
    if (source[position] !== ":") fail();
    position += 1;
    skip();
    const match = source.slice(position).match(/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/);
    if (!match) fail();
    position += match[0].length;
    const length = Number(match[0]);
    if (!Number.isFinite(length) || length < 0) fail();
    return length;
  };
  const readNode = () => {
    skip();
    const children = [];
    if (source[position] === "(") {
      position += 1;
      children.push(readNode());
      while (true) {
        skip();
        if (source[position] !== ",") break;
        position += 1;
        children.push(readNode());
      }
      skip();
      if (source[position] !== ")" || children.length < 2) fail();
      position += 1;
    }
    const name = readLabel();
    return { name, length: readLength(), children };
  };
  const root = readNode();
  skip();
  if (source[position] !== ";") fail();
  position += 1;
  skip();
  if (position !== source.length) fail();
  const nodeNames = [];
  const leaves = [];
  const leafDistances = [];
  let totalLength = 0;
  const visit = (node, distance) => {
    nodeNames.push(node.name);
    totalLength += node.length;
    const nextDistance = distance + node.length;
    if (node.children.length) node.children.forEach((child) => visit(child, nextDistance));
    else {
      leaves.push(node.name);
      leafDistances.push(nextDistance);
    }
  };
  visit(root, 0);
  if (new Set(nodeNames).size !== nodeNames.length || leaves.length < 2) fail();
  return { nodeNames, leaves, totalLength, height: Math.max(...leafDistances) };
}
