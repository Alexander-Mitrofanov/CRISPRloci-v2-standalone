"""Generated Newick grammar and topology validation."""

from __future__ import annotations

import math
import re

from .validation_common import _nonempty_string


_NEWICK_LABEL = re.compile(r"[A-Za-z0-9_.+|-]+")
_NEWICK_NUMBER = re.compile(
    r"[+-]?(?:[0-9]+(?:[.][0-9]*)?|[.][0-9]+)(?:[eE][+-]?[0-9]+)?"
)


def _generated_newick_details(
    value: object,
) -> tuple[tuple[str, ...], float, int] | None:
    """Parse generated Newick and return leaves plus its explicit branch sum."""

    if not _nonempty_string(value):
        return None
    tree = str(value).strip()
    position = 0
    leaf_names: list[str] = []
    branch_total = 0.0
    branch_count = 0

    def skip_whitespace() -> None:
        nonlocal position
        while position < len(tree) and tree[position].isspace():
            position += 1

    def parse_label(*, required: bool, leaf: bool = False) -> bool:
        nonlocal position
        skip_whitespace()
        match = _NEWICK_LABEL.match(tree, position)
        if match is None:
            return not required
        position = match.end()
        if leaf:
            leaf_names.append(match.group(0))
        return True

    def parse_branch_length() -> tuple[bool, bool]:
        nonlocal position, branch_count, branch_total
        skip_whitespace()
        if position >= len(tree) or tree[position] != ":":
            return True, False
        position += 1
        skip_whitespace()
        match = _NEWICK_NUMBER.match(tree, position)
        if match is None:
            return False, False
        try:
            branch_length = float(match.group(0))
        except ValueError:
            return False, False
        if not math.isfinite(branch_length) or branch_length < 0:
            return False, False
        position = match.end()
        branch_count += 1
        branch_total += branch_length
        return True, True

    def parse_subtree() -> tuple[int, bool] | None:
        nonlocal position
        skip_whitespace()
        if position >= len(tree):
            return None
        if tree[position] == "(":
            position += 1
            first_child = parse_subtree()
            if first_child is None or not first_child[1]:
                return None
            leaves = first_child[0]
            child_count = 1
            while True:
                skip_whitespace()
                if position < len(tree) and tree[position] == ",":
                    position += 1
                    nested = parse_subtree()
                    if nested is None or not nested[1]:
                        return None
                    leaves += nested[0]
                    child_count += 1
                    continue
                break
            if child_count < 2 or position >= len(tree) or tree[position] != ")":
                return None
            position += 1
            if not parse_label(required=False):
                return None
        else:
            if not parse_label(required=True, leaf=True):
                return None
            leaves = 1
        valid_length, has_length = parse_branch_length()
        return (leaves, has_length) if valid_length else None

    root = parse_subtree()
    skip_whitespace()
    if root is None or root[0] < 2 or root[0] != len(leaf_names):
        return None
    if position >= len(tree) or tree[position] != ";":
        return None
    position += 1
    skip_whitespace()
    if position != len(tree) or branch_count < len(leaf_names):
        return None
    return tuple(leaf_names), branch_total, branch_count


def _generated_newick_leaves(value: object) -> tuple[str, ...] | None:
    details = _generated_newick_details(value)
    return None if details is None else details[0]


def _valid_generated_newick(value: object) -> bool:
    return _generated_newick_details(value) is not None


def _tree_has_exact_leaves(value: object, expected: frozenset[str]) -> bool:
    details = _generated_newick_details(value)
    if details is None:
        return False
    leaves = details[0]
    return (
        len(leaves) == len(expected)
        and len(set(leaves)) == len(leaves)
        and set(leaves) == expected
    )


def _tree_matches_result(
    value: object, expected_leaves: frozenset[str], expected_length: float
) -> bool:
    details = _generated_newick_details(value)
    if details is None:
        return False
    leaves, branch_total, branch_count = details
    tolerance = branch_count * 5e-6 + 1e-9
    return (
        len(leaves) == len(expected_leaves)
        and len(set(leaves)) == len(leaves)
        and set(leaves) == expected_leaves
        and math.isclose(
            branch_total,
            expected_length,
            rel_tol=1e-9,
            abs_tol=tolerance,
        )
    )


