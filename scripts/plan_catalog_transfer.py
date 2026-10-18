#!/usr/bin/env python3
"""Turn reconciled source inventory into a reproducible compact transfer plan."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path


def plan_transfer(inventory):
    nodes = defaultdict(list)
    raw = defaultdict(list)
    shard_id = 0
    selected_counts = Counter()
    for source in sorted(inventory["sources"], key=lambda s: (int(s["node"]), s["precedence"], s["path"])):
        node = str(source["node"])
        selected_counts.update(entry["status"] for entry in source["entries"])
        if source["kind"] != "bundle":
            raw[node].append(source)
            continue
        shard_id += 1
        nodes[node].append({"id": shard_id, "kind": "bundle", "path": source["path"], "sha256": source["sha256"], "selected": [entry["accession"] for entry in source["entries"]]})
    for node, sources in sorted(raw.items()):
        for start in range(0, len(sources), 50):
            shard_id += 1
            nodes[node].append({"id": shard_id, "kind": "raw_group", "items": sources[start:start + 50]})
    return {"schema": "crisprloci.transfer-plan.v1", "profile_sha256": inventory["profile_sha256"], "counts": dict(selected_counts), "nodes": dict(nodes)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with gzip.open(args.inventory, "rt") as stream:
        inventory = json.load(stream)
    plan = plan_transfer(inventory)
    with args.inventory.open("rb") as stream:
        plan["inventory_sha256"] = hashlib.file_digest(stream, "sha256").hexdigest()
    with gzip.open(args.output, "wt") as stream:
        json.dump(plan, stream, separators=(",", ":"))
    print(json.dumps({"counts": plan["counts"], "shards": {node: len(jobs) for node, jobs in plan["nodes"].items()}}))


if __name__ == "__main__":
    main()
