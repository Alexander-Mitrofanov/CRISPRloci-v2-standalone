# Building and using a modern reference

Run data-intensive commands on the machines that hold the frozen source data.
The build never obtains credentials, connects to machines, updates runtimes or
launches source scientific tools. Operator orchestration is separate.

Install in a separate environment:

```bash
python -m pip install '.[scale,structure,test]'
```

Create a private JSON source plan with the pinned scientific profile and selected
sources. Each source has `kind`, a local `path`, `sha256`, and explicit `selected`
accession versions. Raw JSON inputs additionally require `entries` containing
accession, terminal status and document checksum. Never commit the source plan
when it contains infrastructure paths. A minimal bundle-only plan is:

```json
{
  "profile_sha256": "<64-character scientific profile digest>",
  "sources": [{
    "kind": "bundle",
    "path": "<local immutable source archive>",
    "sha256": "<64-character archive digest>",
    "selected": ["GCA_000000000.1"]
  }]
}
```

Distribute disjoint source selections into independent plans. Each process owns
one output directory. No database is shared between concurrent writers.

```bash
crisprmap catalog-build --plan private-plan.json.gz \
  --profile-sha256 <digest> --output data/part-01
crisprmap catalog-merge --parts data/part-01 data/part-02 --output data/catalog
crisprmap catalog-inspect --catalog data/catalog
crisprmap catalog-index --catalog data/catalog --output data/index
```

The importer commits a source atomically. A failed build can resume with the
identical plan and profile. It reserves free disk space and validates native
record counts. `manifest.json` is published last and is the completion marker.
Do not use `.building.sqlite3` as a reference. Merging rejects overlapping
source identities or accessions. Keep independent copies of source evidence.

Build operational sequence clusters in one or more lanes, then assemble every
lane. The examples use two lanes; larger deployments can distribute independent
lanes across machines with identical catalogue and index hashes.

```bash
crisprmap catalog-cluster-forest --catalog data/catalog --index data/index \
  --lane 0 --lanes 2 --checkpoint-queries 16 --output data/forest-0
crisprmap catalog-cluster-forest --catalog data/catalog --index data/index \
  --lane 1 --lanes 2 --checkpoint-queries 16 --output data/forest-1
crisprmap catalog-cluster-assemble --catalog data/catalog \
  --forests data/forest-0 data/forest-1 --output data/clusters
crisprmap catalog-fold --catalog data/catalog --checkpoint-nodes 100 --output data/structure
crisprmap catalog-search --input repeats.fasta --catalog data/catalog \
  --index data/index --clusters data/clusters --output runs/example
```

Forest workers commit every 16 completed queries in their lane by default;
`--checkpoint-queries` accepts 1–1024. Folding workers commit every 100 newly
folded nodes by default; `--checkpoint-nodes` accepts 1–10000. The final partial
batch is committed before publication. A clean interruption rolls back the
entire uncommitted batch, and resume recomputes it from the previous checkpoint.
Hard termination, such as `SIGKILL` or a machine crash, can leave a hot SQLite
rollback journal. The read-only ownership check cannot recover that journal;
operator recovery is required before resuming, and automatic resume is not
guaranteed in this case.
Both orientations of a folded node and each forest query remain atomic.
Smaller batches reduce repeated work after interruption but require more disk
syncs. Checkpoint sizes may change on resume, including for existing per-node
or per-query checkpoints. Completion manifests record the current invocation's
checkpoint setting separately from scientific source and policy bindings.

Search writes a completion manifest, streamed JSONL and SQLite results. SQLite
retains all best matches even when displayed hits are capped. Requesting cluster
evidence also maps eligible noncore/new queries independently against core
clusters; this does not change cluster membership.

```bash
crisprmap catalog-benchmark --catalog data/catalog --index data/index \
  --legacy-reference data/legacy/2013/reference.json --output runs/benchmark
PYTHON=/path/to/isolated/python bash scripts/check.sh
```

Use the historical `classify` command for original family and motif annotations.
New catalogue sequence groups have a separate vocabulary. See
[scientific policies](MODERN_REFERENCE.md) for interpretation and benchmark limits.
