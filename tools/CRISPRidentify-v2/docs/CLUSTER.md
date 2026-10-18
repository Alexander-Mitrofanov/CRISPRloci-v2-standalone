# Cluster execution

`crispridentify-batch` is the supported restartable interface for cohorts. It
runs one FASTA in one isolated process group, keeps a bounded number of jobs in
flight, validates every `report.json` against its source FASTA, and publishes a
completion marker only after validation. Its memory use is bounded by worker
count rather than manifest length.

The commands below use the `region-screen` preset. This preset targets recall
of high-confidence biological regions; it does not promise identical candidate
representations or boundaries to `balanced` or historical `full` search.

## 1. Install identically on every VM

Check out one immutable release tag on all four machines and run:

```bash
./install.sh
./crispridentify --version
./crispridentify-batch --version
```

Keep the Git revision, environment lock, and model files identical. Prefer a
read-only installation path and node-local storage for inputs, scratch, and
active results.

## 2. Create the frozen manifest

One row represents one genome or assembly FASTA. Uncompressed `.fna`, `.fa`,
`.fasta`, and `.fas` files are supported by the manifest builder.

```bash
./crispridentify-batch manifest \
  --input-root /data/genomes \
  --output /data/manifests/genomes.tsv
```

The output is a deterministic, sample-sorted TSV with these exact columns:

```text
sample_id    fasta    bytes    sha256
```

Paths are relative to the manifest when possible. Before a fresh attempt, the
runner rejects duplicate sample IDs, malformed hashes, changed file sizes, and
changed source hashes. Keep the manifest immutable after launching the study.

## 3. Run four deterministic shards

Use one shard index per VM. Indices are zero-based. The SHA-256 assignment of a
sample ID is stable and independent of manifest order.

On VM 0:

```bash
CRISPRIDENTIFY_ENV=/opt/CRISPRidentify-v2/.envs/crispridentify-core \
./crispridentify-batch run \
  --manifest /data/manifests/genomes.tsv \
  --run-dir /local/results/study-2026/shard-0000 \
  --scratch-root /local/scratch/crispridentify \
  --shard-count 4 \
  --shard-index 0 \
  --preset region-screen \
  --workers 12 \
  --cpus-per-worker 2 \
  --timeout-seconds 1800 \
  --max-attempts 3
```

Repeat with shard/run-directory pairs `1/shard-0001`, `2/shard-0002`, and
`3/shard-0003`. Never point two live runners at the same run directory. A
nonblocking file lock rejects that mistake.

The runner uses at most twice the worker count in scheduled Futures, so a
three-million-row manifest does not create three million Python scheduling
objects. A disk-backed SQLite index holds the selected shard. Each worker is
constrained through `--cpu` and the common OpenMP/BLAS thread variables.

### Starting point for 4 × 32-vCPU / 64-GiB VMs

Use 12 workers × 2 vCPU on each VM initially:

- 24 of 32 vCPU allocated to detector workers;
- 8 vCPU left for the operating system, storage, decompression, and tail
  latency;
- 48 simultaneous genomes and 96 allocated vCPU cluster-wide;
- node-local scratch with at least 100 GiB free is recommended;
- do not enable swap as a substitute for reducing concurrency.

Three million genomes in 14 days requires at least 8,929 completed genomes per
hour cluster-wide (2,233 per VM). Before releasing all shards, run a
10,000–50,000-genome population-representative pilot on the actual VMs. Require
at least 11,200 genomes/hour cluster-wide, zero OOM events, no sustained swap,
stable p95/p99 runtime, and sufficient scratch/output headroom. If RSS or I/O
pressure is high, reduce workers before changing CPUs per worker.

The earlier Diverse-1000 measurement was not performed on these VMs. It is
capacity evidence, not a guarantee for different genome-size distributions,
storage, hypervisors, or Cas annotation.

## Resume and failure contract

Re-run the identical command after interruption. A sample is skipped only when
all of the following still validate:

- manifest-recorded sample ID, byte count, and source-file SHA-256 agree with
  the immutable completion record;
- exact batch configuration SHA-256 and runtime-tree SHA-256;
- every source represented exactly once in the detector output;
- every report schema, source length, source sequence hash, summary count, and
  report-file SHA-256;
- accepted-only categories for `region-screen`.

Each immutable attempt contains `attempt.json`, `stdout.log`, `stderr.log`, and
its isolated result tree. An invalid old completion marker is moved under
`superseded/`; it is never trusted or silently overwritten. Attempts are
bounded by `--max-attempts`. Child processes receive SIGTERM as a process group
at timeout, followed by SIGKILL after the configured grace interval.

Run-level files are:

- `config.json`: effective preset/expansion, arguments, runtime tree hash,
  manifest identity, shard, resource limits, and config hash;
- `state.sqlite`: disk-backed selected-shard state;
- `status.tsv`: deterministic status ledger in manifest order;
- `failures.jsonl`: explicit failed-sample ledger;
- `summary.json`: counts, outcome, and ledger hashes;
- `samples/*/completion.json`: atomic per-sample completion markers.

Exit status is zero only if every selected sample succeeded or was safely
resumed. Any failed, pending, or interrupted sample produces a nonzero exit.
`summary.json` reports elapsed wall time, completion fraction, and
`completed_genomes_per_hour`. The throughput denominator counts only genomes
freshly succeeded in that invocation; safely skipped resume completions are
excluded so a no-op resume cannot inflate the capacity estimate.
Archive the configuration, manifest, ledgers, and completion trees with the
study; do not archive `state.sqlite` as the sole provenance artifact.

Use `--resume-verify-input` for an additional byte-level source recheck on every
resumed sample. It is intentionally optional because rereading millions of
frozen FASTA files can create substantial storage traffic.

Without `--resume-verify-input`, resume trusts the frozen manifest and validated
completion/report hashes and does not read the current FASTA bytes again. This
is appropriate only when the manifest and input store are managed as immutable
study artifacts.

## Scheduler integration

The runner has no dependency on Slurm, PBS, SGE, Kubernetes, or a shared queue.
A site scheduler should launch exactly one command per deterministic shard and
request 32 vCPU plus 64 GiB for each command. Use scheduler-native log capture
and signal forwarding. Do not add another per-genome job array around the
runner: its own worker pool already provides bounded concurrency and retries.

For more than four nodes, increase `--shard-count` before starting the frozen
study and assign each index exactly once. Do not change shard count inside an
existing run directory.

## CasAndra at cohort scale

`region-screen` deliberately disables CasAndra so array-discovery throughput
and Cas annotation are not conflated. For a separately budgeted Cas run, use a
new run directory with `--preset balanced --cas`; never mix it into an existing
region-screen directory:

```bash
./crispridentify-batch run \
  --manifest /data/manifests/pilot.tsv \
  --run-dir /local/results/cas-pilot/shard-0000 \
  --scratch-root /local/scratch/crispridentify-cas \
  --preset balanced \
  --cas \
  --cas-gene-mode auto \
  --cas-timeout-seconds 3600 \
  --workers 4 \
  --cpus-per-worker 6
```

Benchmark CasAndra independently on a representative pilot before selecting
production concurrency. Its gene calling, 341-profile search, and cassette
classification have a different CPU/RAM distribution from array detection.

## Operational cautions

- Use one immutable FASTA per manifest row; duplicate source IDs inside a
  multi-FASTA are rejected.
- Keep active attempts and temporary files on node-local storage. Copy only
  completed artifacts to shared/object storage after checksum verification.
- Set `--min-free-gib` higher than the 10-GiB safety default for production.
- Do not use `--allow-oversubscribe` under a scheduler unless CPU affinity is
  known to be reported incorrectly.
- Extra detector options are accepted only as repeated `--tool-arg=...` tokens
  and become part of the configuration hash. Input, result, CPU, folder-job,
  and preset options cannot be overridden this way.
