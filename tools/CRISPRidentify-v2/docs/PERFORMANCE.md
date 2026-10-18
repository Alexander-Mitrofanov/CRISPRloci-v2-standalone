# Resources and performance

Choose an analysis contract before sizing a run. Presets can change the work
performed and the reported arrays; their runtimes are not interchangeable.

## Choose a preset

| Preset | Work performed | Use |
| --- | --- | --- |
| `balanced` | Fast search, conservative reconstruction and standard reports; optional annotations remain configurable. | General array analysis; the default. |
| `region-screen` | Accepted-only machine output with optional refinement and annotation work disabled. | A first pass for large cohorts. Boundaries, representations and categories may differ from `balanced`. |
| `full` | Historical non-fast search with standard output. | Detailed comparisons that require this search contract. |

`full` and `region-screen` use the legacy reconstruction policy. See the
[README](../README.md#choose-a-preset) for preset compatibility and the
[output contract](OUTPUTS.md) for the scientific meaning of each category.
Cas annotation is optional and adds separate work; measure it explicitly when
it is part of the intended analysis.

## Allocate CPU and memory

For cohorts, use `crispridentify-batch` with a frozen manifest and explicit
`--workers` and `--cpus-per-worker`. Their product must fit the scheduler's
CPU allocation. Each worker processes a sample independently, so increasing
workers also increases concurrent memory and storage demand. Begin with a
small representative shard, including large and repetitive genomes, and size
concurrency from observed aggregate memory use.

The CLI configures native threads before numerical imports: BLAS and NumExpr
use one thread, sequential OpenMP and external work receive the resolved CPU
budget, and detector pool children use one native thread. Both Clustal calls
receive an explicit allocation. This is a thread-allocation policy; use the
scheduler or operating system to enforce CPU and memory limits.

Parallel native matching uses fresh Python workers. Direct Python callers
must use an importable script with an `if __name__ == "__main__":` guard.
For interactive callers, use the CLI or `--parallel False`.

Memory use depends on candidate density as well as sequence length. Streamed
matching and compact retained evidence reduce unnecessary duplication, but
these do not provide a fixed memory bound for every input. Allow headroom and
retain resource failures in the run ledger. A failed sample is not a successful
sample with no arrays.

## Measure the intended workload

Record the source revision, dependency lock, input hashes, preset, resolved
configuration, CPU affinity, memory limit, swap policy and optional outputs.
Compare the same inputs and configurations, reporting failed and completed
samples separately. Keep native analysis, validation, annotation and archiving
costs distinct when they are measured separately. State whether a memory
measurement covers one process or the whole worker process group.

Repeat representative measurements when a throughput or speed claim matters.
Historical runs on different source versions or operating conditions provide
context but cannot isolate the cause of a timing difference. Do not extrapolate
small easy-genome runs to a large cohort without checking its resource tail.

See [cluster execution](CLUSTER.md) for restart and failure handling, and
[CONTRIBUTING.md](../CONTRIBUTING.md) for the performance-change checklist.
Keep campaign data, scientific comparisons and review records outside the tool
checkout; retain their manifests and source identities with the local study.
