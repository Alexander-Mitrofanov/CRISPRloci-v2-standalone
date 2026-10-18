# CRISPR-evOr

CRISPR-evOr predicts the acquisition orientation of a group of CRISPR arrays by running the SpacerPlacer evolutionary
reconstruction in both directions and comparing their block-deletion-model log likelihoods. This directory is an
independent, hardened distribution of the orientation workflow published by the SpacerPlacer authors.

The method is described in Fehrenbach et al., *An evolutionary approach to predict the orientation of CRISPR arrays*,
PLOS Computational Biology (2025), DOI `10.1371/journal.pcbi.1013706`. CRISPR-evOr is most useful for groups of related
arrays; it does not infer transcription direction, leader location, or PAM orientation.

## Installation

This standalone repository contains the fixed server distribution **1.0.4**.
Supported server target: Linux x86-64, Python 3.10–3.12. MAFFT must be on `PATH`;
Graphviz is needed for plots. Install in a dedicated environment because
SpacerPlacer also provides a Python package named `sp_model`.

```bash
git clone https://github.com/Alexander-Mitrofanov/evOr.git
cd evOr
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -c constraints-tested.txt .
crispr-evor --help
```

For a bundled example:

```bash
crispr-evor example_datasets/tomida_datasets/CRISPR1/G1.fa runs/G1 \
  --input_type spacer_fasta \
  --no_plot_reconstruction --no_plot_order_graph
```

Run regression tests separately from SpacerPlacer:

```bash
python -m pip install '.[test]'
python -m pytest -q
```

See [import provenance](PROVENANCE.md). The original authorship and
[GPL-3.0 license](LICENSE) are retained.

## Run

```bash
crispr-evor INPUT OUTPUT \
  --input_type ccdb \
  --tree_path TREES.json \
  --cluster_spacers \
  --orientation_decision_boundary 5 \
  --no_plot_reconstruction \
  --no_plot_order_graph
```

For inputs containing multiple independent groups, add `--alignment_workers N` to align up to `N` groups
concurrently. `--aligned_groups_path FILE` can safely reuse an exact validated alignment cache in a trusted standalone
workflow. The integrated orientation workflow deliberately does not run a provisional input-order SpacerPlacer stage
or forward that stage's estimated tree/cache: without an independent phylogeny, CRISPR-evOr must estimate the forward
and reverse trees separately.

Alignment caches must match the input groups, array names, and ordered spacer IDs after removing alignment gaps.
Group names retain dots in filenames; duplicate group stems are rejected. `-` is reserved for alignment gaps and
cannot be an observed input spacer ID. Supplied trees must have non-negative finite branch lengths and unique node
names. A missing root stem is treated as zero; missing non-root lengths require an explicit positive
`--extend_branches`. Exported tree maps record the actual trees used after internal-node naming and branch extension.

```bash
crispr-evor INPUT OUTPUT ... --alignment_workers 8
```

The dedicated command always enables the forward/reverse comparison. In addition to the detailed SpacerPlacer CSVs,
it writes `crispr_evor_results.json` with one unambiguous `Forward`, `Reverse`, or `ND` result per group. `Forward`
means the input order is preferred; `Reverse` means its reverse is preferred; `ND` means the absolute BDM
log-likelihood difference does not exceed the confidence threshold. `additional_data/orientation_trees.json` preserves
the forward, reverse, and selected trees so every likelihood comparison remains auditable.

At least two related arrays are required per group. A reliable provided phylogeny is preferable. If none is supplied,
the method follows the paper and estimates a separate tree for each orientation; the paper reports lower simulated
accuracy for estimated trees than for true/provided trees.

The alpha and rho bias corrections improve deletion-parameter estimates and remain requested by default. On some
small or sparse groups, the optional unobserved-spacer rho correction has no finite optimum even though the primary IDM
and BDM fits are valid. Strict behavior remains the worker default and fails clearly. A production caller may opt into
`--rho_bias_failure_policy fallback_uncorrected`: each group and orientation first attempts the correction, and only an
affected fit is retried without rho correction. This does not alter the BDM log likelihood used for the orientation
decision; the affected rho/deletion-rate estimates remain uncorrected. The optimizer likelihoods and IDM/BDM likelihood-ratio test are unchanged by rho correction. The warning and effective forward/reverse status are recorded in `crispr_evor_results.json` and the
detailed CSV, as well as emitted to the worker log. `--no_rho_bias_correction` remains the explicit global opt-out.

## Web-server safety

Run every job in a separate subprocess with a new empty output directory and an argument allowlist. Public jobs should
normally use both no-plot flags. Do not expose `--allow_unsafe_pickle`; uploaded pickle files can execute arbitrary
code. Enforce upload/decompressed size, group/array/spacer counts, CPU, memory, concurrency, and a whole-process timeout
in the service layer. Each MAFFT subprocess uses one thread, while `--alignment_workers` controls how many independent
group alignments run concurrently. The MAFFT subprocess timeout is 300 seconds, but reconstruction needs
an outer deadline too.

Simple likelihood is the recommended production mode. ODE mode is supported, but on Python versions incompatible with
the bundled Dill bytecode it regenerates callables from the portable SymPy matrix, which has a substantial cold-start
cost. Exact tested direct dependency versions are in `constraints-tested.txt`.

The underlying SpacerPlacer documentation is retained in [README_SPACERPLACER.md](README_SPACERPLACER.md).
The original server audit records are outside this standalone repository; source provenance is recorded in
[PROVENANCE.md](PROVENANCE.md).
