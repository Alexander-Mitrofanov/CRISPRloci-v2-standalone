# evOr code map

## Execution path

`crispr_evor.cli` → `crispr_evor.main` forces the orientation flag →
`spacerplacer.main` validates CLI paths and dispatches inputs →
`sp_model.run_experiments` aligns/loads groups and reconstructs forward then reverse →
`sp_model.orientation.orientation_tools` compares BDM log likelihoods →
`crispr_evor.write_orientation_json` checks CSV/tree/correction consistency and atomically writes the public result.

| Module | Responsibility |
| --- | --- |
| `crispr_evor.py` | Dedicated CLI and validated `crispr_evor_results.json` contract. |
| `spacerplacer.py` | Legacy CLI options, trusted-pickle gate, input/tree dispatch. |
| `sp_model/input_parser.py`, `crisprfinder_to_spacerfasta.py` | CRISPRCasDB/CRISPRCasFinder adapters and sequence/array grouping. |
| `sp_model/experiments_tools.py` | Spacer FASTA, alignment cache identity, tree validation, pooled estimates. |
| `sp_model/mafft_integration.py` | Spacer encoding, MAFFT subprocesses, ordered alignment conversion. |
| `sp_model/run_experiments.py` | Group orchestration, shared orientation exports and ODE alias resolution; legacy pickle/simulation workflows. |
| `sp_model/run_reconstruction.py` | Per-group fitted-tree provenance, ancestral/event reconstruction, model fits and bias policy. |
| `sp_model/reconstruction_tree.py` | Ancestral states, classified gains/losses, deletion likelihood fitting and event diagnostics. |
| `sp_model/model_tools.py`, `evolution.py` | Likelihood arithmetic, optimization checks, alignment/order transforms, gain/loss transition model. |
| `sp_model/distance_tree_constructor.py`, `data_classes/advanced_tree.py` | Pairwise evolutionary distances, NJ/UPGMA and tree helpers. |
| `sp_model/data_classes/crisprdata.py` | Legacy array/group records; public names also occur in trusted pickles. |
| `sp_model/helpers/sequence.py` | IUPAC complements and deterministic edit-distance/overlap clustering. |
| `sp_model/run_simulation.py`, `simulation_tree.py`, `coalescent_tree_sim.py` | Scientific simulation workflows, separate from normal orientation jobs. |
| `sp_model/spacer_visualization/`, `summary.py` | Optional plots/event exports and human-readable summaries. |
| `sp_model/bdm_likelihood_computation/` | ODE/symbolic likelihood generation and version-sensitive bundled callable cache. |

## Scientific invariants

The signed comparison is forward minus reverse BDM log likelihood. Reverse wins strictly below
`-confidence_threshold`; forward wins strictly above it; equality and all interior values remain `ND`.
`ND` retains the forward fit/tree. Comparison differences and ratios continue to describe input order versus
reverse even when the selected protocol uses reverse fit metrics. Trivial-group filtering uses the selected
orientation's classification. Preserve these distinctions when changing dataframe operations.

The primary orientation BDM likelihood does not change when optional rho correction falls back; rate estimates
change, while optimizer likelihoods and IDM/BDM tests remain unchanged. The five bias provenance columns must remain attached to their corresponding fit.

`get_rec_duplication_type_counts` reports five classified duplicate-event totals followed by repeated gains
of spacers outside the duplicate-candidate set. That last value feeds the diagnostic CSV column
`nb of reconstructed ind. acquisitions, not duplicates` and summaries, after reconstruction and fitting.
It does not feed likelihood fitting or the orientation decision.

## Tests and verification

- `test_orientation_and_numerics.py`: fit selection, threshold equality, underflow, mixed/trivial groups,
  likelihood recurrence equivalence, optimization and rho correction policy.
- `test_orientation_input_identity.py`: cache/input identity, dotted group names and actual fitted trees.
- `test_crispr_evor_cli.py`: JSON schema, correction provenance, scientific consistency, atomic result export.
- `test_pickled_orientation.py`: explicitly authorized legacy pickle CLI, supplied/estimated trees and common/legacy exports.
- `test_reconstructed_event_export.py`: event export without plotting and independent-gain diagnostic accounting.
- Other suites cover MAFFT/CLI safety, adapters, sequence helpers, simulation parameters and packaging.

Run the suites in evOr's own environment. A no-plot smoke command is:

```bash
python crispr_evor.py example_datasets/tomida_datasets/CRISPR1/G1.fa /tmp/evor-smoke \
  --input_type spacer_fasta --no_plot_reconstruction --no_plot_order_graph
```

For refactors compare all fitted numeric columns (excluding runtime measurements), orientation JSON and
forward/reverse/selected tree maps against a frozen pre-change run. Preserve upstream GPL authorship and
historical files under `provenance/`; record new validation separately.
