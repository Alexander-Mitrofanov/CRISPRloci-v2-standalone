# Native-amendment clarification before build/results

The engineering reviewer specifies that repeatability perturbations change process
launch order, output paths and wall-clock start time while preserving exact ordered
FASTA bytes and the declared single-thread setting. Reordering FASTA records or
changing thread count creates a different execution/input and is not a promised
invariance. This clarifies the amendment's “path/order perturbations” wording.
