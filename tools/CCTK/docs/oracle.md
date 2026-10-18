# Reproducing the original scientific oracle

Frozen fixtures were generated before refactoring from upstream 1.0.3 commit
`fe94fde10a2855cbcc3bf118a253ec097d3a5e46` on Python 3.12. Use the exact
`requirements-test.lock` environment (NumPy 2.2.6, Matplotlib 3.10.8,
DendroPy 5.0.8, Pillow 12.3.0). Download or check out that exact upstream revision
as an operator; the regeneration script never downloads anything:

```sh
python scripts/regenerate_oracle.py /path/to/upstream /tmp/new-oracle
```

The script first verifies every upstream CLI and module hash against
`UPSTREAM.json`, then regenerates the pairwise JSON, network, PNG, tree,
constrained tree and simulation outputs. Pair generation uses seed 193, colors
and diffplot use seed 2, tree and simulation use seed 3. It writes interpreter,
package and output hashes in `receipt.json`. Compare the named fixtures with
`tests/fixtures/`; `crisprtree.stdout` must contain
`((B:0,C:11):1,A:52):0;`, and `constrain.stdout` must contain
`((genA.A:5,genB.B:5)Anc_a:5,genC.C:10)Anc_b:0;`.

The modern network intentionally includes singleton clusters that the original
omitted, and applies the threshold before connected components. Native diffplot
pixel parity uses the same pinned plotting stack; the new observed JSON plot
reserves additional viewport space for labels and deletion-only slots and is
checked separately. The original algorithms can have scientific limitations;
parity demonstrates behavior preservation, not biological ground truth.
