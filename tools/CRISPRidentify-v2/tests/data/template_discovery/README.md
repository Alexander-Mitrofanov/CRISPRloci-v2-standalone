# Supplied repeat discovery engineering controls

These small FASTA files are frozen copies of the 2026-10-07 detection council's
synthetic probes. `manifest.json` preserves the audited source-file SHA-256
values, baseline source commit, injection reports and their exact feature
vectors, scores and source intervals. Do not regenerate the fixtures to change
expected truth. The construction used fixed random source flanks and spacers;
the three-copy motif has three polymorphic columns, and the four-copy motif
differs slightly from every observed copy. These are source-derived query
controls, not an independent biological evaluation or a held-out recall estimate.

`three_copy_panel.fa` supplies `AGGACAAAAAAATAGAGAATCTTCCACAGA`.
Its expected zero-based half-open repeats are `[81,111)`, `[143,173)`,
`[205,235)`, with two 32-base spacers. It passes structure but remains Low score
at `0.2027101025001436`, with no accepted call.

`four_copy_panel.fa` supplies `GTTTTTATCGTACCTATGAGGAATTGAAAC`.
Its repeats are `[80,110)`, `[142,172)`, `[204,234)`, `[266,296)`, with three
32-base spacers. Balanced template discovery produces Possible at
`0.7036515235434023`. The balanced no-panel control has no arrays. Sensitive
no-panel discovery already finds this locus with 21-base repeats and score
`0.6204423621476783`; template discovery recovers the full borders.

Additional controls freeze the eleven-base exact-core example, the two-repeat
classification example, identical repeats, and 110-base spacers. Structural
average settings and discovery reachability are tested separately. Expected
report tables are one-based inclusive; the explicit interval lists above and
`expected_half_open_intervals` in the manifest use zero-based half-open notation.
