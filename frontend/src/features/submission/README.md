# Submission feature

`index.js` is the stable component-facing boundary. `useAnalysisForm.js` owns
form state, FASTA inspection, bounded example loading, submission transport,
and credential normalization. The SFCs in `components/submission/` only render
that state and emit user intent.

Keep backend calls and example fetches in this feature composable. Add pure
payload rules to `submission.js` or `fasta.js`, then cover the narrowest owning
module with a focused test.

The application starts with `ModeSelector` links. `useAnalysisNavigation` owns
query-based method pages and browser history; it preserves the recovery hash.
`AnalysisForm` receives a fixed `initialMode`, and its cached draft belongs to
that method. Methods with shipped fixtures offer an input-bound precomputed example; RNA repeats
and Type V tracrRNA have separate examples selected by their input options.
Loading never switches the chosen method. All submitted scientific options must
match before viewing a stored result; edits use the normal live submission path.
The five context/family additions require explicit service advertisement and do
not offer example links until their independently verified assets are shipped.
Locus comparison accepts 2–8 neighborhoods of at most 50,000 nt. RNA-family
analysis accepts 2–50 oriented DNA/RNA sequences of 30–500 nt, with at least two
distinct unambiguous sequences; its sequence JSON port is `rna`.

An explicit result handoff supplies `initialInput` to a fresh form instance. It
preserves the exact JSON and molecule, resets analysis options, and keeps the
previous completed job accessible until submission succeeds. `App` bounds the
form cache and only relaxes the open-job lock for that completed source job.
