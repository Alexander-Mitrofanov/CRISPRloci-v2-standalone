# Scientific-value follow-up — 27 September 2026

The owner requested scientific value and performance beyond comparable tools.
That is an objective to test; it is not established by the current implementation.
The literature review and independent methodology review found that the small
prerelease gate cannot support a superiority claim.

The standalone repository now carries:

- `benchmarks/COMPARATORS.md` and `comparators.json`: a task-matched candidate
  inventory covering original CRISPRleader, orientation methods, deployment
  pipelines, motif annotation and evOr's separate acquisition-orientation task.
- `benchmarks/SUPERIORITY_PROTOCOL.md` and `study-design.json`: a stronger
  comparison design requiring independent labels, homology-aware grouping,
  paired uncertainty, complete error/abstention/failure accounting, practical
  margins, power and multiplicity control. Unresolved fields remain null.
- A generic orientation-reporting command and synthetic engineering tests.
  It evaluates supplied labels/calls; it does not construct scientific labels,
  run competitors or establish a performance improvement.

The proposed primary scientific endpoint is conserved-core localization.
Transcription orientation is secondary. Neither v2 leader-side evidence nor
evOr acquisition orientation is automatically transcription-direction truth.
Conserved cores, full leaders, motifs and experimental functional intervals
remain different annotation targets. See the linked primary papers in the
repository's comparator audit and protocol.

The reference bundle, biological predictor, curated independent truth and actual
competitor runs remain unfinished after the earlier research restrictions.
The targeted methodology/tool-comparison literature review did not resume
biological reference acquisition or route blocked work through another worker.
No numerical biological performance result or universal superiority is claimed.

Existing Identify, CasAndra, CRISPRtracrRNA, evOr and SpacerPlacer responsibilities
remain unchanged. Server integration and deployment are outside this work.

Implementation snapshot pushed to the private repository:
`6ed81eb9026b4739c78216102d582d4495b96b53`.
Local checks: 293 software tests, lint, type checks, sdist/wheel build and an
installed-wheel offline CLI test passed. The checkout matched `origin/main` and
was clean after the push. These remain engineering results, not biological accuracy.
Remote engineering CI passed for this snapshot:
https://github.com/Alexander-Mitrofanov/CRISPRleader-v2/actions/runs/36329895889
