# CRISPRmap-v2

This is an independent 2026 prokaryotic repeat reference mapper, distinct from
historical CRISPRmap v2 (2014). Read docs/COUNCIL.md and docs/SCIENCE.md before
changing scientific behavior.

- Historical families, structural motifs, Cas subtypes and transcription strand
  are separate concepts. Similarity is not calibrated confidence.
- Preserve all equal-best label conflicts before truncating displayed hits.
- Missing annotations remain missing. No host or Cas subtype prediction.
- Runtime analysis is offline. Reference import is an explicit operator action.
- Use bounded input, work and output; fail without publishing partial success.
- Never commit data/, output/, .ops-private/, environments or credentials.
- Install .githooks in real checkouts and run scripts/check.sh before publication.
- Benchmarks must identify original published output comparisons separately
  from new held-out annotation agreement and engineering controls.
