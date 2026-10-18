# Independent biological evaluation

This directory is maintained by worker W3, independently of predictor development
and tuning. Software fixtures are not biological accuracy evidence. No final
predictions have been inspected at the initial split freeze.

`partition.json` freezes the first acquisition panel before model fitting.
Every accession version, locus geometry and label still needs verification
against the source record and primary study. An accession nomination is not
evidence that every array in that strain belongs to the desired family.

`GATES.md` distinguishes release-blocking invariants, quantitative outcomes and
claims that cannot be supported by this panel. `curation.json` will contain
factual locus and assay evidence, independent source hashes and label scope.
Final labels must not be imported by model building or threshold calibration.

This is a prerelease evaluation, not a population-representative accuracy study.
Experimental adaptation labels support strand/junction interpretation only;
they do not by themselves label conserved-core endpoints or full leaders.
Natural alternative flanks are uncertain controls. Sequence shuffles measure
false calls under a declared null and are not proof of biological absence.

## Current status: incomplete

Automated safety restrictions stopped biological research in two worker turns.
The coordinator halted further acquisition, reference construction and final
biological evaluation; it was not rerouted through another worker or source.
The owner's implementation approval remains recorded, but these research tasks
were not completed. Unaffected engineering work continued.

The source audit completed before the stop showed that the initial ED1a
nomination did not establish an I-E leader label. The coordinator therefore
approved an amendment **before any model fitting or biological predictions**:
K-12 training, Sakai calibration, BL21 candidate known-family final and ED1a
candidate novel-family final. The previous assignments and reason are retained
in `partition.json`. BL21(DE3) genomic material is not automatically equivalent
to the BL21-AI material used in the nominated primary study; exact locus
correspondence still requires review.

No curated biological fixture, default biological reference bundle or completed
benchmark report is present here. No scientific gate has been evaluated or
claimed passed. This directory is a pre-output evaluation protocol and an
honest status record; it cannot establish scientific release readiness.

## Comparative-science extension

The owner's subsequent request to establish scientific value and superior
performance is addressed in [COMPARATORS.md](COMPARATORS.md) and
[SUPERIORITY_PROTOCOL.md](SUPERIORITY_PROTOCOL.md). Machine-readable candidate
inventory and unresolved study choices are recorded in `comparators.json` and
`study-design.json`. This does not transform the small source-audit panel into
a confirmatory benchmark.

The generic `benchmark` command evaluates supplied transcription-orientation
labels/calls with synthetic software fixtures. It does not construct biological
references, run a leader predictor or execute competitors. Its intervals are
descriptive and unadjusted for multiple comparisons; it cannot certify superiority.
