"""Independent contract regressions from the frozen module council evidence."""
import copy
import hashlib
import json
from types import SimpleNamespace

import pytest

from crispr_integration.repeat_validation import _pairs as adapter_pairs
from crispr_web.execution.validation_repeats import _pairs as backend_pairs
from crispr_integration.repeat_type import MODEL, validate_report
from crispr_web.execution.validation_repeat_type import validate_repeat_type_pipeline


@pytest.mark.parametrize("sequence,rows", [("GGAAAA", [[0, 1, .005]]), ("GAAAAG", [[0, 5, .005]]), ("CGAAAA", [[0, 1, .005]])])
def test_impossible_repeat_pairs_fail_both_gates(sequence, rows):
    with pytest.raises(ValueError):
        adapter_pairs(rows, sequence)
    with pytest.raises(ValueError):
        backend_pairs(rows, sequence, None)


@pytest.mark.parametrize("sequence,rows", [("GAAAAU", [[0, 5, .5]]), ("UAAAAG", [[0, 5, .5]]), ("GGAAACCA", [[0, 5, .3], [1, 6, .3]])])
def test_wobble_and_crossing_ensemble_pairs_remain_valid(sequence, rows):
    assert adapter_pairs(rows, sequence)
    assert backend_pairs(rows, sequence, None)[0]


def test_repeat_type_reverse_complement_binding_keeps_literal_counts(tmp_path):
    def sha(value):
        return hashlib.sha256(value.encode()).hexdigest()
    records = {"forward": "AACCGG", "reverse": "CCGGTT", "copy": "AACCGG", "ambiguous": "NAAA"}
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for index, (key, sequence) in enumerate(records.items()):
        (inputs / f"{index}.fasta").write_text(f">{key}\n{sequence}\n")
    digest = sha("".join(f">{k}\n{v}\n" for k,v in records.items()))
    counts = {"records": 4, "distinct_classified_sequences": 2, "predicted": 3, "uncertain": 0, "unsupported": 1}
    doc = {"schema": {"name": "repeattyper.classification", "version": "1.0.0"}, "tool": {"name": "RepeatTyper", "version": "2.0.1"}, "execution_completed": True, "model": MODEL,
           "input": {"sha256": digest, "record_count": 4}, "policy": {"threshold": .75, "kmer": 4, "iteration_range": [0,41], "score_semantics": "closed_set_model_score_not_calibrated_confidence", "orientation": "canonical_kmers_no_strand_inference"}, "counts": counts, "results": []}
    for key, sequence in records.items():
        unsupported = "N" in sequence
        doc["results"].append({"id": key, "sequence": sequence, "sequence_sha256": sha(sequence), "native_score": None if unsupported else .9, "winner": None if unsupported else "I-A", "prediction": None if unsupported else "I-A", "status": "unsupported" if unsupported else "predicted", "reason": "ambiguous_iupac" if unsupported else None})
    claim = SimpleNamespace(molecule="DNA", tool_options={})
    integration = {"workflow": "repeat_type", "stages": {"repeattyper": {"status": "completed", "counts": counts, "model_id": MODEL["id"], "model_sha256": MODEL["files"]["xgb_repeats.json"]}}}
    expected = {k: (len(v), sha(v)) for k,v in records.items()}
    output = tmp_path / "repeat-type"
    output.mkdir()
    def backend(document):
        (output / "report.json").write_text(json.dumps(document))
        return validate_repeat_type_pipeline(claim, tmp_path, expected, integration, inputs)
    validate_report(doc, records, digest)
    assert backend(doc)
    altered = copy.deepcopy(doc)
    altered["results"][1]["native_score"] = .8
    with pytest.raises(ValueError, match="Reverse-complement"):
        validate_report(altered, records, digest)
    assert not backend(altered)


@pytest.mark.parametrize("field,value", [("alpha_bdm", "NaN"), ("ln_lh_bdm", "-inf"), ("nb of reconstructed deletions", "-1"), ("nb of reconstructed insertions", "1.5"), ("insertion_rate_bdm", "-0.1")])
def test_published_evolution_metrics_rejected_at_both_gates(field, value):
    from crispr_integration.evolution_validation import validate_metrics
    from crispr_web.execution.validation_evolution import valid_metrics
    with pytest.raises(ValueError):
        validate_metrics({field: value})
    assert not valid_metrics({field: value})


def test_evolution_optional_extreme_signed_and_unpublished_metrics_survive():
    from crispr_integration.evolution_validation import validate_metrics
    from crispr_web.execution.validation_evolution import valid_metrics
    row = {"alpha_bdm": "1e300", "ln_lh_bdm": "-1e200", "test_statistic (-2*ln_lh_ratio)": "-0.01", "deletion_rate_idm": "", "insertion_rate_idm": None, "nb of reconstructed deletions": "0", "nb of existent spacers": "[1, 3, 5]", "alternative_estimator": "NaN"}
    validate_metrics(row)
    assert valid_metrics(row)


@pytest.mark.parametrize("alignment", ["-AAAA:-AAAA", "AA-AA:AA-AA", "AAAA-:AAAA-"])
def test_native_forbidden_double_gap_tir_rejected_independently(alignment):
    import csv
    import io
    from crispr_integration.insertion_evidence import document as adapter_document
    from crispr_web.execution.validation_insertion_evidence import document as backend_document, FIELDS
    from crispr_web.execution.validation_insertion_context import MODELS
    source = b">f000001.s000001\nAAAACCCCTTTT\n"
    row = ["seq000001", "IS1", "IS1_1", "1", "12", "12", "2", "1", "4", "9", "12", "8", "4", "4", "0", "0", "0", "", "1", "1e-20", "-1", "p", "0", "AAAA:AAAA"]
    def tsv(values):
        stream = io.StringIO()
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(FIELDS)
        writer.writerow(values)
        return stream.getvalue().encode()
    for build in (adapter_document, backend_document):
        assert build(source, tsv(row), {"threads": 1, "include_partial": True}, MODELS)["elements"]
    changed = dict(zip(FIELDS, row, strict=True))
    changed.update(tir=alignment, irLen="5", irId="4", nGaps="1")
    assert {"irLen", "irId", "nGaps"} <= set(FIELDS)
    for build in (adapter_document, backend_document):
        with pytest.raises(ValueError, match="gap"):
            build(source, tsv([changed[key] for key in FIELDS]), {"threads": 1, "include_partial": True}, MODELS)


@pytest.mark.parametrize("value", ["+", "-", "both"])
def test_typed_strands_are_lossless_but_identifier_formulas_stay_shielded(value):
    from crispr_web.publication.tabular import scientific_cell
    assert scientific_cell(value, strand=True) == value
    if value != "both":
        assert scientific_cell(value).startswith("'")
    for hostile in ("=1+1", "@cmd", " +1", "-protein", "+SUM(A1:A2)"):
        assert scientific_cell(hostile, strand=True, protein_alignment=True) == "'" + hostile
    assert scientific_cell("--MKWVTF", protein_alignment=True) == "--MKWVTF"
