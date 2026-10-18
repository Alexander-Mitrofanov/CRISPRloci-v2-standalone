import json
from importlib.resources import files
from pathlib import Path

from jsonschema import Draft202012Validator
import pytest

from crisprleader.models import AnalysisConfig
from crisprleader.pipeline import extract


def test_packaged_schemas_accept_real_cli_fixture_and_extracted_context():
    directory = Path(__file__).resolve().parents[1] / 'examples'
    annotations = directory / 'synthetic.arrays.json'
    report = extract(directory / 'synthetic.fa', annotations, config=AnalysisConfig(flank_length=12))
    for filename, document in [('annotations.schema.json', json.loads(annotations.read_text())),
                               ('context-result.schema.json', report)]:
        schema = json.loads(files('crisprleader').joinpath('schemas', filename).read_text())
        Draft202012Validator.check_schema(schema)
        Draft202012Validator(schema).validate(document)
        damaged = dict(document, coordinate_system='one_based_inclusive')
        assert list(Draft202012Validator(schema).iter_errors(damaged))


@pytest.mark.parametrize('filename', ['annotations.schema.json', 'context-result.schema.json'])
def test_structural_schemas_reject_boolean_coordinates(filename):
    schema = json.loads(files('crisprleader').joinpath('schemas', filename).read_text())
    segment_schema = schema['$defs']['segment']
    assert list(Draft202012Validator(segment_schema).iter_errors({'start': True, 'end': 4}))
