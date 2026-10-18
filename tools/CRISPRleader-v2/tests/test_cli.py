"""CLI engineering checks; fixtures carry no biological validity claim."""
import json
from pathlib import Path

import pytest

from crisprleader.cli import main

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'


def inputs():
    return ['--fasta', str(EXAMPLES / 'synthetic.fa'), '--arrays', str(EXAMPLES / 'synthetic.arrays.json')]


def test_extract_validate_and_existing_output(tmp_path, capsys):
    output = tmp_path / 'analysis'
    assert main(['extract', *inputs(), '--flank-length', '12', '--output', str(output)]) == 0
    message = json.loads(capsys.readouterr().out)
    assert message['prediction_performed'] is False
    assert message['analysis_mode'] == 'context_extraction'
    result = json.loads((output / 'result.json').read_text())
    assert result['predictions'] == [] and result['reference'] is None
    assert main(['validate', '--fasta', str(EXAMPLES / 'synthetic.fa'), '--result', str(output / 'result.json')]) == 0
    capsys.readouterr()
    before = (output / 'manifest.json').read_bytes()
    assert main(['extract', *inputs(), '--output', str(output)]) == 2
    message = json.loads(capsys.readouterr().err)
    assert message['execution_status'] == 'failed'
    assert (output / 'manifest.json').read_bytes() == before


def test_native_annotation_validation(capsys):
    assert main(['validate', *inputs()]) == 0
    assert json.loads(capsys.readouterr().out)['arrays'] == 1


def test_prediction_is_explicitly_unavailable_and_writes_nothing(tmp_path, capsys):
    output = tmp_path / 'prediction'
    assert main(['predict', *inputs(), '--output', str(output)]) == 2
    message = json.loads(capsys.readouterr().err)
    assert message['error'] == 'ExecutionError'
    assert 'unavailable' in message['message']
    assert not output.exists()


def test_missing_default_references_is_error_not_no_call(capsys):
    assert main(['references', 'verify']) == 2
    message = json.loads(capsys.readouterr().err)
    assert message['error'] == 'ReferenceError'
    assert message['execution_status'] == 'failed'


@pytest.mark.parametrize('length', ['0', '-1', '5001'])
def test_invalid_window_fails_without_artifacts(length, tmp_path, capsys):
    output = tmp_path / 'result'
    assert main(['extract', *inputs(), '--output', str(output), '--flank-length', length]) == 2
    assert json.loads(capsys.readouterr().err)['error'] == 'InputError'
    assert not output.exists()


@pytest.mark.parametrize('rows', [None, 'bad', [None], [{'id': []}], [{'id': {}, 'topology': []}]])
def test_malformed_result_returns_json_error_not_traceback(rows, tmp_path, capsys):
    document = {'sources': rows, 'schema_version': '1.0.0', 'coordinate_system': 'zero_based_half_open',
                'execution_status': 'completed', 'analysis_mode': 'context_extraction'}
    result = tmp_path / 'bad.json'
    result.write_text(json.dumps(document))
    assert main(['validate', '--fasta', str(EXAMPLES / 'synthetic.fa'), '--result', str(result)]) == 2
    assert json.loads(capsys.readouterr().err)['error'] == 'InputError'


def test_validation_rejects_falsified_context_boundary_claim(tmp_path, capsys):
    output = tmp_path / 'context'
    assert main(['extract', *inputs(), '--output', str(output)]) == 0
    capsys.readouterr()
    result = output / 'result.json'
    document = json.loads(result.read_text())
    document['candidates'][0]['metadata']['window_is_biological_boundary'] = True
    result.write_text(json.dumps(document))
    assert main(['validate', '--fasta', str(EXAMPLES / 'synthetic.fa'), '--result', str(result)]) == 2
    assert json.loads(capsys.readouterr().err)['error'] == 'InputError'
