"""Synthetic reporting examples must never be presented as real performance."""
import hashlib
import json
from pathlib import Path

import pytest

from crisprleader.cli import main
from crisprleader.errors import InputError
from crisprleader.export import write_json_report
import crisprleader.export as exports

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'


def arguments(output):
    return ['benchmark', '--dataset', str(EXAMPLES / 'benchmark.synthetic.dataset.json'),
            '--predictions', str(EXAMPLES / 'benchmark.synthetic.predictions.json'),
            '--output', str(output), '--bootstrap-replicates', '100', '--seed', '8']


def test_benchmark_cli_carries_no_superiority_claim_and_receipt(tmp_path, capsys):
    output = tmp_path / 'report.json'
    assert main(arguments(output)) == 0
    message = json.loads(capsys.readouterr().out)
    assert message['scientific_superiority_established'] is False
    report = json.loads(output.read_text())
    assert report['scientific_superiority_established'] is False
    for field, name in [('dataset', 'benchmark.synthetic.dataset.json'),
                        ('predictions', 'benchmark.synthetic.predictions.json')]:
        document = json.loads((EXAMPLES / name).read_text())
        digest = hashlib.sha256(json.dumps(document, sort_keys=True, separators=(',', ':'),
                                           ensure_ascii=True, allow_nan=False).encode()).hexdigest()
        assert report['input_document_receipts'][field]['canonical_json_sha256'] == digest
    before = output.read_bytes()
    assert main(arguments(output)) == 2
    assert json.loads(capsys.readouterr().err)['execution_status'] == 'failed'
    assert output.read_bytes() == before


def test_benchmark_rejects_mismatched_dataset_without_artifact(tmp_path, capsys):
    predictions = json.loads((EXAMPLES / 'benchmark.synthetic.predictions.json').read_text())
    predictions['dataset_id'] = 'different-dataset'
    source = tmp_path / 'wrong.json'
    source.write_text(json.dumps(predictions))
    output = tmp_path / 'report.json'
    args = arguments(output)
    args[args.index('--predictions') + 1] = str(source)
    assert main(args) == 2
    assert json.loads(capsys.readouterr().err)['error'] == 'InputError'
    assert not output.exists()


def test_atomic_json_publication_never_replaces_racing_file(tmp_path, monkeypatch):
    output = tmp_path / 'report.json'
    original = exports._publish_no_replace

    def competitor_writes_first(source, target):
        target.write_text('preexisting owner data')
        return original(source, target)

    monkeypatch.setattr(exports, '_publish_no_replace', competitor_writes_first)
    with pytest.raises(FileExistsError):
        write_json_report({'valid': True}, output)
    assert output.read_text() == 'preexisting owner data'
    assert list(tmp_path.iterdir()) == [output]


def test_invalid_json_is_never_partially_published(tmp_path):
    output = tmp_path / 'report.json'
    with pytest.raises(ValueError):
        write_json_report({'bad': float('nan')}, output)
    assert not list(tmp_path.iterdir())


def test_existing_symlink_report_is_not_followed(tmp_path):
    output = tmp_path / 'report.json'
    output.symlink_to(tmp_path / 'missing-owner-file')
    with pytest.raises(InputError):
        write_json_report({}, output)
    assert output.is_symlink() and not (tmp_path / 'missing-owner-file').exists()


def test_prediction_remains_unavailable_after_adding_benchmark(tmp_path, capsys):
    output = tmp_path / 'prediction'
    assert main(['predict', '--fasta', str(EXAMPLES / 'synthetic.fa'), '--arrays',
                 str(EXAMPLES / 'synthetic.arrays.json'), '--output', str(output)]) == 2
    assert 'unavailable' in json.loads(capsys.readouterr().err)['message']
    assert not output.exists()
