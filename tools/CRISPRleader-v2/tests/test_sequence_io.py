import hashlib
import json
from dataclasses import replace

import pytest

from crisprleader.errors import InputError
from crisprleader.models import AnalysisConfig
from crisprleader.sequence_io import load_fasta, load_inputs, read_json


@pytest.fixture
def native_files(tmp_path):
    sequence = 'AAAACGTAGACGTTTT'
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>source description\n' + sequence.lower() + '\n')
    doc = {'schema_version': '1.0.0', 'coordinate_system': 'zero_based_half_open',
           'sources': [{'id': 'source', 'length': len(sequence), 'sha256': hashlib.sha256(sequence.encode()).hexdigest(), 'topology': 'linear'}],
           'arrays': [{'id': 'a', 'source_id': 'source', 'segments': [{'start': 3, 'end': 12}],
                       'terminal_repeats': [{'segments': [{'start': 3, 'end': 6}], 'sequence': sequence[3:6]},
                                            {'segments': [{'start': 9, 'end': 12}], 'sequence': sequence[9:12]}],
                       'consensus_repeat': sequence[3:6]}]}
    return fasta, tmp_path / 'arrays.json', doc


def save_load(files, config=None):
    fasta, path, doc = files
    path.write_text(json.dumps(doc))
    return load_inputs(fasta, path, config)


def test_valid_native_bindings(native_files):
    sources, arrays = save_load(native_files)
    assert sources['source'].sequence == 'AAAACGTAGACGTTTT'
    assert sources['source'].header == 'source description'
    assert sources['source'].topology == 'linear'
    assert arrays[0].strand == '?'
    assert arrays[0].strand_constraint is False


@pytest.mark.parametrize('text', [
    '', '>a\n', 'ACGT\n', '>a\nACGT\n>a\nACGT\n', '>a\nACGU\n', '>a\nAC-G\n',
    '>a\nACGT\n>b\n', '>\nACGT\n', '> \nACGT\n', '>a\x00\nACGT\n', '>a\nACGT1\n',
])
def test_invalid_fasta(tmp_path, text):
    path = tmp_path / 'input.fa'
    path.write_text(text)
    with pytest.raises(InputError):
        load_fasta(path)


def test_normalized_multifasta_and_total_base_limit(tmp_path):
    path = tmp_path / 'input.fa'
    path.write_bytes(b'>one first\r\nac gt\r\n\n>two\nryswkmbdhvn\n')
    sources = load_fasta(path, max_bases=15)
    assert list(sources) == ['one', 'two']
    assert sources['two'].sequence == 'RYSWKMBDHVN'
    with pytest.raises(InputError, match='max_bases'):
        load_fasta(path, max_bases=14)
    with pytest.raises(InputError):
        load_fasta(path, max_bases=True)


def test_header_and_json_byte_bounds(tmp_path):
    path = tmp_path / 'input.fa'
    path.write_text('>' + 'a' * 4096 + '\nACGT\n')
    with pytest.raises(InputError, match='header'):
        load_fasta(path)
    path.write_text('{"long":"' + 'a' * 32 + '"}')
    with pytest.raises(InputError, match='exceeds'):
        read_json(path, max_bytes=20)


def test_special_files_rejected_without_blocking(tmp_path):
    import os
    path = tmp_path / 'pipe.fa'
    os.mkfifo(path)
    with pytest.raises(InputError, match='regular file'):
        load_fasta(path)


@pytest.mark.parametrize('text', ['{"a":1,"a":2}', '{"a":NaN}', '{"a":1e999}', '[' * 40 + '0' + ']' * 40, '{broken'])
def test_json_duplicate_nonfinite_deep_and_syntax_rejected(tmp_path, text):
    path = tmp_path / 'input.json'
    path.write_text(text)
    with pytest.raises(InputError):
        read_json(path)


@pytest.mark.parametrize('mutate', [
    lambda d: d.update(schema_version='2.0.0'),
    lambda d: d.update(coordinate_system='one_based'),
    lambda d: d.update(unknown='value'),
    lambda d: d['sources'][0].update(sha256='0' * 64),
    lambda d: d['sources'][0].update(length=True),
    lambda d: d['sources'][0].update(length=99),
    lambda d: d['sources'][0].update(topology='circular'),
    lambda d: d['sources'][0].update(topology='round'),
    lambda d: d['sources'][0].update(unknown='value'),
    lambda d: d['sources'].append(d['sources'][0]),
    lambda d: d['arrays'].append(d['arrays'][0]),
    lambda d: d['arrays'][0].update(source_id='other'),
    lambda d: d['arrays'][0].update(unknown='value'),
    lambda d: d['arrays'][0].update(strand_constraint=True),
    lambda d: d['arrays'][0].update(strand_constraint=1),
    lambda d: d['arrays'][0].update(strand='forward'),
    lambda d: d['arrays'][0].update(consensus_repeat='A-C'),
    lambda d: d['arrays'][0].update(metadata=[]),
    lambda d: d['arrays'][0]['segments'][0].update(start=True),
    lambda d: d['arrays'][0]['segments'][0].update(end=99),
    lambda d: d['arrays'][0]['terminal_repeats'][0].update(sequence='AAA'),
    lambda d: d['arrays'][0]['terminal_repeats'][0]['segments'][0].update(start=2),
])
def test_native_mutations_fail_closed(native_files, mutate):
    mutate(native_files[2])
    with pytest.raises(InputError):
        save_load(native_files)


def test_source_list_cannot_silently_omit_fasta_record(native_files):
    with native_files[0].open('a') as stream:
        stream.write('>extra\nACGT\n')
    with pytest.raises(InputError, match='every FASTA'):
        save_load(native_files)


def test_metadata_optional_fields_and_constraints_survive(native_files):
    doc = native_files[2]
    doc['sources'][0]['metadata'] = {'organism': 'example', 'missing': None}
    doc['arrays'][0].update(strand='-', strand_constraint=True, metadata={'confidence': None})
    sources, arrays = save_load(native_files)
    assert sources['source'].metadata['missing'] is None
    assert arrays[0].strand == '-'
    assert arrays[0].strand_constraint is True
    assert arrays[0].metadata['confidence'] is None


def test_max_array_and_candidate_budgets(native_files):
    for config in (replace(AnalysisConfig(), max_candidates=1), replace(AnalysisConfig(), max_arrays=0)):
        with pytest.raises(InputError):
            save_load(native_files, config)


def test_circular_native_array_and_origin_spanning_anchor(native_files):
    doc = native_files[2]
    doc['sources'][0].update(topology='circular', topology_provenance='INSDC accession test.1')
    seq = 'AAAACGTAGACGTTTT'
    doc['arrays'][0].update(segments=[{'start': 14, 'end': 16}, {'start': 0, 'end': 6}],
        terminal_repeats=[{'segments': [{'start': 14, 'end': 16}, {'start': 0, 'end': 1}], 'sequence': seq[14:] + seq[:1]},
                          {'segments': [{'start': 3, 'end': 6}], 'sequence': seq[3:6]}])
    sources, arrays = save_load(native_files)
    assert sources['source'].topology == 'circular'
    assert len(arrays[0].terminal_repeats[0].segments) == 2
