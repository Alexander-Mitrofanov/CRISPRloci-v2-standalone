"""Independent canonical-report fixtures; never import a detector checkout."""
from copy import deepcopy
import hashlib
import json

import pytest

import crisprleader.export as export_module
from crisprleader.adapters.identify import adapt_identify
from crisprleader.cli import main
from crisprleader.errors import InputError


def interval(start, end):
    return {'start': start + 1, 'end': end, 'length': end - start}


def make_report(sequence='AAAACGTAGACGTTTT', source_id='source', category='Bona-fide', strand='-'):
    array = {
        'id': 'CRISPR-1', 'category': category, 'category_key': 'bona_fide', 'category_position': 0,
        'certainty_score': 0.98, 'source_interval': interval(3, 12), 'source_sequence': sequence[3:12],
        'consensus_repeat': 'ACG', 'repeat_count': 2, 'spacer_count': 1,
        'repeats': [{'index': 1, 'source_interval': interval(3, 6), 'sequence': sequence[3:6]},
                    {'index': 2, 'source_interval': interval(9, 12), 'sequence': sequence[9:12]}],
        'spacers': [{'index': 1, 'source_interval': interval(6, 9), 'sequence': sequence[6:9]}],
        'orientation': {'status': 'computed', 'label': 'Reversed' if strand == '-' else 'Forward', 'strand': strand},
        'validation': {'status': 'valid', 'method': 'independent_original_sequence_slicing', 'source_reconstructed': True, 'issues': []},
        'future_optional_field': {'preserve': 17},
    }
    return {
        'schema': {'name': 'CRISPRidentify-report', 'version': '1.2.0'},
        'coordinate_system': {'indexing': 1, 'end_inclusive': True, 'reference': 'input_forward_strand', 'orientation_is_separate': True,
                              'zero_width_boundaries': {'indexing': 0, 'representation': 'interbase', 'reference': 'input_forward_strand'}},
        'source': {'id': source_id, 'length': len(sequence), 'sha256': hashlib.sha256(sequence.encode()).hexdigest(), 'sha256_scope': 'uppercase_sequence',
                   'header': source_id, 'file_name': 'input.fa'},
        'run': {'parameters': {}, 'analyses': {}}, 'cas': {},
        'summary': {'array_count': 1, 'accepted_array_count': int(category in ('Bona-fide', 'Possible')), 'category_counts': {category: 1}},
        'validation': {'status': 'valid', 'method': 'independent_original_sequence_slicing', 'validated_array_count': 1, 'invalid_array_count': 0},
        'arrays': [array],
    }


@pytest.fixture
def files(tmp_path):
    fasta = tmp_path / 'source.fa'
    fasta.write_text('>source\nAAAACGTAGACGTTTT\n')
    return fasta, tmp_path / 'report.json', make_report()


def run(files):
    fasta, path, report = files
    path.write_text(json.dumps(report))
    return adapt_identify(fasta, [path])


def test_source_forward_consensus_does_not_follow_supplied_reverse_strand(files):
    document = run(files)
    array = document['arrays'][0]
    assert array['segments'] == [{'start': 3, 'end': 12}]
    assert array['consensus_repeat'] == 'ACG'
    assert array['terminal_repeats'][0]['sequence'] == 'ACG'
    assert array['strand'] == '-'
    assert array['strand_constraint'] is False
    assert array['metadata']['identify_original']['future_optional_field'] == {'preserve': 17}
    assert document['sources'][0]['topology'] == 'unknown'


def test_adapter_cli_publishes_complete_native_annotations(files, tmp_path, capsys):
    files[1].write_text(json.dumps(files[2]))
    output = tmp_path / 'native.json'
    command = ['adapt', 'identify', '--fasta', str(files[0]), '--reports', str(files[1]),
               '--output', str(output)]
    assert main(command) == 0
    assert json.loads(capsys.readouterr().out) == {'status': 'valid', 'arrays': 1}
    assert json.loads(output.read_text()) == adapt_identify(files[0], [files[1]])
    previous = output.read_bytes()
    assert main(command) == 2
    capsys.readouterr()
    assert output.read_bytes() == previous


def test_adapter_write_failure_leaves_no_partial_annotations(files, tmp_path, capsys, monkeypatch):
    files[1].write_text(json.dumps(files[2]))
    output = tmp_path / 'native.json'

    def fail_after_partial_write(document, handle, **kwargs):
        handle.write('{"incomplete":')
        raise OSError('simulated disk-full write failure')

    monkeypatch.setattr(export_module.json, 'dump', fail_after_partial_write)
    assert main(['adapt', 'identify', '--fasta', str(files[0]), '--reports', str(files[1]),
                 '--output', str(output)]) == 2
    assert json.loads(capsys.readouterr().err)['execution_status'] == 'failed'
    assert not output.exists()
    assert list(tmp_path.glob('.native.json.*')) == []


def test_unknown_orientation_remains_unknown(files):
    files[2]['arrays'][0]['orientation'] = {'status': 'not_requested', 'strand': None, 'label': 'not computed'}
    assert run(files)['arrays'][0]['strand'] == '?'


@pytest.mark.parametrize('category', ['Alternative', 'Possible Discarded', 'Low score'])
def test_nonaccepted_arrays_explicitly_recorded(files, category):
    files[2].update(make_report(category=category))
    document = run(files)
    assert document['arrays'] == []
    assert document['sources'][0]['metadata']['excluded_arrays'] == [{'id': 'CRISPR-1', 'category': category}]


def test_multiple_sources_do_not_collide_display_ids(files, tmp_path):
    with files[0].open('a') as stream:
        stream.write('>other\nAAAACGTAGACGTTTT\n>unreported\nACGT\n')
    second = tmp_path / 'second.json'
    second.write_text(json.dumps(make_report(source_id='other')))
    files[1].write_text(json.dumps(files[2]))
    document = adapt_identify(files[0], [files[1], second])
    assert len({a['id'] for a in document['arrays']}) == 2
    assert [s['id'] for s in document['sources']] == ['source', 'other', 'unreported']
    assert document['sources'][2]['metadata']['identify_report_status'] == 'not_supplied'


@pytest.mark.parametrize('mutate', [
    lambda r: r['schema'].update(name='legacy'),
    lambda r: r['schema'].update(version='2.0.0'),
    lambda r: r['coordinate_system'].update(indexing=0),
    lambda r: r['coordinate_system'].update(indexing=True),
    lambda r: r['coordinate_system'].update(reference='reverse'),
    lambda r: r['coordinate_system']['zero_width_boundaries'].update(indexing=1),
    lambda r: r['source'].update(id='different'),
    lambda r: r['source'].update(sha256='0' * 64),
    lambda r: r['source'].update(length=True),
    lambda r: r['source'].update(sha256_scope='file_bytes'),
    lambda r: r['arrays'][0].update(category='NewCategory'),
    lambda r: r['arrays'][0].update(consensus_repeat='AC-G'),
    lambda r: r['arrays'][0].update(source_sequence='AAAAAAAAA'),
    lambda r: r['arrays'][0]['repeats'][0].update(sequence='AAA'),
    lambda r: r['arrays'][0]['repeats'][0]['source_interval'].update(start=5),
    lambda r: r['arrays'][0]['repeats'][0].update(index=True),
    lambda r: r['arrays'][0]['spacers'][0].update(sequence='AAA'),
    lambda r: r['arrays'][0]['spacers'][0]['source_interval'].update(start=6),
    lambda r: r['arrays'][0].update(spacers=[]),
    lambda r: r['arrays'][0]['orientation'].update(label='Forward'),
    lambda r: r['arrays'][0]['orientation'].update(status='unknown'),
    lambda r: r['arrays'][0]['validation'].update(source_reconstructed=False),
    lambda r: r['summary'].update(accepted_array_count=0),
    lambda r: r['summary']['category_counts'].update(Unknown=0),
    lambda r: r['validation'].update(validated_array_count=2),
])
def test_report_mutations_rejected(files, mutate):
    mutate(files[2])
    with pytest.raises(InputError):
        run(files)


def topology(report):
    source = report['source']
    return {
        'schema_name': 'CRISPRidentify-sequence-topology-evidence', 'schema_version': '1.0.0',
        'source': {'source_id': source['id'], 'source_length': source['length'], 'source_sha256': source['sha256'],
                   'hash_normalization': 'uppercase_ascii_iupac_without_reformatting'},
        'topology': 'circular', 'inference_used': False, 'evidence_complete': True, 'evidence_ambiguous': False,
        'provenance': {'method': 'frozen_insdc_flatfile', 'accession': 'source', 'locator': 'record.gb', 'date': '2026-09-27', 'evidence_sha256': 'a' * 64},
    }


def test_source_bound_topology_is_retained(files):
    files[2]['extensions'] = {'sequence_topology': topology(files[2])}
    document = run(files)
    source = document['sources'][0]
    assert source['topology'] == 'circular'
    assert 'record.gb' in source['topology_provenance']
    assert source['metadata']['identify_report_metadata']['extensions'] == files[2]['extensions']


@pytest.mark.parametrize('mutate', [
    lambda t: t['source'].update(source_sha256='0' * 64),
    lambda t: t.update(evidence_complete=False),
    lambda t: t.update(inference_used=True),
    lambda t: t.update(schema_version='2.0.0'),
    lambda t: t['provenance'].update(method='guess'),
])
def test_unbound_or_unknown_topology_contract_rejected(files, mutate):
    evidence = topology(files[2])
    mutate(evidence)
    files[2]['extensions'] = {'sequence_topology': evidence}
    with pytest.raises(InputError):
        run(files)


def test_joined_origin_arrays_require_explicit_native_input(files):
    files[2]['extensions'] = {'circular_origin_loci': {
        'schema': {'name': 'CRISPRidentify-circular-origin-loci', 'version': '1.0.0'},
        'source_sha256': files[2]['source']['sha256'], 'joined_loci': [{'id': 'joined'}], 'joined_locus_count': 1,
    }}
    with pytest.raises(InputError, match='provide native annotations'):
        run(files)


def test_duplicate_report_source_rejected(files, tmp_path):
    with files[0].open('a') as stream:
        stream.write('>other\nACGT\n')
    files[1].write_text(json.dumps(files[2]))
    with pytest.raises(InputError, match='duplicate Identify report'):
        adapt_identify(files[0], [files[1], files[1]])


def test_adjacent_repeat_gap_is_zero_width_not_fake_interval(files):
    report = deepcopy(files[2])
    array = report['arrays'][0]
    array['source_interval'] = interval(3, 9)
    array['source_sequence'] = 'ACGTAG'
    array['repeats'][1].update(source_interval=interval(6, 9), sequence='TAG')
    array['spacers'] = []  # Canonical v1.2 stores deletions separately.
    files[2].update(report)
    assert run(files)['arrays'][0]['segments'] == [{'start': 3, 'end': 9}]
