import json

import pytest

from sp_model.crisprfinder_to_spacerfasta import CRISPRFinderParser
from sp_model.experiments_tools import read_fasta
from sp_model.input_parser import InputParser


def test_read_fasta_supports_wrapped_records_and_blank_lines(tmp_path):
    fasta = tmp_path / 'group.fa'
    fasta.write_text('>array-1, metadata\n1, 2,\n3\n\n>array-2\n1, 3\n', encoding='utf-8')
    assert read_fasta(fasta) == {
        'array-1, metadata': '1, 2,3',
        'array-2': '1, 3',
    }


@pytest.mark.parametrize(
    ('content', 'message'),
    [
        ('1,2\n>a\n1\n', 'before the first'),
        ('>\n1\n', 'Empty FASTA header'),
        ('>a\n\n', 'no spacer IDs'),
        ('>a\n1\n>a\n2\n', 'Duplicate FASTA header'),
    ],
)
def test_read_fasta_rejects_malformed_records(tmp_path, content, message):
    fasta = tmp_path / 'bad.fa'
    fasta.write_text(content, encoding='utf-8')
    with pytest.raises(ValueError, match=message):
        read_fasta(fasta)


def test_ccdb_adapter_does_not_create_missing_orientation_directory(tmp_path):
    group = tmp_path / 'group'
    positive = group / 'pos_strand'
    positive.mkdir(parents=True)
    (positive / 'array.fa').write_text('>ignored\nACGT\nARYN\n', encoding='utf-8')
    output = tmp_path / 'group.fa'
    parser = InputParser(group, output)
    assert not (group / 'neg_strand').exists()
    assert parser.file_name_to_spacer_indices == {'array': [0, 1]}
    assert output.read_text(encoding='utf-8') == '>array\n0, 1\n'


def test_ccdb_adapter_reverses_negative_arrays_with_iupac(tmp_path):
    group = tmp_path / 'group'
    negative = group / 'neg_strand'
    negative.mkdir(parents=True)
    (negative / 'reverse.fa').write_text('ACGT\nARYN\n', encoding='utf-8')
    parser = InputParser(group, tmp_path / 'out.fa')
    assert list(parser.spacer_seq_to_number) == ['NRYT', 'ACGT']


def test_ccdb_adapter_resolves_identical_cross_orientation_duplicate_like_upstream(tmp_path):
    group = tmp_path / 'group'
    positive = group / 'pos_strand'
    negative = group / 'neg_strand'
    positive.mkdir(parents=True)
    negative.mkdir()
    for folder in (positive, negative):
        (folder / 'same.fa').write_text('>spacer\nAAGT\n', encoding='utf-8')
    parser = InputParser(group, tmp_path / 'out.fa')
    assert list(parser.spacer_seq_to_number) == ['ACTT']
    assert parser.file_name_to_spacer_indices == {'same': [0]}


def test_ccdb_adapter_rejects_conflicting_duplicate_names(tmp_path):
    group = tmp_path / 'group'
    positive = group / 'pos_strand'
    negative = group / 'neg_strand'
    positive.mkdir(parents=True)
    negative.mkdir()
    (positive / 'same.fa').write_text('AAAA\n', encoding='utf-8')
    (negative / 'same.fa').write_text('CCCC\n', encoding='utf-8')
    with pytest.raises(ValueError, match='Conflicting duplicate'):
        InputParser(group, tmp_path / 'out.fa')


def _crispr(name, spacers, sequence_id='sequence', direction='+'):
    regions = []
    for index, spacer in enumerate(spacers):
        regions.extend([
            {'Type': 'DR', 'Sequence': 'ACGT'},
            {'Type': 'Spacer', 'Sequence': spacer},
        ])
    return {
        'Evidence_Level': 4,
        'Name': name,
        'Start': 1,
        'End': 100,
        'DR_Consensus': 'ACGT',
        'CRISPRDirection': direction,
        'Regions': regions,
    }


def test_crisprcasfinder_keeps_exactly_two_arrays_and_disambiguates_names(tmp_path):
    first = tmp_path / 'first.json'
    second = tmp_path / 'second.json'
    first.write_text(json.dumps({'Sequences': [{'Id': 's1', 'Crisprs': [_crispr('same', ['AAAA', 'AAAT'])]}]}))
    second.write_text(json.dumps({'Sequences': [{'Id': 's2', 'Crisprs': [_crispr('same', ['AAAT', 'AATT'])]}]}))
    output = tmp_path / 'groups'
    small = tmp_path / 'small'
    parser = CRISPRFinderParser(
        [first, second],
        output_folder=output,
        output_folder_for_too_small_groups=small,
        filter_more_than_x_in_group=2,
    )
    group_files = list(output.glob('*.fa'))
    assert len(group_files) == 1
    assert not list(small.glob('*.fa'))
    text = group_files[0].read_text(encoding='utf-8')
    assert '>same ' in text
    assert '>s2::same ' in text
    assert len(parser.get_cluster_spacer_id_to_spacer_seq()['g_0']) == 3


def test_crisprcasfinder_edit_and_overlap_clustering_are_transitive(tmp_path):
    source = tmp_path / 'result.json'
    source.write_text(json.dumps({'Sequences': [{'Id': 's', 'Crisprs': [
        _crispr('a', ['AAAA']),
        _crispr('b', ['AATT']),
        _crispr('bridge', ['AAAT']),
    ]}]}))
    parser = CRISPRFinderParser(source, max_levenshtein_distance_spacers=1)
    assert len(parser.clustered_groups) == 1
    assert len(parser.clustered_groups['g_0']) == 3


def test_crisprcasfinder_rejects_invalid_schema(tmp_path):
    source = tmp_path / 'bad.json'
    source.write_text(json.dumps({'not_sequences': []}))
    with pytest.raises(ValueError, match="'Sequences' list"):
        CRISPRFinderParser(source)
