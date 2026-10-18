"""Geometry checks compare projections with independently indexed source bases."""
from dataclasses import replace
import random

import pytest

from crisprleader.coordinates import extract_candidates, extract_segments, project_interval, reverse_complement
from crisprleader.errors import InputError
from crisprleader.models import AnalysisConfig, ArrayAnnotation, RepeatAnchor, Segment, SourceRecord


def arc(start, length, n):
    if start + length <= n:
        return (Segment(start, start + length),)
    return (Segment(start, n), Segment(0, (start + length) % n))


def annotation(source, start=8, length=10, repeat_length=3):
    n = source.length
    segments = arc(start, length, n)
    first = arc(start, repeat_length, n)
    last = arc((start + length - repeat_length) % n, repeat_length, n)
    return ArrayAnnotation('array', source.id, segments,
                           (RepeatAnchor(first, extract_segments(source.sequence, first)),
                            RepeatAnchor(last, extract_segments(source.sequence, last))),
                           extract_segments(source.sequence, first))


def test_iupac_reverse_complement():
    assert reverse_complement('acgtryswkmbdhvn') == 'NBDHVKMWSRYACGT'
    assert reverse_complement(reverse_complement('ACGTRYSWKMBDHVN')) == 'ACGTRYSWKMBDHVN'
    for invalid in ('A-U', 'A U', 'ACGTX', 15, None):
        with pytest.raises(InputError):
            reverse_complement(invalid)


def test_exact_linear_both_flanks_and_observed_repeat():
    source = SourceRecord('s', 'ACTGACGAAGTCCGGTACGATACTGACGACT', 'linear')
    array = annotation(source)
    left, right = extract_candidates(source, array, AnalysisConfig(flank_length=6))
    assert left.sequence == source.sequence[2:8]
    assert right.sequence == reverse_complement(source.sequence[18:24])
    assert left.segments == (Segment(2, 8),)
    assert right.segments == (Segment(18, 24),)
    assert left.observed_repeat == source.sequence[8:11]
    assert right.observed_repeat == reverse_complement(source.sequence[15:18])
    # Non-identical observed terminal repeats must not be replaced by consensus.
    assert right.observed_repeat != right.consensus_repeat
    assert left.truncated_reason is right.truncated_reason is None
    assert left.metadata['window_is_biological_boundary'] is False


@pytest.mark.parametrize('topology,reason', [('linear', 'contig_boundary'), ('unknown', 'unknown_topology_boundary')])
def test_edges_are_clipped_not_wrapped(topology, reason):
    source = SourceRecord('s', 'ACGTACGTACGTACGTACGT', topology)
    array = annotation(source, 0, 20)
    left, right = extract_candidates(source, array)
    assert left.sequence == right.sequence == ''
    assert left.segments == right.segments == ()
    assert left.truncated_reason == right.truncated_reason == reason
    for candidate in (left, right):
        mapped = project_interval(candidate, 0, 3, include_repeat=True)
        assert extract_segments(source.sequence, mapped, candidate.strand) == candidate.observed_repeat


def test_circular_flanks_stop_before_reentering_array():
    source = SourceRecord('s', 'ACTGACGAAGTCCGGTACGATACTGACGACT', 'circular', 'explicit record')
    array = annotation(source, 25, 10)
    left, right = extract_candidates(source, array)
    assert left.available_length == right.available_length == source.length - 10
    assert left.truncated_reason == right.truncated_reason == 'circular_array_boundary'
    for candidate in (left, right):
        coordinates = [i for s in candidate.segments for i in range(s.start, s.end)]
        array_coordinates = {i for s in array.segments for i in range(s.start, s.end)}
        assert len(coordinates) == len(set(coordinates))
        assert set(coordinates).isdisjoint(array_coordinates)


def test_negative_junction_projection_crosses_origin_correctly():
    source = SourceRecord('s', 'AACCGGTTACGATCCGGTAA', 'circular', 'record')
    array = annotation(source, 10, 9)
    candidate = extract_candidates(source, array, AnalysisConfig(flank_length=5))[1]
    assert candidate.segments == (Segment(19, 20), Segment(0, 4))
    assert project_interval(candidate, 2, 7, include_repeat=True) == (Segment(17, 20), Segment(0, 2))
    mapped = project_interval(candidate, 2, 7, include_repeat=True)
    assert extract_segments(source.sequence, mapped, '-') == candidate.junction_sequence[2:7]


def test_random_all_oriented_subintervals_reconstruct_with_origin_anchors():
    rng = random.Random(4781)
    for _ in range(75):
        n = rng.randrange(8, 60)
        source = SourceRecord('s', ''.join(rng.choices('ACGTRYN', k=n)), 'circular', 'record')
        length = rng.randrange(6, n + 1)
        array = annotation(source, rng.randrange(n), length)
        for candidate in extract_candidates(source, array, AnalysisConfig(flank_length=rng.randrange(1, 70))):
            for include_repeat in (False, True):
                sequence = candidate.junction_sequence if include_repeat else candidate.sequence
                for start in range(len(sequence)):
                    for end in range(start + 1, len(sequence) + 1):
                        mapped = project_interval(candidate, start, end, include_repeat)
                        assert len(mapped) <= 2
                        assert extract_segments(source.sequence, mapped, candidate.strand) == sequence[start:end]
                        positions = [i for s in mapped for i in range(s.start, s.end)]
                        assert len(positions) == len(set(positions))


@pytest.mark.parametrize('mutate', [
    lambda s, a: (replace(s, topology='unknown'), a),
    lambda s, a: (replace(s, topology_provenance=None), a),
    lambda s, a: (s, replace(a, source_id='other')),
    lambda s, a: (s, replace(a, terminal_repeats=a.terminal_repeats[::-1])),
    lambda s, a: (s, replace(a, strand_constraint=True)),
    lambda s, a: (s, replace(a, segments=(Segment(18, 20), Segment(1, 8)))),
    lambda s, a: (s, replace(a, segments=(Segment(3, 20), Segment(0, 8)))),
    lambda s, a: (s, replace(a, terminal_repeats=(replace(a.terminal_repeats[0], sequence='AAA'), a.terminal_repeats[1]))),
])
def test_bad_array_geometry_fails_closed(mutate):
    source = SourceRecord('s', 'ACGTACGTACGTACGTACGT', 'circular', 'record')
    array = annotation(source, 18, 10)
    with pytest.raises(InputError):
        extract_candidates(*mutate(source, array))


@pytest.mark.parametrize('interval', [(0, 0), (-1, 1), (1, 100), (True, 2), (0, 1.5)])
def test_projection_rejects_invalid_offsets(interval):
    source = SourceRecord('s', 'ACGT' * 10)
    candidate = extract_candidates(source, annotation(source))[0]
    with pytest.raises(InputError):
        project_interval(candidate, *interval)
