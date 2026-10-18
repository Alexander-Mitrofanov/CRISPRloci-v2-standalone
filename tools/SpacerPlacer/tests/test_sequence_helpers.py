import random

import pytest

from sp_model.helpers.sequence import (
    cluster_by_edit_distance,
    cluster_by_overlap,
    levenshtein_distance,
    reverse_complement,
)


def reference_distance(first, second):
    matrix = [[0] * (len(second) + 1) for _ in range(len(first) + 1)]
    for row in range(len(first) + 1):
        matrix[row][0] = row
    for column in range(len(second) + 1):
        matrix[0][column] = column
    for row, first_character in enumerate(first, start=1):
        for column, second_character in enumerate(second, start=1):
            matrix[row][column] = min(
                matrix[row - 1][column] + 1,
                matrix[row][column - 1] + 1,
                matrix[row - 1][column - 1] + (first_character != second_character),
            )
    return matrix[-1][-1]


def test_reverse_complement_supports_iupac_and_case():
    assert reverse_complement('ACGTURYSWKMBDHVN') == 'NBDHVKMWSRYAACGT'
    assert reverse_complement('acgtn') == 'nacgt'


def test_bounded_levenshtein_matches_reference_for_random_inputs():
    generator = random.Random(2357)
    alphabet = 'ACGT'
    for _ in range(250):
        first = ''.join(generator.choice(alphabet) for _ in range(generator.randrange(12)))
        second = ''.join(generator.choice(alphabet) for _ in range(generator.randrange(12)))
        expected = reference_distance(first, second)
        assert levenshtein_distance(first, second) == expected
        for cutoff in range(4):
            assert levenshtein_distance(first, second, cutoff) == min(expected, cutoff + 1)


def test_bounded_levenshtein_handles_long_adversarial_inputs_exactly():
    first = 'A' * 2000
    second = 'A' * 999 + 'C' + 'A' * 1000
    assert levenshtein_distance(first, second, 0) == 1
    assert levenshtein_distance(first, second, 1) == 1
    assert levenshtein_distance(first, second + 'G' * 3, 1) == 2


def test_edit_distance_clustering_is_transitive_and_stable():
    assert cluster_by_edit_distance(['AAA', 'ATT', 'AAT', 'AAA'], 1) == [['AAA', 'ATT', 'AAT']]
    assert cluster_by_edit_distance(['CCC', 'AAA'], 0) == [['CCC'], ['AAA']]


def test_overlap_clustering_merges_bridge_arrays_transitively():
    arrays = [['a'], ['b'], ['a', 'b'], ['isolated']]
    assert cluster_by_overlap(arrays) == [[0, 1, 2], [3]]


def test_negative_cutoffs_are_rejected():
    with pytest.raises(ValueError, match='non-negative'):
        levenshtein_distance('A', 'B', -1)
    with pytest.raises(ValueError, match='non-negative'):
        cluster_by_edit_distance(['A'], -1)
