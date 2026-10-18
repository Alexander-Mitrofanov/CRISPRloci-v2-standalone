import random
import unittest
from unittest.mock import mock_open, patch

from components.components_detection import VmatchRun
from components.legacy_hash import py37_hash_int_tuple
from components.module_detection import NATIVE_PATTERN_BATCH_SIZE
from components.sequence_path_vectorizer import BITMASK_16
from components.sequence_path_vectorizer import NATIVE_ASCII_PATH_FEATURES_AVAILABLE
from components.sequence_path_vectorizer import NATIVE_PATH_FEATURE_COUNTS_AVAILABLE
from components.sequence_path_vectorizer import NATIVE_PACKED_PATH_FEATURES_AVAILABLE
from components.sequence_path_vectorizer import SparseFeatureMatrix
from components.sequence_path_vectorizer import _copy_packed_feature_row
from components.sequence_path_vectorizer import _direct_feature_dict
from components.sequence_path_vectorizer import _native_packed_feature_row
from components.sequence_path_vectorizer import vectorize_sequence_paths
from components.sequence_path_vectorizer import _feature_hash3
from components.sequence_path_vectorizer import _feature_hash4


class ExactWidthFeatureHashTests(unittest.TestCase):
    def test_specialized_feature_hashes_match_python37_recurrence(self):
        generator = random.Random(20260821)
        values = [0, 1, 2, 65535, 65536, (1 << 32) - 1]
        values.extend(generator.randrange(1 << 32) for _ in range(2_000))
        for offset in range(0, len(values) - 4, 4):
            first, second, third, fourth = values[offset:offset + 4]
            self.assertEqual(
                (py37_hash_int_tuple(first, second, third) & BITMASK_16) + 1,
                _feature_hash3(first, second, third),
            )
            self.assertEqual(
                (
                    py37_hash_int_tuple(first, second, third, fourth)
                    & BITMASK_16
                ) + 1,
                _feature_hash4(first, second, third, fourth),
            )

    @unittest.skipUnless(
        NATIVE_PATH_FEATURE_COUNTS_AVAILABLE,
        "native path-feature count accelerator is not built",
    )
    def test_native_path_counts_match_events_and_python_fallback(self):
        generator = random.Random(20260823)
        sequences = ["A", "AC", "ACGT", "-", "N", "AAAAA"]
        sequences.extend(
            "".join(
                generator.choice("ACGTN-xyz")
                for _ in range(generator.randint(1, 100))
            )
            for _ in range(128)
        )
        native_counts = [
            _direct_feature_dict(sequence) for sequence in sequences
        ]
        with patch(
                "components.sequence_path_vectorizer."
                "_NATIVE_PATH_FEATURE_COUNTS", None):
            native_events = [
                _direct_feature_dict(sequence) for sequence in sequences
            ]
        with patch(
                "components.sequence_path_vectorizer."
                "NATIVE_PATH_FEATURES_AVAILABLE", False):
            fallback = [
                _direct_feature_dict(sequence) for sequence in sequences
            ]
        self.assertEqual(fallback, native_events)
        self.assertEqual(native_events, native_counts)

    @unittest.skipUnless(
        NATIVE_PACKED_PATH_FEATURES_AVAILABLE,
        "native packed path-feature accelerator is not built",
    )
    def test_native_packed_rows_and_dots_are_bit_exact(self):
        generator = random.Random(20260824)
        sequences = ["A", "AC", "ACGT", "-", "N", "AAAAA"]
        sequences.extend(
            "".join(
                generator.choice("ACGTN-xyz")
                for _ in range(generator.randint(1, 100))
            )
            for _ in range(128)
        )
        expected_rows = [
            _direct_feature_dict(sequence) for sequence in sequences
        ]
        packed_rows = [
            _native_packed_feature_row(sequence) for sequence in sequences
        ]
        self.assertEqual(
            expected_rows,
            [row.as_dict() for row in packed_rows],
        )
        expected = SparseFeatureMatrix(expected_rows)
        packed = vectorize_sequence_paths(sequences)
        self.assertTrue(all(
            expected_value.hex() == observed_value.hex()
            for expected_row, observed_row in zip(
                expected.dot(expected.T).toarray(),
                packed.dot(packed.T).toarray(),
            )
            for expected_value, observed_value in zip(
                expected_row, observed_row
            )
        ))

    @unittest.skipUnless(
        NATIVE_ASCII_PATH_FEATURES_AVAILABLE,
        "native ASCII path-feature accelerator is not built",
    )
    def test_native_ascii_rows_bypass_python_neighborhood_hashing(self):
        sequences = [
            "A",
            "ACGTN-xyz",
            "".join(chr(value) for value in range(128)),
        ]
        expected = [
            _direct_feature_dict(sequence) for sequence in sequences
        ]
        with patch(
                "components.sequence_path_vectorizer."
                "_neighborhood_radius_hashes",
                side_effect=AssertionError("Python neighborhood fallback used")):
            observed = [
                _native_packed_feature_row(sequence).as_dict()
                for sequence in sequences
            ]
        self.assertEqual(expected, observed)

    @unittest.skipUnless(
        NATIVE_ASCII_PATH_FEATURES_AVAILABLE,
        "native ASCII path-feature accelerator is not built",
    )
    def test_native_ascii_sentinel_uses_exact_abi3_fallback(self):
        sequence = "ACGTN-xyz"
        expected = _direct_feature_dict(sequence)
        with patch(
                "components.sequence_path_vectorizer."
                "_NATIVE_PATH_FEATURE_ASCII_ROW",
                return_value=-2) as ascii_row:
            observed = _native_packed_feature_row(sequence).as_dict()
        ascii_row.assert_called_once()
        self.assertEqual(expected, observed)

    def test_native_ascii_rejects_unsupported_inputs_before_ffi(self):
        with patch(
                "components.sequence_path_vectorizer."
                "_NATIVE_PATH_FEATURE_ASCII_ROW") as ascii_row:
            with self.assertRaisesRegex(ValueError, "requires ASCII"):
                _native_packed_feature_row("é")
            self.assertIsNone(_native_packed_feature_row(
                "A" * (BITMASK_16 + 1)
            ))
        ascii_row.assert_not_called()

    def test_bulk_copied_packed_row_is_exact_sized_and_independent(self):
        import ctypes

        source_ids = (ctypes.c_uint32 * 5)(1, 7, 65536, 9, 10)
        source_values = (ctypes.c_double * 5)(0.0, -0.0, 0.25, 1.0, 2.0)
        row = _copy_packed_feature_row(source_ids, source_values, 3)
        source_ids[0] = 99
        source_values[2] = 0.5
        self.assertEqual(3, row.length)
        self.assertEqual([1, 7, 65536], list(row.feature_ids))
        self.assertEqual(
            ["0x0.0p+0", "-0x0.0p+0", "0x1.0000000000000p-2"],
            [value.hex() for value in row.feature_values],
        )

    def test_packed_abi_absence_uses_exact_dictionary_fallback(self):
        sequences = ["A", "ACGT", "N-", "AAAAA"]
        expected = [
            _direct_feature_dict(sequence) for sequence in sequences
        ]
        with patch(
                "components.sequence_path_vectorizer."
                "_NATIVE_PATH_FEATURE_ROW", None), patch(
                "components.sequence_path_vectorizer."
                "_NATIVE_PATH_FEATURE_DOT", None):
            observed = vectorize_sequence_paths(sequences)
        self.assertEqual(expected, list(observed._rows))

    def test_native_batch_cap_preserves_bounded_task_geometry(self):
        self.assertEqual(128, NATIVE_PATTERN_BATCH_SIZE)

    def test_vmatch_index_omits_unused_tables(self):
        runner = VmatchRun.__new__(VmatchRun)
        with patch(
                "components.components_detection.subprocess.run") as run:
            runner._mkvtree_command()
        # Index the call tuple for compatibility with the older mock runtime
        # in the frozen CRISPRidentify test environment.
        command = run.call_args[0][0]
        self.assertEqual(
            [
                "mkvtree", "-db", "new_input.fa", "-dna", "-pl", "-lcp",
                "-suf", "-tis", "-ois", "-bwt",
            ],
            command,
        )

    def test_vmatch_overflow_fails_before_candidate_materialization(self):
        runner = VmatchRun.__new__(VmatchRun)
        runner.dna = "A" * 100
        runner.list_repeat_candidates = []
        raw_output = ">1 0 D 1 2\nA\n>1 3 D 1 5\nA\n"

        with patch.object(VmatchRun, "RESULT_LIMIT", 1):
            with patch("builtins.open", mock_open(read_data=raw_output)):
                with self.assertRaisesRegex(
                        RuntimeError, "refusing silent -best truncation"):
                    runner._take_vmatch_results_new()
        self.assertEqual([], runner.list_repeat_candidates)


if __name__ == "__main__":
    unittest.main()
