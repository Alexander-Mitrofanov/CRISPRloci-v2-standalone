"""NumPy-only inference for the bundled CRISPRstrand repeat model.

The upstream CRISPRstrand model is stored as a Keras HDF5 file.  Loading it at
runtime pulls in TensorFlow, Keras and h5py even though CRISPRidentify only
needs inference for one small, fixed convolutional network.  This module
implements that network directly with NumPy and loads weights from a narrow,
versioned, checksummed ``.npz`` artifact.

The preprocessing and decision rule intentionally preserve CRISPRstrand v2:

* rows are ordered C, G, T, A;
* every non-canonical character is represented like N (all zeroes);
* repeats are padded with zeroes on the right or truncated to 38 bases;
* the full repeat is reverse-complemented *before* truncation; and
* equal forward/reverse probabilities select the reverse orientation.

The portable artifact is inference-only.  The legacy HDF5 model and exporter
remain the authoritative training/export inputs.
"""

from __future__ import division

import functools
import hashlib
import io
import json
import os

import numpy as np

from components.sequence_utils import reverse_complement


FORMAT_NAME = "crispridentify.portable-crisprstrand"
FORMAT_VERSION = 1
SEQUENCE_LENGTH = 38
ALPHABET_ROWS = ("C", "G", "T", "A")
KERNEL_WIDTHS = (4, 5, 6, 8, 10)
FILTERS_PER_BRANCH = 32
BATCH_NORMALIZATION_EPSILON = np.float32(0.001)
FORWARD_LABEL = "Forward"
REVERSE_LABEL = "Reversed"

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MODEL_DIRECTORY = os.path.join(
    PROJECT_ROOT, "trained_models", "strand"
)
DEFAULT_MANIFEST_PATH = os.path.join(
    DEFAULT_MODEL_DIRECTORY, "crisprstrand_model_r_v1.json"
)

_CANONICAL_BASES = frozenset(ALPHABET_ROWS)
_ROW_BY_BASE = dict((base, index) for index, base in enumerate(ALPHABET_ROWS))


def _canonical_json_bytes(value):
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


def _sha256(payload):
    return hashlib.sha256(payload).hexdigest()


def normalize_repeat(sequence):
    """Map non-canonical characters to ``N`` like legacy CRISPRidentify.

    Lowercase bases are deliberately not uppercased: the previous
    ``remove_non_canonical_char_from_string`` implementation treated them as
    non-canonical too.  Production consensus repeats are already uppercase.
    """

    if not isinstance(sequence, str):
        raise TypeError("repeat sequence must be a string")
    if not sequence:
        raise ValueError("repeat sequence must not be empty")
    return "".join(
        base if base in _CANONICAL_BASES else "N" for base in sequence
    )


def encode_repeats(sequences):
    """Return legacy one-hot tensors with shape ``(n, 4, 38, 1)``."""

    if isinstance(sequences, str):
        sequences = [sequences]
    else:
        sequences = list(sequences)
    if not sequences:
        raise ValueError("at least one repeat sequence is required")

    encoded = np.zeros(
        (len(sequences), len(ALPHABET_ROWS), SEQUENCE_LENGTH, 1),
        dtype=np.float32,
    )
    normalized = []
    for sequence_index, sequence in enumerate(sequences):
        repeat = normalize_repeat(sequence)
        normalized.append(repeat)
        for position, base in enumerate(repeat[:SEQUENCE_LENGTH]):
            row = _ROW_BY_BASE.get(base)
            if row is not None:
                encoded[sequence_index, row, position, 0] = np.float32(1.0)
    return normalized, encoded


def _expected_tensor_shapes():
    shapes = {}
    for branch_index, kernel_width in enumerate(KERNEL_WIDTHS, 1):
        shapes["conv_{}_kernel".format(branch_index)] = (
            len(ALPHABET_ROWS),
            kernel_width,
            1,
            FILTERS_PER_BRANCH,
        )
        for parameter in ("gamma", "beta", "moving_mean", "moving_variance"):
            shapes["bn_{}_{}".format(branch_index, parameter)] = (
                FILTERS_PER_BRANCH,
            )
    shapes.update(
        {
            "dense_1_kernel": (FILTERS_PER_BRANCH * len(KERNEL_WIDTHS), 256),
            "dense_1_bias": (256,),
            "dense_2_kernel": (256, 32),
            "dense_2_bias": (32,),
            "dense_3_kernel": (32, 1),
            "dense_3_bias": (1,),
        }
    )
    return shapes


EXPECTED_TENSOR_SHAPES = _expected_tensor_shapes()


def _require_exact_keys(mapping, required, description):
    if not isinstance(mapping, dict):
        raise ValueError("{} must be a JSON object".format(description))
    actual = set(mapping)
    expected = set(required)
    if actual != expected:
        raise ValueError(
            "{} has invalid fields (missing={!r}, unexpected={!r})".format(
                description,
                sorted(expected - actual),
                sorted(actual - expected),
            )
        )


def _validate_manifest(document):
    _require_exact_keys(
        document,
        (
            "architecture",
            "artifact",
            "content_sha256",
            "format",
            "format_version",
            "metadata",
            "preprocessing",
            "tensors",
        ),
        "portable strand manifest",
    )
    if document["format"] != FORMAT_NAME:
        raise ValueError("unsupported portable strand format")
    if document["format_version"] != FORMAT_VERSION:
        raise ValueError(
            "unsupported portable strand format version: {!r}".format(
                document["format_version"]
            )
        )

    expected_content_hash = document["content_sha256"]
    if (
            not isinstance(expected_content_hash, str)
            or len(expected_content_hash) != 64):
        raise ValueError("portable strand content_sha256 is invalid")
    content = dict(document)
    del content["content_sha256"]
    if _sha256(_canonical_json_bytes(content)) != expected_content_hash:
        raise ValueError("portable strand manifest content checksum mismatch")

    architecture = document["architecture"]
    _require_exact_keys(
        architecture,
        (
            "batch_normalization_epsilon",
            "dense_units",
            "filters_per_branch",
            "input_shape",
            "kernel_widths",
            "operation_order",
        ),
        "portable strand architecture",
    )
    expected_architecture = {
        "batch_normalization_epsilon": float(BATCH_NORMALIZATION_EPSILON),
        "dense_units": [256, 32, 1],
        "filters_per_branch": FILTERS_PER_BRANCH,
        "input_shape": [len(ALPHABET_ROWS), SEQUENCE_LENGTH, 1],
        "kernel_widths": list(KERNEL_WIDTHS),
        "operation_order": [
            "conv2d_valid_no_bias",
            "relu",
            "batch_normalization_inference",
            "global_max_pool",
            "concatenate",
            "dense_relu_256",
            "dense_relu_32",
            "dense_sigmoid_1",
        ],
    }
    if architecture != expected_architecture:
        raise ValueError("portable strand architecture does not match runtime")

    preprocessing = document["preprocessing"]
    _require_exact_keys(
        preprocessing,
        (
            "alphabet_rows",
            "decision",
            "noncanonical",
            "padding",
            "reverse_complement_before_truncation",
            "sequence_length",
            "truncation",
        ),
        "portable strand preprocessing",
    )
    expected_preprocessing = {
        "alphabet_rows": list(ALPHABET_ROWS),
        "decision": "forward_probability > reverse_probability; ties reverse",
        "noncanonical": "all-zero N encoding",
        "padding": "right-zero",
        "reverse_complement_before_truncation": True,
        "sequence_length": SEQUENCE_LENGTH,
        "truncation": "keep-leftmost",
    }
    if preprocessing != expected_preprocessing:
        raise ValueError("portable strand preprocessing does not match runtime")

    artifact = document["artifact"]
    _require_exact_keys(
        artifact,
        ("file", "file_sha256", "source_file", "source_sha256"),
        "portable strand artifact",
    )
    filename = artifact["file"]
    if (
            not isinstance(filename, str)
            or not filename
            or os.path.basename(filename) != filename):
        raise ValueError("portable strand artifact filename is unsafe")
    for field in ("file_sha256", "source_sha256"):
        digest = artifact[field]
        if not isinstance(digest, str) or len(digest) != 64:
            raise ValueError(
                "portable strand artifact {} is invalid".format(field)
            )

    tensors = document["tensors"]
    _require_exact_keys(
        tensors, EXPECTED_TENSOR_SHAPES, "portable strand tensor manifest"
    )
    for tensor_name, expected_shape in EXPECTED_TENSOR_SHAPES.items():
        tensor = tensors[tensor_name]
        _require_exact_keys(
            tensor, ("dtype", "sha256", "shape"), tensor_name
        )
        if tensor["dtype"] != "float32":
            raise ValueError("{} must use float32".format(tensor_name))
        if tensor["shape"] != list(expected_shape):
            raise ValueError("{} has an invalid shape".format(tensor_name))
        digest = tensor["sha256"]
        if not isinstance(digest, str) or len(digest) != 64:
            raise ValueError("{} has an invalid checksum".format(tensor_name))

    if not isinstance(document["metadata"], dict):
        raise ValueError("portable strand metadata must be a JSON object")
    return document


class PortableStrandModel(object):
    """Inference-only NumPy implementation of CRISPRstrand ``model_r``."""

    def __init__(self, tensors, metadata=None):
        self._tensors = dict(tensors)
        self.metadata = dict(metadata or {})

    def _branch_output(self, inputs, branch_index, kernel_width):
        kernel = self._tensors[
            "conv_{}_kernel".format(branch_index)
        ][:, :, 0, :]
        position_count = SEQUENCE_LENGTH - kernel_width + 1
        convolution = np.zeros(
            (inputs.shape[0], position_count, FILTERS_PER_BRANCH),
            dtype=np.float32,
        )

        # Keras Conv2D performs cross-correlation, not a flipped convolution.
        # Accumulating one kernel column at a time also keeps the numerical
        # path within float32 and closely matches TensorFlow's CPU inference.
        for offset in range(kernel_width):
            columns = inputs[
                :, :, offset:offset + position_count
            ].transpose(0, 2, 1)
            convolution += np.einsum(
                "bpc,cf->bpf",
                columns,
                kernel[:, offset, :],
                optimize=False,
            )

        np.maximum(convolution, np.float32(0.0), out=convolution)
        gamma = self._tensors["bn_{}_gamma".format(branch_index)]
        beta = self._tensors["bn_{}_beta".format(branch_index)]
        mean = self._tensors["bn_{}_moving_mean".format(branch_index)]
        variance = self._tensors[
            "bn_{}_moving_variance".format(branch_index)
        ]
        convolution -= mean
        convolution *= gamma
        convolution /= np.sqrt(variance + BATCH_NORMALIZATION_EPSILON)
        convolution += beta
        return np.max(convolution, axis=1)

    def predict_encoded(self, encoded):
        """Predict orientation-class probabilities for encoded repeats."""

        matrix = np.asarray(encoded, dtype=np.float32)
        expected_tail = (len(ALPHABET_ROWS), SEQUENCE_LENGTH, 1)
        if matrix.ndim != 4 or matrix.shape[1:] != expected_tail:
            raise ValueError(
                "encoded repeats must have shape (n, 4, 38, 1)"
            )
        if matrix.shape[0] < 1:
            raise ValueError("at least one encoded repeat is required")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("encoded repeats must contain finite values")
        inputs = np.ascontiguousarray(matrix[:, :, :, 0])

        pooled = [
            self._branch_output(inputs, branch_index, kernel_width)
            for branch_index, kernel_width in enumerate(KERNEL_WIDTHS, 1)
        ]
        hidden = np.concatenate(pooled, axis=1)
        hidden = np.dot(
            hidden, self._tensors["dense_1_kernel"]
        ) + self._tensors["dense_1_bias"]
        np.maximum(hidden, np.float32(0.0), out=hidden)
        hidden = np.dot(
            hidden, self._tensors["dense_2_kernel"]
        ) + self._tensors["dense_2_bias"]
        np.maximum(hidden, np.float32(0.0), out=hidden)
        logits = np.dot(
            hidden, self._tensors["dense_3_kernel"]
        ) + self._tensors["dense_3_bias"]
        probabilities = np.float32(1.0) / (
            np.float32(1.0) + np.exp(-logits)
        )
        return probabilities[:, 0]

    def predict_probabilities(self, sequences):
        """Predict probabilities for repeats without comparing orientations."""

        _, encoded = encode_repeats(sequences)
        return self.predict_encoded(encoded)

    def predict_orientation_probabilities(self, sequences):
        """Return ``(forward, reverse)`` probability vectors.

        Reverse complements are formed from the complete normalized repeats,
        then both orientations are independently left-truncated to 38 bases.
        """

        if isinstance(sequences, str):
            sequences = [sequences]
        else:
            sequences = list(sequences)
        normalized = [normalize_repeat(sequence) for sequence in sequences]
        if not normalized:
            raise ValueError("at least one repeat sequence is required")
        reversed_repeats = [
            reverse_complement(sequence) for sequence in normalized
        ]
        probabilities = self.predict_probabilities(
            normalized + reversed_repeats
        )
        count = len(normalized)
        return probabilities[:count], probabilities[count:]

    def predict_strands(self, sequences):
        """Return legacy ``Forward``/``Reversed`` labels for repeats."""

        forward, reverse = self.predict_orientation_probabilities(sequences)
        return [
            FORWARD_LABEL if forward_probability > reverse_probability
            else REVERSE_LABEL
            for forward_probability, reverse_probability in zip(
                forward, reverse
            )
        ]


def load_portable_strand_model(manifest_path=DEFAULT_MANIFEST_PATH):
    """Load and fully validate a portable CRISPRstrand artifact."""

    manifest_path = os.path.abspath(os.fspath(manifest_path))
    with open(manifest_path, "rb") as handle:
        raw_manifest = handle.read()
    try:
        document = json.loads(raw_manifest.decode("ascii"))
    except (UnicodeDecodeError, ValueError) as error:
        raise ValueError("portable strand manifest is not valid JSON") from error
    _validate_manifest(document)

    artifact_path = os.path.join(
        os.path.dirname(manifest_path), document["artifact"]["file"]
    )
    with open(artifact_path, "rb") as handle:
        raw_artifact = handle.read()
    if _sha256(raw_artifact) != document["artifact"]["file_sha256"]:
        raise ValueError("portable strand artifact file checksum mismatch")

    tensors = {}
    try:
        archive_context = np.load(io.BytesIO(raw_artifact), allow_pickle=False)
        with archive_context as archive:
            if set(archive.files) != set(EXPECTED_TENSOR_SHAPES):
                raise ValueError(
                    "portable strand artifact contains unexpected tensors"
                )
            for tensor_name, expected_shape in EXPECTED_TENSOR_SHAPES.items():
                tensor = np.asarray(archive[tensor_name])
                if tensor.dtype != np.dtype("float32"):
                    raise ValueError(
                        "{} must use float32".format(tensor_name)
                    )
                if tensor.shape != expected_shape:
                    raise ValueError(
                        "{} has an invalid shape".format(tensor_name)
                    )
                if not np.all(np.isfinite(tensor)):
                    raise ValueError(
                        "{} contains a non-finite value".format(tensor_name)
                    )
                if (
                        _sha256(tensor.tobytes(order="C"))
                        != document["tensors"][tensor_name]["sha256"]):
                    raise ValueError(
                        "{} checksum mismatch".format(tensor_name)
                    )
                tensor = np.array(tensor, dtype=np.float32, copy=True)
                tensor.setflags(write=False)
                tensors[tensor_name] = tensor
    except (OSError, ValueError, KeyError) as error:
        if isinstance(error, ValueError) and str(error).startswith(
                "portable strand"):
            raise
        raise ValueError("portable strand artifact is invalid") from error

    metadata = dict(document["metadata"])
    metadata.update(
        {
            "artifact_file": artifact_path,
            "artifact_sha256": document["artifact"]["file_sha256"],
            "source_sha256": document["artifact"]["source_sha256"],
        }
    )
    return PortableStrandModel(tensors, metadata=metadata)


@functools.lru_cache(maxsize=1)
def load_default_portable_strand_model():
    """Load the bundled model once per Python process."""

    return load_portable_strand_model(DEFAULT_MANIFEST_PATH)

