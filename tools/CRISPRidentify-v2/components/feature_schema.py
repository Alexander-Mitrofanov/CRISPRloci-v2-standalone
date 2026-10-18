"""Public feature contract for the compact CRISPR classifier."""

from types import MappingProxyType


COMPACT_FEATURE_SCHEMA = (
    "repeat_len",
    "number_repeats",
    "repeat_similarity",
    "at_richness",
    "avg_spacer_len",
    "spacer_similarity",
    "number_mismatches",
    "spacer_evenness",
    "mfe_score",
    "blast_score_1",
)

FEATURE_DISPLAY_NAMES = MappingProxyType({
    "repeat_len": "Maximum repeat length (nt)",
    "number_repeats": "Repeat count",
    "repeat_similarity": "Repeat sequence similarity",
    "at_richness": "Repeat A/T fraction",
    "avg_spacer_len": "Mean spacer length (nt)",
    "spacer_similarity": "Spacer sequence similarity",
    "number_mismatches": "Repeat-consensus mismatch count",
    "spacer_evenness": "Spacer-length evenness",
    "mfe_score": "Consensus-repeat RNA folding energy (kcal/mol)",
    "blast_score_1": "Similarity to known CRISPR repeats (BLAST bit score)",
})


def feature_display_name(feature_name):
    return FEATURE_DISPLAY_NAMES.get(feature_name, feature_name)


def validate_compact_feature_schema(feature_names):
    supplied = tuple(feature_names)
    if supplied != COMPACT_FEATURE_SCHEMA:
        raise ValueError(
            "compact feature schema must exactly match the trained model"
        )
    return supplied


