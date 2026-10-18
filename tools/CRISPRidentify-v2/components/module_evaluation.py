"""Compact-model candidate evaluation and category assignment."""

from components.compact_ml_runtime import evaluate_compact_model
from components.components_detection_refinement import AdvancedFuzzySearchFilter
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.feature_schema import validate_compact_feature_schema
from components.fast_run_rescue import candidate_source_interval
from components.initial_competition import (
    clear_competition_evidence,
    initial_competition_groups,
    with_competition_evidence,
)


def resolve_competitive_component_split(
        *, fast_run, requested, fast_run_rescue):
    """Resolve the source-component competition policy fail closed."""

    values = {
        "fast_run": fast_run,
        "requested": requested,
        "fast_run_rescue": fast_run_rescue,
    }
    if any(not isinstance(value, bool) for value in values.values()):
        raise TypeError(
            "competitive-component policy inputs must be boolean"
        )
    if fast_run and fast_run_rescue and not requested:
        raise ValueError(
            "fast-run rescue requires competitive-component splitting"
        )
    return fast_run and requested


def build_advanced_fuzzy_search_filter(parameters):
    """Build the canonical final structural filter from runtime parameters."""

    required = (
        "param_min_avg_repeat_length",
        "param_max_avg_repeat_length",
        "param_min_avg_spacer_length",
        "param_max_avg_spacer_length",
        "param_min_repeats",
        "param_max_identical_spacers",
        "param_max_identical_cluster_spacers",
    )
    missing = tuple(name for name in required if name not in parameters)
    if missing:
        raise ValueError(
            "structural filter parameters are incomplete: {}".format(
                ", ".join(missing)
            )
        )
    return AdvancedFuzzySearchFilter(
        min_column_dominance_repeat=0.6,
        max_spacer_length=140,
        max_column_dominance_spacer=0.8,
        max_allowed_consecutive_spacers=parameters[
            "param_max_identical_cluster_spacers"
        ],
        max_allowed_same_spacers=parameters[
            "param_max_identical_spacers"
        ],
        max_inconsistent_columns=5,
        min_avg_repeat_length=parameters[
            "param_min_avg_repeat_length"
        ],
        max_avg_repeat_length=parameters[
            "param_max_avg_repeat_length"
        ],
        min_avg_spacer_length=parameters[
            "param_min_avg_spacer_length"
        ],
        max_avg_spacer_length=parameters[
            "param_max_avg_spacer_length"
        ],
        min_repeats=parameters["param_min_repeats"],
    )


class ArrayEvaluation:
    """Score candidates with the one bundled model and assign score bands."""

    def __init__(
            self,
            dict_crispr_array_candidates,
            parameters,
            flag_dev_mode=False,
            precomputed_scored_rows=None):
        self.dict_crispr_array_candidates = dict_crispr_array_candidates
        self.parameters = parameters
        self.flag_dev_mode = flag_dev_mode
        self.feature_schema = validate_compact_feature_schema(
            COMPACT_FEATURE_SCHEMA
        )
        self.split_competitive_components = bool(
            parameters.get("param_split_competitive_components", False)
        )
        self.initial_competition = (
            parameters.get("param_reconstruction_policy", "legacy")
            == "conservative"
        )

        self.dict_scored_result = {}
        self.dict_scored_result_with_all_vectors = {}
        self.dict_bona_fide = {}
        self.dict_alternative = {}
        self.dict_possible = {}
        self.dict_possible_discarded = {}
        self.dict_low_score = {}
        self._load_filter()

        if precomputed_scored_rows is not None:
            self.dict_scored_result = precomputed_scored_rows
        else:
            evaluate_compact_model(self)
        self._split_into_categories()

    @classmethod
    def categorize_scored_rows(
            cls, scored_rows, parameters, flag_dev_mode=False):
        """Apply structural filtering to already scored candidate rows."""

        return cls(
            dict_crispr_array_candidates={},
            parameters=parameters,
            flag_dev_mode=flag_dev_mode,
            precomputed_scored_rows=scored_rows,
        ).output()

    def scored_rows(self):
        return self.dict_scored_result

    def _load_filter(self):
        self.afsf = build_advanced_fuzzy_search_filter(self.parameters)

    @staticmethod
    def _classification_eligible(candidate):
        eligible = getattr(
            candidate, "source_path_classification_eligible", True
        )
        if not isinstance(eligible, bool):
            raise ValueError(
                "source path classification eligibility must be boolean"
            )
        return eligible

    @staticmethod
    def _strict_overlap_components(rows):
        indexed = sorted(
            (
                candidate_source_interval(row[1]),
                row_index,
                row,
            )
            for row_index, row in enumerate(rows)
        )
        if not indexed:
            return []
        components = []
        current = [indexed[0]]
        current_end = indexed[0][0][1]
        for record in indexed[1:]:
            if record[0][0] < current_end:
                current.append(record)
                current_end = max(current_end, record[0][1])
            else:
                components.append(current)
                current = [record]
                current_end = record[0][1]
        components.append(current)
        return [
            [record[2] for record in sorted(
                component, key=lambda item: item[1]
            )]
            for component in components
        ]

    def _retain_competitive_component(self, key, rows):
        high = [row for row in rows if row[0] >= 0.75]
        possible = [row for row in rows if row[0] < 0.75]
        high.sort(key=lambda row: row[0], reverse=True)
        possible.sort(key=lambda row: row[0], reverse=True)

        if high:
            self.dict_bona_fide.setdefault(key, []).append(high[0])
            alternatives = high[1:] + [
                row for row in possible if row[0] >= 0.6
            ]
            alternatives.sort(key=lambda row: row[0], reverse=True)
            if alternatives:
                self.dict_alternative.setdefault(key, []).extend(alternatives)
            discarded = [row for row in possible if row[0] < 0.6]
            if discarded:
                self.dict_possible_discarded.setdefault(key, []).extend(
                    discarded
                )
        elif possible:
            self.dict_possible.setdefault(key, []).append(possible[0])
            if len(possible) > 1:
                self.dict_possible_discarded.setdefault(key, []).extend(
                    possible[1:]
                )

    def _split_into_categories(self):
        for key, original_data in self.dict_scored_result.items():
            if self.initial_competition:
                # A previous classification snapshot cannot describe a fresh
                # ranking. Copy only annotated row containers, never candidates.
                original_data = [clear_competition_evidence(row)
                                 for row in original_data]
            lineage_only = [
                row for row in original_data
                if not self._classification_eligible(row[1])
            ]
            data = [
                row for row in original_data
                if self._classification_eligible(row[1])
            ]
            high = [row for row in data if row[0] >= 0.75]
            possible = [row for row in data if 0.5 <= row[0] < 0.75]
            low = [row for row in data if row[0] < 0.5] + lineage_only

            valid_high = []
            for row in high:
                if self.afsf(row[1]):
                    valid_high.append(row)
                else:
                    low.append(row)

            valid_possible = []
            for row in possible:
                if self.afsf(row[1]):
                    valid_possible.append(row)
                else:
                    low.append(row)

            competitive = valid_high + valid_possible
            conflicts = ()
            if self.initial_competition:
                components = []
                groups, conflicts = initial_competition_groups(
                    competitive, self.split_competitive_components)
                high_owners = set()
                for group in groups:
                    rows = list(group.rows)
                    winner = max(range(len(rows)), key=lambda index: rows[index][0])
                    rows[winner] = with_competition_evidence(rows[winner], group.evidence)
                    if rows[winner][0] >= .75:
                        high_owners.update(owner["component_id"] for owner in group.evidence["ordinary_components"])
                    components.append(rows)
            elif self.split_competitive_components:
                components = self._strict_overlap_components(competitive)
            else:
                components = [competitive] if competitive else []
            for component in components:
                self._retain_competitive_component(key, component)
            for row in conflicts:
                evidence = row[3]["initial_competition"]
                owners = {owner["component_id"] for owner in evidence["ordinary_components"]}
                destination = (self.dict_alternative if row[0] >= .75 or
                               (row[0] >= .6 and owners & high_owners)
                               else self.dict_possible_discarded)
                destination.setdefault(key, []).append(row)

            if low:
                self.dict_low_score[key] = low

    def output(self):
        return [
            self.dict_bona_fide,
            self.dict_alternative,
            self.dict_possible,
            self.dict_possible_discarded,
            self.dict_low_score,
        ]
