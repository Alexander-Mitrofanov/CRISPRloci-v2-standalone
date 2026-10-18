import logging

import numpy as np
import pandas as pd
import pytest
from scipy.stats import poisson
from scipy.optimize import OptimizeResult

from sp_model import model_tools
import sp_model.reconstruction_tree as reconstruction_tree_module
from sp_model.data_classes.crisprdata import CRISPRArray, CRISPRGroup
from sp_model.distance_tree_constructor import LikelihoodDistance
from sp_model.experiments_tools import integer_group_labels, pooling_for_parameter_estimation
from sp_model.helpers.misc import RunTimer, create_logger, remove_completely_same_arrays
from sp_model.helpers.stats import np_rate_exponential, test_significance_ratio_chi2 as chi_square_significance
from sp_model.orientation.orientation_tools import compare_likelihoods_for_orientation, decision_fct
from sp_model.run_experiments import _require_reconstructed_protocols
from sp_model.reconstruction_tree import ReconstructionTree, RhoBiasCorrectionError
from sp_model.run_reconstruction import _compute_lh_ratio_with_bias_policy


def _protocol(lh_idm, ln_lh_idm, lh_bdm, ln_lh_bdm):
    return pd.DataFrame(
        {"lh_idm": [lh_idm], "ln_lh_idm": [ln_lh_idm], "lh_bdm": [lh_bdm], "ln_lh_bdm": [ln_lh_bdm]},
        index=["group"],
    )


def test_reverse_orientation_uses_reverse_metrics_in_selected_protocol():
    forward = _protocol(0.001, -7.0, 0.0001, -10.0)
    reverse = _protocol(0.01, -4.0, 0.006, -5.0)
    comparison, oriented, _, _ = compare_likelihoods_for_orientation(forward, reverse, decision_boundary=1)
    assert comparison.loc["group", "recommend reversing array"]
    assert oriented.loc["group", "predicted orientation"] == "Reverse"
    assert oriented.loc["group", "ln_lh_bdm"] == -5.0
    assert oriented.loc["group", "lh_bdm"] == 0.006


def test_orientation_merge_preserves_forward_reverse_and_selected_bias_provenance():
    forward = _protocol(0.001, -7.0, 0.0001, -10.0)
    reverse = _protocol(0.01, -4.0, 0.006, -5.0)
    for field, forward_value, reverse_value in (
        ("alpha_bias_correction_status", "applied", "applied"),
        ("rho_bias_correction_status", "fallback_uncorrected", "applied"),
        ("rho_bias_failure_policy", "fallback_uncorrected", "fallback_uncorrected"),
        (
            "rho_bias_correction_warning_code",
            "rho_bias_fit_unavailable_used_uncorrected",
            "",
        ),
        ("rho_bias_correction_warning", "forward warning", ""),
    ):
        forward[field] = forward_value
        reverse[field] = reverse_value

    comparison, oriented, _, _ = compare_likelihoods_for_orientation(
        forward, reverse, decision_boundary=1
    )

    assert comparison.loc["group", "rho_bias_correction_status"] == "fallback_uncorrected"
    assert comparison.loc["group", "reversed_rho_bias_correction_status"] == "applied"
    assert oriented.loc["group", "predicted orientation"] == "Reverse"
    assert oriented.loc["group", "rho_bias_correction_status"] == "applied"
    assert oriented.loc["group", "rho_bias_correction_warning"] == ""


def test_orientation_validation_rejects_negative_boundary_and_index_mismatch():
    protocol = _protocol(1, 0, 1, 0)
    with pytest.raises(ValueError, match="non-negative"):
        compare_likelihoods_for_orientation(protocol, protocol, decision_boundary=-1)
    other = protocol.rename(index={"group": "other"})
    with pytest.raises(ValueError, match="same groups"):
        compare_likelihoods_for_orientation(protocol, other)


def test_orientation_aligns_reverse_rows_by_group_and_validates_tree_maps():
    forward = pd.concat(
        [
            _protocol(0.1, -2.0, 0.2, -1.0),
            _protocol(0.01, -4.0, 0.02, -3.0).rename(index={"group": "other"}),
        ]
    )
    reverse = pd.concat(
        [
            _protocol(0.005, -5.0, 0.01, -4.0).rename(index={"group": "other"}),
            _protocol(0.09, -2.1, 0.19, -1.1),
        ]
    )
    comparison, *_ = compare_likelihoods_for_orientation(forward, reverse)
    assert comparison.index.tolist() == ["group", "other"]
    assert comparison.loc["group", "reversed_ln_lh_bdm"] == -1.1
    assert comparison.loc["other", "reversed_ln_lh_bdm"] == -4.0

    with pytest.raises(ValueError, match="forward and reversed tree maps"):
        compare_likelihoods_for_orientation(forward, reverse, dict_trees={})
    with pytest.raises(ValueError, match="every reconstructed group"):
        compare_likelihoods_for_orientation(
            forward,
            reverse,
            dict_trees={"forward": {"group": "tree"}, "reversed": {"group": "tree"}},
        )


def test_orientation_ratio_is_derived_from_logs_after_raw_likelihood_underflow():
    forward = _protocol(0.0, -1000.0, 0.0, -1000.0)
    reverse = _protocol(0.0, -1001.0, 0.0, -1001.0)
    comparison, _, _, _ = compare_likelihoods_for_orientation(forward, reverse, decision_boundary=5)
    assert comparison.loc["group", "lh_idm / reversed_lh_idm"] == pytest.approx(np.e)
    assert comparison.loc["group", "lh_bdm / reversed_lh_bdm"] == pytest.approx(np.e)


@pytest.mark.parametrize(("difference", "boundary"), [(np.nan, 1), (0, np.nan), (0, -1)])
def test_orientation_decision_rejects_non_finite_or_negative_inputs(difference, boundary):
    with pytest.raises(ValueError):
        decision_fct(difference, boundary)


def test_not_determined_orientation_is_not_treated_as_truthy_reverse():
    forward = _protocol(0.1, -2.0, 0.2, -1.0)
    reverse = _protocol(0.09, -2.1, 0.19, -1.1)
    _, oriented, _, _ = compare_likelihoods_for_orientation(forward, reverse, decision_boundary=1)
    assert oriented.loc["group", "predicted orientation"] == "ND"
    assert oriented.loc["group", "ln_lh_bdm"] == -1.0


def test_direct_poisson_log_likelihood_matches_scipy_distribution():
    calculator = LikelihoodDistance(gain_rate=2.5, loss_rate=0.7, alpha=1.4)
    times = [0.01, 0.5, 2.0]
    counts = [0, 2, 8]
    assert np.allclose(
        calculator._gain_lh(times, counts),
        [poisson.logpmf(count, time * 2.5 / (0.7 * 1.4 + 1e-40)) for time, count in zip(times, counts)],
    )


def test_iterative_block_loss_matches_reference_recurrence_and_handles_long_inputs():
    def reference(loss_probability, geometric):
        if not geometric:
            return []
        previous = reference(loss_probability, geometric[:-1])
        convolution = [geometric[index] * previous[-(index + 1)] for index in range(len(geometric) - 1)]
        return previous + [loss_probability * (geometric[-1] + sum(convolution))]

    geometric = [model_tools.geom(1.7, index) for index in range(1, 30)]
    assert np.allclose(model_tools.mult_loss_lh(0.2, 1.7, geometric), reference(0.2, geometric))
    assert len(model_tools.mult_loss_lh(0.1, 1.2, [0.01] * 1200)) == 1200


@pytest.mark.parametrize('loss_probability', [0.0, 1e-8, 0.2, 0.9, 1.0])
@pytest.mark.parametrize('alpha', [1.0, 1.01, 1.7, 10.0])
def test_closed_form_block_loss_matches_original_recurrence(loss_probability, alpha):
    geometric = [model_tools.geom(alpha, index) for index in range(1, 80)]
    expected = np.maximum(model_tools.mult_loss_lh(loss_probability, alpha, geometric), model_tools.EPS)

    actual = model_tools.eff_mult_loss_lh(loss_probability, alpha, len(geometric), log=False)
    actual_log = model_tools.eff_mult_loss_lh(loss_probability, alpha, len(geometric), log=True)

    assert np.allclose(actual, expected, rtol=1e-12, atol=model_tools.EPS)
    assert np.allclose(actual_log, np.log(expected), rtol=1e-12, atol=1e-12)

    block_lengths = [1, 2, 5, 29, 79]
    assert model_tools.sum_log_block_likelihoods(
        loss_probability,
        alpha,
        block_lengths,
    ) == pytest.approx(sum(np.log(expected[length - 1]) for length in block_lengths), abs=1e-12)


def test_closed_form_block_loss_validates_length_and_handles_empty_input():
    assert model_tools.eff_mult_loss_lh(0.2, 1.5, 0).size == 0
    with pytest.raises(ValueError, match='non-negative'):
        model_tools.eff_mult_loss_lh(0.2, 1.5, -1)


def test_long_block_likelihood_adds_full_blocks_and_remainder():
    functions = [
        lambda branch, loss, alpha: 0.9,
        lambda branch, loss, alpha: 0.5,
        lambda branch, loss, alpha: 0.25,
    ]
    actual = model_tools.compute_tree_lh_for_given_lh_fct(
        0.1,
        1.5,
        [1.0],
        [0],
        [[5]],
        lambdifyed_lh_fct=functions,
    )
    assert actual == pytest.approx(2 * np.log(0.25) + np.log(0.5))


def test_failed_or_non_finite_optimizer_results_are_rejected():
    failed = OptimizeResult(success=False, fun=np.nan, x=np.array([np.inf]), message="bad fit")
    with pytest.raises(RuntimeError, match="bad fit"):
        model_tools.validate_optimization_result(failed, "test-method")


def test_rho_bias_fallback_is_explicit_per_fit_and_strict_mode_still_fails(caplog):
    class SparseFit:
        def __init__(self):
            self.calls = []

        def compute_lh_ratio(self, **kwargs):
            self.calls.append(kwargs)
            if kwargs["rho_bias_correction"]:
                raise RhoBiasCorrectionError("non-finite corrected objective")
            return ("lh_idm", "lh_bdm", "idm_fit", "bdm_fit")

    reconstruction = SparseFit()
    result, provenance = _compute_lh_ratio_with_bias_policy(
        reconstruction,
        alpha_bias_correction=True,
        rho_bias_correction=True,
        rho_bias_failure_policy="fallback_uncorrected",
        group_name="sparse_group",
        context="primary",
        logger=logging.getLogger("rho-fallback-test"),
    )

    assert result == ("lh_idm", "lh_bdm", "idm_fit", "bdm_fit")
    assert [call["rho_bias_correction"] for call in reconstruction.calls] == [True, False]
    assert all(call["alpha_bias_correction"] for call in reconstruction.calls)
    assert provenance["rho_bias_correction_status"] == "fallback_uncorrected"
    assert provenance["rho_bias_correction_warning_code"] == (
        "rho_bias_fit_unavailable_used_uncorrected"
    )
    assert "BDM likelihood used for orientation is unchanged" in provenance["rho_bias_correction_warning"]
    assert "uncorrected rho estimates" in caplog.text

    strict = SparseFit()
    with pytest.raises(RhoBiasCorrectionError, match="non-finite corrected objective"):
        _compute_lh_ratio_with_bias_policy(
            strict,
            alpha_bias_correction=True,
            rho_bias_correction=True,
            rho_bias_failure_policy="fail",
            group_name="sparse_group",
            context="primary",
            logger=logging.getLogger("rho-strict-test"),
        )
    assert len(strict.calls) == 1


def test_rho_bias_fallback_does_not_mask_primary_optimizer_failures():
    class PrimaryFailure:
        def compute_lh_ratio(self, **kwargs):
            raise RuntimeError("primary BDM optimizer failed")

    with pytest.raises(RuntimeError, match="primary BDM optimizer failed") as error:
        _compute_lh_ratio_with_bias_policy(
            PrimaryFailure(),
            alpha_bias_correction=True,
            rho_bias_correction=True,
            rho_bias_failure_policy="fallback_uncorrected",
            group_name="group",
            context="primary",
            logger=logging.getLogger("primary-failure-test"),
        )
    assert not isinstance(error.value, RhoBiasCorrectionError)


def test_rho_bias_exception_boundary_wraps_only_the_optional_fit(monkeypatch):
    reconstruction = object.__new__(ReconstructionTree)
    reconstruction.get_data_for_lh_ratio = lambda filter_by=False: ([], [], [])
    reconstruction.gains_as_ls_for_lh_computation = lambda: []
    primary_idm = OptimizeResult(success=True, x=1.0, fun=2.0)
    primary_bdm = OptimizeResult(success=True, x=np.array([1.0, 2.0]), fun=1.0)

    scalar_calls = 0

    def fail_optional_fit(*args, **kwargs):
        nonlocal scalar_calls
        scalar_calls += 1
        if scalar_calls == 1:
            return primary_idm
        raise RuntimeError("non-finite optional corrected objective")

    monkeypatch.setattr(reconstruction_tree_module, "minimize_scalar_fct", fail_optional_fit)
    monkeypatch.setattr(reconstruction_tree_module, "minimize_fct", lambda *args, **kwargs: primary_bdm)
    with pytest.raises(RhoBiasCorrectionError, match="non-finite optional corrected objective"):
        reconstruction.compute_lh_ratio(rho_bias_correction=True)

    def fail_primary_fit(*args, **kwargs):
        raise RuntimeError("primary IDM optimizer failed")

    monkeypatch.setattr(reconstruction_tree_module, "minimize_scalar_fct", fail_primary_fit)
    with pytest.raises(RuntimeError, match="primary IDM optimizer failed") as error:
        reconstruction.compute_lh_ratio(rho_bias_correction=True)
    assert not isinstance(error.value, RhoBiasCorrectionError)


def test_numpy_two_compatibility_and_logger_reconfiguration(tmp_path):
    assert np.isinf(np_rate_exponential(0))
    assert model_tools.log1mexp(0) == -np.inf
    with pytest.raises(ValueError, match="less than or equal to zero"):
        model_tools.log1mexp(0.1)
    first = create_logger("regression-logger", logging.INFO, tmp_path / "first.log")
    assert len(first.handlers) == 2
    second = create_logger("regression-logger", logging.WARNING, tmp_path / "second.log")
    assert second is first
    assert len(second.handlers) == 2
    assert second.propagate is False


@pytest.mark.parametrize(
    ("degrees_of_freedom", "expected_quantile"),
    [(1, 3.8414588206941285), (2, 5.991464547107983)],
)
def test_chi_square_survival_quantile_without_scipy_stats(degrees_of_freedom, expected_quantile):
    _, quantile = chi_square_significance(0, 0.05, n=degrees_of_freedom)
    assert quantile == pytest.approx(expected_quantile)


@pytest.mark.parametrize(("significance", "degrees_of_freedom"), [(0, 1), (1, 1), (0.05, 0)])
def test_chi_square_quantile_rejects_invalid_domain(significance, degrees_of_freedom):
    with pytest.raises(ValueError):
        chi_square_significance(0, significance, n=degrees_of_freedom)


def test_duplicate_array_removal_is_stable_and_validates_lengths():
    arrays, names, duplicates = remove_completely_same_arrays([[1], [2], [1], [1]], ["a", "b", "c", "d"])
    assert arrays == [[1], [2]]
    assert names == ["a", "b"]
    assert duplicates == {"a": ["c", "d"]}
    with pytest.raises(ValueError, match="same length"):
        remove_completely_same_arrays([[1]], [])


def test_crispr_group_cas_info_no_longer_calls_shadowed_attribute():
    array = CRISPRArray(
        "a",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        ["1"],
        "I-A",
        all_cas_types=["I-A"],
        species_fc="species",
    )
    group = CRISPRGroup(None, [array], name="group")
    assert group.get_cas_gene_info() == {"a": ["I-A"]}


def test_mismatch_resolution_supports_non_default_integer_indexes_and_is_deterministic():
    alignment = pd.DataFrame({0: ["a", "b"], 1: ["a", "-"]}, index=[10, 20])

    first = model_tools.get_fix_mismatches_df_alignment(alignment)
    second = model_tools.get_fix_mismatches_df_alignment(alignment)

    pd.testing.assert_frame_equal(first[0], second[0])
    assert first[0].index.tolist() == [10, 20]
    assert first[1] == ["b", "a", "a"]
    assert first[3][1:3] == (1, [("a", "b")])


def test_alignment_helpers_reject_empty_or_inconsistent_inputs():
    with pytest.raises(ValueError, match="empty"):
        model_tools.get_fix_mismatches_df_alignment(pd.DataFrame())
    with pytest.raises(ValueError, match="empty"):
        model_tools.create_df_alignment([], [])
    with pytest.raises(ValueError, match="same length"):
        model_tools.create_df_alignment([["a"]], [])


def test_integer_group_labels_handles_small_remainders_and_validates_domain():
    assert integer_group_labels(0, 5) == []
    assert integer_group_labels(2, 5) == [0, 0]
    assert integer_group_labels(5, 2) == [0, 0, 1, 1, 1]
    assert integer_group_labels(4, np.int64(2)) == [0, 0, 1, 1]
    for invalid_size in (0, -1, True, 1.5):
        with pytest.raises(ValueError, match="positive integer"):
            integer_group_labels(2, invalid_size)


def test_pooling_rejects_missing_or_misaligned_group_labels():
    with pytest.raises(ValueError, match="required"):
        pooling_for_parameter_estimation([], None)
    with pytest.raises(ValueError, match="same length"):
        pooling_for_parameter_estimation([object()], [])


def test_spacer_order_uses_first_occurrence_instead_of_hash_order():
    _, topological_order, (_, labels) = model_tools.determine_order_dict(
        [["a", "b"], ["c", "d"]],
        plot_order=False,
    )
    assert labels == {0: "a", 1: "b", 2: "c", 3: "d"}
    assert topological_order == ["c", "d", "a", "b"]


def test_unobserved_event_filter_does_not_skip_adjacent_events():
    class Tree:
        @staticmethod
        def get_terminals():
            return [type("Leaf", (), {"name": "leaf"})()]

    alignment = pd.DataFrame(
        {0: ["keep", "-"], 1: ["-", "gone"]},
        index=["leaf", "internal"],
    )
    events = {
        "node": [
            ["-", 0, ["gone"], 0, ["gone"]],
            ["-", 0, ["gone"], 0, ["gone"]],
            ["+", 0, ["keep"], 0, ["gone", "keep"]],
        ]
    }
    model_tools.remove_unobserved_spacer(
        alignment,
        Tree(),
        {"leaf": [], "internal": []},
        {"leaf": [], "internal": []},
        events,
        logger=logging.getLogger("event-filter-test"),
    )
    assert len(events["node"]) == 1
    assert events["node"][0][2] == ["keep"]
    assert events["node"][0][4] == ["keep"]


def test_pooling_is_stable_linear_grouping_and_rejects_missing_labels(monkeypatch):
    calls = []

    def fake_likelihood(group_data, **_kwargs):
        calls.append(group_data)
        return 1, 2, object(), object()

    monkeypatch.setattr(model_tools, "compute_lh_ratio_of_multiple_trees", fake_likelihood)
    labels, *_ = pooling_for_parameter_estimation(["a", "b", "c"], ["x", 1, "x"])
    assert labels.tolist() == ["x", 1]
    assert calls == [["a", "c"], ["b"]]
    with pytest.raises(ValueError, match="missing"):
        pooling_for_parameter_estimation(["a"], [np.nan])


def test_likelihood_distance_handles_unequal_arrays_and_rejects_invalid_parameters():
    calculator = LikelihoodDistance(2.5, 0.7, 1.4, skip_letters={"-"})
    fes_position, fes = calculator._find_FES(["a"], ["a", "b"])
    assert (fes_position, fes) == (0, "a")
    assert calculator._find_gains_deletions(["a"], ["a", "b"], fes_position) == (
        [0, 0],
        [1, 2],
        [[1], []],
    )
    for parameters in ((-1, 1, 1), (1, -1, 1), (1, 1, 0.5), (np.nan, 1, 1), (1, 1, None)):
        with pytest.raises(ValueError):
            LikelihoodDistance(*parameters)


def test_logger_creates_nested_parent_directories(tmp_path):
    output = tmp_path / "nested" / "run.log"
    logger = create_logger("nested-logger", logging.INFO, output)
    logger.info("ready")
    assert "ready" in output.read_text(encoding="utf-8")


def test_all_skipped_groups_fail_with_the_skipped_report_path(tmp_path):
    with pytest.raises(ValueError, match="0_protocol_skipped.csv"):
        _require_reconstructed_protocols(pd.DataFrame(), pd.DataFrame(), tmp_path)

    _require_reconstructed_protocols(
        pd.DataFrame(),
        pd.DataFrame({"name": ["trivial"]}),
        tmp_path,
    )


def test_crispr_group_rejects_data_loss_and_handles_partial_metadata():
    def array(name, cas_type, kingdom):
        return CRISPRArray(
            name,
            None,
            None,
            None,
            None,
            None,
            kingdom,
            None,
            ["1"],
            cas_type,
            all_cas_types=[] if cas_type is None else [cas_type],
            species_fc="species",
        )

    with pytest.raises(ValueError, match="at least one array"):
        CRISPRGroup(None, [], name="empty")
    first = array("a", None, None)
    second = array("b", "I-A", "Archaea")
    group = CRISPRGroup(None, [first, second], name="partial")
    assert group.determine_group_cas_type() == ("I-A", 1)
    assert group.kingdom == "Archaea"
    with pytest.raises(ValueError, match="names must be unique"):
        CRISPRGroup(None, [first, array("a", "I-B", "Bacteria")], name="duplicate")


def test_run_timer_reports_non_negative_checkpoint_offsets():
    timer = RunTimer()
    timer.start_time = 10
    timer.checkpoints = [10, 12]
    assert "Checkpoint times (wrt. start): [0, 2]" in str(timer)
