import json
import logging

import numpy as np
import pytest
from Bio.Phylo.Newick import Clade

from sp_model.reconstruction_tree import ReconstructionTree


def _reconstructed_tree(tmp_path):
    leaf_a = Clade(name="leaf_a")
    leaf_b = Clade(name="leaf_b")
    root = Clade(name="root", clades=[leaf_a, leaf_b])
    reconstruction = ReconstructionTree(
        {"model": "gtr", "gain_rate": 0.1, "loss_rate": 1.0},
        str(tmp_path),
        root,
        True,
        save_reconstructed_events=True,
        logger=logging.getLogger("reconstructed-event-export-test"),
    )

    root.rec_spacers = np.array([1, 1])
    leaf_a.rec_spacers = np.array([1, 0])
    leaf_b.rec_spacers = np.array([0, 1])
    reconstruction.top_order = ["s1", "s2"]
    reconstruction.spacer_names_to_numbers = {"s1": 1, "s2": 2}
    reconstruction.spacer_numbers_to_names = {1: "s1", 2: "s2"}
    reconstruction.rec_gain_dict = {"root": [1, 2], "leaf_a": [], "leaf_b": []}
    reconstruction.rec_loss_dict = {"root": [], "leaf_a": [2], "leaf_b": [1]}
    reconstruction.rec_joint_losses_dict = dict(reconstruction.rec_loss_dict)
    return reconstruction


def test_reconstructed_events_are_exported_without_plotting(tmp_path, monkeypatch):
    reconstruction = _reconstructed_tree(tmp_path)

    def fail_if_order_plot_is_attempted(*args, **kwargs):
        raise AssertionError("plotting must not run when plot_tree is false")

    monkeypatch.setattr(reconstruction, "visualize_order", fail_if_order_plot_is_attempted)

    result = reconstruction.visualize_tree(
        name="group_rec",
        do_show=False,
        path=str(tmp_path),
        sort_top_order="top_order",
        plot_tree=False,
    )

    export_path = tmp_path / "reconstructed_events" / "group_rec"
    assert result is None
    assert {path.name for path in export_path.iterdir()} == {
        "group_rec_other_events.json",
        "group_rec_rec_gains_losses.json",
        "group_rec_rec_spacers.json",
        "group_rec_spacer_names_to_numbers.json",
        "group_rec_top_order.json",
        "group_rec_tree.nwk",
    }
    gains_losses = json.loads((export_path / "group_rec_rec_gains_losses.json").read_text())
    assert gains_losses["rec_gains"]["root"] == [1, 2]
    assert gains_losses["rec_losses"]["leaf_a"] == [2]
    assert json.loads((export_path / "group_rec_top_order.json").read_text()) == ["s1", "s2"]
    assert not list(tmp_path.rglob("*.pdf"))
    assert not list(tmp_path.rglob("*.dot"))


def test_duplicate_candidates_are_excluded_from_independent_acquisition_count(tmp_path):
    reconstruction = _reconstructed_tree(tmp_path)
    reconstruction.duplicated_spacers = {"duplicate": {1, 2}}
    reconstruction.rec_gain_dict = {
        "root": [1, 3], "leaf_a": [1, 2, 3], "leaf_b": [2, 3],
    }
    reconstruction.rec_duplications_dict = {"root": [1]}
    reconstruction.rec_rearrangements_dict = {"leaf_a": [2]}
    reconstruction.rec_reacquisition_dict = {"leaf_b": [1]}
    reconstruction.rec_indep_gain_dict = {"leaf_b": [2]}
    reconstruction.rec_other_dup_events_dict = {"root": [2]}

    # Spacer 3 was acquired at three nodes: two independent reacquisitions.
    # Repeated gains of duplicate candidates 1 and 2 belong to separate counts.
    assert reconstruction.get_rec_duplication_type_counts() == (1, 1, 1, 1, 1, 2)
    assert reconstruction.rec_gain_dict["root"] == [1, 3]


def test_duplicate_event_accounting_preserves_unknown_candidate_error(tmp_path):
    reconstruction = _reconstructed_tree(tmp_path)
    reconstruction.duplicated_spacers = {"duplicate": {1}}
    reconstruction.rec_duplications_dict = {"root": [2]}
    with pytest.raises(KeyError, match="2"):
        reconstruction.get_rec_duplication_type_counts()
