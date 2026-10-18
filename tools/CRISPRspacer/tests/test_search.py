import random

import pytest

from crisprspacer.inputs import from_identify, from_spacers
from crisprspacer.io import InputError, LimitError, digest
from crisprspacer.search import Budget, Policy, analyze, find_matches
from conftest import record


def oracle(query, target, k):
    # Independent naive enumeration: no production reverse-complement/window/seed helpers.
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    reverse = "".join(complement[x] for x in reversed(query))
    choices = [(query, "both")] if reverse == query else [(query, "+"), (reverse, "-")]
    text = target["sequence"]
    if len(query) > len(text):
        return []
    count = len(text) if target["topology"] == "circular" else len(text)-len(query)+1
    found = []
    for start in range(count):
        window = "".join(text[(start+i) % len(text)] for i in range(len(query)))
        if set(window) - set("ACGT"):
            continue
        for pattern, strand in choices:
            differences = len([1 for a, b in zip(window, pattern) if a != b])
            if differences <= k:
                found.append((start, strand, differences))
    return sorted(found)


@pytest.mark.parametrize("k", [0, 1, 2])
@pytest.mark.parametrize("topology", ["linear", "circular"])
def test_seed_search_agrees_with_independent_exhaustive_oracle(k, topology):
    rng = random.Random(7249)
    for _ in range(350):
        query = "".join(rng.choice("ACGT") for _ in range(rng.randrange(1, 13)))
        target = record("".join(rng.choice("ACGTN") for _ in range(rng.randrange(1, 50))), topology=topology)
        assert find_matches(query, target, k, Budget(Policy())) == oracle(query, target, k)


def queries(*sequences):
    return from_spacers([record(seq, f"q{i}") for i, seq in enumerate(sequences)], digest(b"fixture"))


def test_palindrome_is_one_both_strand_location():
    result = analyze(queries("ACGT"), [record("ACGTACGT")], Policy())
    assert [(r["target_start"], r["relative_strand"]) for r in result["matches"]] == [(1, "both"), (5, "both")]


def test_duplicate_queries_retain_occurrences_and_deduplicate_computation():
    result = analyze(queries("ACGT", "ACGT"), [record("ACGT")], Policy())
    assert result["counts"]["unique_eligible_sequences"] == 1
    assert result["counts"]["candidate_matches"] == 2
    assert result["coverage"][0]["query_occurrence_count"] == 2


def test_ambiguity_is_not_a_match_or_a_silent_negative():
    result = analyze(queries("ACN", "ACG"), [record("NNN")], Policy())
    assert result["queries"][0]["status"] == "ambiguous_query"
    assert result["outcome"] == "no_assessable_windows"
    assert result["counts"]["ambiguous_windows"] == 1


def test_no_eligible_queries_differs_from_no_matches():
    assert analyze(queries("NN"), [record("ACGT")], Policy())["outcome"] == "no_eligible_queries"
    assert analyze(queries("AAAA"), [record("ACGT")], Policy())["outcome"] == "no_candidate_matches"


def test_circular_segments_are_source_reconstructible():
    result = analyze(queries("TTAA"), [record("AACCTT", topology="circular")], Policy())
    hit = result["matches"][0]
    assert hit["target_segments"] == [{"start": 5, "end": 6}, {"start": 1, "end": 2}]
    assert hit["target_start"] == 5 and hit["target_end"] == 2
    assert hit["wraps_origin"] and hit["relative_strand"] == "both"
    assert hit["target_sequence_forward"] == "TTAA"


def test_linear_flanks_are_clipped_without_invented_bases():
    result = analyze(queries("ACG"), [record("ACGAAA")], Policy(flank=10))
    hit = result["matches"][0]
    assert hit["left_flank_forward"] is None
    assert hit["right_flank_forward"] == "AAA"
    assert hit["pam_status"] == "not_evaluated"


def test_self_excludes_originating_array_keeps_external_same_record(native_fixture):
    source, _, path, _ = native_fixture
    inventory = from_identify([source], [path])
    result = analyze(inventory, [source], Policy(scope="same-record"), True)
    assert result["counts"]["candidate_matches"] == 1
    assert result["counts"]["excluded_array_matches"] == 1
    assert {r["array_overlap"] for r in result["matches"]} == {"source_array", "outside_reported_arrays"}
    assert all(r["same_source_record"] is True for r in result["matches"])
    assert all(r["same_declared_assembly"] is None for r in result["matches"])


def test_all_reported_arrays_excluded_even_unselected_category(native_fixture):
    source, _, path, _ = native_fixture
    inventory = from_identify([source], [path])
    inventory["arrays"].append({"id": "other", "native_id": "low", "source_id": source["id"],
                               "start": 15, "end": 18, "category": "Low score", "orientation": None})
    result = analyze(inventory, [source], Policy(scope="same-record"), True)
    assert result["counts"]["candidate_matches"] == 0
    assert {r["array_overlap"] for r in result["matches"]} == {"source_array", "other_array"}


def test_declared_assembly_scope_and_unknown_annotation(native_fixture):
    source, _, path, _ = native_fixture
    source["assembly_id"] = "assembly-A"
    second = record("GATT", "second", assembly="assembly-A")
    other = record("GATT", "other", assembly="assembly-B")
    inventory = from_identify([source, second, other], [path])
    result = analyze(inventory, inventory["sources"], Policy(scope="same-assembly"), True)
    second_hit = next(x for x in result["matches"] if x["target_id"] == "second")
    assert second_hit["same_source_record"] is False
    assert second_hit["same_declared_assembly"] is True
    assert second_hit["array_overlap"] == "unknown"
    assert not any(x["target_id"] == "other" for x in result["matches"])


def test_same_upload_does_not_imply_assembly(native_fixture):
    source, _, path, _ = native_fixture
    inventory = from_identify([source], [path])
    with pytest.raises(InputError, match="assembly_id"):
        analyze(inventory, [source], Policy(scope="same-assembly"), True)


def test_external_matching_id_and_hash_do_not_imply_source_identity(native_fixture):
    source, _, path, _ = native_fixture
    result = analyze(from_identify([source], [path]), [source], Policy())
    assert all(x["same_source_record"] is None and x["array_overlap"] == "unknown" for x in result["matches"])


@pytest.mark.parametrize("policy", [Policy(max_hits=1), Policy(max_candidates=1), Policy(max_work=1)])
def test_limits_abort_instead_of_truncating(policy):
    with pytest.raises(LimitError):
        analyze(queries("AAA"), [record("AAAAAAAAAA")], policy)


def test_deadline_checked_during_discovery():
    budget = Budget(Policy())
    budget.deadline = 0
    with pytest.raises(LimitError, match="deadline"):
        find_matches("AAA", record("AAAA"), 0, budget)


def test_deterministic_replay():
    inventory, targets, policy = queries("GATT", "ACGT"), [record("GATTACGT")], Policy(substitutions=1)
    assert analyze(inventory, targets, policy) == analyze(inventory, targets, policy)


@pytest.mark.parametrize("values", [{"substitutions": True}, {"max_hits": 0}, {"timeout": 0},
                                     {"min_length": 10, "max_length": 2}, {"scope": "upload"}])
def test_invalid_policies(values):
    with pytest.raises(InputError):
        Policy(**values).validate()


@pytest.mark.parametrize("policy", [Policy(max_pairs=2), Policy(max_work=2), Policy(max_evidence_bytes=1)])
def test_zero_window_pairs_are_bounded(policy):
    with pytest.raises(LimitError):
        analyze(queries("A"*20), [record("A", f"t{i}") for i in range(2000)], policy)


def test_zero_window_deadline_checked(monkeypatch):
    import crisprspacer.search as search
    calls = 0
    def clock():
        nonlocal calls
        calls += 1
        return 0 if calls == 1 else 1000
    monkeypatch.setattr(search.time, "monotonic", clock)
    with pytest.raises(LimitError, match="deadline"):
        analyze(queries("A"*20), [record("A")], Policy(timeout=1))


def test_circular_slicing_does_not_copy_full_target():
    from crisprspacer.search import forward_slice
    class GuardedSequence(str):
        def __add__(self, other):
            raise AssertionError("Entire circular target was duplicated")
    target = record("AACCTT", topology="circular")
    target["sequence"] = GuardedSequence(target["sequence"])
    assert forward_slice(target, 4, 4) == "TTAA"


@pytest.mark.parametrize("topology", ["linear", "circular"])
@pytest.mark.parametrize("query", ["AAA", "ACGT", "AGTC"])
def test_exact_search_preserves_overlaps_strands_and_budget(topology, query):
    target = record("AAAAAACGTAGTCGACTACGTNAAA", topology=topology)
    budget = Budget(Policy())
    found = find_matches(query, target, 0, budget)
    assert found == oracle(query, target, 0)
    orientations = 1 if query == "ACGT" else 2
    text_size = target["length"] + (len(query) - 1 if topology == "circular" else 0)
    assert budget.candidates == len(found)
    assert budget.work == orientations * text_size + len(found) * len(query)


def test_exact_discovery_needs_no_target_window_copies():
    class NoWindowCopies(str):
        def __getitem__(self, index):
            if isinstance(index, slice):
                raise AssertionError("A full-query exact seed was redundantly sliced")
            return super().__getitem__(index)

    target = record("AAAANAAAA")
    target["sequence"] = NoWindowCopies(target["sequence"])
    assert find_matches("AAA", target, 0, Budget(Policy())) == oracle("AAA", record("AAAANAAAA"), 0)
