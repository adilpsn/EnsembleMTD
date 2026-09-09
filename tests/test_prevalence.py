"""The prevalence metric and the reaction tallies."""

from __future__ import annotations

from collections import Counter, defaultdict

import pytest

from ensemblemtd.prevalence import (
    add_reactions,
    build_reaction_rows,
    key_to_str,
    reaction_event_weight,
    side_key,
    species_observation_stats,
)

pytestmark = pytest.mark.unit


def test_prevalence_is_the_fraction_of_runs_a_species_appeared_in():
    observed, pct = species_observation_stats(["run1", "run2", "run3"], n_runs=20)
    assert observed == 3
    assert pct == pytest.approx(15.0)


def test_a_run_counts_once_however_often_the_species_appeared():
    observed, pct = species_observation_stats(
        ["run1", "run1", "run1", "run2"], n_runs=4
    )
    assert observed == 2
    assert pct == pytest.approx(50.0)


def test_prevalence_of_zero_runs_is_zero_rather_than_a_division_error():
    assert species_observation_stats([], n_runs=0) == (0, 0.0)


@pytest.mark.parametrize(
    "record, expected",
    [({"n": 7}, 7), ({"n": 0}, 1), ({}, 1), ({"n": None}, 1), ({"n": -3}, 1)],
)
def test_a_reaction_record_always_weighs_at_least_one_event(record, expected):
    assert reaction_event_weight(record) == expected


def _fold(records, run_name, events=None, support=None, species=None):
    events = Counter() if events is None else events
    support = defaultdict(set) if support is None else support
    species = Counter() if species is None else species
    add_reactions(records, run_name, events, support, species)
    return events, support, species


def test_reaction_sides_are_sorted_so_member_order_does_not_split_a_key():
    one, _s, _c = _fold([{"l": ["B", "A"], "r": ["C"], "n": 1}], "run1")
    two, _s, _c = _fold([{"l": ["A", "B"], "r": ["C"], "n": 1}], "run2")
    assert list(one) == list(two)


def test_a_record_with_identical_sides_is_not_a_reaction():
    events, support, species = _fold([{"l": ["A"], "r": ["A"], "n": 5}], "run1")
    assert events == Counter()
    assert support == {}
    assert species == Counter()


def test_a_record_with_an_empty_side_is_dropped():
    events, _s, _c = _fold([{"l": [], "r": ["A"], "n": 2}], "run1")
    assert events == Counter()


def test_non_string_species_entries_are_ignored():
    events, _s, _c = _fold([{"l": ["A", 7, None], "r": ["B"], "n": 1}], "run1")
    assert list(events) == [(("A",), ("B",))]


def test_events_accumulate_while_run_support_stays_binary():
    events, support, species = _fold([{"l": ["A"], "r": ["B"], "n": 3}], "run1")
    _fold([{"l": ["A"], "r": ["B"], "n": 4}], "run1", events, support, species)
    _fold([{"l": ["A"], "r": ["B"], "n": 1}], "run2", events, support, species)

    key = (("A",), ("B",))
    assert events[key] == 8
    assert support[key] == {"run1", "run2"}
    # every species on either side is credited with the record's events
    assert species["A"] == 8
    assert species["B"] == 8


def test_reaction_rows_are_ordered_by_run_support_then_events():
    events = Counter(
        {
            (("A",), ("B",)): 2,
            (("A",), ("C",)): 99,
            (("A",), ("D",)): 5,
        }
    )
    support = defaultdict(
        set,
        {
            (("A",), ("B",)): {"r1", "r2", "r3"},
            (("A",), ("C",)): {"r1"},
            (("A",), ("D",)): {"r1", "r2", "r3"},
        },
    )

    rows = build_reaction_rows(events, support, n_runs=4)

    # D before B: same support (3), more events (5 > 2).  C last despite 99
    # events, because it was seen in only one run.
    assert [key[1][0] for key, _rs, _pct, _n in rows] == ["D", "B", "C"]
    assert rows[0][1:] == (3, 75.0, 5)


def test_reaction_percentages_are_zero_when_no_run_succeeded():
    rows = build_reaction_rows(
        Counter({(("A",), ("B",)): 1}),
        defaultdict(set, {(("A",), ("B",)): {"r1"}}),
        n_runs=0,
    )
    assert rows[0][2] == 0.0


@pytest.mark.parametrize(
    "side, expected",
    [
        (["A"], "A"),
        (["A", "B"], "A + B"),
        (["A", "A"], "2*A"),
        (["B", "A", "A"], "2*A + B"),
        ([], "∅"),
    ],
)
def test_a_reaction_side_shows_stoichiometry(side, expected):
    assert side_key(side) == expected


def test_a_reaction_key_renders_as_lhs_to_rhs():
    assert key_to_str((("A", "A"), ("B", "C"))) == "2*A => B + C"
