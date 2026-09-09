"""The species table built on top of the collapsed ensemble."""

from __future__ import annotations

from collections import Counter, defaultdict

import pytest

from ensemblemtd.aggregate import (
    CollapsedEnsemble,
    RunHarvest,
    build_species_table,
    collapse_harvest,
    support_audit_warnings,
)

pytestmark = pytest.mark.unit


def ensemble(**kw) -> CollapsedEnsemble:
    defaults = dict(
        species_run_support=defaultdict(set),
        species_run_support_json=defaultdict(set),
        species_run_support_timeline=defaultdict(set),
        species_event_counter=Counter(),
        initial_species_counts=Counter(),
        initial_first_seen_order={},
        first_frame_support={},
        reaction_rows_full=[],
        reaction_rows_abcd_full=[],
        reaction_run_support=defaultdict(set),
        rep_to_inchi={},
        rep_to_members={},
        n_raw_species=0,
        n_collapsed_species=0,
        obabel_version="",
    )
    defaults.update(kw)
    return CollapsedEnsemble(**defaults)


def test_prevalence_reaches_the_species_table():
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(
                set, {"A": {"r1", "r2", "r3"}, "B": {"r1"}}
            )
        ),
        n_runs=4,
        support_source="timeline",
    )
    assert table.stats["A"]["pct_runs_observed"] == pytest.approx(75.0)
    assert table.stats["B"]["observed_run_support"] == 1


def test_a_reactant_consumed_before_the_first_dump_keeps_its_support():
    # Under a strong bias a reactant can be gone by the time the first frame is
    # written, so displayed support falls back to the first-frame count.
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(set, {"EC": {"r1"}}),
            initial_species_counts=Counter({"EC": 5}),
            first_frame_support={"EC": 5},
        ),
        n_runs=5,
        support_source="timeline",
    )
    assert table.stats["EC"]["observed_run_support"] == 1
    assert table.stats["EC"]["first_frame_run_support"] == 5
    assert table.stats["EC"]["display_run_support"] == 5
    assert table.stats["EC"]["is_initial"] is True


def test_display_support_never_falls_below_what_was_observed():
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(set, {"EC": {"r1", "r2", "r3"}}),
            initial_species_counts=Counter({"EC": 1}),
            first_frame_support={"EC": 1},
        ),
        n_runs=3,
        support_source="timeline",
    )
    assert table.stats["EC"]["display_run_support"] == 3


def test_reactants_come_first_in_the_order_they_were_first_seen():
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(
                set, {"EC1": {"r1"}, "EC2": {"r1"}, "P": {"r1", "r2"}}
            ),
            initial_species_counts=Counter({"EC1": 1, "EC2": 1}),
            initial_first_seen_order={"EC1": (1, 1), "EC2": (1, 0)},
        ),
        n_runs=2,
        support_source="timeline",
    )
    # EC2 appeared at position 0 of run 1, so it is S1 even though EC1 sorts
    # first alphabetically; the product follows the reactants
    assert table.id_order == ["EC2", "EC1", "P"]
    assert table.index == {"EC2": 1, "EC1": 2, "P": 3}


def test_products_are_ordered_by_prevalence_then_events_then_name():
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(
                set, {"low": {"r1"}, "high": {"r1", "r2"}, "mid": {"r1"}}
            ),
            species_event_counter=Counter({"mid": 10, "low": 1}),
        ),
        n_runs=2,
        support_source="timeline",
    )
    assert table.id_order == ["high", "mid", "low"]


def test_the_rank_order_leads_with_displayed_support():
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(set, {"P": {"r1", "r2"}, "EC": {"r1"}}),
            initial_species_counts=Counter({"EC": 3}),
            first_frame_support={"EC": 3},
        ),
        n_runs=3,
        support_source="timeline",
    )
    # EC shows 3 (first-frame) against the product's 2
    assert table.rank_order[0] == "EC"


def test_a_species_seen_only_in_the_timeline_is_still_in_the_table():
    # ReacNet's JSON species list is capped by --maxspecies; the timeline is not
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(set, {"A": {"r1", "r2"}}),
            species_run_support_timeline=defaultdict(set, {"A": {"r1", "r2"}}),
            species_run_support_json=defaultdict(set, {"A": {"r1"}}),
        ),
        n_runs=2,
        support_source="timeline",
    )
    assert table.stats["A"]["observed_run_support"] == 2
    assert table.stats["A"]["observed_run_support_json"] == 1
    assert table.stats["A"]["support_delta_json"] == 1
    assert table.discrepancies == [("A", 2, 1, 1)]


def test_the_widest_support_gap_is_reported_first():
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(
                set, {"small": {"r1", "r2"}, "big": {"r1", "r2", "r3", "r4"}}
            ),
            species_run_support_json=defaultdict(set, {"small": {"r1"}, "big": set()}),
        ),
        n_runs=4,
        support_source="timeline",
    )
    assert [row[0] for row in table.discrepancies] == ["big", "small"]


def test_a_widespread_support_gap_raises_a_warning():
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(set, {"A": {"r1"}, "B": {"r1"}}),
            species_run_support_json=defaultdict(set),
        ),
        n_runs=1,
        support_source="timeline",
    )
    warnings = dict(support_audit_warnings(table, universe_size=2))
    assert "support_audit:aggregate" in warnings
    assert "2/2 species (100.0%)" in warnings["support_audit:aggregate"]
    assert "support_audit:top_species" in warnings


def test_no_gap_means_no_warning():
    table = build_species_table(
        ensemble(
            species_run_support=defaultdict(set, {"A": {"r1"}}),
            species_run_support_json=defaultdict(set, {"A": {"r1"}}),
        ),
        n_runs=1,
        support_source="timeline",
    )
    assert support_audit_warnings(table, universe_size=1) == []


def test_an_empty_harvest_collapses_to_an_empty_ensemble():
    collapsed = collapse_harvest(
        RunHarvest(), n_runs=0, support_source="timeline", obabel_bin="obabel"
    )
    assert collapsed.n_raw_species == 0
    assert collapsed.n_collapsed_species == 0
    assert collapsed.collapse_ratio == 1.0
    assert collapsed.universe() == set()
    assert collapsed.obabel_version == ""  # Open Babel was never called


def test_the_harvest_universe_spans_every_tally():
    harvest = RunHarvest()
    harvest.species_run_support_timeline["A"].add("r1")
    harvest.species_run_support_json["B"].add("r1")
    harvest.species_event_counter["C"] += 1
    harvest.initial_species_counts["D"] += 1
    assert harvest.species_universe() == {"A", "B", "C", "D"}


def test_collapsing_rekeys_every_tally_onto_representatives(monkeypatch):
    monkeypatch.setattr(
        "ensemblemtd.aggregate.get_obabel_version", lambda binary: "Open Babel 3.1.1"
    )
    monkeypatch.setattr(
        "ensemblemtd.collapse.smiles_to_inchi_fixedh",
        lambda smiles, obabel_bin: "InChI=1/X" if smiles in ("a", "b") else "InChI=1/Y",
    )

    harvest = RunHarvest()
    harvest.species_run_support_timeline["a"].add("r1")
    harvest.species_run_support_timeline["b"].add("r2")
    harvest.species_run_support_timeline["z"].add("r1")
    harvest.species_event_counter["a"] += 2
    harvest.species_event_counter["b"] += 3
    harvest.initial_species_counts["a"] += 1
    harvest.initial_run_support["a"].add("r1")
    harvest.initial_first_seen_order["a"] = (1, 0)
    harvest.reaction_events[(("a",), ("z",))] = 4
    harvest.reaction_run_support[(("a",), ("z",))] = {"r1"}

    collapsed = collapse_harvest(
        harvest, n_runs=2, support_source="timeline", obabel_bin="obabel"
    )

    assert collapsed.n_raw_species == 3
    assert collapsed.n_collapsed_species == 2
    assert collapsed.obabel_version == "Open Babel 3.1.1"
    assert collapsed.species_run_support["a"] == {"r1", "r2"}
    assert collapsed.species_event_counter["a"] == 5
    assert collapsed.first_frame_support == {"a": 1}
    assert collapsed.reaction_rows_full[0][0] == (("a",), ("z",))
    assert collapsed.rep_to_members["a"] == ["a", "b"]
