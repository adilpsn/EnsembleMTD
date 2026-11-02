"""Edge projection, thresholds and the cleaned display."""

from __future__ import annotations

from collections import defaultdict

import pytest

from ensemblemtd.network import (
    build_adjacency,
    build_linkreac,
    compute_node_flow_stats,
    filter_edge_candidates,
    net_bidirectional_edges,
    project_reactions_to_edges,
    select_display_subgraph,
    sort_candidates,
)

pytestmark = pytest.mark.unit


def rows(*records):
    """Build (reaction rows, run support) from (lhs, rhs, events, runs) tuples."""
    reaction_rows = []
    support = defaultdict(set)
    for lhs, rhs, events, runs in records:
        key = (tuple(lhs), tuple(rhs))
        reaction_rows.append((key, len(runs), 0.0, events))
        support[key] = set(runs)
    return reaction_rows, support


def test_a_one_to_one_record_puts_all_its_events_on_one_edge():
    events, support = project_reactions_to_edges(*rows((["A"], ["B"], 6, {"r1"})))
    assert events == {("A", "B"): 6.0}
    assert support == {("A", "B"): {"r1"}}


def test_a_records_events_are_shared_over_the_pairs_it_implies():
    # one lhs, two rhs -> denominator 1*2, so 6 events give 3 per edge
    events, _s = project_reactions_to_edges(*rows((["A"], ["B", "C"], 6, {"r1"})))
    assert events == {("A", "B"): 3.0, ("A", "C"): 3.0}
    assert sum(events.values()) == 6.0


def test_the_denominator_counts_self_pairs_before_they_are_dropped():
    # A + B -> A + C has denominator 2*2 = 4, and A->A is then omitted, so the
    # surviving three edges carry 3/4 of the record rather than all of it
    events, _s = project_reactions_to_edges(*rows((["A", "B"], ["A", "C"], 8, {"r1"})))
    assert set(events) == {("A", "C"), ("B", "A"), ("B", "C")}
    assert all(v == pytest.approx(2.0) for v in events.values())
    assert sum(events.values()) == pytest.approx(6.0)


def test_repeated_species_on_a_side_count_once_in_the_denominator():
    events, _s = project_reactions_to_edges(*rows((["A", "A"], ["B"], 4, {"r1"})))
    assert events == {("A", "B"): 4.0}


def test_edge_run_support_is_the_union_over_contributing_records():
    events, support = project_reactions_to_edges(
        *rows(
            (["A"], ["B"], 2, {"r1", "r2"}),
            (["A", "X"], ["B"], 2, {"r2", "r3"}),
        )
    )
    assert support[("A", "B")] == {"r1", "r2", "r3"}
    assert events[("A", "B")] == pytest.approx(3.0)  # 2 + 2/2


def test_edges_below_either_threshold_are_dropped():
    events = {("A", "B"): 5.0, ("A", "C"): 0.4, ("A", "D"): 5.0}
    support = defaultdict(
        set, {("A", "B"): {"r1", "r2", "r3"}, ("A", "C"): {"r1", "r2", "r3"}, ("A", "D"): {"r1"}}
    )

    kept = filter_edge_candidates(events, support, min_edge_runs=3, min_edge_events=1.0)

    assert [(u, v) for u, v, _e, _r in kept] == [("A", "B")]


def test_candidates_are_ordered_by_run_support_then_event_weight():
    candidates = [
        ("A", "B", 1.0, 2),
        ("A", "C", 9.0, 5),
        ("A", "D", 1.0, 5),
    ]
    ordered = sort_candidates(candidates, {}, {"A": 1, "B": 2, "C": 3, "D": 4})
    assert [v for _u, v, _e, _r in ordered] == ["C", "D", "B"]


def test_node_flow_separates_incoming_from_outgoing_weight():
    in_events, out_events, in_runs, out_runs = compute_node_flow_stats(
        [("A", "B", 3.0, 2), ("C", "B", 1.0, 1)]
    )
    assert out_events == {"A": 3.0, "C": 1.0}
    assert in_events == {"B": 4.0}
    assert in_runs == {"B": 3}
    assert out_runs == {"A": 2, "C": 1}


def test_adjacency_reports_degrees_and_the_nodes_the_edges_span():
    adjacency, outdeg, indeg, nodes = build_adjacency(
        [("A", "B", 1.0, 1), ("A", "C", 1.0, 1)], {"A": 1, "B": 2, "C": 3}
    )
    assert nodes == {"A", "B", "C"}
    assert outdeg["A"] == 2
    assert indeg["B"] == 1
    assert [v for v, _r, _e in adjacency["A"]] == ["B", "C"]


class TestCleanedDisplay:
    def test_an_opposing_pair_becomes_one_net_edge(self):
        edges, ties, stats = net_bidirectional_edges(
            [("A", "B", 7.0, 4), ("B", "A", 2.0, 3)], mandatory_edge_keys=set()
        )
        assert edges == [("A", "B", 5.0, 4)]
        assert ties == set()
        assert stats == {"pairs_collapsed": 1, "tie_pairs_kept": 0}

    def test_the_net_edge_points_the_way_the_larger_weight_points(self):
        edges, _t, _s = net_bidirectional_edges(
            [("A", "B", 2.0, 3), ("B", "A", 7.0, 4)], mandatory_edge_keys=set()
        )
        assert edges == [("B", "A", 5.0, 4)]

    def test_an_exact_tie_is_kept_as_a_zero_weight_bidirectional_edge(self):
        edges, ties, stats = net_bidirectional_edges(
            [("A", "B", 4.0, 3), ("B", "A", 4.0, 2)], mandatory_edge_keys=set()
        )
        assert edges == [("A", "B", 0.0, 3)]
        assert ties == {("A", "B")}
        assert stats["tie_pairs_kept"] == 1

    def test_an_edge_on_a_highlighted_path_is_never_netted_away(self):
        edges, _t, stats = net_bidirectional_edges(
            [("A", "B", 7.0, 4), ("B", "A", 2.0, 3)],
            mandatory_edge_keys={("A", "B")},
        )
        assert ("A", "B", 7.0, 4) in edges
        assert ("B", "A", 2.0, 3) in edges
        assert stats["pairs_collapsed"] == 0

    def test_a_one_way_edge_passes_through_untouched(self):
        edges, _t, stats = net_bidirectional_edges(
            [("A", "B", 7.0, 4)], mandatory_edge_keys=set()
        )
        assert edges == [("A", "B", 7.0, 4)]
        assert stats["pairs_collapsed"] == 0

    def test_a_self_loop_passes_through_untouched(self):
        edges, _t, _s = net_bidirectional_edges(
            [("A", "A", 3.0, 2)], mandatory_edge_keys=set()
        )
        assert edges == [("A", "A", 3.0, 2)]


class TestDisplaySubgraph:
    def test_mandatory_nodes_survive_a_limit_that_would_exclude_them(self):
        edges, nodes = select_display_subgraph(
            candidates=[],
            mandatory_nodes={"A", "B", "C"},
            mandatory_edges=[],
            species_id_order=["A", "B", "C"],
            node_limit=1,
            edge_limit=1,
            include_isolates=False,
        )
        assert nodes == {"A", "B", "C"}
        assert edges == []

    def test_strong_edges_are_taken_first_until_the_edge_limit(self):
        candidates = [("A", "B", 9.0, 5), ("C", "D", 8.0, 4), ("E", "F", 7.0, 3)]
        edges, _nodes = select_display_subgraph(
            candidates=candidates,
            mandatory_nodes=set(),
            mandatory_edges=[],
            species_id_order=[],
            node_limit=99,
            edge_limit=2,
            include_isolates=False,
        )
        assert [(u, v) for u, v, _e, _r in edges] == [("A", "B"), ("C", "D")]

    def test_an_edge_that_would_exceed_the_node_limit_is_skipped_not_stopped_on(self):
        candidates = [("A", "B", 9.0, 5), ("C", "D", 8.0, 4), ("A", "B2", 7.0, 3)]
        edges, nodes = select_display_subgraph(
            candidates=candidates,
            mandatory_nodes={"A", "B"},
            mandatory_edges=[],
            species_id_order=[],
            node_limit=3,
            edge_limit=99,
            include_isolates=False,
        )
        # C->D needs two new nodes and does not fit; A->B2 needs one and does
        assert [(u, v) for u, v, _e, _r in edges] == [("A", "B"), ("A", "B2")]
        assert nodes == {"A", "B", "B2"}

    def test_isolates_are_added_only_when_asked_for(self):
        _edges, nodes = select_display_subgraph(
            candidates=[],
            mandatory_nodes={"A"},
            mandatory_edges=[],
            species_id_order=["A", "B", "C"],
            node_limit=3,
            edge_limit=3,
            include_isolates=True,
        )
        assert nodes == {"A", "B", "C"}


def test_linkreac_is_undirected_and_lists_every_species():
    reaction_rows, _support = rows((["A"], ["B", "C"], 1, {"r1"}))
    links = build_linkreac(reaction_rows, [], species_order=["A", "B", "C", "Z"])

    assert links["A"] == ["B", "C"]
    assert links["B"] == ["A"]
    assert links["Z"] == []  # present but unconnected


def test_linkreac_ignores_a_species_reacting_with_itself():
    reaction_rows, _support = rows((["A", "B"], ["A"], 1, {"r1"}))
    links = build_linkreac(reaction_rows, [], species_order=["A", "B"])
    assert "A" not in links["A"]
