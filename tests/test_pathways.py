"""Route selection over the consensus graph."""

from __future__ import annotations

import pytest

from ensemblemtd.network import build_adjacency, compute_node_flow_stats
from ensemblemtd.pathways import (
    PathContext,
    best_paths_for_threshold,
    choose_path,
    choose_probabilistic,
    choose_support_first,
    dijkstra_cost,
    pseudo_sink_scores,
    reconstruct_path,
    score_better,
    select_target_candidates,
)

pytestmark = pytest.mark.unit


def context(candidates, initial=("A",), n_runs=10, target_mode="strict-sink", **kw):
    species = []
    for u, v, _e, _r in candidates:
        for s in (u, v):
            if s not in species:
                species.append(s)
    species_idx = {s: i + 1 for i, s in enumerate(sorted(species))}
    _adj, outdeg, indeg, nodes = build_adjacency(candidates, species_idx)
    in_events, out_events, _ir, _or = compute_node_flow_stats(candidates)
    return PathContext(
        candidates=tuple(candidates),
        nodes=frozenset(nodes),
        species_idx=species_idx,
        species_support_obs=kw.get("support", {s: 1 for s in species}),
        outdeg=outdeg,
        indeg=indeg,
        in_events=in_events,
        out_events=out_events,
        initial_start_candidates=tuple(s for s in initial if s in nodes),
        n_runs=n_runs,
        target_mode=target_mode,
        path_min_steps=kw.get("path_min_steps", 1),
        top_targets=kw.get("top_targets", 1),
    )


class TestGraphSearch:
    def test_dijkstra_finds_the_cheapest_route_not_the_shortest(self):
        adjacency = {
            "A": [("B", 1.0), ("C", 10.0)],
            "B": [("C", 1.0)],
        }
        dist, prev = dijkstra_cost(adjacency, "A")
        assert dist["C"] == pytest.approx(2.0)
        assert reconstruct_path(prev, "A", "C") == ["A", "B", "C"]

    def test_an_unreachable_node_has_no_route(self):
        dist, prev = dijkstra_cost({"A": [("B", 1.0)]}, "A")
        assert "Z" not in dist
        assert reconstruct_path(prev, "A", "Z") is None

    def test_the_route_to_the_start_is_the_start(self):
        assert reconstruct_path({}, "A", "A") == ["A"]

    def test_a_broken_predecessor_chain_yields_no_route(self):
        assert reconstruct_path({"C": "B"}, "A", "C") is None

    def test_layered_search_prefers_the_best_supported_route_into_a_layer(self):
        # both routes reach C in two steps; the one through B2 has the stronger
        # weakest edge (4 vs 1) and must win
        adjacency = {
            "A": [("B1", 1, 5.0), ("B2", 4, 5.0)],
            "B1": [("C", 9, 5.0)],
            "B2": [("C", 4, 5.0)],
        }
        dist, best, prev = best_paths_for_threshold(adjacency, "A")
        assert dist["C"] == 2
        assert reconstruct_path(prev, "A", "C") == ["A", "B2", "C"]
        assert best["C"][0] == 4  # the weakest edge along the chosen route

    def test_a_higher_score_tuple_wins(self):
        assert score_better((3, 0, 0.0, 0), None)
        assert score_better((3, 0, 0.0, 0), (2, 99, 99.0, 0))
        assert not score_better((2, 0, 0.0, 0), (3, 0, 0.0, 0))


class TestTargetSelection:
    NODES = ["A", "B", "C"]

    def test_a_strict_sink_is_a_node_nothing_leaves(self):
        targets = select_target_candidates(
            self.NODES, "A", {"A": 1, "B": 1}, {"B": 1, "C": 1},
            "strict-sink", {}, {},
        )
        assert targets == ["C"]

    def test_with_no_sink_at_all_every_other_node_is_allowed(self):
        # reverse edges commonly leave the graph without a true sink
        targets = select_target_candidates(
            self.NODES, "A", {"A": 1, "B": 1, "C": 1}, {"A": 1, "B": 1, "C": 1},
            "strict-sink", {}, {},
        )
        assert targets == ["B", "C"]

    def test_a_pseudo_sink_only_needs_incoming_weight(self):
        targets = select_target_candidates(
            self.NODES, "A", {}, {}, "pseudo-sink", {"B": 2.0, "C": 0.0}, {}
        )
        assert targets == ["B"]

    def test_target_mode_all_accepts_anything_but_the_start(self):
        targets = select_target_candidates(
            self.NODES, "A", {}, {}, "all", {}, {}
        )
        assert targets == ["B", "C"]

    def test_a_product_like_node_scores_high(self):
        scores = pseudo_sink_scores(["A", "B"], {"B": 4.0}, {"A": 4.0})
        assert scores["B"] > scores["A"]


class TestSupportFirst:
    def test_the_route_with_the_strongest_weakest_link_is_chosen(self):
        candidates = [
            ("A", "B1", 5.0, 1),
            ("B1", "C", 5.0, 9),
            ("A", "B2", 5.0, 4),
            ("B2", "C", 5.0, 4),
        ]
        result = choose_support_first(context(candidates, target_mode="all"))
        assert result["path"] == ["A", "B2", "C"]

    def test_the_route_starts_from_an_initial_fragment(self):
        candidates = [("A", "B", 5.0, 3), ("X", "B", 5.0, 3)]
        result = choose_support_first(context(candidates, initial=("X",), target_mode="all"))
        assert result["start"] == "X"

    def test_with_no_initial_fragment_left_the_dominant_node_is_used(self, capsys):
        candidates = [("A", "B", 5.0, 3)]
        result = choose_support_first(
            context(candidates, initial=("gone",), support={"A": 9, "B": 1}, target_mode="all")
        )
        assert result["start"] == "A"
        assert "No initial fragments survived" in capsys.readouterr().err

    def test_a_minimum_step_count_rules_out_shorter_routes(self):
        candidates = [("A", "B", 5.0, 3), ("B", "C", 5.0, 3)]
        result = choose_support_first(
            context(candidates, target_mode="all", path_min_steps=2)
        )
        assert result["path"] == ["A", "B", "C"]

    def test_an_exact_step_count_can_be_requested(self):
        candidates = [("A", "B", 5.0, 3), ("B", "C", 5.0, 3)]
        assert choose_support_first(
            context(candidates, target_mode="all"), exact_steps=1
        )["path"] == ["A", "B"]
        assert choose_support_first(
            context(candidates, target_mode="all"), exact_steps=2
        )["path"] == ["A", "B", "C"]

    def test_an_impossible_step_count_yields_no_route(self):
        candidates = [("A", "B", 5.0, 3)]
        assert choose_support_first(
            context(candidates, target_mode="all"), exact_steps=7
        ) is None

    def test_an_empty_graph_yields_no_route(self):
        assert choose_support_first(context([])) is None


class TestProbabilistic:
    def test_the_likelier_branch_is_chosen(self):
        # A branches 9:1 towards B; both branches then reach the sink D
        candidates = [
            ("A", "B", 9.0, 5),
            ("A", "C", 1.0, 5),
            ("B", "D", 5.0, 5),
            ("C", "D", 5.0, 5),
        ]
        result = choose_probabilistic(context(candidates, n_runs=5))
        assert result["target"] == "D"
        assert result["path"] == ["A", "B", "D"]

    def test_run_support_is_folded_into_the_cost(self):
        # equal branching weight, but the route through C is seen in more runs
        candidates = [
            ("A", "B", 5.0, 1),
            ("A", "C", 5.0, 5),
            ("B", "D", 5.0, 1),
            ("C", "D", 5.0, 5),
        ]
        result = choose_probabilistic(context(candidates, n_runs=5))
        assert result["path"] == ["A", "C", "D"]

    def test_the_start_is_the_most_prevalent_node(self):
        candidates = [("A", "B", 5.0, 3), ("B", "C", 5.0, 3)]
        result = choose_probabilistic(
            context(candidates, target_mode="all", support={"A": 9, "B": 2, "C": 1})
        )
        assert result["start"] == "A"

    def test_the_cost_is_reported_alongside_the_route(self):
        # one edge carrying all the weight in every run: P_eff = 1, cost = 0
        result = choose_probabilistic(context([("A", "B", 5.0, 5)], n_runs=5))
        assert result["cost"] == pytest.approx(0.0, abs=1e-9)

    def test_an_empty_graph_yields_no_route(self):
        assert choose_probabilistic(context([])) is None


class TestDispatch:
    CANDIDATES = [("A", "B", 5.0, 3), ("B", "C", 5.0, 3)]

    @pytest.mark.parametrize("mode", ["support-first", "probabilistic"])
    def test_both_modes_return_a_route(self, mode):
        result = choose_path(context(self.CANDIDATES, target_mode="all"), mode)
        assert result is not None and len(result["path"]) >= 2

    def test_mode_none_draws_nothing(self):
        assert choose_path(context(self.CANDIDATES), "none") is None
