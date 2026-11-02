"""Consensus-network assembly and its DOT rendering."""

from __future__ import annotations

from collections import defaultdict

import pytest

from ensemblemtd.consensus import build_consensus_network
from ensemblemtd.dotgraph import (
    build_dot,
    escape_dot,
    node_id_map,
    scan_overlay,
)
from ensemblemtd.settings import Settings

pytestmark = pytest.mark.unit


def species_stats(*names, initial=("A",)):
    return {
        name: {
            "display_run_support": 5,
            "observed_run_support": 5,
            "pct_runs_observed": 100.0,
            "is_initial": name in initial,
            "initial_frame_count": 1 if name in initial else 0,
        }
        for name in names
    }


def reaction_rows(*records):
    rows = []
    support = defaultdict(set)
    for lhs, rhs, events, runs in records:
        key = (tuple(lhs), tuple(rhs))
        rows.append((key, len(runs), 0.0, events))
        support[key] = set(runs)
    return rows, support


class TestConsensusNetwork:
    def test_it_writes_the_dot_graph_and_its_two_side_tables(self, tmp_path):
        rows, support = reaction_rows((["A"], ["B"], 6, {"r1", "r2", "r3"}))
        html, obs, meta = build_consensus_network(
            reaction_rows_full=rows,
            reaction_run_support=support,
            species_stats=species_stats("A", "B"),
            species_id_order=["A", "B"],
            n_runs=3,
            settings=Settings(inputs=("x",), outdir=str(tmp_path)),
            outdir=tmp_path,
        )
        assert (tmp_path / "aggregate_reacnet_network.dot").exists()
        assert (tmp_path / "aggregate_reacnet_network_nodes.tsv").exists()
        assert (tmp_path / "aggregate_reacnet_mechanistic_path.tsv").exists()
        assert "Aggregated Network" in html
        assert obs == {"A": 5, "B": 5}
        assert meta["candidate_edges"] == 1
        assert meta["shown_edges"] == 1

    def test_an_edge_below_the_run_support_floor_never_becomes_a_candidate(
        self, tmp_path
    ):
        rows, support = reaction_rows((["A"], ["B"], 6, {"r1"}))
        _html, _obs, meta = build_consensus_network(
            reaction_rows_full=rows,
            reaction_run_support=support,
            species_stats=species_stats("A", "B"),
            species_id_order=["A", "B"],
            n_runs=3,
            settings=Settings(inputs=("x",), outdir=str(tmp_path)),
            outdir=tmp_path,
        )
        assert meta["candidate_edges"] == 0
        assert meta["min_edge_runs"] == 3

    def test_no_route_is_drawn_by_default(self, tmp_path):
        rows, support = reaction_rows(
            (["A"], ["B"], 6, {"r1", "r2", "r3"}),
            (["B"], ["C"], 6, {"r1", "r2", "r3"}),
        )
        _html, _obs, meta = build_consensus_network(
            reaction_rows_full=rows,
            reaction_run_support=support,
            species_stats=species_stats("A", "B", "C"),
            species_id_order=["A", "B", "C"],
            n_runs=3,
            settings=Settings(inputs=("x",), outdir=str(tmp_path)),
            outdir=tmp_path,
        )
        assert meta["pathway_mode"] == "none"
        assert meta["mechanistic_path"] == "none"
        assert (tmp_path / "aggregate_reacnet_mechanistic_path.tsv").read_text() == (
            "step\tnode_id\tsmiles\n"
        )

    def test_a_requested_route_is_reported_and_highlighted(self, tmp_path):
        rows, support = reaction_rows(
            (["A"], ["B"], 6, {"r1", "r2", "r3"}),
            (["B"], ["C"], 6, {"r1", "r2", "r3"}),
        )
        _html, _obs, meta = build_consensus_network(
            reaction_rows_full=rows,
            reaction_run_support=support,
            species_stats=species_stats("A", "B", "C"),
            species_id_order=["A", "B", "C"],
            n_runs=3,
            settings=Settings(
                inputs=("x",), outdir=str(tmp_path), pathway_mode="support-first"
            ),
            outdir=tmp_path,
        )
        assert meta["mechanistic_path"] == "S1 -> S2 -> S3"
        assert meta["mechanistic_path_min_edge_runs"] == 3
        assert meta["mechanistic_path_sum_edge_runs"] == 6
        assert meta["mechanistic_path_sum_edge_events"] == pytest.approx(12.0)
        path_tsv = (tmp_path / "aggregate_reacnet_mechanistic_path.tsv").read_text()
        assert path_tsv.splitlines()[1:] == ["1\tS1\tA", "2\tS2\tB", "3\tS3\tC"]
        # highlighted nodes are drawn in the path colour
        assert "#ef4444" in (tmp_path / "aggregate_reacnet_network.dot").read_text()

    def test_a_path_length_scan_is_reported_per_length(self, tmp_path):
        rows, support = reaction_rows(
            (["A"], ["B"], 6, {"r1", "r2", "r3"}),
            (["B"], ["C"], 6, {"r1", "r2", "r3"}),
        )
        _html, _obs, meta = build_consensus_network(
            reaction_rows_full=rows,
            reaction_run_support=support,
            species_stats=species_stats("A", "B", "C"),
            species_id_order=["A", "B", "C"],
            n_runs=3,
            settings=Settings(
                inputs=("x",),
                outdir=str(tmp_path),
                pathway_mode="support-first",
                target_mode="all",
                scan_pathlength="1-2",
            ).validated(),
            outdir=tmp_path,
        )
        scan = {row["path_length"]: row for row in meta["pathlength_scan"]}
        assert scan[1]["path"] == ["S1", "S2"]
        assert scan[2]["path"] == ["S1", "S2", "S3"]
        assert scan[1]["mode"] == "support-first"
        assert scan[1]["color"] != scan[2]["color"]
        assert "cluster_scan_legend" in (
            tmp_path / "aggregate_reacnet_network.dot"
        ).read_text()

    def test_cleaning_nets_the_opposing_pair(self, tmp_path):
        rows, support = reaction_rows(
            (["A"], ["B"], 9, {"r1", "r2", "r3"}),
            (["B"], ["A"], 3, {"r1", "r2", "r3"}),
        )
        _html, _obs, meta = build_consensus_network(
            reaction_rows_full=rows,
            reaction_run_support=support,
            species_stats=species_stats("A", "B"),
            species_id_order=["A", "B"],
            n_runs=3,
            settings=Settings(inputs=("x",), outdir=str(tmp_path), clean=True),
            outdir=tmp_path,
        )
        assert meta["clean_enabled"] is True
        assert meta["clean_pairs_collapsed"] == 1
        assert meta["shown_edges"] == 1
        assert "e=6.0" in (tmp_path / "aggregate_reacnet_network.dot").read_text()

    def test_an_exact_tie_survives_cleaning_as_a_dashed_edge(self, tmp_path):
        rows, support = reaction_rows(
            (["A"], ["B"], 6, {"r1", "r2", "r3"}),
            (["B"], ["A"], 6, {"r1", "r2", "r3"}),
        )
        _html, _obs, meta = build_consensus_network(
            reaction_rows_full=rows,
            reaction_run_support=support,
            species_stats=species_stats("A", "B"),
            species_id_order=["A", "B"],
            n_runs=3,
            settings=Settings(inputs=("x",), outdir=str(tmp_path), clean=True),
            outdir=tmp_path,
        )
        assert meta["clean_tie_pairs_kept"] == 1
        dot = (tmp_path / "aggregate_reacnet_network.dot").read_text()
        assert 'style="dashed", dir="both"' in dot
        assert "e=0.0" in dot

    def test_an_empty_ensemble_produces_an_empty_graph_not_a_crash(self, tmp_path):
        html, obs, meta = build_consensus_network(
            reaction_rows_full=[],
            reaction_run_support=defaultdict(set),
            species_stats={},
            species_id_order=[],
            n_runs=0,
            settings=Settings(inputs=("x",), outdir=str(tmp_path)),
            outdir=tmp_path,
        )
        assert obs == {}
        assert meta["candidate_edges"] == 0
        assert meta["shown_nodes"] == 0
        assert "Aggregated Network" in html


class TestDotRendering:
    def test_node_ids_follow_the_species_order(self):
        node_list, ids = node_id_map({"B", "A", "C"}, {"A": 1, "B": 2, "C": 3})
        assert node_list == ["A", "B", "C"]
        assert ids == {"A": "N0", "B": "N1", "C": "N2"}

    def test_a_species_without_an_index_sorts_last(self):
        node_list, _ids = node_id_map({"A", "ghost"}, {"A": 1})
        assert node_list == ["A", "ghost"]

    @pytest.mark.parametrize(
        "text, expected",
        [
            ('say "hi"', 'say \\"hi\\"'),
            ("back\\slash", "back\\\\slash"),
            ("two\nlines", "two lines"),
        ],
    )
    def test_dot_strings_are_escaped(self, text, expected):
        assert escape_dot(text) == expected

    def test_each_scanned_length_gets_its_own_colour(self):
        node_colors, edge_colors, entries = scan_overlay(
            {2: {"path": ["A", "B"]}, 3: {"path": ["A", "C", "D"]}}
        )
        assert [length for length, _c in entries] == [2, 3]
        assert entries[0][1] != entries[1][1]
        # the shorter path claims shared nodes first
        assert node_colors["A"] == entries[0][1]
        assert edge_colors[("A", "C")] == entries[1][1]

    def test_pen_width_grows_with_event_weight(self):
        stats = species_stats("A", "B", "C", initial=())
        node_list, ids = node_id_map({"A", "B", "C"}, {"A": 1, "B": 2, "C": 3})
        dot = build_dot(
            node_list=node_list,
            node_ids=ids,
            display_edges=[("A", "B", 100.0, 5), ("A", "C", 1.0, 5)],
            species_stats=stats,
            species_idx={"A": 1, "B": 2, "C": 3},
            species_support={"A": 5, "B": 5, "C": 5},
            species_support_obs={"A": 5, "B": 5, "C": 5},
            species_pct_obs={"A": 100.0, "B": 100.0, "C": 100.0},
            path_nodes=[],
            path_edges=set(),
            scan_node_colors={},
            scan_edge_colors={},
            scan_entries=[],
            tie_edge_keys=set(),
        )
        widths = {
            line.split("penwidth=")[1].split(",")[0]
            for line in dot.splitlines()
            if "penwidth=" in line
        }
        assert widths == {"5.00", "1.40"}

    def test_an_edge_to_a_node_outside_the_subgraph_is_skipped(self):
        node_list, ids = node_id_map({"A"}, {"A": 1})
        dot = build_dot(
            node_list=node_list,
            node_ids=ids,
            display_edges=[("A", "offscreen", 1.0, 5)],
            species_stats=species_stats("A"),
            species_idx={"A": 1},
            species_support={"A": 5},
            species_support_obs={"A": 5},
            species_pct_obs={"A": 100.0},
            path_nodes=[],
            path_edges=set(),
            scan_node_colors={},
            scan_edge_colors={},
            scan_entries=[],
            tie_edge_keys=set(),
        )
        assert "->" not in dot
