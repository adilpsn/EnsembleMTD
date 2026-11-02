"""Assembling the consensus network from the collapsed ensemble tallies.

This ties the pieces together: project reactions onto directed species pairs,
threshold them, optionally pick a route to highlight, choose the subgraph to
show, and hand the result to the DOT renderer.
"""

from __future__ import annotations

import html
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

from . import dotgraph
from .network import (
    Edge,
    EdgeKey,
    ReactionKey,
    build_adjacency,
    compute_node_flow_stats,
    filter_edge_candidates,
    net_bidirectional_edges,
    project_reactions_to_edges,
    select_display_subgraph,
    sort_candidates,
)
from .pathways import PathContext, choose_path
from .settings import (
    DEFAULT_GRAPH_INCLUDE_ISOLATES,
    DEFAULT_GRAPH_MAX_EDGES,
    DEFAULT_GRAPH_MAX_NODES,
    Settings,
)

EDGE_LOGIC_NOTE = (
    "Edge logic: reactions are decomposed into directed species pairs; edge "
    "event weight is distributed by lhs×rhs pair count; pair-edge support "
    "uses union-of-runs across contributing reactions."
)

PATH_LOGIC = {
    "probabilistic": (
        "min cumulative -log(P_eff)",
        "Path logic: minimize cumulative -log(P_eff), where "
        "P_eff(u->v)=P(u->v)*(r_uv/N) and P(u->v)=w_uv/sum_k(w_uk).",
    ),
    "support-first": (
        "support-first",
        "Path logic: maximize minimum edge run support first, then choose the "
        "strongest shortest path within that support threshold.",
    ),
    "none": (
        "disabled",
        "Path logic: pathway-mode=none (no pathfinding applied).",
    ),
}


def _species_views(
    species_stats: Dict[str, Dict[str, object]], species_id_order: Sequence[str]
) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, float], Dict[str, int]]:
    return (
        {s: int(v["display_run_support"]) for s, v in species_stats.items()},
        {s: int(v["observed_run_support"]) for s, v in species_stats.items()},
        {s: float(v["pct_runs_observed"]) for s, v in species_stats.items()},
        {s: i + 1 for i, s in enumerate(species_id_order)},
    )


def _highlighted_path(
    ctx: PathContext,
    settings: Settings,
    candidates: Sequence[Edge],
    species_idx: Dict[str, int],
) -> Tuple[List[str], Set[EdgeKey], Dict[str, object], Dict[int, Dict[str, object]]]:
    """Run the path selector, plus the path-length scan if one was requested."""
    main_result = choose_path(ctx, settings.pathway_mode)

    scan_results: Dict[int, Dict[str, object]] = {}
    if settings.scan_pathlength_values and settings.pathway_mode != "none":
        for length in sorted(settings.scan_pathlength_values):
            result = choose_path(ctx, settings.pathway_mode, exact_steps=length)
            if not result:
                continue
            scan_results[length] = {
                "target": str(result["target"]),
                "path": [str(x) for x in result["path"]],
                "score": result["score"],
                "mode": settings.pathway_mode,
            }

    summary: Dict[str, object] = {
        "description": "none",
        "min_edge_runs": None,
        "sum_edge_runs": None,
        "sum_edge_events": None,
        "start": None,
        "target": None,
    }
    path_nodes: List[str] = []
    path_edges: Set[EdgeKey] = set()

    if main_result is not None:
        path_nodes = [str(x) for x in main_result["path"]]
        path_edges = {
            (path_nodes[i], path_nodes[i + 1]) for i in range(len(path_nodes) - 1)
        }
        edge_map = {(u, v): (float(e), int(r)) for u, v, e, r in candidates}
        values = [edge_map[key] for key in path_edges if key in edge_map]
        summary["description"] = " -> ".join(f"S{species_idx[s]}" for s in path_nodes)
        summary["start"] = str(main_result["start"])
        summary["target"] = str(main_result["target"])
        if values:
            summary["min_edge_runs"] = min(r for _e, r in values)
            summary["sum_edge_runs"] = sum(r for _e, r in values)
            summary["sum_edge_events"] = round(sum(e for e, _r in values), 6)

    return path_nodes, path_edges, summary, scan_results


def _report_block(
    node_count: int,
    species_count: int,
    shown_edges: int,
    candidate_edges: int,
    n_runs: int,
    path_summary: Dict[str, object],
    settings: Settings,
    scan_results: Dict[int, Dict[str, object]],
) -> str:
    title, logic = PATH_LOGIC[settings.pathway_mode]
    if settings.pathway_mode != "none":
        logic += (
            f" Target mode={settings.target_mode}, "
            f"min steps={settings.path_min_steps}."
        )
    scan_note = ""
    if scan_results:
        lengths = ", ".join(f"L={x}" for x in sorted(scan_results))
        scan_note = (
            f"<p><b>Path-length scan:</b> {html.escape(lengths)} (exact edges).</p>"
        )
    return (
        "<div class='container py-3'>"
        "<h3>Aggregated Network</h3>"
        f"<p>Nodes shown: {node_count} / {species_count}; edges shown: "
        f"{shown_edges} / {candidate_edges}; runs: {n_runs}.</p>"
        f"<p>{EDGE_LOGIC_NOTE}</p>"
        f"<p><b>Highlighted mechanistic path ({title}):</b> "
        f"{html.escape(str(path_summary['description']))}</p>"
        f"<p>{html.escape(logic)}</p>"
        f"{scan_note}"
        "<p>Node labels (S#) map directly to entries in the Species table.</p>"
        "</div>"
    )


def _scan_metadata(
    scan_results: Dict[int, Dict[str, object]],
    scan_entries: Sequence[Tuple[int, str]],
    species_idx: Dict[str, int],
    pathway_mode: str,
) -> List[dict]:
    colors = dict(scan_entries)
    rows = []
    for length in sorted(scan_results):
        row = scan_results[length]
        score = row["score"]
        if pathway_mode == "probabilistic":
            score_text = f"cost={float(score[0]):.6f}"
        else:
            score_text = "rank=" + ",".join(f"{float(x):.6g}" for x in score)
        path = row["path"]
        rows.append(
            {
                "path_length": int(length),
                "target": str(row["target"]),
                "target_id": f"S{species_idx.get(str(row['target']), 0)}",
                "path": [f"S{species_idx.get(s, 0)}" for s in path],
                "path_smiles": path,
                "score": score_text,
                "mode": pathway_mode,
                "color": colors.get(length, dotgraph.SCAN_COLORS[0]),
            }
        )
    return rows


def build_consensus_network(
    reaction_rows_full: List[Tuple[ReactionKey, int, float, int]],
    reaction_run_support: Dict[ReactionKey, Set[str]],
    species_stats: Dict[str, Dict[str, object]],
    species_id_order: Sequence[str],
    n_runs: int,
    settings: Settings,
    outdir: Path,
) -> Tuple[str, Dict[str, int], Dict[str, object]]:
    """Build, write and describe the consensus network.

    Returns the HTML fragment for the report, the per-species observed run
    support, and a metadata dict that goes into the run summary.
    """
    species_support, species_support_obs, species_pct_obs, species_idx = _species_views(
        species_stats, species_id_order
    )

    edge_events, edge_run_names = project_reactions_to_edges(
        reaction_rows_full, reaction_run_support
    )
    candidates = sort_candidates(
        filter_edge_candidates(
            edge_events,
            edge_run_names,
            settings.graph_min_edge_runs,
            settings.graph_min_edge_events,
        ),
        species_support_obs,
        species_idx,
    )

    ctx = _path_context(
        candidates=candidates,
        species_stats=species_stats,
        species_id_order=species_id_order,
        species_idx=species_idx,
        species_support_obs=species_support_obs,
        n_runs=n_runs,
        settings=settings,
    )
    path_nodes, path_edges, path_summary, scan_results = _highlighted_path(
        ctx, settings, candidates, species_idx
    )

    display_edges, selected_nodes, tie_edge_keys, clean_stats = _display_subgraph(
        candidates=candidates,
        species_stats=species_stats,
        species_id_order=species_id_order,
        species_idx=species_idx,
        path_nodes=path_nodes,
        path_edges=path_edges,
        scan_results=scan_results,
        settings=settings,
    )
    node_list, node_ids = dotgraph.node_id_map(selected_nodes, species_idx)
    scan_node_colors, scan_edge_colors, scan_entries = dotgraph.scan_overlay(
        scan_results
    )

    dot_text = dotgraph.build_dot(
        node_list=node_list,
        node_ids=node_ids,
        display_edges=display_edges,
        species_stats=species_stats,
        species_idx=species_idx,
        species_support=species_support,
        species_support_obs=species_support_obs,
        species_pct_obs=species_pct_obs,
        path_nodes=path_nodes,
        path_edges=path_edges,
        scan_node_colors=scan_node_colors,
        scan_edge_colors=scan_edge_colors,
        scan_entries=scan_entries,
        tie_edge_keys=tie_edge_keys,
    )

    dot_path = outdir / "aggregate_reacnet_network.dot"
    dot_path.write_text(dot_text, encoding="utf-8")
    dotgraph.write_node_legend(
        outdir / "aggregate_reacnet_network_nodes.tsv",
        node_list,
        species_idx,
        species_stats,
        species_support,
        species_support_obs,
        species_pct_obs,
    )
    dotgraph.write_mechanistic_path(
        outdir / "aggregate_reacnet_mechanistic_path.tsv", path_nodes, species_idx
    )
    graphviz_ok, svg_body = dotgraph.render_svg(
        dot_path, outdir / "aggregate_reacnet_network.svg"
    )

    report_block = _report_block(
        node_count=len(node_list),
        species_count=len(species_id_order),
        shown_edges=len(display_edges),
        candidate_edges=len(candidates),
        n_runs=n_runs,
        path_summary=path_summary,
        settings=settings,
        scan_results=scan_results,
    )
    if graphviz_ok:
        network_html = report_block + svg_body
    else:
        network_html = report_block + (
            "<div class='container py-3'><p>Graphviz render unavailable. See "
            "aggregate_reacnet_network.dot for the graph specification.</p></div>"
        )

    meta = _graph_meta(
        graphviz_ok=graphviz_ok,
        candidates=candidates,
        display_edges=display_edges,
        node_list=node_list,
        clean_stats=clean_stats,
        path_summary=path_summary,
        scan_results=scan_results,
        scan_entries=scan_entries,
        species_idx=species_idx,
        settings=settings,
    )
    return network_html, species_support_obs, meta


def _path_context(
    candidates: Sequence[Edge],
    species_stats: Dict[str, Dict[str, object]],
    species_id_order: Sequence[str],
    species_idx: Dict[str, int],
    species_support_obs: Dict[str, int],
    n_runs: int,
    settings: Settings,
) -> PathContext:
    _adjacency, outdeg, indeg, graph_nodes = build_adjacency(candidates, species_idx)
    in_events, out_events, _in_runs, _out_runs = compute_node_flow_stats(candidates)
    return PathContext(
        candidates=tuple(candidates),
        nodes=frozenset(graph_nodes),
        species_idx=species_idx,
        species_support_obs=species_support_obs,
        outdeg=outdeg,
        indeg=indeg,
        in_events=in_events,
        out_events=out_events,
        initial_start_candidates=tuple(
            s
            for s in species_id_order
            if bool(species_stats[s]["is_initial"]) and s in graph_nodes
        ),
        n_runs=n_runs,
        target_mode=settings.target_mode,
        path_min_steps=max(settings.path_min_steps, 0),
        top_targets=max(settings.top_targets, 1),
    )


def _display_subgraph(
    candidates: Sequence[Edge],
    species_stats: Dict[str, Dict[str, object]],
    species_id_order: Sequence[str],
    species_idx: Dict[str, int],
    path_nodes: Sequence[str],
    path_edges: Set[EdgeKey],
    scan_results: Dict[int, Dict[str, object]],
    settings: Settings,
) -> Tuple[List[Edge], Set[str], Set[EdgeKey], Dict[str, int]]:
    """Decide what is drawn, and net opposing pairs if asked to."""
    scan_nodes: Set[str] = set()
    scan_edges: Set[EdgeKey] = set()
    for row in scan_results.values():
        path = [str(x) for x in row["path"]]
        scan_nodes.update(path)
        scan_edges.update((path[i], path[i + 1]) for i in range(len(path) - 1))

    # The initial fragments and any highlighted route are always drawn, whatever
    # the display limits say; everything else competes for the remaining slots.
    mandatory_nodes = {
        s for s in species_id_order if bool(species_stats[s]["is_initial"])
    }
    mandatory_nodes.update(path_nodes)
    mandatory_nodes.update(scan_nodes)
    mandatory_edges = [
        edge for edge in candidates if (edge[0], edge[1]) in path_edges | scan_edges
    ]

    node_limit = max(DEFAULT_GRAPH_MAX_NODES, len(mandatory_nodes))
    edge_limit = max(DEFAULT_GRAPH_MAX_EDGES, len(mandatory_edges))
    if settings.pathway_mode == "none":
        # Nothing is being highlighted, so there is no reason to trim: show the
        # whole thresholded network.
        node_limit = max(node_limit, len(species_id_order))
        edge_limit = max(edge_limit, len(candidates))

    selected_edges, selected_nodes = select_display_subgraph(
        candidates=candidates,
        mandatory_nodes=mandatory_nodes,
        mandatory_edges=mandatory_edges,
        species_id_order=species_id_order,
        node_limit=node_limit,
        edge_limit=edge_limit,
        include_isolates=DEFAULT_GRAPH_INCLUDE_ISOLATES,
    )

    if not settings.clean:
        return list(selected_edges), selected_nodes, set(), {
            "pairs_collapsed": 0,
            "tie_pairs_kept": 0,
        }

    display_edges, tie_edge_keys, clean_stats = net_bidirectional_edges(
        selected_edges=selected_edges,
        mandatory_edge_keys={(u, v) for u, v, _e, _r in mandatory_edges},
        species_rank=species_idx,
    )
    return display_edges, selected_nodes, tie_edge_keys, clean_stats


def _graph_meta(
    graphviz_ok: bool,
    candidates: Sequence[Edge],
    display_edges: Sequence[Edge],
    node_list: Sequence[str],
    clean_stats: Dict[str, int],
    path_summary: Dict[str, object],
    scan_results: Dict[int, Dict[str, object]],
    scan_entries: Sequence[Tuple[int, str]],
    species_idx: Dict[str, int],
    settings: Settings,
) -> Dict[str, object]:
    return {
        "graphviz_ok": graphviz_ok,
        "candidate_edges": len(candidates),
        "shown_edges": len(display_edges),
        "clean_enabled": bool(settings.clean),
        "clean_pairs_collapsed": int(clean_stats["pairs_collapsed"]),
        "clean_zero_pairs_dropped": 0,
        "clean_tie_pairs_kept": int(clean_stats["tie_pairs_kept"]),
        "display_edges_after_clean": len(display_edges),
        "shown_nodes": len(node_list),
        "max_nodes": DEFAULT_GRAPH_MAX_NODES,
        "max_edges": DEFAULT_GRAPH_MAX_EDGES,
        "min_edge_runs": int(settings.graph_min_edge_runs),
        "min_edge_events": float(settings.graph_min_edge_events),
        "mechanistic_path": path_summary["description"],
        "mechanistic_path_min_edge_runs": path_summary["min_edge_runs"],
        "mechanistic_path_sum_edge_runs": path_summary["sum_edge_runs"],
        "mechanistic_path_sum_edge_events": path_summary["sum_edge_events"],
        "mechanistic_path_start_species": path_summary["start"],
        "mechanistic_path_target_species": path_summary["target"],
        "pathway_mode": settings.pathway_mode,
        "pathlength_scan": _scan_metadata(
            scan_results, scan_entries, species_idx, settings.pathway_mode
        ),
    }
