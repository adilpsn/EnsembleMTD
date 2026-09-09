"""Rendering the consensus network as Graphviz DOT, and to SVG if dot is there.

Edge labels carry both numbers the network is built on, ``e`` (event weight) and
``r`` (run support), because they say different things and a reader needs both.
Pen width follows sqrt(e/e_max) so a dominant edge does not flatten the rest.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

from .network import Edge, EdgeKey, UNRANKED

# Colours for the optional path-length scan overlay, cycled if the scan is long.
SCAN_COLORS = ["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"]

INITIAL_FILL, INITIAL_EDGE = "#eef2ff", "#4f46e5"
PLAIN_FILL, PLAIN_EDGE = "#f8fafc", "#94a3b8"
PATH_FILL, PATH_EDGE = "#fee2e2", "#ef4444"
PATH_ARROW = "#dc2626"
TIE_ARROW = "#6b7280"
DEFAULT_ARROW = "#64748b"

MIN_PEN, PEN_SPAN = 1.0, 4.0
TIE_MIN_PEN, SCAN_MIN_PEN, PATH_MIN_PEN = 2.2, 3.8, 4.5

DOT_PREAMBLE = [
    "digraph G {",
    "  rankdir=LR;",
    '  graph [bgcolor="white", labelloc="t", labeljust="l", fontsize=16, fontname="Helvetica"];',
    '  node [shape=circle, style="filled", fillcolor="#f8fafc", color="#94a3b8", fontname="Helvetica", fontsize=10, fixedsize=true, width=0.45];',
    '  edge [color="#64748b", fontname="Helvetica", fontsize=9, arrowsize=0.7];',
]


def escape_dot(text: str) -> str:
    return text.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ")


def scan_overlay(
    path_scan_results: Dict[int, Dict[str, object]],
) -> Tuple[Dict[str, str], Dict[EdgeKey, str], List[Tuple[int, str]]]:
    """Assign one colour per scanned path length to its nodes and edges."""
    node_colors: Dict[str, str] = {}
    edge_colors: Dict[EdgeKey, str] = {}
    entries: List[Tuple[int, str]] = []
    for index, length in enumerate(sorted(path_scan_results)):
        color = SCAN_COLORS[index % len(SCAN_COLORS)]
        path = path_scan_results[length]["path"]
        entries.append((length, color))
        for node in path:
            node_colors.setdefault(node, color)
        for i in range(len(path) - 1):
            edge_colors.setdefault((path[i], path[i + 1]), color)
    return node_colors, edge_colors, entries


def node_lines(
    node_list: Sequence[str],
    node_ids: Dict[str, str],
    species_stats: Dict[str, Dict[str, object]],
    species_idx: Dict[str, int],
    species_support: Dict[str, int],
    species_support_obs: Dict[str, int],
    species_pct_obs: Dict[str, float],
    path_nodes: Sequence[str],
    scan_node_colors: Dict[str, str],
) -> List[str]:
    lines = []
    for species in node_list:
        sid = species_idx.get(species, 0)
        is_initial = bool(species_stats[species]["is_initial"])
        tooltip = (
            f"S{sid} | runs={species_support.get(species, 0)} "
            f"| obs={species_support_obs.get(species, 0)} "
            f"| pct={species_pct_obs.get(species, 0.0):.1f}% "
            f"| initial={'yes' if is_initial else 'no'} | {species}"
        )
        if species in path_nodes:
            fill, color = PATH_FILL, PATH_EDGE
        elif species in scan_node_colors:
            fill, color = PLAIN_FILL, scan_node_colors[species]
        elif is_initial:
            fill, color = INITIAL_FILL, INITIAL_EDGE
        else:
            fill, color = PLAIN_FILL, PLAIN_EDGE
        lines.append(
            f'  {node_ids[species]} [label="{escape_dot(f"S{sid}")}", '
            f'tooltip="{escape_dot(tooltip)}", fillcolor="{fill}", color="{color}"];'
        )
    return lines


def edge_lines(
    display_edges: Sequence[Edge],
    node_ids: Dict[str, str],
    path_edges: Set[EdgeKey],
    scan_edge_colors: Dict[EdgeKey, str],
    tie_edge_keys: Set[EdgeKey],
) -> List[str]:
    max_events = max((e[2] for e in display_edges), default=0.0) or 1.0
    lines = []
    for u, v, events, support in display_edges:
        if u not in node_ids or v not in node_ids:
            continue
        on_path = (u, v) in path_edges
        on_scan = (u, v) in scan_edge_colors
        is_tie = (u, v) in tie_edge_keys

        pen = MIN_PEN + PEN_SPAN * (events / max_events) ** 0.5
        if is_tie:
            pen = max(pen, TIE_MIN_PEN)
        if on_path:
            pen = max(pen, PATH_MIN_PEN)
        elif on_scan:
            pen = max(pen, SCAN_MIN_PEN)

        if on_path:
            color = PATH_ARROW
        elif on_scan:
            color = scan_edge_colors[(u, v)]
        elif is_tie:
            color = TIE_ARROW
        else:
            color = DEFAULT_ARROW

        # A tie is drawn as a dashed double arrow: the two directions carried
        # the same event weight, so no net direction can be claimed.
        extra = (
            ' style="dashed", dir="both", arrowhead="normal", arrowtail="normal"'
            if is_tie
            else ""
        )
        label = f"e={events:.1f}, r={support}"
        tooltip = f"events={events:.2f}; run_support={support}; {u} -> {v}"
        lines.append(
            f'  {node_ids[u]} -> {node_ids[v]} [label="{escape_dot(label)}", '
            f"penwidth={pen:.2f}, color=\"{color}\", fontcolor=\"{color}\", "
            f'tooltip="{escape_dot(tooltip)}"{extra}];'
        )
    return lines


def scan_legend_lines(scan_entries: Sequence[Tuple[int, str]]) -> List[str]:
    if not scan_entries:
        return []
    lines = [
        "  subgraph cluster_scan_legend {",
        '    label="Path-length scan";',
        '    color="#e2e8f0"; fontsize=10;',
    ]
    for length, color in scan_entries:
        lines.append(
            f'    legend_L{length} [shape=plaintext, label="L={length}", '
            f'fontcolor="{color}"];'
        )
    for i in range(len(scan_entries) - 1):
        lines.append(
            f"    legend_L{scan_entries[i][0]} -> "
            f"legend_L{scan_entries[i + 1][0]} [style=invis];"
        )
    lines.append("  }")
    return lines


def render_svg(dot_path: Path, svg_path: Path) -> Tuple[bool, str]:
    """Run ``dot -Tsvg``.  Returns (succeeded, inline SVG body).

    Graphviz is optional; the .dot file is the artefact that matters and can be
    rendered later anywhere.
    """
    try:
        proc = subprocess.run(
            ["dot", "-Tsvg", str(dot_path), "-o", str(svg_path)],
            capture_output=True,
            text=True,
        )
    except OSError:
        return False, ""
    if proc.returncode != 0 or not svg_path.exists():
        return False, ""
    raw = svg_path.read_text(encoding="utf-8", errors="ignore")
    start = raw.find("<svg")
    body = raw[start:] if start >= 0 else raw
    return ("<svg" in body), body


def write_node_legend(
    path: Path,
    node_list: Sequence[str],
    species_idx: Dict[str, int],
    species_stats: Dict[str, Dict[str, object]],
    species_support: Dict[str, int],
    species_support_obs: Dict[str, int],
    species_pct_obs: Dict[str, float],
) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(
            "node_id\tsmiles\tdisplay_run_support\tobserved_run_support\t"
            "pct_runs_observed\tis_initial\tinitial_frame_count\n"
        )
        for species in node_list:
            stats = species_stats[species]
            f.write(
                f"S{species_idx.get(species, 0)}\t{species}\t"
                f"{species_support.get(species, 0)}\t"
                f"{species_support_obs.get(species, 0)}\t"
                f"{species_pct_obs.get(species, 0.0):.3f}\t"
                f"{int(bool(stats['is_initial']))}\t"
                f"{int(stats['initial_frame_count'])}\n"
            )


def write_mechanistic_path(
    path: Path, path_nodes: Sequence[str], species_idx: Dict[str, int]
) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("step\tnode_id\tsmiles\n")
        for step, species in enumerate(path_nodes, start=1):
            f.write(f"{step}\tS{species_idx.get(species, 0)}\t{species}\n")


def node_id_map(
    selected_nodes: Set[str], species_idx: Dict[str, int]
) -> Tuple[List[str], Dict[str, str]]:
    node_list = sorted(selected_nodes, key=lambda s: species_idx.get(s, UNRANKED))
    return node_list, {s: f"N{i}" for i, s in enumerate(node_list)}


def build_dot(
    node_list: Sequence[str],
    node_ids: Dict[str, str],
    display_edges: Sequence[Edge],
    species_stats: Dict[str, Dict[str, object]],
    species_idx: Dict[str, int],
    species_support: Dict[str, int],
    species_support_obs: Dict[str, int],
    species_pct_obs: Dict[str, float],
    path_nodes: Sequence[str],
    path_edges: Set[EdgeKey],
    scan_node_colors: Dict[str, str],
    scan_edge_colors: Dict[EdgeKey, str],
    scan_entries: Sequence[Tuple[int, str]],
    tie_edge_keys: Set[EdgeKey],
) -> str:
    lines = list(DOT_PREAMBLE)
    lines += node_lines(
        node_list,
        node_ids,
        species_stats,
        species_idx,
        species_support,
        species_support_obs,
        species_pct_obs,
        path_nodes,
        scan_node_colors,
    )
    lines += edge_lines(
        display_edges, node_ids, path_edges, scan_edge_colors, tie_edge_keys
    )
    lines += scan_legend_lines(scan_entries)
    lines.append("}")
    return "\n".join(lines) + "\n"
