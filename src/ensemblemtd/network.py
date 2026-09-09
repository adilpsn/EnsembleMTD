"""Projecting collapsed reaction records onto a directed species graph.

A ReacNet record can have several species on each side, but the consensus figure
is a graph over single species, so each record is split over the ordered pairs
it implies:

    e_uv = sum over records q with u in L_q, v in R_q, u != v  of  n_q/(|L_q||R_q|)

The denominator counts unique species on both sides *before* self-pairs are
dropped, so a record's weight is conserved across the pairs it is shared over
and e_uv is generally not an integer.  The companion quantity r_uv is the number
of distinct successful runs contributing at least one record to u -> v.  Both are
edge properties: a species with high prevalence Pi_i can still have low run
support on any single formation edge.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Set, Tuple

ReactionKey = Tuple[Tuple[str, ...], Tuple[str, ...]]
EdgeKey = Tuple[str, str]

# (u, v, event weight, run support)
Edge = Tuple[str, str, float, int]

# Sorting sentinel for a species that carries no index.
UNRANKED = 10 ** 9

TIE_TOLERANCE = 1e-15


def project_reactions_to_edges(
    reaction_rows_full: Sequence[Tuple[ReactionKey, int, float, int]],
    reaction_run_support: Dict[ReactionKey, Set[str]],
) -> Tuple[Dict[EdgeKey, float], Dict[EdgeKey, Set[str]]]:
    """Split every reaction record over its directed species pairs."""
    edge_events: Dict[EdgeKey, float] = defaultdict(float)
    edge_run_names: Dict[EdgeKey, Set[str]] = defaultdict(set)

    for key, _support, _pct, n_events in reaction_rows_full:
        lhs, rhs = key
        left = sorted(set(lhs))
        right = sorted(set(rhs))
        if not left or not right:
            continue
        share = float(n_events) / float(max(len(left) * len(right), 1))
        support_names = reaction_run_support.get(key, set())
        for u in left:
            for v in right:
                if u == v:
                    continue
                edge_events[(u, v)] += share
                edge_run_names[(u, v)].update(support_names)

    return edge_events, edge_run_names


def filter_edge_candidates(
    edge_events: Dict[EdgeKey, float],
    edge_run_names: Dict[EdgeKey, Set[str]],
    min_edge_runs: int,
    min_edge_events: float,
) -> List[Edge]:
    """Keep the edges that clear both the run-support and event-weight floors."""
    candidates: List[Edge] = []
    for (u, v), event_weight in edge_events.items():
        run_support = len(edge_run_names[(u, v)])
        if run_support < min_edge_runs or event_weight < min_edge_events:
            continue
        candidates.append((u, v, event_weight, run_support))
    return candidates


def sort_candidates(
    candidates: Sequence[Edge],
    species_support_obs: Dict[str, int],
    species_idx: Dict[str, int],
) -> List[Edge]:
    """Strongest first: run support, event weight, then species prevalence."""
    return sorted(
        candidates,
        key=lambda edge: (
            -edge[3],
            -edge[2],
            -species_support_obs.get(edge[0], 0),
            -species_support_obs.get(edge[1], 0),
            species_idx.get(edge[0], UNRANKED),
            species_idx.get(edge[1], UNRANKED),
        ),
    )


def compute_node_flow_stats(
    candidates: Sequence[Edge],
) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, int], Dict[str, int]]:
    """In/out event weight and in/out summed run support per node."""
    in_events: Dict[str, float] = {}
    out_events: Dict[str, float] = {}
    in_runs: Dict[str, int] = {}
    out_runs: Dict[str, int] = {}
    for u, v, events, support in candidates:
        out_events[u] = out_events.get(u, 0.0) + float(events)
        in_events[v] = in_events.get(v, 0.0) + float(events)
        out_runs[u] = out_runs.get(u, 0) + int(support)
        in_runs[v] = in_runs.get(v, 0) + int(support)
    return in_events, out_events, in_runs, out_runs


def build_adjacency(
    candidates: Sequence[Edge], species_idx: Dict[str, int]
) -> Tuple[Dict[str, List[Tuple[str, int, float]]], Dict[str, int], Dict[str, int], Set[str]]:
    """Adjacency list plus degrees and the node set the edges span."""
    adjacency: Dict[str, List[Tuple[str, int, float]]] = defaultdict(list)
    outdeg: Dict[str, int] = defaultdict(int)
    indeg: Dict[str, int] = defaultdict(int)
    for u, v, events, support in candidates:
        adjacency[u].append((v, support, float(events)))
        outdeg[u] += 1
        indeg[v] += 1
    for u in adjacency:
        adjacency[u].sort(key=lambda item: species_idx.get(item[0], UNRANKED))

    nodes = set(adjacency)
    for neighbours in adjacency.values():
        for v, _support, _events in neighbours:
            nodes.add(v)
    return adjacency, outdeg, indeg, nodes


def _collapse_opposing_pair(edges: Sequence[Edge]) -> Optional[Tuple[Edge, bool]]:
    """Net one u->v / v->u pair, or return None if this is not such a pair.

    The second element of the result says whether the two directions tied.
    """
    if len(edges) != 2:
        return None
    forward, reverse = edges
    if forward[0] != reverse[1] or forward[1] != reverse[0]:
        return None

    support = max(int(forward[3]), int(reverse[3]))
    if abs(forward[2] - reverse[2]) <= TIE_TOLERANCE:
        u, v = sorted((forward[0], forward[1]))
        return (u, v, 0.0, support), True
    stronger, weaker = (
        (forward, reverse) if forward[2] > reverse[2] else (reverse, forward)
    )
    return (stronger[0], stronger[1], float(stronger[2] - weaker[2]), support), False


def net_bidirectional_edges(
    selected_edges: Sequence[Edge],
    mandatory_edge_keys: Set[EdgeKey],
    species_rank: Optional[Dict[str, int]] = None,
) -> Tuple[List[Edge], Set[EdgeKey], Dict[str, int]]:
    """Collapse opposing edge pairs for the cleaned display.

    An opposing pair u->v and v->u becomes one edge with e_net = |e_uv - e_vu|
    pointing the way the larger event weight points, and r = max(r_uv, r_vu).
    An exact tie is kept as a bidirectional edge with e = 0, which is a real
    statement about the ensemble rather than a rendering artefact: the two
    directions were seen equally often.

    Edges on a highlighted path are exempt so the path stays readable.
    """
    mandatory = [e for e in selected_edges if (e[0], e[1]) in mandatory_edge_keys]
    background = [e for e in selected_edges if (e[0], e[1]) not in mandatory_edge_keys]

    pair_groups: Dict[Tuple[str, str], List[Edge]] = defaultdict(list)
    self_loops: List[Edge] = []
    for u, v, events, support in background:
        if u == v:
            self_loops.append((u, v, events, support))
            continue
        pair_groups[(u, v) if u < v else (v, u)].append((u, v, events, support))

    cleaned: List[Edge] = []
    tie_edge_keys: Set[EdgeKey] = set()
    pairs_collapsed = 0
    tie_pairs_kept = 0

    for edges in pair_groups.values():
        collapsed = _collapse_opposing_pair(edges)
        if collapsed is None:
            cleaned.extend(edges)
            continue
        edge, is_tie = collapsed
        cleaned.append(edge)
        pairs_collapsed += 1
        if is_tie:
            tie_edge_keys.add((edge[0], edge[1]))
            tie_pairs_kept += 1

    cleaned.extend(self_loops)
    rank = species_rank or {}
    cleaned.sort(
        key=lambda edge: (
            -edge[3],
            -edge[2],
            rank.get(edge[0], UNRANKED),
            rank.get(edge[1], UNRANKED),
            edge[0],
            edge[1],
        )
    )
    stats = {"pairs_collapsed": pairs_collapsed, "tie_pairs_kept": tie_pairs_kept}
    return mandatory + cleaned, tie_edge_keys, stats


def select_display_subgraph(
    candidates: Sequence[Edge],
    mandatory_nodes: Set[str],
    mandatory_edges: Sequence[Edge],
    species_id_order: Sequence[str],
    node_limit: int,
    edge_limit: int,
    include_isolates: bool,
) -> Tuple[List[Edge], Set[str]]:
    """Grow the shown subgraph from the mandatory core, strongest edges first."""
    mandatory_keys = {(u, v) for u, v, _e, _r in mandatory_edges}
    selected_edges: List[Edge] = list(mandatory_edges)
    selected_nodes: Set[str] = set(mandatory_nodes)

    for u, v, events, support in candidates:
        if (u, v) in mandatory_keys:
            continue
        if len(selected_edges) >= edge_limit:
            break
        grown = selected_nodes | {u, v}
        if len(grown) > node_limit:
            continue
        selected_edges.append((u, v, events, support))
        selected_nodes = grown

    if include_isolates:
        for species in species_id_order:
            if len(selected_nodes) >= node_limit:
                break
            selected_nodes.add(species)

    if not selected_nodes:
        selected_nodes.update(species_id_order[:node_limit])

    return selected_edges, selected_nodes


def build_linkreac(
    reaction_rows_full: Sequence[Tuple[ReactionKey, int, float, int]],
    reaction_rows_abcd_full: Sequence[Tuple[ReactionKey, int, float, int]],
    species_order: Sequence[str],
) -> Dict[str, List[str]]:
    """Undirected neighbour lists for the interactive ReacNet-style report."""
    adjacency: Dict[str, set] = defaultdict(set)
    for key, _support, _pct, _events in list(reaction_rows_full) + list(
        reaction_rows_abcd_full
    ):
        lhs, rhs = key
        for u in set(lhs):
            for v in set(rhs):
                if u == v:
                    continue
                adjacency[u].add(v)
                adjacency[v].add(u)
    for species in species_order:
        adjacency.setdefault(species, set())
    return {k: sorted(v) for k, v in sorted(adjacency.items())}
