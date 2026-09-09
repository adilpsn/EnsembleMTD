"""Optional pathway highlighting on the consensus graph.

This is a reading aid, not a result.  The consensus network already carries all
the evidence (e_uv and r_uv per edge, Pi_i per node); the highlighter picks one
route through it so a figure has something to follow, under one of two rules:

``support-first``
    Maximise the weakest edge's run support along the route, then take the
    strongest shortest path that survives at that support threshold.  This
    answers "which chain of steps is reproduced in the most trajectories".

``probabilistic``
    Minimise the cumulative -log P_eff with
    P_eff(u->v) = P(u->v) * r_uv/N and P(u->v) = e_uv / sum_k e_uk.
    Branching weights and run support are folded into one cost.

Neither is a kinetic path and neither yields barriers or rates.  The default is
``none``: no route is drawn.
"""

from __future__ import annotations

import heapq
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .network import UNRANKED, Edge, compute_node_flow_stats

DIJKSTRA_TOLERANCE = 1e-15

# Floor on an effective probability, so a zero-weight edge stays traversable at
# a large but finite cost instead of producing -log(0).
MIN_EFFECTIVE_PROBABILITY = 1e-12

# Guard against dividing by a node with no outgoing weight.
EPSILON = 1e-12


@dataclass(frozen=True)
class PathContext:
    """Everything the two selection rules read, fixed up front."""

    candidates: Tuple[Edge, ...]
    nodes: frozenset
    species_idx: Dict[str, int]
    species_support_obs: Dict[str, int]
    outdeg: Dict[str, int]
    indeg: Dict[str, int]
    in_events: Dict[str, float]
    out_events: Dict[str, float]
    initial_start_candidates: Tuple[str, ...]
    n_runs: int
    target_mode: str
    path_min_steps: int
    top_targets: int

    def rank(self, species: str) -> int:
        return self.species_idx.get(species, UNRANKED)

    def by_index(self, species: Sequence[str]) -> List[str]:
        return sorted(species, key=self.rank)

    def dominant_node(self) -> str:
        """Node with the highest prevalence, earliest species id breaking ties."""
        return max(
            self.nodes,
            key=lambda s: (self.species_support_obs.get(s, 0), -self.rank(s)),
        )


def score_better(
    lhs: Tuple[int, int, float, int], rhs: Optional[Tuple[int, int, float, int]]
) -> bool:
    return rhs is None or lhs > rhs


def select_target_candidates(
    nodes: Sequence[str],
    start: str,
    outdeg: Dict[str, int],
    indeg: Dict[str, int],
    target_mode: str,
    in_events: Dict[str, float],
    out_events: Dict[str, float],
) -> List[str]:
    """Which nodes may end a highlighted route.

    ``strict-sink`` wants a true terminal product (something consumes nothing
    further) and falls back to every node when the graph has no sink, which is
    common once reverse edges survive the thresholds.  ``pseudo-sink`` accepts
    anything with incoming weight.
    """
    if target_mode == "strict-sink":
        sinks = [
            s
            for s in nodes
            if s != start and outdeg.get(s, 0) == 0 and indeg.get(s, 0) > 0
        ]
        return sinks if sinks else [s for s in nodes if s != start]
    if target_mode == "pseudo-sink":
        return [s for s in nodes if s != start and in_events.get(s, 0.0) > 0.0]
    return [s for s in nodes if s != start]


def reconstruct_path(prev: Dict[str, str], src: str, dst: str) -> Optional[List[str]]:
    if src == dst:
        return [src]
    if dst not in prev:
        return None
    path = [dst]
    current = dst
    while current != src:
        current = prev.get(current)
        if current is None:
            return None
        path.append(current)
    path.reverse()
    return path


def dijkstra_cost(
    weighted_adj: Dict[str, List[Tuple[str, float]]], src: str
) -> Tuple[Dict[str, float], Dict[str, str]]:
    dist: Dict[str, float] = {src: 0.0}
    prev: Dict[str, str] = {}
    heap: List[Tuple[float, str]] = [(0.0, src)]
    while heap:
        d, u = heapq.heappop(heap)
        if d > dist.get(u, float("inf")) + DIJKSTRA_TOLERANCE:
            continue
        for v, w in weighted_adj.get(u, []):
            nd = d + w
            if nd + DIJKSTRA_TOLERANCE < dist.get(v, float("inf")):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(heap, (nd, v))
    return dist, prev


def best_paths_for_threshold(
    adjacency: Dict[str, List[Tuple[str, int, float]]], src: str
) -> Tuple[Dict[str, int], Dict[str, Tuple[int, int, float, int]], Dict[str, str]]:
    """Breadth-first layering, then the best-supported route into each layer.

    The score compared per node is (weakest edge support, summed support,
    summed event weight, -depth), maximised lexicographically.
    """
    dist: Dict[str, int] = {src: 0}
    queue: List[str] = [src]
    head = 0
    while head < len(queue):
        u = queue[head]
        head += 1
        for v, _support, _events in adjacency.get(u, []):
            if v not in dist:
                dist[v] = dist[u] + 1
                queue.append(v)

    best_score: Dict[str, Tuple[int, int, float, int]] = {src: (0, 0, 0.0, 0)}
    prev: Dict[str, str] = {}
    for u in sorted(dist, key=lambda s: (dist[s], s)):
        base = best_score.get(u)
        if base is None:
            continue
        for v, support, events in adjacency.get(u, []):
            if dist.get(v) != dist[u] + 1:
                continue
            candidate = (
                support if u == src else min(base[0], support),
                base[1] + support,
                base[2] + events,
                -(dist[u] + 1),
            )
            if score_better(candidate, best_score.get(v)):
                best_score[v] = candidate
                prev[v] = u
    return dist, best_score, prev


def pseudo_sink_scores(
    nodes, in_events: Dict[str, float], out_events: Dict[str, float]
) -> Dict[str, float]:
    """Ratio of incoming to outgoing event weight; high means product-like."""
    return {
        s: in_events.get(s, 0.0) / (out_events.get(s, 0.0) + EPSILON)
        for s in nodes
    }


def choose_probabilistic(
    ctx: PathContext, exact_steps: Optional[int] = None
) -> Optional[Dict[str, object]]:
    if not ctx.candidates or not ctx.nodes:
        return None

    start = ctx.dominant_node()
    weighted_adj: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
    for u, v, events, support in ctx.candidates:
        branching = (
            events / ctx.out_events[u] if ctx.out_events.get(u, 0.0) > 0.0 else 0.0
        )
        support_frac = support / max(float(ctx.n_runs), 1.0)
        p_eff = max(
            branching * max(support_frac, MIN_EFFECTIVE_PROBABILITY),
            MIN_EFFECTIVE_PROBABILITY,
        )
        weighted_adj[u].append((v, -math.log(p_eff)))
    for u in weighted_adj:
        weighted_adj[u].sort(key=lambda item: ctx.rank(item[0]))

    targets = select_target_candidates(
        ctx.by_index(ctx.nodes),
        start,
        ctx.outdeg,
        ctx.indeg,
        ctx.target_mode,
        ctx.in_events,
        ctx.out_events,
    )
    sink_scores = pseudo_sink_scores(ctx.nodes, ctx.in_events, ctx.out_events)

    dist, prev = dijkstra_cost(weighted_adj, start)
    ranked: List[Tuple[Tuple[float, ...], str, List[str]]] = []
    for target in targets:
        if target not in dist:
            continue
        path = reconstruct_path(prev, start, target)
        if not path or len(path) < 2:
            continue
        steps = max(len(path) - 1, 0)
        if steps < ctx.path_min_steps:
            continue
        if exact_steps is not None and steps != exact_steps:
            continue
        if ctx.target_mode == "pseudo-sink":
            score = (
                dist[target],
                -sink_scores.get(target, 0.0),
                -ctx.species_support_obs.get(target, 0),
                len(path),
                ctx.rank(target),
            )
        else:
            score = (
                dist[target],
                -ctx.species_support_obs.get(target, 0),
                len(path),
                ctx.rank(target),
            )
        ranked.append((score, target, path))

    if not ranked:
        return None
    ranked.sort(key=lambda item: item[0])
    best_score, best_target, best_path = ranked[0]
    return {
        "start": start,
        "target": best_target,
        "path": best_path,
        "score": best_score,
        "cost": float(dist[best_target]),
        "top": ranked[: ctx.top_targets],
    }


def choose_support_first(
    ctx: PathContext, exact_steps: Optional[int] = None
) -> Optional[Dict[str, object]]:
    start_candidates = list(ctx.initial_start_candidates)
    if not start_candidates and ctx.nodes:
        fallback = ctx.dominant_node()
        print(
            "[warn] No initial fragments survived path-graph filtering; falling "
            f"back to dominant node S{ctx.rank(fallback)}.",
            file=sys.stderr,
        )
        start_candidates = [fallback]
    if not start_candidates or not ctx.candidates:
        return None

    # Walk the distinct run-support values downwards and stop at the first one
    # that still connects a start to a target: that is the "support-first" part.
    thresholds = sorted({support for _u, _v, _e, support in ctx.candidates}, reverse=True)
    for threshold in thresholds:
        adjacency: Dict[str, List[Tuple[str, int, float]]] = defaultdict(list)
        edges: List[Edge] = []
        nodes: Set[str] = set()
        outdeg: Dict[str, int] = defaultdict(int)
        indeg: Dict[str, int] = defaultdict(int)
        for u, v, events, support in ctx.candidates:
            if support < threshold:
                continue
            adjacency[u].append((v, support, float(events)))
            edges.append((u, v, float(events), int(support)))
            nodes.update((u, v))
            outdeg[u] += 1
            indeg[v] += 1
        for u in adjacency:
            adjacency[u].sort(key=lambda item: ctx.rank(item[0]))

        in_events, out_events, _in_runs, _out_runs = compute_node_flow_stats(edges)
        sink_scores = pseudo_sink_scores(nodes, in_events, out_events)

        ranked: List[Tuple[Tuple[float, ...], str, str, List[str]]] = []
        for start in start_candidates:
            if start not in nodes:
                continue
            dist, best_score, prev = best_paths_for_threshold(adjacency, start)
            targets = select_target_candidates(
                ctx.by_index(nodes),
                start,
                outdeg,
                indeg,
                ctx.target_mode,
                in_events,
                out_events,
            )
            for target in targets:
                if target not in dist or target not in best_score:
                    continue
                path = reconstruct_path(prev, start, target)
                if not path or len(path) < 2:
                    continue
                steps = max(len(path) - 1, 0)
                if steps < ctx.path_min_steps:
                    continue
                if exact_steps is not None and steps != exact_steps:
                    continue
                score = best_score[target]
                if ctx.target_mode == "pseudo-sink":
                    rank = (
                        float(threshold),
                        float(score[1]),
                        float(score[2]),
                        sink_scores.get(target, 0.0),
                        float(ctx.species_support_obs.get(target, 0)),
                        float(score[3]),
                        -float(ctx.rank(target)),
                    )
                else:
                    is_sink = outdeg.get(target, 0) == 0 and indeg.get(target, 0) > 0
                    rank = (
                        float(threshold),
                        float(score[1]),
                        float(score[2]),
                        float(ctx.species_support_obs.get(target, 0)),
                        1.0 if is_sink else 0.0,
                        float(score[3]),
                        -float(ctx.rank(target)),
                    )
                ranked.append((rank, start, target, path))

        if ranked:
            ranked.sort(key=lambda item: item[0], reverse=True)
            rank, start, target, path = ranked[0]
            return {"start": start, "target": target, "path": path, "score": rank}

    return None


def choose_path(
    ctx: PathContext, pathway_mode: str, exact_steps: Optional[int] = None
) -> Optional[Dict[str, object]]:
    if pathway_mode == "probabilistic":
        return choose_probabilistic(ctx, exact_steps)
    if pathway_mode == "support-first":
        return choose_support_first(ctx, exact_steps)
    return None
