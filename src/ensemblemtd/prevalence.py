"""Species prevalence and reaction run support across the ensemble.

The metric the paper reports is

    Pi_i = 100% * N_i^runs / N_total^runs,   N_i^runs = sum_k I_ik,  I_ik in {0,1}

with I_ik = 1 when collapsed species i appears in at least one frame of
successful run k.  It is a node property and it is binary per run: a species
that persists for a whole trajectory and one that flickers for a single frame
both contribute exactly 1.  Pi_i therefore measures how reproducible a species
is across independent trajectories -- not its abundance, residence time, yield,
or branching ratio.
"""

from __future__ import annotations

from collections import Counter
from typing import Dict, Iterable, List, Sequence, Set, Tuple

ReactionKey = Tuple[Tuple[str, ...], Tuple[str, ...]]

EMPTY_SIDE = "∅"


def species_observation_stats(
    supporting_runs: Iterable[str], n_runs: int
) -> Tuple[int, float]:
    """Run support and prevalence for one collapsed species."""
    observed = len(set(supporting_runs))
    prevalence = 100.0 * observed / n_runs if n_runs else 0.0
    return observed, prevalence


def reaction_event_weight(rxn: dict) -> int:
    """Event count of a ReacNet reaction record, floored at one.

    Some records arrive without an ``n`` field; they still happened once.
    """
    n = int(rxn.get("n", 0) or 0)
    return n if n > 0 else 1


def add_reactions(
    reactions: Sequence[dict],
    run_name: str,
    reaction_events: Counter,
    reaction_run_support: Dict[ReactionKey, Set[str]],
    species_event_counter: Counter,
) -> None:
    """Fold one run's reaction records into the ensemble tallies.

    Sides are sorted so that the same reaction written in either member order
    lands on one key, and records with identical sides are dropped.
    """
    for rxn in reactions:
        lhs = tuple(sorted(x for x in rxn.get("l", []) if isinstance(x, str)))
        rhs = tuple(sorted(x for x in rxn.get("r", []) if isinstance(x, str)))
        if not lhs or not rhs or lhs == rhs:
            continue
        n = reaction_event_weight(rxn)
        key: ReactionKey = (lhs, rhs)
        reaction_events[key] += n
        reaction_run_support[key].add(run_name)
        for species in set(lhs + rhs):
            species_event_counter[species] += n


def build_reaction_rows(
    reaction_events: Counter,
    reaction_run_support: Dict[ReactionKey, Set[str]],
    n_runs: int,
) -> List[Tuple[ReactionKey, int, float, int]]:
    """Reaction table ordered by run support, then by total events."""
    rows: List[Tuple[ReactionKey, int, float, int]] = []
    for key, total_events in reaction_events.items():
        support = len(reaction_run_support[key])
        pct = 100.0 * support / n_runs if n_runs else 0.0
        rows.append((key, support, pct, int(total_events)))
    rows.sort(key=lambda row: (-row[1], -row[3], row[0][0], row[0][1]))
    return rows


def side_key(side: Sequence[str]) -> str:
    """Render one side of a reaction, with stoichiometry as "2*X"."""
    if not side:
        return EMPTY_SIDE
    counts = Counter(side)
    parts = []
    for species in sorted(counts):
        n = counts[species]
        parts.append(f"{n}*{species}" if n > 1 else species)
    return " + ".join(parts)


def key_to_str(key: ReactionKey) -> str:
    lhs, rhs = key
    return f"{side_key(lhs)} => {side_key(rhs)}"
