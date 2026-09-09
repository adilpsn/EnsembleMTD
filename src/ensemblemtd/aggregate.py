"""Collapsing the raw ensemble tallies and deriving the species table.

The per-run loop in :mod:`ensemblemtd.pipeline` accumulates raw ReacNet species
and reaction records; everything here consumes that accumulator without
modifying it and returns fresh mappings keyed on collapsed representatives.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple

from .collapse import (
    build_inchi_collapse_maps,
    collapse_counter,
    collapse_first_seen_order,
    collapse_reaction_maps,
    collapse_ratio,
    collapse_support_map,
    get_obabel_version,
)
from .prevalence import (
    ReactionKey,
    build_reaction_rows,
    species_observation_stats,
)

# Sentinel ordering position for a species that never appeared in a first frame.
NEVER_SEEN = (10 ** 9, 10 ** 9)

# Fraction of species whose JSON support undercounts the timeline before the
# discrepancy is worth a warning rather than a note.
SUPPORT_AUDIT_WARN_FRACTION = 0.10


@dataclass
class RunHarvest:
    """Raw tallies collected while walking the ensemble, before collapsing.

    Accumulated in place by the per-run loop, then read-only from that point on.
    """

    reaction_events: Counter = field(default_factory=Counter)
    reaction_events_abcd: Counter = field(default_factory=Counter)
    reaction_run_support: Dict[ReactionKey, Set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )
    reaction_run_support_abcd: Dict[ReactionKey, Set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )
    species_run_support_timeline: Dict[str, Set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )
    species_run_support_json: Dict[str, Set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )
    species_event_counter: Counter = field(default_factory=Counter)
    initial_species_counts: Counter = field(default_factory=Counter)
    initial_run_support: Dict[str, Set[str]] = field(
        default_factory=lambda: defaultdict(set)
    )
    initial_first_seen_order: Dict[str, Tuple[int, int]] = field(default_factory=dict)

    def species_universe(self) -> Set[str]:
        return (
            set(self.species_run_support_timeline)
            | set(self.species_run_support_json)
            | set(self.species_event_counter)
            | set(self.initial_species_counts)
        )


@dataclass(frozen=True)
class CollapsedEnsemble:
    """The ensemble after InChI collapsing, ready to be tabulated."""

    species_run_support: Dict[str, Set[str]]
    species_run_support_json: Dict[str, Set[str]]
    species_run_support_timeline: Dict[str, Set[str]]
    species_event_counter: Counter
    initial_species_counts: Counter
    initial_first_seen_order: Dict[str, Tuple[int, int]]
    first_frame_support: Dict[str, int]
    reaction_rows_full: List[Tuple[ReactionKey, int, float, int]]
    reaction_rows_abcd_full: List[Tuple[ReactionKey, int, float, int]]
    reaction_run_support: Dict[ReactionKey, Set[str]]
    rep_to_inchi: Dict[str, str]
    rep_to_members: Dict[str, List[str]]
    n_raw_species: int
    n_collapsed_species: int
    obabel_version: str

    @property
    def collapse_ratio(self) -> float:
        return collapse_ratio(self.n_collapsed_species, self.n_raw_species)

    def universe(self) -> Set[str]:
        """Every species named by any tally, whichever support source is active."""
        return (
            set(self.species_run_support)
            | set(self.species_run_support_json)
            | set(self.species_run_support_timeline)
            | set(self.species_event_counter)
            | set(self.initial_species_counts)
        )


def collapse_harvest(
    harvest: RunHarvest,
    n_runs: int,
    support_source: str,
    obabel_bin: str,
) -> CollapsedEnsemble:
    """Group raw species by fixed-H InChI and re-key every tally onto groups."""
    raw_universe = harvest.species_universe()
    n_raw = len(raw_universe)

    obabel_version = ""
    raw_to_rep: Dict[str, str] = {}
    rep_to_inchi: Dict[str, str] = {}
    rep_to_members: Dict[str, List[str]] = {}

    if raw_universe:
        obabel_version = get_obabel_version(obabel_bin)
        reference = (
            harvest.species_run_support_json
            if support_source == "json"
            else harvest.species_run_support_timeline
        )
        raw_to_rep, _raw_to_inchi, rep_to_inchi, rep_to_members = (
            build_inchi_collapse_maps(
                raw_species_universe=raw_universe,
                raw_support_for_rep={
                    s: len(reference.get(s, set())) for s in raw_universe
                },
                obabel_bin=obabel_bin,
            )
        )

    if raw_universe:
        timeline = collapse_support_map(harvest.species_run_support_timeline, raw_to_rep)
        json_support = collapse_support_map(harvest.species_run_support_json, raw_to_rep)
        initial_support = collapse_support_map(harvest.initial_run_support, raw_to_rep)
        event_counter = collapse_counter(harvest.species_event_counter, raw_to_rep)
        initial_counts = collapse_counter(harvest.initial_species_counts, raw_to_rep)
        first_seen = collapse_first_seen_order(
            harvest.initial_first_seen_order, raw_to_rep
        )
        reaction_events, reaction_run_support = collapse_reaction_maps(
            harvest.reaction_events, harvest.reaction_run_support, raw_to_rep
        )
        reaction_events_abcd, reaction_run_support_abcd = collapse_reaction_maps(
            harvest.reaction_events_abcd, harvest.reaction_run_support_abcd, raw_to_rep
        )
    else:
        timeline = harvest.species_run_support_timeline
        json_support = harvest.species_run_support_json
        initial_support = defaultdict(set)
        event_counter = harvest.species_event_counter
        initial_counts = harvest.initial_species_counts
        first_seen = harvest.initial_first_seen_order
        reaction_events = harvest.reaction_events
        reaction_run_support = harvest.reaction_run_support
        reaction_events_abcd = harvest.reaction_events_abcd
        reaction_run_support_abcd = harvest.reaction_run_support_abcd

    return CollapsedEnsemble(
        species_run_support=(json_support if support_source == "json" else timeline),
        species_run_support_json=json_support,
        species_run_support_timeline=timeline,
        species_event_counter=event_counter,
        initial_species_counts=initial_counts,
        initial_first_seen_order=first_seen,
        first_frame_support={s: len(runs) for s, runs in initial_support.items()},
        reaction_rows_full=build_reaction_rows(
            reaction_events, reaction_run_support, n_runs
        ),
        reaction_rows_abcd_full=build_reaction_rows(
            reaction_events_abcd, reaction_run_support_abcd, n_runs
        ),
        reaction_run_support=reaction_run_support,
        rep_to_inchi=rep_to_inchi,
        rep_to_members=rep_to_members,
        n_raw_species=n_raw,
        n_collapsed_species=len(rep_to_members) if rep_to_members else n_raw,
        obabel_version=obabel_version,
    )


@dataclass(frozen=True)
class SpeciesTable:
    """Per-species statistics plus the two orderings the reports use."""

    stats: Dict[str, Dict[str, object]]
    id_order: List[str]
    rank_order: List[str]
    discrepancies: List[Tuple[str, int, int, int]]

    @property
    def index(self) -> Dict[str, int]:
        return {s: i + 1 for i, s in enumerate(self.id_order)}


def build_species_table(
    ensemble: CollapsedEnsemble, n_runs: int, support_source: str
) -> SpeciesTable:
    """Tabulate prevalence and provenance for every collapsed species.

    ``is_initial`` marks a species present in the first analysed frame.  Under
    an RMSD bias a reactant can be consumed before the first dump, so the
    displayed support takes ``max(observed, first_frame)`` -- otherwise a
    reactant that reacted immediately in every run would look rare.
    """
    universe = ensemble.universe()
    stats: Dict[str, Dict[str, object]] = {}
    discrepancies: List[Tuple[str, int, int, int]] = []

    for species in universe:
        observed, pct = species_observation_stats(
            ensemble.species_run_support[species], n_runs
        )
        observed_json = len(ensemble.species_run_support_json[species])
        first_frame = int(ensemble.first_frame_support.get(species, 0))
        stats[species] = {
            "smiles": species,
            "collapsed_inchi": ensemble.rep_to_inchi.get(species, ""),
            "representative_smiles": species,
            "member_raw_species_count": len(
                ensemble.rep_to_members.get(species, [species])
            ),
            "observed_run_support": observed,
            "observed_run_support_json": observed_json,
            "support_delta_json": observed - observed_json,
            "support_source": support_source,
            "display_run_support": max(observed, first_frame),
            "first_frame_run_support": first_frame,
            "pct_runs_observed": pct,
            "involvement_events": int(ensemble.species_event_counter[species]),
            "is_initial": species in ensemble.initial_species_counts,
            "initial_frame_count": int(ensemble.initial_species_counts[species]),
            "initial_first_seen_order": ensemble.initial_first_seen_order.get(
                species, NEVER_SEEN
            ),
        }
        if observed - observed_json > 0:
            discrepancies.append((species, observed, observed_json, observed - observed_json))

    discrepancies.sort(key=lambda row: (-row[3], -row[1], row[0]))

    # Reactants first, in the order they were first seen; then everything else
    # by prevalence.  This is what makes the S-numbers stable and readable.
    initial = sorted(
        (s for s in universe if bool(stats[s]["is_initial"])),
        key=lambda s: (stats[s]["initial_first_seen_order"], s),
    )
    remaining = sorted(
        (s for s in universe if not bool(stats[s]["is_initial"])),
        key=lambda s: (
            -int(stats[s]["observed_run_support"]),
            -int(stats[s]["involvement_events"]),
            s,
        ),
    )
    rank_order = sorted(
        universe,
        key=lambda s: (
            -int(stats[s]["display_run_support"]),
            -int(stats[s]["observed_run_support"]),
            -int(stats[s]["involvement_events"]),
            s,
        ),
    )
    return SpeciesTable(
        stats=stats,
        id_order=initial + remaining,
        rank_order=rank_order,
        discrepancies=discrepancies,
    )


def support_audit_warnings(
    table: SpeciesTable, universe_size: int
) -> List[Tuple[str, str]]:
    """Flag species whose ReacNet JSON support undercounts the timeline.

    The JSON species list is truncated by ``--maxspecies``, so it can miss
    species the per-frame timeline does record.  Support is counted from the
    timeline for exactly this reason; these warnings quantify the gap.
    """
    warnings: List[Tuple[str, str]] = []
    if universe_size:
        fraction = len(table.discrepancies) / float(universe_size)
        if fraction > SUPPORT_AUDIT_WARN_FRACTION:
            warnings.append(
                (
                    "support_audit:aggregate",
                    f"JSON species support undercounts timeline support for "
                    f"{len(table.discrepancies)}/{universe_size} species "
                    f"({100.0 * fraction:.1f}%).",
                )
            )
    if table.discrepancies:
        species, observed, observed_json, delta = table.discrepancies[0]
        warnings.append(
            (
                "support_audit:top_species",
                f"Largest JSON-vs-timeline support gap: {species} "
                f"(timeline={observed}, json={observed_json}, delta={delta}).",
            )
        )
    return warnings
