"""Collapsing raw ReacNet SMILES onto fixed-H InChI groups.

Stripping Li leaves fragments whose written bond order and radical placement
depend on where the Li happened to sit in that frame, so the same backbone comes
back from ReacNetGenerator under several different SMILES.  Mapping each raw
SMILES to a fixed-H InChI (Open Babel ``-oinchi -xF``) groups the ones with the
same heavy-atom skeleton *and* the same hydrogen inventory, and one member is
promoted to represent the group.

The consequence is deliberate and worth stating: the collapsed labels describe
the Li-stripped organic framework, not charge-resolved speciation.
"""

from __future__ import annotations

import subprocess
from collections import Counter, defaultdict
from typing import Dict, List, Sequence, Set, Tuple

ReactionKey = Tuple[Tuple[str, ...], Tuple[str, ...]]


class CollapseError(RuntimeError):
    """Open Babel could not produce a usable InChI."""


def get_obabel_version(obabel_bin: str) -> str:
    """Version banner of the Open Babel binary, recorded in the run summary.

    The InChI grouping depends on the Open Babel build, so which one was used
    is part of the provenance of every collapsed species label.
    """
    cmd = [obabel_bin, "-V"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise CollapseError(f"Failed to run {' '.join(cmd)}: {proc.stderr.strip()}")
    for line in proc.stdout.splitlines() + proc.stderr.splitlines():
        line = line.strip()
        if line:
            return line
    raise CollapseError(f"Could not read Open Babel version from {' '.join(cmd)}")


def smiles_to_inchi_fixedh(smiles: str, obabel_bin: str) -> str:
    """Fixed-H InChI for one SMILES.  ``-xF`` is what fixes the H layer."""
    cmd = [obabel_bin, f"-:{smiles}", "-oinchi", "-xF"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        err = proc.stderr.strip() or proc.stdout.strip() or "unknown Open Babel error"
        raise CollapseError(f"InChI conversion failed for species {smiles}: {err}")
    for line in proc.stdout.splitlines():
        line = line.strip()
        if line.startswith("InChI="):
            return line
    raise CollapseError(f"InChI conversion returned no InChI line for species {smiles}")


def build_inchi_collapse_maps(
    raw_species_universe: Set[str],
    raw_support_for_rep: Dict[str, int],
    obabel_bin: str,
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str], Dict[str, List[str]]]:
    """Group every raw species by fixed-H InChI and pick a representative.

    The representative is the group member with the largest run support, with
    the SMILES string itself as a tie-break so the choice is reproducible.

    Returns ``(raw -> representative, raw -> InChI, representative -> InChI,
    representative -> sorted members)``.
    """
    raw_to_inchi: Dict[str, str] = {}
    inchi_to_members: Dict[str, List[str]] = defaultdict(list)
    for raw in sorted(raw_species_universe):
        inchi = smiles_to_inchi_fixedh(raw, obabel_bin)
        raw_to_inchi[raw] = inchi
        inchi_to_members[inchi].append(raw)

    inchi_to_rep: Dict[str, str] = {}
    for inchi, members in inchi_to_members.items():
        inchi_to_rep[inchi] = sorted(
            members, key=lambda s: (-raw_support_for_rep.get(s, 0), s)
        )[0]

    # Defensive: each raw SMILES belongs to exactly one InChI group, so two
    # groups cannot share a representative.  The species tables are keyed on the
    # representative, and a collision would silently merge two chemically
    # distinct entries, so it is worth checking rather than assuming.
    seen_rep: Dict[str, str] = {}
    for inchi, rep in inchi_to_rep.items():
        previous = seen_rep.get(rep)
        if previous is not None and previous != inchi:
            raise CollapseError(
                f"Representative SMILES collision across InChI groups: {rep} "
                f"(InChI {previous} and {inchi})"
            )
        seen_rep[rep] = inchi

    raw_to_rep = {raw: inchi_to_rep[raw_to_inchi[raw]] for raw in raw_species_universe}
    rep_to_inchi = {rep: inchi for inchi, rep in inchi_to_rep.items()}
    rep_to_members = {
        inchi_to_rep[inchi]: sorted(members)
        for inchi, members in inchi_to_members.items()
    }
    return raw_to_rep, raw_to_inchi, rep_to_inchi, rep_to_members


def collapse_support_map(
    raw_map: Dict[str, Set[str]], raw_to_rep: Dict[str, str]
) -> Dict[str, Set[str]]:
    """Union the run sets of every member of a collapsed group.

    Union, not sum: run support is binary per run, so a run that saw two members
    of the same group still counts once.
    """
    out: Dict[str, Set[str]] = defaultdict(set)
    for raw, runs in raw_map.items():
        rep = raw_to_rep.get(raw)
        if rep is None:
            continue
        out[rep].update(runs)
    return out


def collapse_counter(counts: Counter, raw_to_rep: Dict[str, str]) -> Counter:
    """Re-key an event counter onto representatives, summing the members."""
    out: Counter = Counter()
    for raw, n in counts.items():
        rep = raw_to_rep.get(raw)
        if rep is not None:
            out[rep] += int(n)
    return out


def collapse_first_seen_order(
    order: Dict[str, Tuple[int, int]], raw_to_rep: Dict[str, str]
) -> Dict[str, Tuple[int, int]]:
    """Keep the earliest (run index, position) over the members of a group."""
    out: Dict[str, Tuple[int, int]] = {}
    for raw, position in order.items():
        rep = raw_to_rep.get(raw)
        if rep is None:
            continue
        current = out.get(rep)
        if current is None or position < current:
            out[rep] = position
    return out


def collapse_reaction_maps(
    raw_events: Counter,
    raw_support: Dict[ReactionKey, Set[str]],
    raw_to_rep: Dict[str, str],
) -> Tuple[Counter, Dict[ReactionKey, Set[str]]]:
    """Re-key reaction records onto representatives.

    Records whose two sides collapse onto the same multiset are dropped: after
    collapsing they describe no change.
    """
    out_events: Counter = Counter()
    out_support: Dict[ReactionKey, Set[str]] = defaultdict(set)
    for key, total_events in raw_events.items():
        lhs_raw, rhs_raw = key
        lhs = tuple(sorted(raw_to_rep[s] for s in lhs_raw if s in raw_to_rep))
        rhs = tuple(sorted(raw_to_rep[s] for s in rhs_raw if s in raw_to_rep))
        if not lhs or not rhs or lhs == rhs:
            continue
        collapsed: ReactionKey = (lhs, rhs)
        out_events[collapsed] += int(total_events)
        out_support[collapsed].update(raw_support.get(key, set()))
    return out_events, out_support


def collapse_ratio(n_collapsed: int, n_raw: int) -> float:
    return (float(n_collapsed) / float(n_raw)) if n_raw else 1.0


def representative_members(
    rep_to_members: Dict[str, List[str]], rep: str
) -> Sequence[str]:
    return rep_to_members.get(rep, [rep])
