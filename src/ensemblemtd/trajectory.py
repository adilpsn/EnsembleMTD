"""Reading the raw xtb trajectories and stripping lithium out of them.

Li is removed frame by frame before any connectivity is inferred.  A Li+ that
hops between two carbonyl oxygens changes the Open Babel bond graph without any
covalent chemistry happening, and ReacNetGenerator would faithfully report that
as a reaction.  Removing Li first costs the coordination information but leaves
the organic framework, which is what the consensus network describes.
"""

from __future__ import annotations

import glob
import os
import re
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

# How many frames to look at when guessing the element list for ReacNet's -a
# option.  A handful is enough: the composition cannot change during an MD run.
AUTO_ELEMENT_SCAN_FRAMES = 20

# setup_reactor.sh writes "energy: <E>" into the xyz comment line.  xtb emits
# exactly 0.0 there when an SCF has collapsed, which is how a diverged run is
# recognised downstream.
ENERGY_COMMENT_RE = re.compile(
    r"\benergy\s*:\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
)
ENERGY_ZERO_TOL = 1e-12

LITHIUM = "Li"


def natural_key(s: str):
    """Sort key that orders run10 after run9 rather than after run1."""
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", s)]


def expand_inputs(patterns: Sequence[str]) -> List[str]:
    """Resolve globs and plain paths to a de-duplicated, run-ordered list."""
    files: List[str] = []
    for pat in patterns:
        hits = glob.glob(pat)
        if hits:
            files.extend(hits)
        elif os.path.exists(pat):
            files.append(pat)
    return sorted({os.path.abspath(x) for x in files}, key=natural_key)


def parse_comment_energy(comment_line: str) -> Optional[float]:
    match = ENERGY_COMMENT_RE.search(comment_line)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def strip_li_xyz(src: Path, dst: Path) -> Tuple[int, int, int, List[str], Optional[int]]:
    """Write ``src`` to ``dst`` with every Li atom dropped.

    Returns the frame count, the atom counts before and after stripping, the
    non-Li elements in order of first appearance, and the 1-based index of the
    first frame whose comment reports zero energy (``None`` if there is none).
    """
    frames = 0
    atoms_in = 0
    atoms_out = 0
    seen: set = set()
    ordered_elements: List[str] = []
    first_zero_energy_frame: Optional[int] = None

    with src.open("r", encoding="utf-8", errors="replace") as fin, \
            dst.open("w", encoding="utf-8") as fout:
        while True:
            first = fin.readline()
            if not first:
                break
            first = first.strip()
            if not first:
                continue
            try:
                nat = int(first)
            except ValueError as exc:
                raise RuntimeError(
                    f"{src}: invalid XYZ frame atom-count line: {first!r}"
                ) from exc

            comment = fin.readline()
            if not comment:
                raise RuntimeError(f"{src}: truncated XYZ, missing comment line")
            energy = parse_comment_energy(comment)
            if (
                energy is not None
                and abs(energy) <= ENERGY_ZERO_TOL
                and first_zero_energy_frame is None
            ):
                first_zero_energy_frame = frames + 1

            atoms = [fin.readline() for _ in range(nat)]
            if any(x == "" for x in atoms):
                raise RuntimeError(f"{src}: truncated XYZ atom block")

            kept: List[str] = []
            for line in atoms:
                toks = line.split()
                if not toks:
                    continue
                elem = toks[0]
                if elem == LITHIUM:
                    continue
                if elem not in seen:
                    seen.add(elem)
                    ordered_elements.append(elem)
                kept.append(line)

            fout.write(f"{len(kept)}\n")
            fout.write(comment.rstrip("\n") + " | Li-removed\n")
            for line in kept:
                fout.write(line)

            frames += 1
            atoms_in += nat
            atoms_out += len(kept)

    return frames, atoms_in, atoms_out, ordered_elements, first_zero_energy_frame


def detect_elements_xyz(path: Path, max_frames: Optional[int] = None) -> List[str]:
    """Non-Li elements in ``path``, in order of first appearance.

    Only used as a fallback: the element list normally comes straight out of
    :func:`strip_li_xyz`, which has already read every frame.
    """
    seen: set = set()
    ordered: List[str] = []
    scanned = 0
    frames_to_scan = (
        AUTO_ELEMENT_SCAN_FRAMES if max_frames is None else max(1, int(max_frames))
    )

    with path.open("r", encoding="utf-8", errors="replace") as fin:
        while scanned < frames_to_scan:
            first = fin.readline()
            if not first:
                break
            first = first.strip()
            if not first:
                continue
            try:
                nat = int(first)
            except ValueError as exc:
                raise RuntimeError(
                    f"{path}: invalid XYZ frame atom-count line while detecting "
                    f"elements: {first!r}"
                ) from exc

            if not fin.readline():
                break
            for _ in range(nat):
                line = fin.readline()
                if not line:
                    break
                toks = line.split()
                if not toks:
                    continue
                elem = toks[0]
                if elem == LITHIUM:
                    continue
                if elem not in seen:
                    seen.add(elem)
                    ordered.append(elem)
            scanned += 1

    if not ordered:
        raise RuntimeError(
            f"Could not auto-detect non-Li elements from {path}. "
            "Use --elements explicitly."
        )
    return ordered
