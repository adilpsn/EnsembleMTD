"""Driving ReacNetGenerator on one Li-free trajectory and reading it back.

ReacNetGenerator is run per trajectory, never on a concatenation of them: the
ensemble members are independent, so joining them would invent reactions at the
seams.  HMM filtering is switched off (``--nohmm``) because with its default
settings it removes short-lived intermediates, and under an RMSD bias those are
exactly the radicals worth counting.
"""

from __future__ import annotations

import json
import re
import subprocess
from collections import deque
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .settings import Settings


def tail_lines(path: Path, n: int = 60) -> str:
    """Last ``n`` lines of a file, or "" if it cannot be read."""
    try:
        buf: deque = deque(maxlen=max(int(n), 1))
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                buf.append(line.rstrip("\n"))
        return "\n".join(buf)
    except OSError:
        return ""


def reacnet_command(
    input_xyz: Path, settings: Settings, atom_types: Sequence[str]
) -> List[str]:
    cmd = [
        settings.reacnet_bin,
        "--type", settings.traj_type,
        "-i", input_xyz.name,
        "-a", *list(atom_types),
        "--nopbc",
        "--miso", str(settings.miso),
        "--maxspecies", str(settings.maxspecies),
        "--split", str(settings.split),
        "-n", str(settings.nproc),
        "--stepinterval", str(settings.stepinterval),
    ]
    if settings.nohmm:
        cmd.append("--nohmm")
    return cmd


def run_reacnet(
    input_xyz: Path,
    tmpdir: Path,
    settings: Settings,
    atom_types: Sequence[str],
) -> Tuple[Path, Optional[Path]]:
    """Run ReacNetGenerator in ``tmpdir`` and return its JSON and HTML output.

    The HTML is optional and only used as a template for the aggregate report,
    so a missing one is not an error.  A missing JSON is.
    """
    if not atom_types:
        raise RuntimeError("No atom types available for ReacNet -a list")

    cmd = reacnet_command(input_xyz, settings, atom_types)
    log = tmpdir / f"{input_xyz.name}.reacnet.log"
    with log.open("w", encoding="utf-8") as logf:
        proc = subprocess.run(
            cmd, cwd=str(tmpdir), stdout=logf, stderr=subprocess.STDOUT, text=True
        )
    if proc.returncode != 0:
        raise RuntimeError(
            f"ReacNet failed for {input_xyz.name}\n"
            f"Command: {' '.join(cmd)}\n"
            f"Log: {log}\n"
            f"{tail_lines(log, n=60)}"
        )

    out_json = tmpdir / f"{input_xyz.name}.json"
    out_html = tmpdir / f"{input_xyz.name}.html"
    if not out_json.exists():
        raise RuntimeError(f"Expected output not found: {out_json}")
    return out_json, (out_html if out_html.exists() else None)


def load_payload(out_json: Path) -> Dict:
    return json.loads(out_json.read_text(encoding="utf-8"))


def unwrap_species(raw) -> List[str]:
    """Pull the SMILES list out of ReacNet's JSON species field.

    The field is sometimes a bare list of strings and sometimes a single-element
    list wrapping a list of ``{"s": smiles}`` dicts, depending on version.
    """
    if isinstance(raw, list) and len(raw) == 1 and isinstance(raw[0], list):
        raw = raw[0]
    if not isinstance(raw, list):
        return []
    out = []
    for x in raw:
        if isinstance(x, dict) and "s" in x:
            out.append(str(x["s"]))
        elif isinstance(x, str):
            out.append(x)
    return out


def unwrap_reactions(raw) -> List[dict]:
    """Same single-element-list unwrapping for the reaction records."""
    if isinstance(raw, list) and len(raw) == 1 and isinstance(raw[0], list):
        raw = raw[0]
    if not isinstance(raw, list):
        return []
    return [x for x in raw if isinstance(x, dict)]


TIMESTEP_RE = re.compile(r"^Timestep\s+\S+:\s*(.*)$")


def parse_species_timeline(species_path: Path) -> Tuple[List[Tuple[str, int]], Set[str]]:
    """Read the per-frame ``.species`` file.

    Returns the species of the first frame with their counts, and the set of
    every species seen in any frame.  The set, not the JSON species list, is
    what run support is counted from: the JSON list is truncated by
    ``--maxspecies`` and would silently undercount.
    """
    if not species_path.exists():
        raise RuntimeError(f"Expected species file not found: {species_path}")

    initial_species: List[Tuple[str, int]] = []
    species_present: Set[str] = set()
    saw_timestep = False
    captured_initial = False

    with species_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            match = TIMESTEP_RE.match(line)
            if not match:
                raise RuntimeError(
                    f"Unexpected species line format in {species_path}: {line!r}"
                )
            saw_timestep = True
            payload = match.group(1).strip()
            frame_species: List[Tuple[str, int]] = []
            if payload:
                toks = payload.split()
                if len(toks) % 2 != 0:
                    raise RuntimeError(
                        f"Malformed species line in {species_path}: {line!r}"
                    )
                for i in range(0, len(toks), 2):
                    species = toks[i]
                    try:
                        count = int(toks[i + 1])
                    except ValueError as exc:
                        raise RuntimeError(
                            f"Malformed species count in {species_path}: {line!r}"
                        ) from exc
                    if count > 0:
                        frame_species.append((species, count))
                        species_present.add(species)
            if not captured_initial:
                initial_species = frame_species
                captured_initial = True

    if not saw_timestep:
        raise RuntimeError(f"Species file is empty: {species_path}")
    return initial_species, species_present
