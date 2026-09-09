"""Frozen configuration for one aggregation run.

Everything the pipeline needs is fixed before the first trajectory is read, so
the settings travel as one immutable value rather than as an argparse Namespace
threaded through a dozen call sites.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional, Sequence, Tuple

# Report defaults.  These are display limits only and do not affect any of the
# numbers written to the TSV files.
DEFAULT_SPECIES_SHOW_NUM = 30
DEFAULT_REACTIONS_SHOW_NUM = 20
DEFAULT_GRAPH_MAX_NODES = 70
DEFAULT_GRAPH_MAX_EDGES = 160
DEFAULT_GRAPH_INCLUDE_ISOLATES = False

TRAJ_TYPES = ("xyz", "extxyz", "dump", "bond")
PATHWAY_MODES = ("support-first", "probabilistic", "none")
TARGET_MODES = ("strict-sink", "pseudo-sink", "all")
SUPPORT_SOURCES = ("timeline", "json", "audit")


class SettingsError(ValueError):
    """Raised for a combination of options that cannot be honoured."""


@dataclass(frozen=True)
class Settings:
    inputs: Tuple[str, ...] = ()
    outdir: Optional[str] = None
    rebuild_from: Optional[str] = None

    # Passed through to ReacNetGenerator.
    nproc: int = 36
    stepinterval: int = 1
    miso: int = 1
    nohmm: bool = False
    maxspecies: int = 20
    split: int = 1
    traj_type: str = "xyz"
    reacnet_bin: str = "reacnetgenerator"
    obabel_bin: str = "obabel"
    elements: Optional[Tuple[str, ...]] = None

    # Aggregation and reporting.
    top_n: int = 100
    species_support_source: str = "timeline"

    # Consensus-graph thresholds.  Three runs is the value used in the paper:
    # an edge seen in one or two trajectories out of twenty is not a recurrent
    # channel.
    graph_min_edge_runs: int = 3
    graph_min_edge_events: float = 1.0
    clean: bool = False

    # Optional path highlighting on top of the consensus graph.
    pathway_mode: str = "none"
    target_mode: str = "strict-sink"
    path_min_steps: int = 1
    top_targets: int = 1
    scan_pathlength: str = ""
    scan_pathlength_values: Tuple[int, ...] = ()

    # Bookkeeping.
    keep_temp: bool = False
    tmp_root: Optional[str] = None
    save_noli_dir: Optional[str] = None
    reacnet_template_html: Optional[str] = None
    dry_run: bool = False

    def validated(self) -> "Settings":
        """Return a copy with the path-length scan parsed, or raise."""
        if self.pathway_mode == "none" and self.scan_pathlength.strip():
            raise SettingsError(
                "--pathway-mode none cannot be used with --scan-pathlength"
            )
        if self.inputs and self.rebuild_from:
            raise SettingsError("Use either --inputs or --rebuild-from, not both")
        if not self.inputs and not self.rebuild_from:
            raise SettingsError("Provide either --inputs or --rebuild-from")
        if self.inputs and not self.outdir:
            raise SettingsError("--outdir is required when using --inputs")
        if self.traj_type not in TRAJ_TYPES:
            raise SettingsError(f"Unknown --traj-type: {self.traj_type}")
        if self.pathway_mode not in PATHWAY_MODES:
            raise SettingsError(f"Unknown --pathway-mode: {self.pathway_mode}")
        if self.target_mode not in TARGET_MODES:
            raise SettingsError(f"Unknown --target-mode: {self.target_mode}")
        if self.species_support_source not in SUPPORT_SOURCES:
            raise SettingsError(
                f"Unknown --species-support-source: {self.species_support_source}"
            )
        if self.graph_min_edge_runs < 1:
            raise SettingsError("--graph-min-edge-runs must be at least 1")
        if self.nproc < 1:
            raise SettingsError("--nproc must be at least 1")

        return replace(
            self,
            scan_pathlength_values=parse_pathlength_range(self.scan_pathlength),
        )


def parse_pathlength_range(text: str) -> Tuple[int, ...]:
    """Parse "2-4" into (2, 3, 4).  An empty string means no scan."""
    if not text.strip():
        return ()
    try:
        start_text, end_text = text.split("-")
        start, end = int(start_text), int(end_text)
    except ValueError as exc:
        raise SettingsError(
            f"Invalid --scan-pathlength {text!r}: expected a range like 2-4"
        ) from exc
    if start < 1 or end < start:
        raise SettingsError(
            "--scan-pathlength must look like 2-4 with 1 <= start <= end"
        )
    return tuple(range(start, end + 1))


def dedup(values: Optional[Sequence[str]]) -> Optional[Tuple[str, ...]]:
    """Drop repeats while keeping order (used for the element list)."""
    if values is None:
        return None
    return tuple(dict.fromkeys(values))
