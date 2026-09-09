"""Command-line entry point for the ensemble aggregation."""

from __future__ import annotations

import argparse
import sys
from typing import Optional, Sequence

from .pipeline import EXIT_BAD_INPUT, run_full_pipeline
from .rebuild import run_rebuild_pipeline
from .settings import (
    PATHWAY_MODES,
    SUPPORT_SOURCES,
    Settings,
    SettingsError,
    TARGET_MODES,
    TRAJ_TYPES,
    dedup,
)

EXAMPLE = (
    "ensemble-mtd-aggregate --inputs 'k0.5_a0.6run*.trj' --outdir rcng --nohmm"
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ensemble-mtd-aggregate",
        description=(
            "Aggregate an ensemble of xtb RMSD-metadynamics trajectories into "
            "species prevalence and a consensus reaction network."
        ),
        epilog=f"Example: {EXAMPLE}",
    )

    source = p.add_argument_group("input")
    source.add_argument(
        "--inputs", nargs="+", default=None, help="Trajectory files or globs"
    )
    source.add_argument("--outdir", default=None, help="Where to write the aggregate")
    source.add_argument(
        "--rebuild-from",
        default=None,
        help="Re-render reports from an existing aggregate directory instead of "
        "re-analysing trajectories",
    )

    rn = p.add_argument_group("ReacNetGenerator")
    rn.add_argument("--nproc", type=int, default=36)
    rn.add_argument("--stepinterval", type=int, default=1)
    rn.add_argument("--miso", type=int, default=1, choices=[0, 1, 2])
    rn.add_argument(
        "--nohmm",
        action="store_true",
        help="Disable HMM filtering. Recommended: its defaults suppress the "
        "short-lived intermediates that matter here.",
    )
    rn.add_argument(
        "--maxspecies",
        type=int,
        default=20,
        help="Per-run species cap passed to ReacNet. Affects its JSON species "
        "list only; run support is counted from the per-frame timeline.",
    )
    rn.add_argument("--split", type=int, default=1, help="Per-run ReacNet time split")
    rn.add_argument("--traj-type", default="xyz", choices=list(TRAJ_TYPES))
    rn.add_argument("--reacnet-bin", default="reacnetgenerator")
    rn.add_argument(
        "--obabel-bin", default="obabel", help="Open Babel used for InChI conversion"
    )
    rn.add_argument(
        "--elements",
        nargs="+",
        default=None,
        help="Atom symbols for ReacNet's -a list. Auto-detected from the "
        "Li-removed trajectory when omitted.",
    )

    agg = p.add_argument_group("aggregation")
    agg.add_argument("--top-n", type=int, default=100, help="Rows kept in the reports")
    agg.add_argument(
        "--species-support-source",
        choices=list(SUPPORT_SOURCES),
        default="timeline",
        help="Where run support comes from: the per-frame timeline (default), "
        "ReacNet's JSON species list, or timeline with a JSON audit",
    )

    graph = p.add_argument_group("consensus graph")
    graph.add_argument(
        "--graph-min-edge-runs",
        type=int,
        default=3,
        help="Minimum number of independent runs supporting a shown edge",
    )
    graph.add_argument(
        "--graph-min-edge-events",
        type=float,
        default=1.0,
        help="Minimum aggregated event weight for a shown edge",
    )
    graph.add_argument(
        "--clean",
        action="store_true",
        help="Net opposing edge pairs for display: e_net = |e_uv - e_vu| in the "
        "direction of the larger weight, exact ties kept as dashed both-ways",
    )

    path = p.add_argument_group("path highlighting (optional)")
    path.add_argument(
        "--pathway-mode", choices=list(PATHWAY_MODES), default="none"
    )
    path.add_argument(
        "--target-mode",
        choices=list(TARGET_MODES),
        default="strict-sink",
        help="Which nodes may end a highlighted route",
    )
    path.add_argument("--path-min-steps", type=int, default=1)
    path.add_argument("--top-targets", type=int, default=1)
    path.add_argument(
        "--scan-pathlength",
        default="",
        help="Also highlight the best route of each exact length in a range, e.g. 2-4",
    )

    misc = p.add_argument_group("bookkeeping")
    misc.add_argument(
        "--keep-temp", action="store_true", help="Keep the per-run temp directories"
    )
    misc.add_argument("--tmp-root", default=None, help="Parent for the temp directory")
    misc.add_argument(
        "--save-noli-dir",
        default=None,
        help="Keep the Li-removed trajectories in this directory",
    )
    misc.add_argument("--reacnet-template-html", default=None)
    misc.add_argument(
        "--dry-run", action="store_true", help="List the resolved inputs and stop"
    )
    return p


def settings_from_args(args: argparse.Namespace) -> Settings:
    return Settings(
        inputs=tuple(args.inputs or ()),
        outdir=args.outdir,
        rebuild_from=args.rebuild_from,
        nproc=args.nproc,
        stepinterval=args.stepinterval,
        miso=args.miso,
        nohmm=args.nohmm,
        maxspecies=args.maxspecies,
        split=args.split,
        traj_type=args.traj_type,
        reacnet_bin=args.reacnet_bin,
        obabel_bin=args.obabel_bin,
        elements=dedup(args.elements),
        top_n=args.top_n,
        species_support_source=args.species_support_source,
        graph_min_edge_runs=args.graph_min_edge_runs,
        graph_min_edge_events=args.graph_min_edge_events,
        clean=args.clean,
        pathway_mode=args.pathway_mode,
        target_mode=args.target_mode,
        path_min_steps=args.path_min_steps,
        top_targets=args.top_targets,
        scan_pathlength=args.scan_pathlength,
        keep_temp=args.keep_temp,
        tmp_root=args.tmp_root,
        save_noli_dir=args.save_noli_dir,
        reacnet_template_html=args.reacnet_template_html,
        dry_run=args.dry_run,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        settings = settings_from_args(args).validated()
    except SettingsError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_BAD_INPUT

    if settings.rebuild_from:
        return run_rebuild_pipeline(settings)
    return run_full_pipeline(settings)


if __name__ == "__main__":
    raise SystemExit(main())
