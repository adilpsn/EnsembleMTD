"""Writing the aggregate artefacts to the output directory.

Filenames are kept exactly as earlier versions of this pipeline wrote them, so
existing analysis scripts and the figures in the paper keep working.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .aggregate import CollapsedEnsemble, SpeciesTable
from .prevalence import ReactionKey, key_to_str
from .report import write_tsv
from .settings import Settings

SUMMARY_JSON = "aggregate_reacnet_summary.json"
REACTIONS_TSV = "aggregate_reacnet_reactions.tsv"
SPECIES_TSV = "aggregate_reacnet_species.tsv"
MEMBERSHIP_TSV = "aggregate_reacnet_species_run_membership.tsv"
ALIAS_MAP_TSV = "aggregate_reacnet_inchi_alias_map.tsv"
PATHLENGTH_SCAN_TSV = "aggregate_reacnet_pathlength_scan.tsv"
FAILURES_TSV = "aggregate_reacnet_failures.tsv"
WARNINGS_TSV = "aggregate_reacnet_warnings.tsv"
REPORT_HTML = "aggregate_reacnet_report.html"
RNGDATA_JSON = "aggregate_reacnet_rngdata.json"
REACNETSTYLE_HTML = "aggregate_reacnet_reacnetstyle.html"
NETWORK_DOT = "aggregate_reacnet_network.dot"
NETWORK_SVG = "aggregate_reacnet_network.svg"

SPECIES_COLUMNS = (
    "species_id",
    "species",
    "collapsed_inchi",
    "representative_smiles",
    "member_raw_species_count",
    "display_run_support",
    "observed_run_support",
    "observed_run_support_json",
    "support_delta_json",
    "support_source",
    "pct_runs_observed",
    "is_initial",
    "first_frame_run_support",
    "initial_frame_count",
    "involvement_events",
)


def simple_reaction_rows(
    reaction_rows_full: Sequence[Tuple[ReactionKey, int, float, int]]
) -> List[Tuple[str, int, float, int]]:
    return [
        (key_to_str(key), support, pct, events)
        for key, support, pct, events in reaction_rows_full
    ]


def write_reactions(outdir: Path, rows: Sequence[Tuple[str, int, float, int]]) -> None:
    write_tsv(
        outdir / REACTIONS_TSV,
        ["reaction", "run_support", "pct_runs", "total_events"],
        rows,
    )


def write_species(outdir: Path, table: SpeciesTable) -> None:
    index = table.index
    write_tsv(
        outdir / SPECIES_TSV,
        SPECIES_COLUMNS,
        [
            (
                f"S{index[s]}",
                s,
                str(table.stats[s]["collapsed_inchi"]),
                str(table.stats[s]["representative_smiles"]),
                int(table.stats[s]["member_raw_species_count"]),
                int(table.stats[s]["display_run_support"]),
                int(table.stats[s]["observed_run_support"]),
                int(table.stats[s]["observed_run_support_json"]),
                int(table.stats[s]["support_delta_json"]),
                str(table.stats[s]["support_source"]),
                f"{float(table.stats[s]['pct_runs_observed']):.3f}",
                int(bool(table.stats[s]["is_initial"])),
                int(table.stats[s]["first_frame_run_support"]),
                int(table.stats[s]["initial_frame_count"]),
                int(table.stats[s]["involvement_events"]),
            )
            for s in table.id_order
        ],
    )


def write_run_membership(
    outdir: Path, table: SpeciesTable, ensemble: CollapsedEnsemble
) -> None:
    """Per-run presence matrix, one row per species with the run names.

    The other tables store only support counts, which is not enough to redo the
    Pi_i-versus-N convergence analysis or any per-run resampling.  The run sets
    here are the same collapsed, source-selected sets ``pct_runs_observed`` is
    computed from, so the counts reproduce exactly.
    """
    index = table.index
    write_tsv(
        outdir / MEMBERSHIP_TSV,
        ["species_id", "species", "observed_run_support", "pct_runs_observed", "runs"],
        [
            (
                f"S{index[s]}",
                s,
                int(table.stats[s]["observed_run_support"]),
                f"{float(table.stats[s]['pct_runs_observed']):.3f}",
                ";".join(sorted(ensemble.species_run_support[s])),
            )
            for s in table.id_order
        ],
    )


def write_alias_map(
    outdir: Path, table: SpeciesTable, ensemble: CollapsedEnsemble
) -> None:
    """Every raw ReacNet SMILES and the collapsed group it was folded into."""
    index = table.index
    rows: List[Tuple[str, str, str, str, int]] = []
    for species in table.id_order:
        for raw in ensemble.rep_to_members.get(species, [species]):
            rows.append(
                (
                    raw,
                    str(ensemble.rep_to_inchi.get(species, "")),
                    f"S{index[species]}",
                    str(table.stats[species]["representative_smiles"]),
                    int(table.stats[species]["observed_run_support"]),
                )
            )
    rows.sort(key=lambda row: (row[2], row[0]))
    write_tsv(
        outdir / ALIAS_MAP_TSV,
        [
            "raw_species_smiles",
            "collapsed_inchi",
            "species_id",
            "representative_smiles",
            "observed_run_support",
        ],
        rows,
    )


def write_pathlength_scan(outdir: Path, scan_rows: Sequence[dict]) -> None:
    if not scan_rows:
        return
    write_tsv(
        outdir / PATHLENGTH_SCAN_TSV,
        ["path_length", "target", "path", "score", "mode"],
        [
            (
                int(row.get("path_length", 0)),
                str(row.get("target_id", row.get("target", ""))),
                " -> ".join(row.get("path", [])),
                str(row.get("score", "")),
                str(row.get("mode", "")),
            )
            for row in scan_rows
        ],
    )


def build_summary(
    settings: Settings,
    ensemble: CollapsedEnsemble,
    table: SpeciesTable,
    graph_meta: Dict[str, object],
    input_files: Sequence[str],
    run_stats: Sequence[dict],
    excluded_runs: Sequence[str],
    failures: Sequence[Tuple[str, str]],
    warnings: Sequence[Tuple[str, str]],
    collapse_failures: Sequence[str],
    n_runs: int,
    reaction_rows_simple: Sequence[Tuple[str, int, float, int]],
) -> dict:
    index = table.index
    return {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "reacnet_settings": {
            "traj_type": settings.traj_type,
            "miso": int(settings.miso),
            "nohmm": bool(settings.nohmm),
            "maxspecies": int(settings.maxspecies),
            "split": int(settings.split),
            "stepinterval": int(settings.stepinterval),
            "nproc": int(settings.nproc),
        },
        "species_support_source": settings.species_support_source,
        "initial_support_mode": "max(observed,first_frame)",
        "collapse_mode": "inchi_fixedh_full_analysis",
        "n_raw_species": ensemble.n_raw_species,
        "n_collapsed_species": ensemble.n_collapsed_species,
        "collapse_ratio": ensemble.collapse_ratio,
        "obabel_version": ensemble.obabel_version,
        "collapse_failures": list(collapse_failures),
        "n_input_runs": len(input_files),
        "n_success_runs": n_runs,
        "n_excluded_runs": len(excluded_runs),
        "excluded_runs": list(excluded_runs),
        "n_failed_runs": len(failures),
        "input_files": list(input_files),
        "run_stats": list(run_stats),
        "graph": graph_meta,
        "failures": [{"run": r, "error": e} for r, e in failures],
        "warnings": [{"run": r, "warning": w} for r, w in warnings],
        "n_species_support_discrepancies": len(table.discrepancies),
        "top_species_support_discrepancies": [
            {
                "species": s,
                "species_id": f"S{index[s]}",
                "observed_run_support": observed,
                "observed_run_support_json": observed_json,
                "support_delta_json": delta,
            }
            for s, observed, observed_json, delta in table.discrepancies[
                : settings.top_n
            ]
        ],
        "top_reactions": [
            {
                "reaction": r,
                "run_support": support,
                "pct_runs": round(pct, 3),
                "total_events": events,
            }
            for r, support, pct, events in reaction_rows_simple[: settings.top_n]
        ],
        "top_species": [
            {
                "species": s,
                "species_id": f"S{index[s]}",
                "collapsed_inchi": str(table.stats[s]["collapsed_inchi"]),
                "representative_smiles": str(table.stats[s]["representative_smiles"]),
                "member_raw_species_count": int(
                    table.stats[s]["member_raw_species_count"]
                ),
                "display_run_support": int(table.stats[s]["display_run_support"]),
                "observed_run_support": int(table.stats[s]["observed_run_support"]),
                "observed_run_support_json": int(
                    table.stats[s]["observed_run_support_json"]
                ),
                "support_delta_json": int(table.stats[s]["support_delta_json"]),
                "pct_runs_observed": round(
                    float(table.stats[s]["pct_runs_observed"]), 3
                ),
                "is_initial": bool(table.stats[s]["is_initial"]),
                "first_frame_run_support": int(
                    table.stats[s]["first_frame_run_support"]
                ),
            }
            for s in table.rank_order[: settings.top_n]
        ],
    }


def write_summary(outdir: Path, summary: dict) -> None:
    (outdir / SUMMARY_JSON).write_text(json.dumps(summary, indent=2), encoding="utf-8")


def top_species_rows(
    table: SpeciesTable, top_n: int
) -> List[Tuple[str, int, float]]:
    return [
        (
            s,
            int(table.stats[s]["display_run_support"]),
            float(table.stats[s]["pct_runs_observed"]),
        )
        for s in table.rank_order[:top_n]
    ]


def report_written_paths(outdir: Path, has_failures: bool, has_warnings: bool) -> List[Path]:
    """The files a completed run should announce, in a stable order."""
    always = [
        SUMMARY_JSON,
        REACTIONS_TSV,
        SPECIES_TSV,
        ALIAS_MAP_TSV,
        REPORT_HTML,
        RNGDATA_JSON,
        NETWORK_DOT,
    ]
    optional = [NETWORK_SVG, REACNETSTYLE_HTML, PATHLENGTH_SCAN_TSV]
    paths = [outdir / name for name in always]
    paths += [outdir / name for name in optional if (outdir / name).exists()]
    if has_failures:
        paths.append(outdir / FAILURES_TSV)
    if has_warnings:
        paths.append(outdir / WARNINGS_TSV)
    return paths
