"""Re-rendering the reports from an existing aggregate directory.

Useful when only the presentation has to change -- a Graphviz version, a
template tweak -- because re-running ReacNetGenerator over twenty trajectories
costs minutes per ensemble and would produce the same numbers.  Nothing here
recomputes prevalence: the TSV files are read back as written.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from . import artifacts
from .report import (
    inject_rngdata_into_template,
    load_top_rows_from_tsv,
    patch_reacnet_templates_for_support_labels,
    write_simple_html,
    write_tsv,
)
from .settings import Settings

EXIT_OK = 0
EXIT_BAD_INPUT = 2
EXIT_PARTIAL = 3

REACTION_FIELDS = ("reaction", "run_support", "pct_runs", "total_events")
SPECIES_FIELDS = ("species", "display_run_support", "pct_runs_observed")


def _as_int(text: str) -> int:
    return int(float(text))


def _load_summary(
    srcdir: Path, warnings: List[Tuple[str, str]]
) -> Dict[str, object]:
    path = srcdir / artifacts.SUMMARY_JSON
    if not path.exists():
        warnings.append(
            (
                "rebuild_summary",
                f"Missing {artifacts.SUMMARY_JSON}; report fallback will be limited.",
            )
        )
        return {}
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        warnings.append(("rebuild_summary", f"Failed to parse summary JSON: {exc}"))
        return {}
    return loaded if isinstance(loaded, dict) else {}


def _rerender_svg(
    srcdir: Path, outdir: Path, failures: List[Tuple[str, str]]
) -> None:
    dot_src = srcdir / artifacts.NETWORK_DOT
    if not dot_src.exists():
        failures.append(("rebuild_svg", f"Missing required file: {dot_src}"))
        return
    dot_dst = outdir / artifacts.NETWORK_DOT
    svg_dst = outdir / artifacts.NETWORK_SVG
    try:
        if dot_src.resolve() != dot_dst.resolve():
            shutil.copy2(dot_src, dot_dst)
        proc = subprocess.run(
            ["dot", "-Tsvg", str(dot_dst), "-o", str(svg_dst)],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            message = (proc.stderr or proc.stdout).strip()
            failures.append(("rebuild_svg", f"Graphviz failed: {message}"))
    except OSError as exc:
        failures.append(("rebuild_svg", str(exc)))


def _reactions_from(
    srcdir: Path, summary: Dict[str, object], top_n: int
) -> List[Tuple[str, int, float, int]]:
    rows = [
        (str(r[0]), int(r[1]), float(r[2]), int(r[3]))
        for r in load_top_rows_from_tsv(
            tsv_path=srcdir / artifacts.REACTIONS_TSV,
            required_fields=REACTION_FIELDS,
            converters={
                "reaction": str,
                "run_support": _as_int,
                "pct_runs": float,
                "total_events": _as_int,
            },
            top_n=top_n,
        )
    ]
    if rows:
        return rows
    for row in summary.get("top_reactions", [])[:top_n]:
        try:
            rows.append(
                (
                    str(row.get("reaction", "")),
                    int(row.get("run_support", 0)),
                    float(row.get("pct_runs", 0.0)),
                    int(row.get("total_events", 0)),
                )
            )
        except (TypeError, ValueError):
            continue
    return rows


def _species_from(
    srcdir: Path, summary: Dict[str, object], top_n: int
) -> List[Tuple[str, int, float]]:
    rows = [
        (str(r[0]), int(r[1]), float(r[2]))
        for r in load_top_rows_from_tsv(
            tsv_path=srcdir / artifacts.SPECIES_TSV,
            required_fields=SPECIES_FIELDS,
            converters={
                "species": str,
                "display_run_support": _as_int,
                "pct_runs_observed": float,
            },
            top_n=top_n,
        )
    ]
    if rows:
        return rows
    for row in summary.get("top_species", [])[:top_n]:
        try:
            rows.append(
                (
                    str(row.get("species", "")),
                    int(row.get("display_run_support", 0)),
                    float(row.get("pct_runs_observed", 0.0)),
                )
            )
        except (TypeError, ValueError):
            continue
    return rows


def _rerender_reacnetstyle(
    srcdir: Path,
    outdir: Path,
    settings: Settings,
    failures: List[Tuple[str, str]],
    warnings: List[Tuple[str, str]],
) -> None:
    rngdata_path = srcdir / artifacts.RNGDATA_JSON
    if not rngdata_path.exists():
        failures.append(
            ("rebuild_reacnetstyle", f"Missing required file: {rngdata_path}")
        )
        return
    try:
        rngdata = json.loads(rngdata_path.read_text(encoding="utf-8"))
        template = _rebuild_template(srcdir, settings)
        if not template:
            failures.append(
                (
                    "rebuild_reacnetstyle",
                    "No template available. Provide --reacnet-template-html or "
                    f"source {artifacts.REACNETSTYLE_HTML}",
                )
            )
            return
        patched, patch_warnings = patch_reacnet_templates_for_support_labels(template)
        for warning in patch_warnings:
            warnings.append(("rebuild_reacnetstyle_template_patch", warning))
        (outdir / artifacts.REACNETSTYLE_HTML).write_text(
            inject_rngdata_into_template(patched, rngdata), encoding="utf-8"
        )
        target = outdir / artifacts.RNGDATA_JSON
        if target.resolve() != rngdata_path.resolve():
            target.write_text(json.dumps(rngdata), encoding="utf-8")
    except (OSError, ValueError, RuntimeError) as exc:
        failures.append(("rebuild_reacnetstyle", str(exc)))


def _rebuild_template(srcdir: Path, settings: Settings) -> Optional[str]:
    if settings.reacnet_template_html:
        return Path(settings.reacnet_template_html).read_text(
            encoding="utf-8", errors="ignore"
        )
    candidate = srcdir / artifacts.REACNETSTYLE_HTML
    if candidate.exists():
        return candidate.read_text(encoding="utf-8", errors="ignore")
    return None


def run_rebuild_pipeline(settings: Settings) -> int:
    srcdir = Path(str(settings.rebuild_from)).resolve()
    if not srcdir.is_dir():
        print(f"Invalid --rebuild-from directory: {srcdir}")
        return EXIT_BAD_INPUT

    outdir = Path(settings.outdir).resolve() if settings.outdir else srcdir
    outdir.mkdir(parents=True, exist_ok=True)

    if settings.dry_run:
        print(f"Rebuild dry-run from: {srcdir}")
        print(f"Output dir: {outdir}")
        return EXIT_OK

    failures: List[Tuple[str, str]] = []
    warnings: List[Tuple[str, str]] = []

    summary = _load_summary(srcdir, warnings)
    _rerender_svg(srcdir, outdir, failures)

    top_n = settings.top_n
    try:
        write_simple_html(
            outdir / artifacts.REPORT_HTML,
            n_runs=int(summary.get("n_success_runs", 0)),
            input_files=list(summary.get("input_files", []) or []),
            top_reactions=_reactions_from(srcdir, summary, top_n),
            top_species=_species_from(srcdir, summary, top_n),
        )
    except (OSError, TypeError, ValueError) as exc:
        failures.append(("rebuild_report", str(exc)))

    _rerender_reacnetstyle(srcdir, outdir, settings, failures, warnings)

    if failures:
        write_tsv(outdir / artifacts.FAILURES_TSV, ["run", "error"], failures)
    if warnings:
        write_tsv(outdir / artifacts.WARNINGS_TSV, ["run", "warning"], warnings)

    print("Wrote:")
    for name in (
        artifacts.REPORT_HTML,
        artifacts.NETWORK_DOT,
        artifacts.NETWORK_SVG,
        artifacts.REACNETSTYLE_HTML,
        artifacts.RNGDATA_JSON,
    ):
        if (outdir / name).exists():
            print(outdir / name)
    if failures:
        print(outdir / artifacts.FAILURES_TSV)
    if warnings:
        print(outdir / artifacts.WARNINGS_TSV)

    return EXIT_OK if not failures else EXIT_PARTIAL
