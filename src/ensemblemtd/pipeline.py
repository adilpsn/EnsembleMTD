"""The aggregation pipeline: one pass over the ensemble, then the tables.

Per trajectory: strip Li, run ReacNetGenerator, read the species timeline and
the reaction records, fold them into the ensemble tallies.  A run whose xyz
comment reports zero energy is excluded rather than analysed -- that is xtb
telling us the SCF collapsed, and its connectivity is meaningless.  A run that
raises for any other reason is recorded as a failure and the rest continue: an
ensemble is worth reporting with 19 of 20 members.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from . import artifacts
from .aggregate import RunHarvest, build_species_table, collapse_harvest, support_audit_warnings
from .collapse import CollapseError
from .consensus import build_consensus_network
from .network import build_linkreac
from .prevalence import add_reactions
from .reacnet import (
    load_payload,
    parse_species_timeline,
    run_reacnet,
    unwrap_reactions,
    unwrap_species,
)
from .report import (
    build_reacnet_rngdata,
    inject_rngdata_into_template,
    link_or_copy,
    patch_reacnet_templates_for_support_labels,
    write_simple_html,
    write_tsv,
)
from .settings import (
    DEFAULT_REACTIONS_SHOW_NUM,
    DEFAULT_SPECIES_SHOW_NUM,
    Settings,
)
from .trajectory import detect_elements_xyz, expand_inputs, strip_li_xyz

EXIT_OK = 0
EXIT_BAD_INPUT = 2
EXIT_NO_SUCCESSFUL_RUNS = 3

MISSING_SPECIES_PREVIEW = 8

# Where to point someone whose environment is not set up yet.
TOOL_HINTS = {
    "reacnetgenerator": "pip install reacnetgenerator  (or: pip install 'ensemblemtd[reacnet]')",
    "obabel": "conda install -c conda-forge openbabel  (or: apt install openbabel)",
}


def missing_tools(settings: Settings) -> List[Tuple[str, str]]:
    """External programs the run needs that are not on PATH.

    Checked before any trajectory is read.  Without this, a missing binary
    surfaces as one opaque OSError per run and an output directory full of
    empty tables, with the real cause buried in the failures file.
    """
    missing = []
    for binary in (settings.reacnet_bin, settings.obabel_bin):
        if shutil.which(binary) is None:
            hint = TOOL_HINTS.get(os.path.basename(binary), "")
            missing.append((binary, hint))
    return missing


def _excluded_run_record(
    run_name: str, frames: int, atoms_in: int, atoms_out: int, frame: int
) -> dict:
    return {
        "run": run_name,
        "frames": frames,
        "atoms_in": atoms_in,
        "atoms_out": atoms_out,
        "atom_types_used": [],
        "species_unique": 0,
        "species_unique_timeline": 0,
        "species_unique_json": 0,
        "species_missing_in_json": 0,
        "reactions_unique": 0,
        "reactionsabcd_unique": 0,
        "initial_species_unique": 0,
        "excluded": 1,
        "excluded_reason": "zero_energy_frame_detected",
        "excluded_frame": int(frame),
    }


def _support_audit_warning(
    run_name: str, missing_in_json: Sequence[str]
) -> List[Tuple[str, str]]:
    """Note species the per-frame timeline saw but ReacNet's JSON list did not.

    The JSON list is capped by ``--maxspecies``; the timeline is not.  Support
    is counted from the timeline for that reason, and this quantifies the gap.
    """
    if not missing_in_json:
        return []
    preview = ", ".join(missing_in_json[:MISSING_SPECIES_PREVIEW])
    if len(missing_in_json) > MISSING_SPECIES_PREVIEW:
        preview += ", ..."
    return [
        (
            f"support_audit:{run_name}",
            f"Timeline species absent from JSON species list "
            f"({len(missing_in_json)}): {preview}",
        )
    ]


def _harvest_run(
    index: int,
    in_file: str,
    run_tmp: Path,
    save_noli_dir: Optional[Path],
    settings: Settings,
    harvest: RunHarvest,
    explicit_atom_types: Optional[Sequence[str]],
) -> Tuple[dict, List[Tuple[str, str]], Optional[str]]:
    """Analyse one trajectory and fold it into ``harvest``.

    Returns its run-stats record, any warnings it raised, and the ReacNet HTML
    shell if this run produced one.
    """
    run_name = os.path.basename(in_file)
    src = Path(in_file)
    noli_file = (
        run_tmp / f"{src.stem}.noli.xyz"
        if save_noli_dir is None
        else save_noli_dir / f"{src.stem}.noli.xyz"
    )

    frames, atoms_in, atoms_out, stripped_types, zero_frame = strip_li_xyz(src, noli_file)
    if zero_frame is not None:
        record = _excluded_run_record(run_name, frames, atoms_in, atoms_out, zero_frame)
        warning = (
            f"excluded:{run_name}",
            f"Run excluded from aggregation: frame {zero_frame} has energy: 0.",
        )
        return record, [warning], None

    local_input = run_tmp / noli_file.name
    if noli_file.resolve() != local_input.resolve():
        link_or_copy(noli_file, local_input)

    if explicit_atom_types is not None:
        atom_types = list(explicit_atom_types)
    else:
        atom_types = stripped_types
        if not atom_types:
            atom_types = detect_elements_xyz(local_input)

    out_json, out_html = run_reacnet(local_input, run_tmp, settings, atom_types)
    payload = load_payload(out_json)
    initial_species, timeline_species = parse_species_timeline(
        run_tmp / f"{local_input.name}.species"
    )

    species = unwrap_species(payload.get("species", []))
    present = (
        set(species)
        if settings.species_support_source == "json"
        else set(timeline_species)
    )
    missing_in_json = sorted(present - set(species))

    reactions = unwrap_reactions(payload.get("reactions", []))
    reactions_abcd = unwrap_reactions(payload.get("reactionsabcd", []))

    for order, (species_name, count) in enumerate(initial_species):
        harvest.initial_species_counts[species_name] += count
        harvest.initial_run_support[species_name].add(run_name)
        harvest.initial_first_seen_order.setdefault(species_name, (index, order))
    for name in timeline_species:
        harvest.species_run_support_timeline[name].add(run_name)
    for name in species:
        harvest.species_run_support_json[name].add(run_name)

    add_reactions(
        reactions=reactions,
        run_name=run_name,
        reaction_events=harvest.reaction_events,
        reaction_run_support=harvest.reaction_run_support,
        species_event_counter=harvest.species_event_counter,
    )
    if reactions_abcd:
        add_reactions(
            reactions=reactions_abcd,
            run_name=run_name,
            reaction_events=harvest.reaction_events_abcd,
            reaction_run_support=harvest.reaction_run_support_abcd,
            species_event_counter=harvest.species_event_counter,
        )

    record = {
        "run": run_name,
        "frames": frames,
        "atoms_in": atoms_in,
        "atoms_out": atoms_out,
        "atom_types_used": atom_types,
        "species_unique": len(present),
        "species_unique_timeline": len(timeline_species),
        "species_unique_json": len(species),
        "species_missing_in_json": len(missing_in_json),
        "reactions_unique": len(reactions),
        "reactionsabcd_unique": len(reactions_abcd),
        "initial_species_unique": len(initial_species),
        "excluded": 0,
        "excluded_reason": "",
    }
    warnings = _support_audit_warning(run_name, missing_in_json)

    template_html = None
    if out_html is not None and out_html.exists():
        template_html = out_html.read_text(encoding="utf-8", errors="ignore")
    return record, warnings, template_html


def _harvest_ensemble(
    files: Sequence[str], settings: Settings, save_noli_dir: Optional[Path]
) -> Tuple[RunHarvest, List[dict], List[Tuple[str, str]], List[Tuple[str, str]], Optional[str]]:
    harvest = RunHarvest()
    run_stats: List[dict] = []
    failures: List[Tuple[str, str]] = []
    warnings: List[Tuple[str, str]] = []
    template_html: Optional[str] = None
    explicit_atom_types = settings.elements

    # mkdtemp, not TemporaryDirectory: the latter registers a finalizer that
    # deletes the tree at interpreter exit, which silently defeated --keep-temp.
    temp_root = Path(tempfile.mkdtemp(prefix="rng_liagg_", dir=settings.tmp_root))
    try:
        for index, in_file in enumerate(files, 1):
            run_name = os.path.basename(in_file)
            print(f"[{index}/{len(files)}] {run_name}")
            run_tmp = temp_root / f"run_{index:03d}"
            run_tmp.mkdir(parents=True, exist_ok=True)
            try:
                record, run_warnings, run_template = _harvest_run(
                    index=index,
                    in_file=in_file,
                    run_tmp=run_tmp,
                    save_noli_dir=save_noli_dir,
                    settings=settings,
                    harvest=harvest,
                    explicit_atom_types=explicit_atom_types,
                )
                run_stats.append(record)
                warnings.extend(run_warnings)
                if template_html is None and run_template is not None:
                    template_html = run_template
            except Exception as exc:  # one bad run must not lose the ensemble
                failures.append((run_name, str(exc)))
            if not settings.keep_temp and run_tmp.exists():
                shutil.rmtree(run_tmp, ignore_errors=True)
    finally:
        if settings.keep_temp:
            print(f"Kept temp root: {temp_root}")
        else:
            shutil.rmtree(temp_root, ignore_errors=True)

    return harvest, run_stats, failures, warnings, template_html


def _write_reacnetstyle_html(
    outdir: Path,
    rngdata: dict,
    settings: Settings,
    harvested_template: Optional[str],
    failures: List[Tuple[str, str]],
) -> None:
    template = None
    if settings.reacnet_template_html:
        template = Path(settings.reacnet_template_html).read_text(
            encoding="utf-8", errors="ignore"
        )
    elif harvested_template is not None:
        template = harvested_template

    target = outdir / artifacts.REACNETSTYLE_HTML
    if template is None:
        message = (
            "No ReacNet HTML template available from successful runs. Pass "
            "--reacnet-template-html, and/or fix per-run failures."
        )
        failures.append(("reacnetstyle_html", message))
        target.write_text(
            "<!doctype html><html><head><meta charset='utf-8'>"
            "<title>ReacNet Aggregate</title></head><body>"
            "<h2>aggregate_reacnet_reacnetstyle.html not rendered</h2>"
            f"<p>{message}</p>"
            "<p>See aggregate_reacnet_failures.tsv for run-level errors.</p>"
            "</body></html>",
            encoding="utf-8",
        )
        return

    try:
        patched, patch_warnings = patch_reacnet_templates_for_support_labels(template)
        for warning in patch_warnings:
            failures.append(("reacnetstyle_template_patch", warning))
        target.write_text(
            inject_rngdata_into_template(patched, rngdata), encoding="utf-8"
        )
    except (RuntimeError, OSError) as exc:
        failures.append(("reacnetstyle_html", str(exc)))


def run_full_pipeline(settings: Settings) -> int:
    files = expand_inputs(settings.inputs)
    if not files:
        print("No input files matched.", file=sys.stderr)
        return EXIT_BAD_INPUT

    outdir = Path(settings.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    save_noli_dir = None
    if settings.save_noli_dir:
        save_noli_dir = Path(settings.save_noli_dir).resolve()
        save_noli_dir.mkdir(parents=True, exist_ok=True)

    if settings.dry_run:
        print("Dry run inputs:")
        for path in files:
            print(path)
        return EXIT_OK

    absent = missing_tools(settings)
    if absent:
        for binary, hint in absent:
            print(f"Required program not found on PATH: {binary}", file=sys.stderr)
            if hint:
                print(f"  {hint}", file=sys.stderr)
        return EXIT_BAD_INPUT

    harvest, run_stats, failures, warnings, harvested_template = _harvest_ensemble(
        files, settings, save_noli_dir
    )

    successful_runs = [r["run"] for r in run_stats if int(r.get("excluded", 0)) == 0]
    excluded_runs = [r["run"] for r in run_stats if int(r.get("excluded", 0)) == 1]
    n_runs = len(successful_runs)

    collapse_failures: List[str] = []
    try:
        ensemble = collapse_harvest(
            harvest=harvest,
            n_runs=n_runs,
            support_source=settings.species_support_source,
            obabel_bin=settings.obabel_bin,
        )
    except CollapseError as exc:
        collapse_failures.append(str(exc))
        print(f"InChI collapse failed: {exc}", file=sys.stderr)
        return EXIT_BAD_INPUT

    table = build_species_table(
        ensemble, n_runs, settings.species_support_source
    )
    warnings.extend(support_audit_warnings(table, len(ensemble.universe())))

    network_html, _support, graph_meta = build_consensus_network(
        reaction_rows_full=ensemble.reaction_rows_full,
        reaction_run_support=ensemble.reaction_run_support,
        species_stats=table.stats,
        species_id_order=table.id_order,
        n_runs=n_runs,
        settings=settings,
        outdir=outdir,
    )
    linkreac = build_linkreac(
        reaction_rows_full=ensemble.reaction_rows_full,
        reaction_rows_abcd_full=ensemble.reaction_rows_abcd_full,
        species_order=table.id_order,
    )

    reaction_rows_simple = artifacts.simple_reaction_rows(ensemble.reaction_rows_full)
    artifacts.write_summary(
        outdir,
        artifacts.build_summary(
            settings=settings,
            ensemble=ensemble,
            table=table,
            graph_meta=graph_meta,
            input_files=files,
            run_stats=run_stats,
            excluded_runs=excluded_runs,
            failures=failures,
            warnings=warnings,
            collapse_failures=collapse_failures,
            n_runs=n_runs,
            reaction_rows_simple=reaction_rows_simple,
        ),
    )
    artifacts.write_reactions(outdir, reaction_rows_simple)
    artifacts.write_species(outdir, table)
    artifacts.write_run_membership(outdir, table, ensemble)
    artifacts.write_alias_map(outdir, table, ensemble)
    artifacts.write_pathlength_scan(outdir, graph_meta.get("pathlength_scan", []))

    write_simple_html(
        outdir / artifacts.REPORT_HTML,
        n_runs=n_runs,
        input_files=successful_runs,
        top_reactions=reaction_rows_simple[: settings.top_n],
        top_species=artifacts.top_species_rows(table, settings.top_n),
    )

    rngdata = build_reacnet_rngdata(
        species_id_order=table.id_order,
        species_stats=table.stats,
        reaction_rows_full=ensemble.reaction_rows_full,
        reaction_rows_abcd_full=ensemble.reaction_rows_abcd_full,
        speciesshownum=DEFAULT_SPECIES_SHOW_NUM,
        reactionsshownum=DEFAULT_REACTIONS_SHOW_NUM,
        network_html=network_html,
        linkreac=linkreac,
    )
    (outdir / artifacts.RNGDATA_JSON).write_text(json.dumps(rngdata), encoding="utf-8")
    _write_reacnetstyle_html(outdir, rngdata, settings, harvested_template, failures)

    if failures:
        write_tsv(outdir / artifacts.FAILURES_TSV, ["run", "error"], failures)
    if warnings:
        write_tsv(outdir / artifacts.WARNINGS_TSV, ["run", "warning"], warnings)

    print("Wrote:")
    for path in artifacts.report_written_paths(outdir, bool(failures), bool(warnings)):
        print(path)

    return EXIT_OK if n_runs > 0 else EXIT_NO_SUCCESSFUL_RUNS
