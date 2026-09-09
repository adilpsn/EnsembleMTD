"""End-to-end runs of the aggregator against a replayed ReacNet.

Three synthetic runs of a 2 Li + 2 EC system: two open the ring and lose CO2,
one only opens the ring.  So the expected prevalences are 100% for the ring and
the ring-opened fragment, and 2/3 for CO2 -- and the two SMILES the stub emits
for each product have to collapse onto single nodes for those numbers to come
out right.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from conftest import DATA, run_cli

EC_RING = "[H][C]1([H])[O][C]([O])[O][C]1([H])[H]"

pytestmark = pytest.mark.integration


def species_by_inchi(outdir: Path) -> dict:
    with (outdir / "aggregate_reacnet_species.tsv").open() as handle:
        return {row["collapsed_inchi"]: row for row in csv.DictReader(handle, delimiter="\t")}


def formula(inchi: str) -> str:
    return inchi.split("/")[1]


@pytest.fixture
def aggregate(tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel):
    outdir = tmp_path / "rcng"
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(outdir),
            "--nohmm",
            "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
        ]
    )
    assert result.returncode == 0, result.stderr
    return outdir


def test_every_expected_artifact_is_written(aggregate):
    expected = {
        "aggregate_reacnet_summary.json",
        "aggregate_reacnet_reactions.tsv",
        "aggregate_reacnet_species.tsv",
        "aggregate_reacnet_species_run_membership.tsv",
        "aggregate_reacnet_inchi_alias_map.tsv",
        "aggregate_reacnet_network.dot",
        "aggregate_reacnet_network_nodes.tsv",
        "aggregate_reacnet_mechanistic_path.tsv",
        "aggregate_reacnet_report.html",
        "aggregate_reacnet_rngdata.json",
        "aggregate_reacnet_reacnetstyle.html",
    }
    assert expected <= {p.name for p in aggregate.iterdir()}


def test_the_two_smiles_per_product_collapse_onto_one_node(aggregate):
    summary = json.loads((aggregate / "aggregate_reacnet_summary.json").read_text())
    # 6 raw SMILES (ring, 2x open, 2x CO2, vinoxy) -> 4 chemical species
    assert summary["n_raw_species"] == 6
    assert summary["n_collapsed_species"] == 4
    assert summary["collapse_mode"] == "inchi_fixedh_full_analysis"
    assert summary["obabel_version"]


def test_prevalence_counts_runs_not_frames(aggregate):
    species = species_by_inchi(aggregate)
    by_formula = {formula(k): v for k, v in species.items()}

    assert by_formula["C3H4O3"]["pct_runs_observed"] == "100.000"  # the EC ring
    assert by_formula["C3H3O3"]["pct_runs_observed"] == "100.000"  # ring-opened
    assert by_formula["CO2"]["pct_runs_observed"] == "66.667"      # 2 of 3 runs
    assert by_formula["C2H3O"]["observed_run_support"] == "0"      # only a product


def test_a_collapsed_node_reports_how_many_raw_smiles_it_absorbed(aggregate):
    by_formula = {formula(k): v for k, v in species_by_inchi(aggregate).items()}
    assert by_formula["C3H3O3"]["member_raw_species_count"] == "2"
    assert by_formula["CO2"]["member_raw_species_count"] == "2"
    assert by_formula["C3H4O3"]["member_raw_species_count"] == "1"


def test_the_reactant_is_flagged_as_initial_and_ordered_first(aggregate):
    with (aggregate / "aggregate_reacnet_species.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["species_id"] == "S1"
    assert rows[0]["species"] == EC_RING
    assert rows[0]["is_initial"] == "1"
    assert all(row["is_initial"] == "0" for row in rows[1:])


def test_the_alias_map_names_every_raw_smiles_and_its_group(aggregate):
    with (aggregate / "aggregate_reacnet_inchi_alias_map.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 6
    groups = {}
    for row in rows:
        groups.setdefault(row["species_id"], set()).add(row["raw_species_smiles"])
    assert sorted(len(v) for v in groups.values()) == [1, 1, 2, 2]


def test_the_membership_table_lists_the_runs_behind_each_count(aggregate):
    with (aggregate / "aggregate_reacnet_species_run_membership.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    for row in rows:
        runs = [r for r in row["runs"].split(";") if r]
        assert len(runs) == int(row["observed_run_support"])
        assert runs == sorted(runs)


def test_reactions_are_reported_with_run_support(aggregate):
    with (aggregate / "aggregate_reacnet_reactions.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    ring_opening = next(r for r in rows if r["reaction"].startswith(EC_RING))
    # all three runs opened the ring, with 3 + 4 + 5 events between them
    assert ring_opening["run_support"] == "3"
    assert ring_opening["pct_runs"] == "100.0"
    assert ring_opening["total_events"] == "12"


def test_the_dot_graph_labels_edges_with_both_e_and_r(aggregate):
    dot = (aggregate / "aggregate_reacnet_network.dot").read_text()
    assert dot.startswith("digraph G {")
    assert dot.rstrip().endswith("}")
    assert "label=\"e=" in dot and ", r=" in dot


def test_an_edge_below_the_run_support_floor_is_not_drawn(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    strict = tmp_path / "strict"
    loose = tmp_path / "loose"
    for outdir, floor in ((strict, "3"), (loose, "1")):
        result = run_cli(
            [
                "--inputs", *ensemble_inputs,
                "--outdir", str(outdir),
                "--nohmm", "--nproc", "2",
                "--reacnet-bin", str(reacnet_stub),
                "--graph-min-edge-runs", floor,
            ]
        )
        assert result.returncode == 0, result.stderr

    def edges(d):
        return json.loads((d / "aggregate_reacnet_summary.json").read_text())["graph"][
            "candidate_edges"
        ]

    assert edges(strict) < edges(loose)


def test_a_zero_energy_run_is_excluded_rather_than_analysed(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    outdir = tmp_path / "rcng"
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            str(DATA / "k0.5_a0.6run4_zeroenergy.trj"),
            "--outdir", str(outdir),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
        ]
    )
    assert result.returncode == 0, result.stderr

    summary = json.loads((outdir / "aggregate_reacnet_summary.json").read_text())
    assert summary["n_input_runs"] == 4
    assert summary["n_success_runs"] == 3
    assert summary["excluded_runs"] == ["k0.5_a0.6run4_zeroenergy.trj"]
    # prevalence still divides by the three runs that were usable
    by_formula = {formula(k): v for k, v in species_by_inchi(outdir).items()}
    assert by_formula["CO2"]["pct_runs_observed"] == "66.667"


def test_a_broken_run_is_recorded_as_a_failure_and_the_ensemble_survives(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    outdir = tmp_path / "rcng"
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            str(DATA / "k0.5_a0.6run5_truncated.trj"),
            "--outdir", str(outdir),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
        ]
    )
    assert result.returncode == 0, result.stderr

    summary = json.loads((outdir / "aggregate_reacnet_summary.json").read_text())
    assert summary["n_failed_runs"] == 1
    assert summary["n_success_runs"] == 3
    assert "truncated" in summary["failures"][0]["error"]
    failures = (outdir / "aggregate_reacnet_failures.tsv").read_text()
    assert "k0.5_a0.6run5_truncated.trj" in failures


def test_the_li_free_trajectories_can_be_kept(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    noli = tmp_path / "noli"
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(tmp_path / "rcng"),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
            "--save-noli-dir", str(noli),
        ]
    )
    assert result.returncode == 0, result.stderr
    written = sorted(p.name for p in noli.glob("*.noli.xyz"))
    assert len(written) == 3
    assert "Li " not in (noli / written[0]).read_text()


def test_keeping_the_temp_directory_actually_keeps_it(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    tmp_root = tmp_path / "scratch"
    tmp_root.mkdir()
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(tmp_path / "rcng"),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
            "--keep-temp",
            "--tmp-root", str(tmp_root),
        ]
    )
    assert result.returncode == 0, result.stderr
    kept = list(tmp_root.glob("rng_liagg_*"))
    assert len(kept) == 1
    assert sorted(p.name for p in kept[0].iterdir()) == ["run_001", "run_002", "run_003"]


def test_the_temp_directory_is_removed_by_default(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    tmp_root = tmp_path / "scratch"
    tmp_root.mkdir()
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(tmp_path / "rcng"),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
            "--tmp-root", str(tmp_root),
        ]
    )
    assert result.returncode == 0, result.stderr
    assert list(tmp_root.iterdir()) == []


def test_the_interactive_report_carries_the_aggregate_and_the_support_labels(aggregate):
    html = (aggregate / "aggregate_reacnet_reacnetstyle.html").read_text()
    assert '{"placeholder": true}' not in html
    assert "ORIGINAL SPECIES ROW" not in html
    assert "runs={{:rs}}" in html
    assert EC_RING in html


def test_the_reports_can_be_rebuilt_without_touching_the_trajectories(
    aggregate, tmp_path
):
    rebuilt = tmp_path / "again"
    result = run_cli(["--rebuild-from", str(aggregate), "--outdir", str(rebuilt)])
    # graphviz may be absent, which is a recorded failure rather than a crash
    assert result.returncode in (0, 3), result.stderr
    assert (rebuilt / "aggregate_reacnet_report.html").exists()
    assert (rebuilt / "aggregate_reacnet_reacnetstyle.html").exists()


def test_a_dry_run_lists_the_inputs_and_writes_nothing(
    tmp_path, ensemble_inputs
):
    outdir = tmp_path / "rcng"
    result = run_cli(
        ["--inputs", *ensemble_inputs, "--outdir", str(outdir), "--dry-run"]
    )
    assert result.returncode == 0
    assert result.stdout.startswith("Dry run inputs:")
    assert len([line for line in result.stdout.splitlines()[1:] if line.strip()]) == 3
    assert list(outdir.iterdir()) == []


def test_inputs_that_match_nothing_are_an_error(tmp_path):
    result = run_cli(
        ["--inputs", str(tmp_path / "nothing-*.trj"), "--outdir", str(tmp_path / "o")]
    )
    assert result.returncode == 2
    assert "No input files matched" in result.stderr


def test_the_installed_console_script_runs(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    from conftest import run_cli_subprocess

    result = run_cli_subprocess(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(tmp_path / "rcng"),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
        ]
    )
    assert result.returncode == 0, result.stderr
    assert "Wrote:" in result.stdout


def test_an_external_reacnet_template_is_preferred_over_the_harvested_one(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    template = tmp_path / "shell.html"
    template.write_text(
        '<html><script id="specTmpl">X</script><script id="rTmpl">Y</script>'
        '<script id="rngdata">{"old": true}</script></html>'
    )
    outdir = tmp_path / "rcng"
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(outdir),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
            "--reacnet-template-html", str(template),
        ]
    )
    assert result.returncode == 0, result.stderr
    page = (outdir / "aggregate_reacnet_reacnetstyle.html").read_text()
    assert '{"old": true}' not in page
    assert "runs={{:rs}}" in page


def test_a_template_without_the_data_block_is_recorded_as_a_failure(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    template = tmp_path / "useless.html"
    template.write_text("<html><body>nothing here</body></html>")
    outdir = tmp_path / "rcng"
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(outdir),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
            "--reacnet-template-html", str(template),
        ]
    )
    assert result.returncode == 0, result.stderr
    failures = (outdir / "aggregate_reacnet_failures.tsv").read_text()
    assert "script id=rngdata" in failures


def test_every_run_failing_is_reported_as_such(
    tmp_path, reacnet_stub, stub_env, obabel
):
    outdir = tmp_path / "rcng"
    result = run_cli(
        [
            "--inputs", str(DATA / "k0.5_a0.6run5_truncated.trj"),
            "--outdir", str(outdir),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
        ]
    )
    # exit 3: nothing was aggregated, but the failure record was still written
    assert result.returncode == 3
    summary = json.loads((outdir / "aggregate_reacnet_summary.json").read_text())
    assert summary["n_success_runs"] == 0
    assert summary["n_failed_runs"] == 1
    fallback = (outdir / "aggregate_reacnet_reacnetstyle.html").read_text()
    assert "not rendered" in fallback


def test_explicit_elements_override_what_was_detected(
    tmp_path, reacnet_stub, stub_env, ensemble_inputs, obabel
):
    outdir = tmp_path / "rcng"
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(outdir),
            "--nohmm", "--nproc", "2",
            "--reacnet-bin", str(reacnet_stub),
            "--elements", "C", "H", "O", "C",
        ]
    )
    assert result.returncode == 0, result.stderr
    summary = json.loads((outdir / "aggregate_reacnet_summary.json").read_text())
    # repeats dropped, order kept
    assert summary["run_stats"][0]["atom_types_used"] == ["C", "H", "O"]


@pytest.mark.unit
def test_a_missing_external_program_is_reported_before_any_work(tmp_path, ensemble_inputs):
    outdir = tmp_path / "rcng"
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(outdir),
            "--reacnet-bin", "definitely-not-a-real-binary",
        ]
    )
    assert result.returncode == 2
    assert "Required program not found on PATH" in result.stderr
    assert "definitely-not-a-real-binary" in result.stderr
    # and nothing half-written to confuse the next reader
    assert list(outdir.iterdir()) == []


@pytest.mark.unit
def test_the_hint_names_how_to_install_reacnetgenerator(tmp_path, ensemble_inputs):
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(tmp_path / "rcng"),
            "--reacnet-bin", "reacnetgenerator-not-installed-here",
        ]
    )
    assert result.returncode == 2
    # the basename does not match a known tool, so no hint is invented
    assert "pip install" not in result.stderr


@pytest.mark.unit
def test_a_dry_run_does_not_require_the_external_programs(tmp_path, ensemble_inputs):
    result = run_cli(
        [
            "--inputs", *ensemble_inputs,
            "--outdir", str(tmp_path / "rcng"),
            "--reacnet-bin", "definitely-not-a-real-binary",
            "--dry-run",
        ]
    )
    assert result.returncode == 0
    assert result.stdout.startswith("Dry run inputs:")
