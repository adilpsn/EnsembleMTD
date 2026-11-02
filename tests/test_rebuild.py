"""Re-rendering reports from an existing aggregate directory."""

from __future__ import annotations

import json

import pytest

from conftest import run_cli
from ensemblemtd.artifacts import (
    NETWORK_DOT,
    REACNETSTYLE_HTML,
    REACTIONS_TSV,
    REPORT_HTML,
    RNGDATA_JSON,
    SPECIES_TSV,
    SUMMARY_JSON,
    WARNINGS_TSV,
)
from ensemblemtd.report import write_tsv

pytestmark = pytest.mark.unit

TEMPLATE = (
    "<html><body>\n"
    '<script type="text/x-jsrender" id="specTmpl">OLD SPECIES</script>\n'
    '<script type="text/x-jsrender" id="rTmpl">OLD REACTION</script>\n'
    '<script type="application/json" id="rngdata">{"old": true}</script>\n'
    "</body></html>"
)


@pytest.fixture
def aggregate_dir(tmp_path):
    """A minimal but complete aggregate directory."""
    src = tmp_path / "rcng"
    src.mkdir()
    (src / SUMMARY_JSON).write_text(
        json.dumps(
            {
                "n_success_runs": 20,
                "input_files": ["/scratch/ec/k0.5_a0.6run1.trj"],
                "top_reactions": [
                    {"reaction": "A => B", "run_support": 9, "pct_runs": 45.0,
                     "total_events": 30}
                ],
                "top_species": [
                    {"species": "A", "display_run_support": 20,
                     "pct_runs_observed": 100.0}
                ],
            }
        )
    )
    write_tsv(
        src / REACTIONS_TSV,
        ["reaction", "run_support", "pct_runs", "total_events"],
        [("A => B", 18, 90.0, 42)],
    )
    write_tsv(
        src / SPECIES_TSV,
        ["species", "display_run_support", "pct_runs_observed"],
        [("A", 20, 100.0)],
    )
    (src / NETWORK_DOT).write_text("digraph G {\n  N0 -> N1;\n}\n")
    (src / RNGDATA_JSON).write_text(json.dumps({"species": [[{"s": "A"}]]}))
    (src / REACNETSTYLE_HTML).write_text(TEMPLATE)
    return src


def test_the_reports_are_rebuilt_from_the_tsv_tables(aggregate_dir, tmp_path):
    out = tmp_path / "again"
    result = run_cli(["--rebuild-from", str(aggregate_dir), "--outdir", str(out)])

    # Graphviz is optional; without it the SVG step is a recorded failure
    assert result.returncode in (0, 3)
    html = (out / REPORT_HTML).read_text()
    assert "<b>Successful runs:</b> 20" in html
    # the TSV numbers win over the summary's copy of them
    assert "<td>90.0</td>" in html and "<td>42</td>" in html


def test_the_interactive_report_is_rebuilt_from_the_stored_payload(
    aggregate_dir, tmp_path
):
    out = tmp_path / "again"
    run_cli(["--rebuild-from", str(aggregate_dir), "--outdir", str(out)])
    page = (out / REACNETSTYLE_HTML).read_text()
    assert '{"old": true}' not in page
    assert '{"s": "A"}' in page
    assert "runs={{:rs}}" in page  # row templates were patched again


def test_an_external_template_can_be_supplied(aggregate_dir, tmp_path):
    (aggregate_dir / REACNETSTYLE_HTML).unlink()
    template = tmp_path / "shell.html"
    template.write_text(TEMPLATE.replace("OLD SPECIES", "EXTERNAL SHELL"))

    out = tmp_path / "again"
    run_cli(
        [
            "--rebuild-from", str(aggregate_dir),
            "--outdir", str(out),
            "--reacnet-template-html", str(template),
        ]
    )
    assert '{"s": "A"}' in (out / REACNETSTYLE_HTML).read_text()


def test_with_no_template_anywhere_the_failure_is_recorded(aggregate_dir, tmp_path):
    (aggregate_dir / REACNETSTYLE_HTML).unlink()
    out = tmp_path / "again"
    result = run_cli(["--rebuild-from", str(aggregate_dir), "--outdir", str(out)])

    assert result.returncode == 3
    failures = (out / "aggregate_reacnet_failures.tsv").read_text()
    assert "No template available" in failures


def test_a_missing_payload_is_recorded_rather_than_raised(aggregate_dir, tmp_path):
    (aggregate_dir / RNGDATA_JSON).unlink()
    out = tmp_path / "again"
    result = run_cli(["--rebuild-from", str(aggregate_dir), "--outdir", str(out)])

    assert result.returncode == 3
    assert "Missing required file" in (out / "aggregate_reacnet_failures.tsv").read_text()


def test_a_missing_dot_file_is_recorded(aggregate_dir, tmp_path):
    (aggregate_dir / NETWORK_DOT).unlink()
    out = tmp_path / "again"
    result = run_cli(["--rebuild-from", str(aggregate_dir), "--outdir", str(out)])

    assert result.returncode == 3
    assert "Missing required file" in (out / "aggregate_reacnet_failures.tsv").read_text()


def test_a_missing_summary_falls_back_to_a_limited_report(aggregate_dir, tmp_path):
    (aggregate_dir / SUMMARY_JSON).unlink()
    out = tmp_path / "again"
    run_cli(["--rebuild-from", str(aggregate_dir), "--outdir", str(out)])

    assert "report fallback will be limited" in (out / WARNINGS_TSV).read_text()
    # the tables are still there, so the reaction rows survive
    assert "A =&gt; B" in (out / REPORT_HTML).read_text()


def test_an_unparseable_summary_is_warned_about_not_fatal(aggregate_dir, tmp_path):
    (aggregate_dir / SUMMARY_JSON).write_text("{not json")
    out = tmp_path / "again"
    run_cli(["--rebuild-from", str(aggregate_dir), "--outdir", str(out)])
    assert "Failed to parse summary JSON" in (out / WARNINGS_TSV).read_text()


def test_the_summary_supplies_the_rows_when_the_tables_are_gone(
    aggregate_dir, tmp_path
):
    (aggregate_dir / REACTIONS_TSV).unlink()
    (aggregate_dir / SPECIES_TSV).unlink()
    out = tmp_path / "again"
    run_cli(["--rebuild-from", str(aggregate_dir), "--outdir", str(out)])

    html = (out / REPORT_HTML).read_text()
    # falls back to the summary's copies: 45.0 % and 30 events, not 90.0 / 42
    assert "<td>45.0</td>" in html and "<td>30</td>" in html


def test_rebuilding_in_place_needs_no_output_directory(aggregate_dir):
    result = run_cli(["--rebuild-from", str(aggregate_dir)])
    assert result.returncode in (0, 3)
    assert (aggregate_dir / REPORT_HTML).exists()


def test_a_directory_that_does_not_exist_is_an_error(tmp_path):
    result = run_cli(["--rebuild-from", str(tmp_path / "nowhere")])
    assert result.returncode == 2
    assert "Invalid --rebuild-from" in result.stdout


def test_a_file_is_not_a_valid_rebuild_source(tmp_path):
    not_a_dir = tmp_path / "file.txt"
    not_a_dir.write_text("x")
    result = run_cli(["--rebuild-from", str(not_a_dir)])
    assert result.returncode == 2


def test_a_dry_run_reports_the_directories_and_writes_nothing(aggregate_dir, tmp_path):
    out = tmp_path / "again"
    result = run_cli(
        ["--rebuild-from", str(aggregate_dir), "--outdir", str(out), "--dry-run"]
    )
    assert result.returncode == 0
    assert "Rebuild dry-run from:" in result.stdout
    assert list(out.iterdir()) == []
