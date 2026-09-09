"""Prevalence against ensemble size."""

from __future__ import annotations

import csv

import pytest

from ensemblemtd.convergence import (
    ConvergenceError,
    Species,
    build_curves,
    main,
    native_curve,
    read_ensemble,
    run_number,
    select_species,
    write_csv,
)
from ensemblemtd.report import write_tsv

pytestmark = pytest.mark.unit

numpy = pytest.importorskip("numpy")


def aggregate(tmp_path, rows):
    """Write a species table and membership matrix from (id, smiles, initial, runs)."""
    d = tmp_path / "rcng"
    d.mkdir(exist_ok=True)
    n_runs = len({r for _i, _s, _init, runs in rows for r in runs})
    write_tsv(
        d / "aggregate_reacnet_species.tsv",
        ["species_id", "species", "collapsed_inchi", "observed_run_support",
         "pct_runs_observed", "is_initial"],
        [
            (sid, smiles, f"InChI=1/{sid}", len(runs),
             f"{100.0 * len(runs) / n_runs:.3f}", int(initial))
            for sid, smiles, initial, runs in rows
        ],
    )
    write_tsv(
        d / "aggregate_reacnet_species_run_membership.tsv",
        ["species_id", "species", "observed_run_support", "pct_runs_observed", "runs"],
        [
            (sid, smiles, len(runs), "0.0",
             ";".join(f"k0.5_a0.6run{n}.trj" for n in sorted(runs)))
            for sid, smiles, _initial, runs in rows
        ],
    )
    return d


FOUR_RUNS = [
    ("S1", "EC", True, [1, 2, 3, 4]),
    ("S2", "P_always", False, [1, 2, 3, 4]),
    ("S3", "P_half", False, [1, 3]),
    ("S4", "P_rare", False, [4]),
]


def test_run_numbers_come_out_of_the_trajectory_names():
    assert run_number("k0.5_a0.6run17.trj") == 17
    with pytest.raises(ConvergenceError, match="Cannot read a run number"):
        run_number("mystery.trj")


def test_the_species_table_and_membership_are_read_together(tmp_path):
    species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
    assert runs == [1, 2, 3, 4]
    by_id = {s.species_id: s for s in species}
    assert by_id["S3"].runs == frozenset({1, 3})
    assert by_id["S3"].prevalence == pytest.approx(50.0)
    assert by_id["S1"].is_initial is True


def test_a_missing_table_is_reported_clearly(tmp_path):
    d = aggregate(tmp_path, FOUR_RUNS)
    (d / "aggregate_reacnet_species_run_membership.tsv").unlink()
    with pytest.raises(ConvergenceError, match="Missing aggregate_reacnet_species_run"):
        read_ensemble(d)


def test_an_aggregate_with_no_runs_is_an_error(tmp_path):
    d = tmp_path / "empty"
    d.mkdir()
    write_tsv(d / "aggregate_reacnet_species.tsv", ["species_id"], [])
    write_tsv(
        d / "aggregate_reacnet_species_run_membership.tsv",
        ["species_id", "species", "observed_run_support", "pct_runs_observed", "runs"],
        [],
    )
    with pytest.raises(ConvergenceError, match="No runs recorded"):
        read_ensemble(d)


def test_a_species_present_in_no_run_still_loads(tmp_path):
    species, _runs = read_ensemble(
        aggregate(tmp_path, FOUR_RUNS + [("S5", "P_never", False, [])])
    )
    assert next(s for s in species if s.species_id == "S5").runs == frozenset()


class TestSelection:
    def test_reactants_are_left_out_of_the_automatic_choice(self, tmp_path):
        species, _runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        chosen = select_species(species, top=4)
        assert [s.species_id for s in chosen] == ["S2", "S3", "S4"]

    def test_reactants_can_be_asked_for(self, tmp_path):
        species, _runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        chosen = select_species(species, top=4, include_initial=True)
        assert "S1" in [s.species_id for s in chosen]

    def test_the_automatic_choice_is_ordered_by_run_support(self, tmp_path):
        species, _runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        chosen = select_species(species, top=2)
        assert [s.species_id for s in chosen] == ["S2", "S3"]

    def test_an_explicit_list_keeps_the_order_given(self, tmp_path):
        species, _runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        chosen = select_species(species, smiles=["P_rare", "P_always"])
        assert [s.smiles for s in chosen] == ["P_rare", "P_always"]

    def test_an_unknown_smiles_is_reported_rather_than_skipped(self, tmp_path):
        species, _runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        with pytest.raises(ConvergenceError, match="Not in the aggregate"):
            select_species(species, smiles=["not-a-species"])


class TestCurves:
    def test_the_native_curve_is_a_running_percentage(self):
        assert native_curve([1, 0, 1, 1]) == pytest.approx(
            (100.0, 50.0, 200 / 3, 75.0)
        )

    def test_a_species_seen_everywhere_stays_at_a_hundred(self):
        assert native_curve([1, 1, 1]) == (100.0, 100.0, 100.0)

    def test_the_curve_ends_at_the_full_ensemble_prevalence(self, tmp_path):
        species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        curves = build_curves(select_species(species, top=3), runs, n_orderings=200)
        for curve in curves:
            assert curve.native[-1] == pytest.approx(curve.species.prevalence)

    def test_the_band_closes_at_the_full_ensemble_size(self, tmp_path):
        # sampling is without replacement, so at N = N_total every ordering
        # gives the same answer and the band has to collapse
        species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        curve = build_curves(select_species(species, top=1), runs, n_orderings=200)[0]
        assert curve.band_lo[-1] == pytest.approx(curve.band_hi[-1])
        assert curve.half_width_at(len(runs)) == pytest.approx(0.0)

    def test_the_band_mean_equals_the_full_prevalence_at_every_size(self, tmp_path):
        # exchangeability, not convergence -- which is why the mean is not the
        # diagnostic the report prints
        species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        curve = build_curves(select_species(species, top=1), runs, n_orderings=500)[0]
        for value in curve.band_mean:
            assert value == pytest.approx(curve.species.prevalence, abs=1e-9)

    def test_the_band_brackets_the_native_curve(self, tmp_path):
        species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        for curve in build_curves(select_species(species, top=3), runs, n_orderings=500):
            for lo, hi in zip(curve.band_lo, curve.band_hi):
                assert lo <= hi

    def test_the_same_seed_gives_the_same_band(self, tmp_path):
        species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        chosen = select_species(species, top=1)
        first = build_curves(chosen, runs, n_orderings=100, seed=7)[0]
        again = build_curves(chosen, runs, n_orderings=100, seed=7)[0]
        assert first.band_lo == again.band_lo

    def test_a_flat_species_settles_immediately(self, tmp_path):
        species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        always = select_species(species, smiles=["P_always"])
        curve = build_curves(always, runs, n_orderings=100)[0]
        assert curve.settles_at(tolerance_pp=5.0) == 1

    def test_a_late_appearing_species_settles_late(self, tmp_path):
        species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        rare = select_species(species, smiles=["P_rare"])  # only run 4 of 4
        curve = build_curves(rare, runs, n_orderings=100)[0]
        assert curve.settles_at(tolerance_pp=5.0) == 4

    def test_a_half_width_outside_the_ensemble_is_undefined(self, tmp_path):
        species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
        curve = build_curves(select_species(species, top=1), runs, n_orderings=50)[0]
        assert curve.half_width_at(99) is None
        assert curve.half_width_at(0) is None


def test_the_csv_carries_one_row_per_species_and_size(tmp_path):
    species, runs = read_ensemble(aggregate(tmp_path, FOUR_RUNS))
    curves = build_curves(select_species(species, top=3), runs, n_orderings=100)
    out = tmp_path / "curves" / "convergence.csv"
    write_csv(out, curves)

    rows = list(csv.DictReader(out.open()))
    assert len(rows) == 3 * len(runs)
    assert rows[0]["N"] == "1"
    assert set(rows[0]) == {
        "species_id", "smiles", "collapsed_inchi", "is_initial", "pi_full",
        "N", "pi_native", "pi_band_mean", "pi_band_lo", "pi_band_hi",
    }


class TestCommandLine:
    def test_it_reports_the_diagnostics_and_writes_the_csv(self, tmp_path, capsys):
        d = aggregate(tmp_path, FOUR_RUNS)
        out = tmp_path / "convergence.csv"
        code = main(
            [str(d), "--out-csv", str(out), "--orderings", "100", "--production-size", "2"]
        )
        assert code == 0
        printed = capsys.readouterr().out
        assert "4 runs, 3 species traced" in printed
        assert "band half-width at N=2" in printed
        assert out.exists()

    def test_a_figure_is_written_when_asked_for(self, tmp_path):
        pytest.importorskip("matplotlib")
        d = aggregate(tmp_path, FOUR_RUNS)
        code = main(
            [
                str(d),
                "--out-csv", str(tmp_path / "c.csv"),
                "--plot", str(tmp_path / "c.pdf"),
                "--orderings", "50",
                "--production-size", "2",
            ]
        )
        assert code == 0
        assert (tmp_path / "c.pdf").stat().st_size > 0

    def test_a_bad_aggregate_directory_exits_with_a_message(self, tmp_path, capsys):
        assert main([str(tmp_path / "nowhere")]) == 2
        assert "Missing" in capsys.readouterr().out

    def test_an_unknown_species_exits_with_a_message(self, tmp_path, capsys):
        d = aggregate(tmp_path, FOUR_RUNS)
        assert main([str(d), "--species", "nope"]) == 2
        assert "Not in the aggregate" in capsys.readouterr().out

    def test_selecting_nothing_exits_with_a_message(self, tmp_path, capsys):
        d = aggregate(tmp_path, FOUR_RUNS)
        assert main([str(d), "--top", "0"]) == 2
        assert "No species selected" in capsys.readouterr().out
