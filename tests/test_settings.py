"""Option validation."""

from __future__ import annotations

import pytest

from ensemblemtd.settings import (
    Settings,
    SettingsError,
    dedup,
    parse_pathlength_range,
)

pytestmark = pytest.mark.unit


def base(**kw) -> Settings:
    return Settings(inputs=("run1.trj",), outdir="out", **kw)


def test_a_valid_configuration_survives_validation():
    assert base().validated().outdir == "out"


def test_settings_cannot_be_mutated_after_construction():
    with pytest.raises(Exception):
        base().nproc = 1


def test_either_inputs_or_a_rebuild_source_is_required():
    with pytest.raises(SettingsError, match="Provide either"):
        Settings().validated()


def test_inputs_and_a_rebuild_source_together_are_refused():
    with pytest.raises(SettingsError, match="not both"):
        Settings(inputs=("a.trj",), rebuild_from="d").validated()


def test_inputs_without_an_output_directory_are_refused():
    with pytest.raises(SettingsError, match="--outdir is required"):
        Settings(inputs=("a.trj",)).validated()


def test_a_rebuild_needs_no_output_directory():
    assert Settings(rebuild_from="d").validated().rebuild_from == "d"


def test_a_path_scan_without_a_pathway_mode_is_refused():
    with pytest.raises(SettingsError, match="cannot be used with --scan-pathlength"):
        base(pathway_mode="none", scan_pathlength="2-4").validated()


def test_a_path_scan_is_parsed_into_its_lengths():
    settings = base(pathway_mode="support-first", scan_pathlength="2-4").validated()
    assert settings.scan_pathlength_values == (2, 3, 4)


@pytest.mark.parametrize(
    "field, value",
    [
        ("traj_type", "pdb"),
        ("pathway_mode", "magic"),
        ("target_mode", "somewhere"),
        ("species_support_source", "guess"),
        ("graph_min_edge_runs", 0),
        ("nproc", 0),
    ],
)
def test_an_out_of_range_option_is_refused(field, value):
    with pytest.raises(SettingsError):
        base(**{field: value}).validated()


@pytest.mark.parametrize(
    "text, expected",
    [("", ()), ("   ", ()), ("2-4", (2, 3, 4)), ("3-3", (3,))],
)
def test_a_path_length_range_is_inclusive(text, expected):
    assert parse_pathlength_range(text) == expected


@pytest.mark.parametrize("text", ["4-2", "0-3", "2", "a-b", "1-2-3", "-"])
def test_a_malformed_path_length_range_is_refused(text):
    with pytest.raises(SettingsError):
        parse_pathlength_range(text)


def test_repeated_elements_are_dropped_but_order_is_kept():
    assert dedup(["C", "H", "C", "O", "H"]) == ("C", "H", "O")
    assert dedup(None) is None
