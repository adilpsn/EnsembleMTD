"""Reading trajectories and stripping Li."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensemblemtd.trajectory import (
    detect_elements_xyz,
    expand_inputs,
    natural_key,
    parse_comment_energy,
    strip_li_xyz,
)

DATA = Path(__file__).parent / "data"

pytestmark = pytest.mark.unit


def test_stripping_removes_every_lithium(tmp_path):
    out = tmp_path / "noli.xyz"
    frames, atoms_in, atoms_out, elements, zero_frame = strip_li_xyz(
        DATA / "k0.5_a0.6run1.trj", out
    )

    assert frames == 5
    assert atoms_in == 5 * 8
    assert atoms_out == 5 * 6  # two Li per frame
    assert "Li" not in elements
    assert zero_frame is None
    # comment lines carry the "Li-removed" marker, so check the atom lines
    atom_lines = [
        line for line in out.read_text().splitlines() if line.startswith(("C", "O", "H", "Li"))
    ]
    assert atom_lines and not any(line.startswith("Li") for line in atom_lines)


def test_stripping_keeps_elements_in_order_of_first_appearance(tmp_path):
    _f, _i, _o, elements, _z = strip_li_xyz(
        DATA / "k0.5_a0.6run1.trj", tmp_path / "noli.xyz"
    )
    assert elements == ["C", "O", "H"]


def test_stripped_frames_keep_the_original_comment_plus_a_marker(tmp_path):
    out = tmp_path / "noli.xyz"
    strip_li_xyz(DATA / "k0.5_a0.6run1.trj", out)
    lines = out.read_text().splitlines()

    assert lines[0] == "6"
    assert lines[1].startswith(" energy: -41.83")
    assert lines[1].endswith("| Li-removed")


def test_a_zero_energy_frame_is_reported_with_its_one_based_index(tmp_path):
    _f, _i, _o, _e, zero_frame = strip_li_xyz(
        DATA / "k0.5_a0.6run4_zeroenergy.trj", tmp_path / "noli.xyz"
    )
    # the third frame's energy was zeroed
    assert zero_frame == 3


def test_a_truncated_atom_block_is_an_error(tmp_path):
    with pytest.raises(RuntimeError, match="truncated XYZ atom block"):
        strip_li_xyz(DATA / "k0.5_a0.6run5_truncated.trj", tmp_path / "noli.xyz")


def test_a_non_numeric_atom_count_is_an_error(tmp_path):
    src = tmp_path / "bad.xyz"
    src.write_text("not-a-number\ncomment\nC 0.0 0.0 0.0\n")
    with pytest.raises(RuntimeError, match="invalid XYZ frame atom-count line"):
        strip_li_xyz(src, tmp_path / "out.xyz")


def test_a_missing_comment_line_is_an_error(tmp_path):
    src = tmp_path / "bad.xyz"
    src.write_text("1\n")
    with pytest.raises(RuntimeError, match="missing comment line"):
        strip_li_xyz(src, tmp_path / "out.xyz")


@pytest.mark.parametrize(
    "comment, expected",
    [
        (" energy: -41.838678746550 gnorm: 0.41 xtb: 6.7.1", -41.83867874655),
        ("energy: 0.000000000000", 0.0),
        ("energy:1.5e-3", 1.5e-3),
        ("MTD params: charge=0 kpush_factor=0.5", None),
        ("", None),
    ],
)
def test_energy_is_read_from_the_comment_line(comment, expected):
    assert parse_comment_energy(comment) == expected


def test_run_numbers_sort_numerically_not_lexically():
    names = ["run10.trj", "run2.trj", "run1.trj"]
    assert sorted(names, key=natural_key) == ["run1.trj", "run2.trj", "run10.trj"]


def test_inputs_are_globbed_deduplicated_and_run_ordered(tmp_path):
    for name in ("k0.5run1.trj", "k0.5run2.trj", "k0.5run10.trj"):
        (tmp_path / name).write_text("1\nc\nC 0 0 0\n")

    resolved = expand_inputs([str(tmp_path / "*.trj"), str(tmp_path / "k0.5run1.trj")])

    assert [Path(p).name for p in resolved] == [
        "k0.5run1.trj",
        "k0.5run2.trj",
        "k0.5run10.trj",
    ]


def test_a_pattern_matching_nothing_contributes_nothing(tmp_path):
    assert expand_inputs([str(tmp_path / "nothing-here-*.trj")]) == []


def test_elements_can_be_detected_without_stripping_first(tmp_path):
    out = tmp_path / "noli.xyz"
    strip_li_xyz(DATA / "k0.5_a0.6run1.trj", out)
    assert detect_elements_xyz(out) == ["C", "O", "H"]


def test_detecting_elements_in_an_empty_file_is_an_error(tmp_path):
    empty = tmp_path / "empty.xyz"
    empty.write_text("")
    with pytest.raises(RuntimeError, match="Could not auto-detect"):
        detect_elements_xyz(empty)


def test_a_lithium_only_frame_leaves_no_elements_to_detect(tmp_path):
    src = tmp_path / "li.xyz"
    src.write_text("2\ncomment\nLi 0 0 0\nLi 1 0 0\n")
    with pytest.raises(RuntimeError, match="Use --elements explicitly"):
        detect_elements_xyz(src)
