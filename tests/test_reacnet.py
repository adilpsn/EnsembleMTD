"""Reading ReacNetGenerator's output and building its command line."""

from __future__ import annotations

from pathlib import Path

import pytest

from ensemblemtd.reacnet import (
    load_payload,
    parse_species_timeline,
    reacnet_command,
    run_reacnet,
    tail_lines,
    unwrap_reactions,
    unwrap_species,
)
from ensemblemtd.settings import Settings

pytestmark = pytest.mark.unit


def timeline(tmp_path: Path, *lines: str) -> Path:
    path = tmp_path / "run.species"
    path.write_text("\n".join(lines) + "\n")
    return path


class TestSpeciesTimeline:
    def test_the_first_frame_gives_the_initial_species_with_counts(self, tmp_path):
        initial, _present = parse_species_timeline(
            timeline(tmp_path, "Timestep 0: A 2 B 1", "Timestep 1: A 1")
        )
        assert initial == [("A", 2), ("B", 1)]

    def test_presence_across_frames_is_a_set_not_a_count(self, tmp_path):
        _initial, present = parse_species_timeline(
            timeline(
                tmp_path,
                "Timestep 0: A 2",
                "Timestep 1: A 1 B 1",
                "Timestep 2: A 3",
            )
        )
        assert present == {"A", "B"}

    def test_a_species_with_a_zero_count_is_absent(self, tmp_path):
        initial, present = parse_species_timeline(
            timeline(tmp_path, "Timestep 0: A 1 B 0")
        )
        assert initial == [("A", 1)]
        assert present == {"A"}

    def test_an_empty_first_frame_leaves_no_initial_species(self, tmp_path):
        initial, present = parse_species_timeline(
            timeline(tmp_path, "Timestep 0:", "Timestep 1: A 1")
        )
        assert initial == []
        assert present == {"A"}

    def test_blank_lines_are_skipped(self, tmp_path):
        path = tmp_path / "run.species"
        path.write_text("Timestep 0: A 1\n\n\nTimestep 1: B 1\n")
        _initial, present = parse_species_timeline(path)
        assert present == {"A", "B"}

    def test_a_missing_file_is_an_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="Expected species file not found"):
            parse_species_timeline(tmp_path / "absent.species")

    def test_an_empty_file_is_an_error(self, tmp_path):
        path = tmp_path / "run.species"
        path.write_text("")
        with pytest.raises(RuntimeError, match="Species file is empty"):
            parse_species_timeline(path)

    def test_a_line_that_is_not_a_timestep_is_an_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="Unexpected species line format"):
            parse_species_timeline(timeline(tmp_path, "garbage"))

    def test_an_odd_number_of_tokens_is_an_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="Malformed species line"):
            parse_species_timeline(timeline(tmp_path, "Timestep 0: A 1 B"))

    def test_a_non_numeric_count_is_an_error(self, tmp_path):
        with pytest.raises(RuntimeError, match="Malformed species count"):
            parse_species_timeline(timeline(tmp_path, "Timestep 0: A many"))


class TestJsonUnwrapping:
    def test_a_list_of_dicts_wrapped_in_a_list_is_unwrapped(self):
        assert unwrap_species([[{"s": "A"}, {"s": "B"}]]) == ["A", "B"]

    def test_a_bare_list_of_strings_is_accepted(self):
        assert unwrap_species(["A", "B"]) == ["A", "B"]

    def test_entries_without_a_species_field_are_ignored(self):
        assert unwrap_species([[{"x": 1}, {"s": "A"}]]) == ["A"]

    @pytest.mark.parametrize("raw", [None, {}, 7, "A"])
    def test_anything_that_is_not_a_list_gives_nothing(self, raw):
        assert unwrap_species(raw) == []
        assert unwrap_reactions(raw) == []

    def test_reaction_records_are_unwrapped_and_filtered_to_dicts(self):
        assert unwrap_reactions([[{"l": ["A"], "r": ["B"]}, "junk"]]) == [
            {"l": ["A"], "r": ["B"]}
        ]

    def test_a_payload_round_trips_from_disk(self, tmp_path):
        path = tmp_path / "run.json"
        path.write_text('{"species": [["A"]], "reactions": [[]]}')
        assert load_payload(path) == {"species": [["A"]], "reactions": [[]]}


class TestCommandLine:
    def test_the_command_carries_the_element_list_and_the_settings(self):
        cmd = reacnet_command(
            Path("/tmp/run.noli.xyz"),
            Settings(nohmm=True, nproc=8, maxspecies=25, miso=2, split=3),
            ["C", "H", "O"],
        )
        assert cmd[0] == "reacnetgenerator"
        assert cmd[cmd.index("-i") + 1] == "run.noli.xyz"
        assert cmd[cmd.index("-a") + 1 : cmd.index("-a") + 4] == ["C", "H", "O"]
        assert cmd[cmd.index("-n") + 1] == "8"
        assert cmd[cmd.index("--maxspecies") + 1] == "25"
        assert "--nopbc" in cmd
        assert "--nohmm" in cmd

    def test_hmm_filtering_is_left_on_unless_asked_for(self):
        cmd = reacnet_command(Path("run.xyz"), Settings(), ["C"])
        assert "--nohmm" not in cmd

    def test_a_custom_binary_is_honoured(self):
        cmd = reacnet_command(
            Path("run.xyz"), Settings(reacnet_bin="/opt/rng"), ["C"]
        )
        assert cmd[0] == "/opt/rng"


class TestRunning:
    def test_an_empty_element_list_is_refused_before_launching_anything(self, tmp_path):
        with pytest.raises(RuntimeError, match="No atom types"):
            run_reacnet(tmp_path / "run.xyz", tmp_path, Settings(), [])

    def test_a_failing_binary_reports_its_log_tail(self, tmp_path):
        failing = tmp_path / "boom.sh"
        failing.write_text("#!/bin/sh\necho 'segmentation fault' >&2\nexit 1\n")
        failing.chmod(0o755)
        with pytest.raises(RuntimeError, match="segmentation fault"):
            run_reacnet(
                tmp_path / "run.xyz",
                tmp_path,
                Settings(reacnet_bin=str(failing)),
                ["C"],
            )

    def test_a_silent_success_without_json_is_still_an_error(self, tmp_path):
        quiet = tmp_path / "quiet.sh"
        quiet.write_text("#!/bin/sh\nexit 0\n")
        quiet.chmod(0o755)
        with pytest.raises(RuntimeError, match="Expected output not found"):
            run_reacnet(
                tmp_path / "run.xyz", tmp_path, Settings(reacnet_bin=str(quiet)), ["C"]
            )


def test_the_log_tail_is_capped_and_unreadable_files_are_tolerated(tmp_path):
    log = tmp_path / "log"
    log.write_text("\n".join(str(i) for i in range(100)) + "\n")
    assert tail_lines(log, n=3).splitlines() == ["97", "98", "99"]
    assert tail_lines(tmp_path / "absent") == ""
