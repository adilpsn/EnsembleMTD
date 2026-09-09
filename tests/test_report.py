"""TSV, HTML and ReacNet-shell output."""

from __future__ import annotations

import json

import pytest

from ensemblemtd.report import (
    build_reacnet_rngdata,
    inject_rngdata_into_template,
    link_or_copy,
    load_top_rows_from_tsv,
    patch_reacnet_templates_for_support_labels,
    write_simple_html,
    write_tsv,
)

pytestmark = pytest.mark.unit

TEMPLATE = (
    "<html><body>\n"
    '<script type="text/x-jsrender" id="specTmpl">OLD SPECIES</script>\n'
    '<script type="text/x-jsrender" id="rTmpl">OLD REACTION</script>\n'
    '<script type="application/json" id="rngdata">{"old": true}</script>\n'
    "</body></html>"
)


def stats(**kw):
    base = {
        "collapsed_inchi": "InChI=1/X",
        "representative_smiles": "A",
        "member_raw_species_count": 2,
        "display_run_support": 5,
        "observed_run_support": 4,
        "observed_run_support_json": 3,
        "support_delta_json": 1,
        "support_source": "timeline",
        "involvement_events": 7,
        "is_initial": True,
        "initial_frame_count": 2,
        "pct_runs_observed": 80.0,
    }
    base.update(kw)
    return base


class TestTsv:
    def test_rows_are_tab_separated_with_a_header(self, tmp_path):
        path = tmp_path / "t.tsv"
        write_tsv(path, ["a", "b"], [(1, "x"), (2, "y")])
        assert path.read_text() == "a\tb\n1\tx\n2\ty\n"

    def test_a_table_with_no_rows_still_gets_its_header(self, tmp_path):
        path = tmp_path / "t.tsv"
        write_tsv(path, ["a", "b"], [])
        assert path.read_text() == "a\tb\n"

    def test_the_first_rows_can_be_read_back_with_conversions(self, tmp_path):
        path = tmp_path / "t.tsv"
        write_tsv(path, ["name", "n"], [("a", 1), ("b", 2), ("c", 3)])
        rows = load_top_rows_from_tsv(
            path, ["name", "n"], {"name": str, "n": int}, top_n=2
        )
        assert rows == [("a", 1), ("b", 2)]

    def test_a_missing_column_yields_nothing_rather_than_a_crash(self, tmp_path):
        path = tmp_path / "t.tsv"
        write_tsv(path, ["name"], [("a",)])
        assert load_top_rows_from_tsv(path, ["name", "n"], {}, top_n=5) == []

    def test_an_absent_file_yields_nothing(self, tmp_path):
        assert load_top_rows_from_tsv(tmp_path / "gone.tsv", ["a"], {}, 5) == []

    def test_a_row_that_will_not_convert_is_skipped(self, tmp_path):
        path = tmp_path / "t.tsv"
        path.write_text("name\tn\na\t1\nb\tnope\nc\t3\n")
        rows = load_top_rows_from_tsv(
            path, ["name", "n"], {"name": str, "n": int}, top_n=5
        )
        assert rows == [("a", 1), ("c", 3)]

    def test_a_short_row_is_skipped(self, tmp_path):
        path = tmp_path / "t.tsv"
        path.write_text("name\tn\na\t1\nb\n")
        rows = load_top_rows_from_tsv(
            path, ["name", "n"], {"name": str, "n": int}, top_n=5
        )
        assert rows == [("a", 1)]


class TestSummaryPage:
    def test_the_page_reports_the_runs_and_both_tables(self, tmp_path):
        path = tmp_path / "r.html"
        write_simple_html(
            path,
            n_runs=20,
            input_files=["/scratch/ec/k0.5_a0.6run1.trj"],
            top_reactions=[("A => B", 18, 90.0, 42)],
            top_species=[("A", 20, 100.0)],
        )
        html = path.read_text()
        assert "<b>Successful runs:</b> 20" in html
        assert "k0.5_a0.6run1.trj" in html
        assert "/scratch/ec/" not in html  # basenames only
        assert "A =&gt; B" in html  # reaction arrows are escaped
        assert "<td>90.0</td>" in html and "<td>42</td>" in html

    def test_species_names_are_escaped(self, tmp_path):
        path = tmp_path / "r.html"
        write_simple_html(path, 1, [], [], [("<script>x</script>", 1, 100.0)])
        assert "<script>x</script>" not in path.read_text()
        assert "&lt;script&gt;" in path.read_text()


class TestReacnetShell:
    def test_the_payload_carries_the_ensemble_fields(self):
        rngdata = build_reacnet_rngdata(
            species_id_order=["A"],
            species_stats={"A": stats()},
            reaction_rows_full=[((("A",), ("B",)), 4, 80.0, 9)],
            reaction_rows_abcd_full=[],
            speciesshownum=30,
            reactionsshownum=20,
            network_html="<svg/>",
            linkreac={"A": ["B"]},
        )
        species = rngdata["species"][0][0]
        assert species["s"] == "A" and species["i"] == 1
        assert species["rs"] == 5 and species["rs_obs"] == 4
        assert species["m"] == 2 and species["init"] == 1
        reaction = rngdata["reactions"][0][0]
        assert reaction == {"i": 1, "l": ["A"], "r": ["B"], "n": 9, "rs": 4}
        assert rngdata["network"] == ["<svg/>"]
        assert rngdata["linkreac"] == {"A": ["B"]}

    def test_a_non_initial_species_is_flagged_as_such(self):
        rngdata = build_reacnet_rngdata(
            ["A"], {"A": stats(is_initial=False)}, [], [], 30, 20, "", {}
        )
        assert rngdata["species"][0][0]["init"] == 0

    def test_the_payload_replaces_the_templates_own_data(self):
        out = inject_rngdata_into_template(TEMPLATE, {"new": True})
        assert '{"old": true}' not in out
        assert '{"new": true}' in out
        assert out.count("</script>") == TEMPLATE.count("</script>")

    def test_a_template_without_a_data_block_is_an_error(self):
        with pytest.raises(RuntimeError, match="script id=rngdata"):
            inject_rngdata_into_template("<html></html>", {})

    def test_an_unterminated_data_block_is_an_error(self):
        with pytest.raises(RuntimeError, match="closing </script>"):
            inject_rngdata_into_template('<script id="rngdata">{}', {})

    def test_the_row_templates_gain_run_support_labels(self):
        out, warnings = patch_reacnet_templates_for_support_labels(TEMPLATE)
        assert warnings == []
        assert "OLD SPECIES" not in out and "OLD REACTION" not in out
        assert "runs={{:rs}}" in out
        assert "obs={{:rs_obs}}" in out

    def test_a_template_missing_a_row_block_warns_but_still_returns_html(self):
        out, warnings = patch_reacnet_templates_for_support_labels(
            '<script id="specTmpl">OLD</script>'
        )
        assert "runs={{:rs}}" in out
        assert warnings == ["template patch skipped: could not find script id=rTmpl"]

    def test_patching_then_injecting_gives_a_usable_page(self):
        patched, _warnings = patch_reacnet_templates_for_support_labels(TEMPLATE)
        page = inject_rngdata_into_template(patched, {"species": [[]]})
        assert json.loads(page.split('id="rngdata">')[1].split("</script>")[0]) == {
            "species": [[]]
        }


class TestLinkOrCopy:
    def test_the_content_arrives_at_the_destination(self, tmp_path):
        src = tmp_path / "a.xyz"
        src.write_text("payload")
        dst = tmp_path / "b.xyz"
        link_or_copy(src, dst)
        assert dst.read_text() == "payload"

    def test_an_existing_destination_is_replaced(self, tmp_path):
        src = tmp_path / "a.xyz"
        src.write_text("new")
        dst = tmp_path / "b.xyz"
        dst.write_text("old")
        link_or_copy(src, dst)
        assert dst.read_text() == "new"

    def test_a_copy_is_used_when_linking_is_not_possible(self, tmp_path, monkeypatch):
        monkeypatch.setattr("os.link", lambda s, d: (_ for _ in ()).throw(OSError))
        monkeypatch.setattr("os.symlink", lambda s, d: (_ for _ in ()).throw(OSError))
        src = tmp_path / "a.xyz"
        src.write_text("payload")
        dst = tmp_path / "b.xyz"
        link_or_copy(src, dst)
        assert dst.read_text() == "payload"
        assert not dst.is_symlink()
