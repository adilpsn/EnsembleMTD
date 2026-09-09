"""Grouping raw ReacNet SMILES by fixed-H InChI."""

from __future__ import annotations

from collections import Counter, defaultdict

import pytest

from ensemblemtd.collapse import (
    CollapseError,
    build_inchi_collapse_maps,
    collapse_counter,
    collapse_first_seen_order,
    collapse_ratio,
    collapse_reaction_maps,
    collapse_support_map,
    get_obabel_version,
    representative_members,
    smiles_to_inchi_fixedh,
)

# Two SMILES ReacNetGenerator really does emit for the same ring-opened EC
# fragment: one writes the C-C as a double bond, the other as two radical
# centres.  Same atoms, same hydrogens, one fixed-H InChI.
EC_OPEN = "[H][C]([H])=[C]([H])[O][C]([O])=[O]"
EC_OPEN_ALT = "[H][C]([H])[C]([H])[O][C]([O])=[O]"
CO2 = "[O]=[C]=[O]"
CO2_ALT = "[O][C][O]"


def fake_inchi(monkeypatch, mapping):
    """Replace the Open Babel call with a lookup table."""
    monkeypatch.setattr(
        "ensemblemtd.collapse.smiles_to_inchi_fixedh",
        lambda smiles, obabel_bin: mapping[smiles],
    )


@pytest.mark.unit
class TestCollapseMaps:
    def test_two_smiles_with_one_inchi_become_one_group(self, monkeypatch):
        fake_inchi(monkeypatch, {EC_OPEN: "InChI=1/X", EC_OPEN_ALT: "InChI=1/X"})

        raw_to_rep, raw_to_inchi, rep_to_inchi, rep_to_members = (
            build_inchi_collapse_maps(
                {EC_OPEN, EC_OPEN_ALT}, {EC_OPEN: 5, EC_OPEN_ALT: 2}, "obabel"
            )
        )

        assert raw_to_rep == {EC_OPEN: EC_OPEN, EC_OPEN_ALT: EC_OPEN}
        assert raw_to_inchi[EC_OPEN_ALT] == "InChI=1/X"
        assert rep_to_inchi == {EC_OPEN: "InChI=1/X"}
        assert rep_to_members == {EC_OPEN: sorted([EC_OPEN, EC_OPEN_ALT])}

    def test_the_group_is_represented_by_its_best_supported_member(self, monkeypatch):
        fake_inchi(monkeypatch, {EC_OPEN: "InChI=1/X", EC_OPEN_ALT: "InChI=1/X"})
        raw_to_rep, _i, _r, _m = build_inchi_collapse_maps(
            {EC_OPEN, EC_OPEN_ALT}, {EC_OPEN: 1, EC_OPEN_ALT: 9}, "obabel"
        )
        assert set(raw_to_rep.values()) == {EC_OPEN_ALT}

    def test_equal_support_is_broken_on_the_smiles_string_for_reproducibility(
        self, monkeypatch
    ):
        fake_inchi(monkeypatch, {"aaa": "InChI=1/X", "bbb": "InChI=1/X"})
        raw_to_rep, _i, _r, _m = build_inchi_collapse_maps(
            {"aaa", "bbb"}, {"aaa": 3, "bbb": 3}, "obabel"
        )
        assert set(raw_to_rep.values()) == {"aaa"}

    def test_distinct_inchis_stay_distinct(self, monkeypatch):
        fake_inchi(monkeypatch, {EC_OPEN: "InChI=1/X", CO2: "InChI=1/CO2"})
        _rep, _i, rep_to_inchi, members = build_inchi_collapse_maps(
            {EC_OPEN, CO2}, {}, "obabel"
        )
        assert len(rep_to_inchi) == 2
        assert all(len(m) == 1 for m in members.values())

    def test_the_collapse_ratio_reports_how_much_was_merged(self):
        assert collapse_ratio(3, 12) == pytest.approx(0.25)
        assert collapse_ratio(0, 0) == 1.0


@pytest.mark.unit
class TestRekeying:
    def test_run_support_is_unioned_not_summed(self):
        collapsed = collapse_support_map(
            {EC_OPEN: {"r1", "r2"}, EC_OPEN_ALT: {"r2", "r3"}},
            {EC_OPEN: EC_OPEN, EC_OPEN_ALT: EC_OPEN},
        )
        # r2 saw both members; it must still count as one run
        assert collapsed[EC_OPEN] == {"r1", "r2", "r3"}

    def test_a_species_with_no_representative_is_dropped(self):
        assert collapse_support_map({"ghost": {"r1"}}, {}) == {}

    def test_event_counts_are_summed_over_members(self):
        collapsed = collapse_counter(
            Counter({EC_OPEN: 4, EC_OPEN_ALT: 6}),
            {EC_OPEN: EC_OPEN, EC_OPEN_ALT: EC_OPEN},
        )
        assert collapsed == Counter({EC_OPEN: 10})

    def test_the_earliest_first_appearance_wins(self):
        collapsed = collapse_first_seen_order(
            {EC_OPEN: (3, 0), EC_OPEN_ALT: (1, 2)},
            {EC_OPEN: EC_OPEN, EC_OPEN_ALT: EC_OPEN},
        )
        assert collapsed == {EC_OPEN: (1, 2)}

    def test_reactions_are_rekeyed_onto_representatives(self):
        events, support = collapse_reaction_maps(
            Counter({(("ring",), (EC_OPEN,)): 3, (("ring",), (EC_OPEN_ALT,)): 4}),
            defaultdict(
                set,
                {
                    (("ring",), (EC_OPEN,)): {"r1"},
                    (("ring",), (EC_OPEN_ALT,)): {"r2"},
                },
            ),
            {"ring": "ring", EC_OPEN: EC_OPEN, EC_OPEN_ALT: EC_OPEN},
        )
        key = (("ring",), (EC_OPEN,))
        assert events == Counter({key: 7})
        assert support[key] == {"r1", "r2"}

    def test_a_reaction_whose_sides_collapse_together_is_dropped(self):
        events, support = collapse_reaction_maps(
            Counter({((EC_OPEN,), (EC_OPEN_ALT,)): 5}),
            defaultdict(set, {((EC_OPEN,), (EC_OPEN_ALT,)): {"r1"}}),
            {EC_OPEN: EC_OPEN, EC_OPEN_ALT: EC_OPEN},
        )
        assert events == Counter()
        assert support == {}

    def test_a_reaction_with_an_unmapped_side_is_dropped(self):
        events, _s = collapse_reaction_maps(
            Counter({(("ghost",), (EC_OPEN,)): 1}),
            defaultdict(set),
            {EC_OPEN: EC_OPEN},
        )
        assert events == Counter()

    def test_members_default_to_the_representative_itself(self):
        assert representative_members({}, "X") == ["X"]


@pytest.mark.integration
class TestOpenBabel:
    def test_open_babel_reports_a_version(self, obabel):
        assert get_obabel_version(obabel)

    def test_the_two_ring_opened_smiles_share_one_fixed_h_inchi(self, obabel):
        assert smiles_to_inchi_fixedh(EC_OPEN, obabel) == smiles_to_inchi_fixedh(
            EC_OPEN_ALT, obabel
        )

    def test_carbon_dioxide_is_not_grouped_with_the_ring_opened_radical(self, obabel):
        assert smiles_to_inchi_fixedh(CO2, obabel) != smiles_to_inchi_fixedh(
            EC_OPEN, obabel
        )

    def test_the_two_carbon_dioxide_forms_also_share_one_inchi(self, obabel):
        assert smiles_to_inchi_fixedh(CO2, obabel) == smiles_to_inchi_fixedh(
            CO2_ALT, obabel
        )

    def test_hydrogen_count_is_part_of_the_identity(self, obabel):
        # the point of -xF: same C3O3 skeleton, one H apart, so not one group
        assert smiles_to_inchi_fixedh(
            "[H][C]([H])[C]([H])([H])[O][C]([O])=[O]", obabel
        ) != smiles_to_inchi_fixedh(
            "[H][C]([H])([H])[C]([H])([H])[O][C]([O])=[O]", obabel
        )

    def test_an_unparseable_smiles_is_an_error(self, obabel):
        with pytest.raises(CollapseError):
            smiles_to_inchi_fixedh("this is not a molecule (((", obabel)

    def test_a_missing_binary_is_an_error(self):
        with pytest.raises((CollapseError, FileNotFoundError, OSError)):
            get_obabel_version("obabel-that-does-not-exist")
