"""Histone strategy picker (STRATEGY_PICKER_2.md section 6)."""

import time

import pandas as pd
import pytest

from proteasescope import histones as H
from proteasescope import strategy as S
from proteasescope.chemistry import Chemistry, Propionylation
from proteasescope.digest import DigestSettings
from proteasescope.fasta import parse_fasta
from proteasescope.presets import PRESETS
from proteasescope.rules import available_enzymes, get_rule

from .conftest import ROOT

GOALS = list(S.GOALS)
PTM = PRESETS["histone_ptm"].baseline


@pytest.fixture(scope="module")
def inputs():
    return H.load_histone_inputs(ROOT / "data" / "histones")


def candidates(mode="psi-ms"):
    return S.default_candidates([get_rule(n, mode) for n in available_enzymes(mode)])


EXPASY_BASE = DigestSettings((get_rule("Arg-C", "expasy"),), missed_cleavages=0, min_length=7,
                             max_length=30, chemistry=Chemistry(carbamidomethyl=False))  # fmt: skip


@pytest.fixture(scope="module")
def ptm_results(inputs):
    """Histone PTM preset: every default entry with each of its region presets."""
    hs, ann, marks = inputs
    out = []
    for acc in H.cooccurrence_entries(ann):
        fam = ann.set_index("accession").loc[acc, "family"]
        for region in S.region_presets(fam, len(hs.sequences[hs.index_of(acc)]) - 1):
            out.append(S.pick_strategies(hs, ann, marks, acc, region, PTM, candidates()))
    return out


@pytest.fixture(scope="module")
def h2b_expasy(inputs):
    hs, ann, marks = inputs
    return S.pick_strategies(hs, ann, marks, "P62807", (1, 35), EXPASY_BASE,
                             candidates("expasy"))  # fmt: skip


def _goal(res, code, goal, n=1):
    s = res.get(code)
    b = res.goal_set(res.base, goal, n)
    a = res.goal_set(res.items[s], goal, n)
    return len(b), len(a), res.names(goal, a - b), res.names(goal, b - a)


# ----- exact cases ------------------------------------------------------------------------------
def test_lysines_seen_h2b(h2b_expasy):
    sc = h2b_expasy.scope
    assert sc.shift == 1
    region_k = "K5 K11 K12 K15 K16 K20 K23 K24 K27 K28 K30 K34".split()
    assert [sc.k_name(k) for k in sc.lysines] == region_k
    assert [(a + 1 - sc.shift, b - sc.shift) for a, b, _ in h2b_expasy.base.peptides] == [(1, 29)]
    assert _goal(h2b_expasy, "parallel + Asp-N", "see") == (10, 12, ["K30", "K34"], [])
    assert _goal(h2b_expasy, "sequential + Asp-N", "see") == (10, 9, ["K34"], ["K27", "K28"])
    b, a, gained, lost = _goal(h2b_expasy, "Lys-N alone (unmodified)", "see")
    assert (b, a, gained, len(lost)) == (10, 1, ["K34"], 10)
    assert _goal(h2b_expasy, "Lys-N alone (unmodified)", "alone", n=1)[2] == ["K34"]


def test_previous_regression_pairs_and_resolution(h2b_expasy):
    """STRATEGY_PICKER_1 regression values that still apply (pairs, <= 2 lysines)."""
    r = h2b_expasy
    assert _goal(r, "parallel + Asp-N", "together")[:2] == (45, 50)
    b, a, gained, lost = _goal(r, "sequential + Asp-N", "together")
    assert (b, a, len(gained), len(lost)) == (45, 28, 0, 17)
    assert _goal(r, "Lys-C alone (unmodified)", "together")[:2] == (45, 0)
    assert _goal(r, "Lys-N alone (unmodified)", "alone", n=2)[:3] == (0, 1, ["K34"])


def test_separable_pairs_h4(inputs):
    hs, ann, marks = inputs
    r = S.pick_strategies(hs, ann, marks, "P62805", (1, 25), EXPASY_BASE, candidates("expasy"))
    s = r.get("parallel + Chymotrypsin low specificity")
    new = r.names("separate", r.goal_set(r.items[s], "separate") - r.goal_set(r.base, "separate"))
    assert new == ["K5–K12", "K5–K16", "K8–K12", "K8–K16"]
    # CUSTOM_STRATEGY.md section 8: all four were tangled (on 4-17) in the standard digest
    assert {(i, j) for i, j in r.items[s].separable} - r.base.together == set()
    sep = r.names("separate", r.items[s].separable)
    assert "K5–K8" not in sep and "K12–K16" not in sep
    # separable needs both lysines seen: K20 is on no identifiable peptide here
    assert "K20" not in r.names("see", r.items[s].seen) and "K5–K20" not in sep
    assert all(i in r.items[s].seen and j in r.items[s].seen for i, j in r.items[s].separable)
    # the baseline peptide 4-17 that holds all four lysines is still there
    assert (3 + r.scope.shift, 17 + r.scope.shift) in r.items[s].key


def test_specific_lysines_scope(inputs):
    """STRATEGY_PICKER_2.md section 6 case, on the v2 candidate set (before the
    literature-defined enzymes)."""
    hs, ann, marks = inputs
    v2 = [r for r in candidates("expasy") if r.mode != "literature"]
    r = S.pick_strategies(hs, ann, marks, "P62805", None, EXPASY_BASE, v2,
                          lysines=[16, 20])  # fmt: skip
    assert [r.scope.k_name(k) for k in r.scope.lysines] == ["K16", "K20"]
    cards = r.cards("together")
    assert [c.option.label for c in cards[:2]] == ["Use Asp-N instead",
                                                   "Add a separate Asp-N digest"]  # fmt: skip
    assert cards[0].headline_items == ["K16–K20"]
    assert all(c.gives_up == "Loses (of the lysines you selected): nothing" for c in cards[:2])
    assert cards[0].extra_work == "same work"
    assert cards[1].extra_work == "+1 sample prep, +1 LC-MS run"
    # every peptide in scope carries a selected lysine
    for s, e, _ in r.items[r.get("Asp-N alone")].peptides:
        assert any(s <= k < e for k in r.scope.lysines)
    with pytest.raises(ValueError):
        S.pick_strategies(hs, ann, marks, "P62805", None, EXPASY_BASE, [], lysines=[17])
    # With the literature-defined enzymes in the set, ProAlanase 16-32 (KRHRKVLRDNIQGITKP) and
    # recombinant chymotrypsin 11-22 (GKGGAKRHRKVL) also read K16 + K20: both tie Asp-N alone
    # (same gain, nothing lost, 1 run) and rank after it by name.
    full = S.pick_strategies(hs, ann, marks, "P62805", None, EXPASY_BASE, candidates("expasy"),
                             lysines=[16, 20])  # fmt: skip
    assert [c.option.label for c in full.cards("together")] == [
        "Use Asp-N instead", "Use ProAlanase instead",
        "Use recombinant chymotrypsin instead"]  # fmt: skip


def test_no_improvement_message(inputs):
    hs, ann, marks = inputs
    r = S.pick_strategies(hs, ann, marks, "P62807", (1, 35), PTM, candidates())
    assert len(r.base.seen) == len(r.scope.lysines) == 12
    assert r.cards("see") == []
    assert r.status("see") == "already complete"
    msg = r.message("see")
    assert msg.startswith("The standard digest already **sees all 12 lysines** in ")
    assert "1–35 with these settings, so no strategy can improve this." in msg


# ----- invariants -------------------------------------------------------------------------------
def test_net_identity_every_goal(ptm_results, h2b_expasy):
    n_checked = 0
    for res in [*ptm_results, h2b_expasy]:
        t = res.advanced_table(n=1)
        for g in S.GOALS.values():
            c = {k: t[f"{g.short}: {k}"] for k in ("standard", "strategy", "gained", "lost", "net")}
            assert (c["net"] == c["strategy"] - c["standard"]).all()
            assert (c["net"] == c["gained"] - c["lost"]).all()
            n_checked += len(t)
        for goal in GOALS:
            for _opt, nums in res.ranked(goal, positive_only=False):
                assert nums["net"] == len(nums["gained"]) - len(nums["lost"])
    assert n_checked > 1000


def test_parallel_never_loses(ptm_results, h2b_expasy):
    for res in [*ptm_results, h2b_expasy]:
        t = res.advanced_table(n=1)
        par = t[t["kind"] == "parallel"]
        assert len(par)
        for g in S.GOALS.values():
            assert (par[f"{g.short}: lost"] == 0).all(), (res.scope.description, g.key)
        assert (par["lysines lost"] == 0).all()


def test_ranking_rules(ptm_results, h2b_expasy):
    for res in [*ptm_results, h2b_expasy]:
        opts = res.options()
        members = [m for o in opts for m in o.members]
        assert len(members) == len(set(members))  # merged duplicates appear once
        assert len({(o.kind, res.items[o.members[0]].key) for o in opts}) == len(opts)
        for o in opts:  # identical to the standard digest: never listed
            assert res.items[o.members[0]].key != res.base.key
        same = [s for s in res.strategies if res.items[s].key == res.base.key]
        assert not set(same) & set(members)
        t = res.advanced_table()
        assert set(t.loc[t["same_as_standard"], "code"]) == {s.code for s in same}
        for goal in GOALS:
            for n in (1, 2):
                ranked = res.ranked(goal, n)
                assert all(nums["net"] > 0 for _, nums in ranked)
                keys = [(-x["net"], len(x["lysines_lost"]), o.costs[2], o.costs[0], o.label)
                        for o, x in ranked]  # fmt: skip
                assert keys == sorted(keys)
                cards = res.cards(goal, n)
                assert len(cards) == min(3, len(ranked))
                assert all(c.net > 0 for c in cards)
                assert (res.message(goal, n) is None) == bool(cards)
                again = [c.option.label for c in res.cards(goal, n)]
                assert [c.option.label for c in cards] == again


def test_merging_and_titles():
    asp = get_rule("Asp-N", "psi-ms")
    both = [S.Strategy("parallel", asp, "propionylated"), S.Strategy("parallel", asp, "unmodified")]
    assert S.option_title(both, PTM)[0] == (
        "Add a separate **Asp-N** digest (with or without propionylation)"
    )
    lysc, lyscp = get_rule("Lys-C", "psi-ms"), get_rule("Lys-C/P", "psi-ms")
    pair = [S.Strategy("alone", lysc, "unmodified"), S.Strategy("alone", lyscp, "unmodified")]
    assert S.option_title(pair, PTM)[0] == (
        "Use **Lys-C** or **Lys-C/P** instead (no propionylation)"
    )
    seq = [S.Strategy("sequential", asp)]
    assert S.option_title(seq, PTM)[0] == "Follow trypsin with **Asp-N** in the same tube"
    assert S.extra_work((1, 2, 1)) == "+1 enzyme step, no extra run"
    assert S.extra_work((2, 2, 2)) == "+1 sample prep, +1 LC-MS run"
    assert S.extra_work((1, 1, 1)) == "same work"
    assert S.short_list(["K27", "K28", "K30", "K34", "K36"]) == "K27, K28, K30 and 2 more"
    assert S.short_list(["K30", "K34"]) == "K30 and K34"
    assert S.standard_line(PTM) == (
        "Compared with: propionylation + trypsin · 4–40 aa · ≤1 missed cleavage"
    )


def test_merge_happens_in_practice(ptm_results):
    merged = [o for res in ptm_results for o in res.options() if len(o.members) > 1]
    assert merged
    for o in merged:
        assert len({m.kind for m in o.members}) == 1
        assert " or " in o.label


# ----- region rules on a synthetic family -----------------------------------------------------
SYNTH = """>sp|T00001|TA_HUMAN Test A OS=Homo sapiens OX=9606 GN=TA PE=1 SV=1
MPAAKAAAKAARGGGGKGGKGGRWWWWWW
>sp|T00002|TB_HUMAN Test B OS=Homo sapiens OX=9606 GN=TB PE=1 SV=1
MPSAKAAAKAARGGGGKGGKGGRWWWWWW
>sp|T00003|TC_HUMAN Test C OS=Homo sapiens OX=9606 GN=TC PE=1 SV=1
MPAAKAAAKAARGGGGKGGKGGRWWAWWW
>sp|T00004|TD_HUMAN Test D OS=Homo sapiens OX=9606 GN=TD PE=1 SV=1
MPAAKAAAKAARGGGGKGGKGGRWWWWWW
"""
# Mature numbering (Met removed): P1 A2 A3 K4 ... K8 ... R11 | G12 ... K16 G17 G18 K19 ... R22 |
# W23-W28. Arg-C peptides: 1-11 (K4, K8), 12-22 (K16, K19), 23-28 (no K).


@pytest.fixture(scope="module")
def synth():
    ps = parse_fasta(SYNTH)
    ann = pd.DataFrame({"accession": ps.accessions, "display_name": ["A", "B", "C", "D"],
                        "family": ["HX"] * 4})  # fmt: skip
    marks = pd.DataFrame({"accession": ["T00001"] * 3, "position_uniprot": [17, 20, 20],
                          "residue": ["K"] * 3, "mod_abbrev": ["ac", "me1", "me2"]})  # fmt: skip
    base = DigestSettings((get_rule("Arg-C", "expasy"),), missed_cleavages=0, min_length=3,
                          max_length=40, chemistry=Chemistry(carbamidomethyl=False))  # fmt: skip
    return ps, ann, marks, base


def test_region_rules(synth):
    ps, ann, marks, base = synth
    res = S.pick_strategies(ps, ann, marks, "T00001", (10, 16), base, [])
    sc, b = res.scope, res.base
    assert [sc.k_name(k) for k in sc.lysines] == ["K16"]
    assert sc.identical == ["D"]
    # in scope: 1-11 (overlaps residues 10-11) and 12-22; 23-28 does not overlap
    assert [(s + 1 - sc.shift, e - sc.shift) for s, e, _ in b.peptides] == [(1, 11), (12, 22)]
    assert res.names("see", res.goal_set(b, "see")) == ["K16"]  # K4, K8, K19 are outside
    assert res.goal_set(b, "alone", n=1) == set()  # peptide 12-22 carries 2 K (K19 counts)
    assert res.names("alone", res.goal_set(b, "alone", n=2)) == ["K16"]
    assert res.names("marks", res.goal_set(b, "marks")) == ["K16"]  # K19's marks: outside
    assert res.goal_set(b, "marks", marks="methylation") == set()  # K16 has only ac
    assert res.goal_set(b, "together") == set()  # K16-K19 share a peptide; K19 is outside
    # B differs at residue 2 (peptide 1-11 overlaps the region); C differs at residue 25
    # (peptide 23-28, outside); D is identical and excluded
    assert b.variants == {"T00002"}
    assert res.goal_total("variants") == 2

    whole = S.pick_strategies(ps, ann, marks, "T00001", "Whole protein", base, [])
    wb = whole.base
    assert whole.names("see", wb.seen) == ["K4", "K8", "K16", "K19"]
    assert whole.names("together", wb.together) == ["K4–K8", "K16–K19"]
    assert whole.names("separate", wb.separable) == ["K4–K16", "K4–K19", "K8–K16", "K8–K19"]
    assert wb.variants == {"T00002", "T00003"}
    assert whole.scope.region == (1, 28)
    assert whole.status("variants") == "already complete"
    assert whole.status("marks", marks="ubiquitination") == "nothing to measure"
    assert whole.message("marks", marks="ubiquitination").startswith("Nothing to measure")


SELECTED = """>sp|T00005|TE_HUMAN Test E OS=Homo sapiens OX=9606 GN=TE PE=1 SV=1
MPAAKAAARGGGGGGRAAKAAR
>sp|T00006|TF_HUMAN Test F OS=Homo sapiens OX=9606 GN=TF PE=1 SV=1
MPAAKAAARGGSGGGRAAKAAR
"""
# Mature: 1-8 (K4), 9-15 (no K; F differs here, G11S), 16-21 (K18).


def test_selected_lysines_use_peptides_that_carry_them(synth):
    _, _, marks, base = synth
    ps = parse_fasta(SELECTED)
    ann = pd.DataFrame({"accession": ps.accessions, "display_name": ["E", "F"],
                        "family": ["HY", "HY"]})  # fmt: skip
    sel = S.pick_strategies(ps, ann, marks, "T00005", None, base, [], lysines=[4, 18])
    sc = sel.scope
    # peptide 9-15 lies between the selected lysines but carries neither: not in scope
    assert [(s + 1 - sc.shift, e - sc.shift) for s, e, _ in sel.base.peptides] == [(1, 8), (16, 21)]
    assert sel.base.variants == set()  # only peptide 9-15 tells E from F
    assert sel.status("variants") == "no improvement"
    span = S.pick_strategies(ps, ann, marks, "T00005", (4, 18), base, [])
    assert span.base.variants == {"T00006"}


def test_region_validation(synth):
    ps, ann, marks, base = synth
    with pytest.raises(ValueError):
        S.pick_strategies(ps, ann, marks, "T00001", (30, 40), base, [])
    assert S.pick_strategies(ps, ann, marks, "T00001", (20, 99), base, []).scope.region == (20, 28)


def test_region_presets():
    assert S.region_presets("H2A", 129) == {"N-terminal tail": (1, 20),
                                            "C-terminal tail": (110, 129),
                                            "Whole protein": None}  # fmt: skip
    assert "C-terminal tail" not in S.region_presets("H2A", 371)  # macroH2A
    assert S.region_presets("H1", 194) == {"Whole protein": None}
    assert S.region_presets("H4", 102)["N-terminal tail"] == (1, 25)


# ----- strategies and chemistry --------------------------------------------------------------
def test_enumeration_and_chemistry():
    rules = candidates()
    assert "Proteinase K" not in [r.name for r in rules]
    strats = S.enumerate_strategies(PTM, rules)
    assert [s.code for s in strats if s.rule.name == "Trypsin"] == [
        "Trypsin alone (unmodified)", "parallel + Trypsin (unmodified)"]  # fmt: skip
    unmod = {s.rule.name for s in strats if s.variant == "unmodified"}
    assert unmod == {"Trypsin", "Trypsin/P", "Lys-C", "Lys-C/P", "Lys-N"}
    k_dep = sum(S.cleaves_at_lysine(r) for r in rules)
    assert len(strats) == 3 * len(rules) + 2 * k_dep - 3  # trypsin: no propionylated three
    assert len({s.code for s in strats}) == len(strats)
    asp = get_rule("Asp-N", "psi-ms")
    (comp,) = S.Strategy("sequential", asp).components(PTM)
    assert comp.chemistry.propionylation is Propionylation.PROTEIN_PEPTIDE_NTERM
    assert comp.rules == (*PTM.rules, asp)
    lysc = get_rule("Lys-C", "psi-ms")
    alone_u = S.Strategy("alone", lysc, "unmodified").components(PTM)
    assert alone_u[0].chemistry.propionylation is Propionylation.NONE
    par = S.Strategy("parallel", lysc).components(PTM)
    assert par[0] == PTM and par[1].chemistry == PTM.chemistry
    assert S.Strategy("sequential", asp).costs == (1, 2, 1)
    with pytest.raises(ValueError):
        S.Strategy("sequential", asp, "unmodified")


def test_cleaves_at_lysine():
    for name, expected in [("Trypsin", True), ("Lys-C", True), ("Lys-N", True),
                           ("Arg-C", False), ("Asp-N", False), ("Glu-C", False),
                           ("Chymotrypsin", False)]:  # fmt: skip
        assert S.cleaves_at_lysine(get_rule(name, "psi-ms")) is expected, name


def test_sequential_k_cleavers_are_never_recommended(ptm_results):
    """Same tube, lysines blocked: sequential Lys-C/Lys-N equals the standard digest."""
    for res in ptm_results:
        listed = {m for o in res.options() for m in o.members}
        for s in res.strategies:
            if s.kind == "sequential" and s.rule.name in ("Lys-C", "Lys-N"):
                assert res.items[s].key == res.base.key and s not in listed


def test_warnings(ptm_results):
    res = ptm_results[0]
    w = {s.code: res.warnings[s] for s in res.strategies}
    assert "unmodified" in w["Lys-C alone (unmodified)"]
    assert "unmodified" in w["parallel + Trypsin (unmodified)"]
    assert "unmodified" not in w["Lys-C alone"]
    assert not any("middle_down" in v for v in w.values())  # max length 40
    assert S.WARNINGS["unmodified"][0] == "Unmodified lysines only"


# ----- all histones and performance -----------------------------------------------------------
def test_compare_all_histones(inputs):
    hs, ann, marks = inputs
    grid = S.compare_all_histones(hs, ann, marks, PTM, candidates())
    rows = S.grid_rows(hs, ann, PTM)
    assert len(rows) == 6  # H2A N + C tail, H2B, H3.1, H3.3, H4
    assert len(grid) == len(rows) * len(GOALS)
    assert set(grid["status"]) <= {"best", "already complete", "no improvement",
                                   "nothing to measure"}  # fmt: skip
    for row in grid.itertuples():
        res = S.pick_strategies(hs, ann, marks, row.accession, row.region_name, PTM,
                                candidates())  # fmt: skip
        b = len(res.goal_set(res.base, row.goal))
        total = res.goal_total(row.goal)
        if row.status == "already complete":
            assert b == total > 0
        elif row.status == "no improvement":
            assert b < total and not res.cards(row.goal)
        elif row.status == "best":
            assert row.best and row.gain.startswith("+")
    h2b_see = grid[(grid["accession"] == "P62807") & (grid["goal"] == "see")]
    assert h2b_see["status"].tolist() == ["already complete"]


def test_ranked_cards_under_one_second(inputs):
    hs, ann, marks = inputs
    for acc in H.cooccurrence_entries(ann):
        t = time.perf_counter()
        res = S.pick_strategies(hs, ann, marks, acc, "Whole protein", PTM, candidates())
        for goal in GOALS:
            res.cards(goal)
        res.advanced_table()
        assert time.perf_counter() - t < 1.0, acc


# ----- ProAlanase -------------------------------------------------------------------------------
PROALANASE_EXPECTED = {
    # accession: (start_mature, end_mature, peptide, lysines) that ProAlanase alone must give
    "P68431": [(1, 7, "ARTKQTA", "K4"), (8, 15, "RKSTGGKA", "K9,K14"), (22, 25, "TKAA", "K23"),
               (31, 38, "ATGGVKKP", "K36,K37")],
    "P62807": [(11, 17, "KKGSKKA", "K11,K12,K15,K16"), (18, 21, "VTKA", "K20")],
}  # fmt: skip


def test_proalanase_is_a_propionylated_candidate():
    rule = get_rule("ProAlanase")
    assert rule in candidates() and rule in candidates("expasy")
    assert not S.cleaves_at_lysine(rule)
    found = S.enumerate_strategies(PTM, [rule], S.KINDS)
    assert {(s.kind, s.variant) for s in found} == {
        ("alone", "propionylated"), ("parallel", "propionylated"),
        ("sequential", "propionylated")}  # fmt: skip


@pytest.mark.parametrize("acc", sorted(PROALANASE_EXPECTED))
def test_proalanase_alone_histone_ptm_peptides(inputs, acc):
    """Histone PTM preset (4-40 aa, <=1 missed cleavage, Met removed, mature numbering)."""
    hs, ann, marks = inputs
    assert (PTM.min_length, PTM.max_length, PTM.missed_cleavages) == (4, 40, 1)
    rule = get_rule("ProAlanase")
    res = S.pick_strategies(hs, ann, marks, acc, S.N_TAIL, PTM, [rule])
    assert res.scope.shift == 1  # Met removed: mature numbering
    pep = res.peptides(S.Strategy("alone", rule))
    got = set(zip(pep.start_mature, pep.end_mature, pep.peptide, pep.K_positions, strict=True))
    missing = set(PROALANASE_EXPECTED[acc]) - got
    assert not missing, missing
