"""Strategy model, per-digest chemistry and custom strategies (CUSTOM_STRATEGY.md section 9)."""

import json

import pytest

from proteasescope import histones as H
from proteasescope import strategy as S
from proteasescope import workflow as W
from proteasescope.chemistry import Chemistry, Propionylation
from proteasescope.digest import DigestSettings
from proteasescope.presets import PRESETS, DigestRecipe
from proteasescope.rules import available_enzymes, get_rule

from .conftest import ROOT

GOALS = list(S.GOALS)
PTM = PRESETS["histone_ptm"].baseline
#: RIPUP regression settings: PSI-MS rules, 7-30 aa, 0 missed cleavages, Met removed (MetAP).
RIPUP_BASE = DigestSettings(
    (get_rule("Trypsin", "psi-ms"),), missed_cleavages=0, min_length=7, max_length=30,
    chemistry=Chemistry(propionylation="before_and_after"),
)  # fmt: skip
RIPUP = W.builder_preset("RIPUP (Turner et al. 2026)")
STANDARD = W.builder_preset("Standard histone digest")


@pytest.fixture(scope="module")
def inputs():
    return H.load_histone_inputs(ROOT / "data" / "histones")


def _candidates():
    return S.default_candidates([get_rule(n) for n in available_enzymes("psi-ms")])


def _goal_sets(res, strategy):
    s = res.get(strategy)
    return {g: res.goal_set(res.items[s], g) for g in GOALS}


# ----- the model ----------------------------------------------------------------------------------
def test_builder_presets():
    assert [(d.enzyme.name, d.chemistry.value) for d in RIPUP.digests] == [
        ("Arg-C Ultra", "none"), ("Recombinant chymotrypsin", "none")]  # fmt: skip
    assert RIPUP.combine == "separate"
    assert RIPUP.work == "2 sample preps, 2 enzymes, 2 LC-MS runs"
    argc_prop = W.builder_preset("Arg-C Ultra + propionylation (Turner et al. 2026)")
    assert [(d.enzyme.name, d.chemistry) for d in argc_prop.digests] == [
        ("Arg-C Ultra", Propionylation.AFTER_ONLY)]  # fmt: skip
    assert [(d.enzyme.name, d.chemistry) for d in STANDARD.digests] == [
        ("Trypsin", Propionylation.BEFORE_AND_AFTER)]  # fmt: skip
    for p in W.BUILDER_PRESETS.values():
        assert p.citation
    assert "Turner" in W.BUILDER_PRESETS["RIPUP (Turner et al. 2026)"].citation
    order = W.builder_enzymes(available_enzymes("psi-ms"))
    assert order[:2] == ["Arg-C Ultra", "Recombinant chymotrypsin"]
    assert sorted(order) == sorted(available_enzymes("psi-ms"))


def test_costs():
    tryp, aspn = get_rule("Trypsin"), get_rule("Asp-N")
    steps = (W.DigestStep(tryp), W.DigestStep(aspn))
    assert W.Workflow("x", steps, "separate").costs == (2, 2, 2)
    assert W.Workflow("x", steps, "same_tube").costs == (1, 2, 1)
    assert W.Workflow("x", steps, "same_tube").work == "1 sample prep, 2 enzymes, 1 LC-MS run"
    # the picker's kinds keep their costs
    for kind in S.KINDS:
        s = S.Strategy(kind, aspn)
        assert s.workflow_for(PTM).costs == S.COSTS[kind]


def test_json_round_trip(tmp_path):
    text = RIPUP.to_json()
    data = json.loads(text)
    assert data == {"name": "RIPUP (Turner et al. 2026)", "combine": "separate", "digests": [
        {"enzyme": "Arg-C Ultra", "source": "literature", "chemistry": "none"},
        {"enzyme": "Recombinant chymotrypsin", "source": "literature",
         "chemistry": "none"}]}  # fmt: skip
    assert W.Workflow.from_json(text) == RIPUP
    # old chemistry names migrate; enzymes resolve in the given mode
    old = {"name": "old", "combine": "same_tube", "digests": [
        {"enzyme": "Arg-C", "chemistry": "protein-level"},
        {"enzyme": "Asp-N", "chemistry": "protein + peptide N-term"}]}  # fmt: skip
    w = W.Workflow.from_dict(old, mode="expasy")
    assert [d.enzyme.label for d in w.digests] == ["Arg-C (ExPASy)", "Asp-N (ExPASy)"]
    assert {d.chemistry for d in w.digests} == {Propionylation.BEFORE_AND_AFTER}
    # a custom enzyme keeps its regex
    from proteasescope.rules import custom_rule

    cust = W.Workflow("c", (W.DigestStep(custom_rule(after="W"), "none"),))
    back = W.Workflow.from_json(cust.to_json())
    assert back.digests[0].enzyme.regex == "(?<=[W])"
    assert back.digests[0].enzyme.mode == "custom"


def test_description_for_methods():
    assert RIPUP.description == (
        "RIPUP (Turner et al. 2026): Arg-C Ultra, no derivatization; Recombinant chymotrypsin, "
        "no derivatization (separate digests, results combined)")  # fmt: skip


# ----- guard rails -----------------------------------------------------------------------------
def test_same_tube_chemistries_are_synchronized():
    steps = (W.DigestStep(get_rule("Trypsin"), "before_and_after"),
             W.DigestStep(get_rule("Asp-N"), "none"))  # fmt: skip
    w = W.Workflow("x", steps, "same_tube")
    assert [d.chemistry for d in w.digests] == [Propionylation.BEFORE_AND_AFTER] * 2
    assert w.synced
    assert (None, W.NOTE_SAME_TUBE, False) in w.notes()
    sep = W.Workflow("x", steps, "separate")
    assert not sep.synced and [d.chemistry.value for d in sep.digests] == ["before_and_after",
                                                                           "none"]  # fmt: skip
    assert W.NOTE_SAME_TUBE not in [n for _, n, _ in sep.notes()]


def test_guard_rail_notes():
    w = W.Workflow("x", (W.DigestStep(get_rule("Lys-C"), "before_and_after"),
                         W.DigestStep(get_rule("Glu-C"), "before_and_after"),
                         W.DigestStep(get_rule("Asp-N"), "none")))  # fmt: skip
    notes = w.notes()
    assert (0, W.NOTE_K_BLOCKED, False) in notes
    assert (1, W.NOTE_K_BLOCKED, False) not in notes  # Glu-C does not cut at K
    assert (2, W.NOTE_NO_DERIVATIZATION, False) in notes
    assert w.can_compare
    for name in ("Trypsin", "Lys-N"):
        n = W.Workflow("x", (W.DigestStep(get_rule(name), "before_and_after"),)).notes()
        assert (0, W.NOTE_K_BLOCKED, False) in n
    assert W.Workflow("x", (W.DigestStep(get_rule("Lys-C"), "after_only"),)).notes() == []


def test_identical_rows_block_compare(inputs):
    hs, ann, marks = inputs
    step = W.DigestStep(get_rule("Asp-N"), "none")
    w = W.Workflow("dup", (step, step, W.DigestStep(get_rule("Asp-N"), "before_and_after")))
    assert w.duplicates() == [1]
    assert (1, W.NOTE_DUPLICATE, True) in w.notes()
    assert not w.can_compare
    with pytest.raises(ValueError, match="already in the strategy"):
        w.validate()
    with pytest.raises(ValueError, match="already in the strategy"):
        S.compare_workflows(hs, ann, marks, "P62805", None, PTM, w)


# ----- equivalence (one code path) -------------------------------------------------------------
def test_standard_only_custom_shows_no_change(inputs):
    hs, ann, marks = inputs
    for acc in ("P62805", "P68431", "P62807", "P0C0S8"):
        r = S.compare_workflows(hs, ann, marks, acc, None, PTM, STANDARD)
        summary = r.all_goals_summary(STANDARD.name)
        assert len(summary) == len(GOALS)
        assert (summary[["Gained", "Lost", "Net"]] == 0).all().all(), acc
        assert r.items[r.get(STANDARD.name)].key == r.base.key


@pytest.mark.parametrize("enzyme", ["Asp-N", "Glu-C", "Chymotrypsin", "ProAlanase",
                                    "Recombinant chymotrypsin", "Lys-C"])  # fmt: skip
def test_custom_separate_equals_parallel(inputs, enzyme):
    hs, ann, marks = inputs
    rule = get_rule(enzyme)
    custom = W.Workflow("mine", (W.DigestStep(get_rule("Trypsin"), "before_and_after"),
                                 W.DigestStep(rule, "before_and_after")), "separate")  # fmt: skip
    custom_u = W.Workflow("mine (u)", (W.DigestStep(get_rule("Trypsin"), "before_and_after"),
                                       W.DigestStep(rule, "none")), "separate")  # fmt: skip
    for acc in ("P62805", "P62807"):
        r = S.pick_strategies(hs, ann, marks, acc, None, PTM, [rule], all_unmodified=True,
                              custom=[custom, custom_u])  # fmt: skip
        for code, mine in ((f"parallel + {enzyme}", "mine"),
                           (f"parallel + {enzyme} (unmodified)", "mine (u)")):  # fmt: skip
            a, b = r.get(code), r.get(mine)
            assert r.items[a].key == r.items[b].key
            assert _goal_sets(r, a) == _goal_sets(r, b)
            assert a.components(PTM) == b.components(PTM)


@pytest.mark.parametrize("enzyme", ["Asp-N", "Glu-C", "Chymotrypsin low specificity",
                                    "ProAlanase", "Recombinant chymotrypsin"])  # fmt: skip
def test_custom_same_tube_equals_sequential(inputs, enzyme):
    hs, ann, marks = inputs
    rule = get_rule(enzyme)
    custom = W.Workflow("mine", (W.DigestStep(get_rule("Trypsin"), "before_and_after"),
                                 W.DigestStep(rule, "before_and_after")), "same_tube")  # fmt: skip
    for acc in ("P62805", "P62807", "P68431"):
        r = S.pick_strategies(hs, ann, marks, acc, None, PTM, [rule], custom=[custom])
        a, b = r.get(f"sequential + {enzyme}"), r.get("mine")
        assert r.items[a].key == r.items[b].key
        assert _goal_sets(r, a) == _goal_sets(r, b)
        assert a.components(PTM) == b.components(PTM)


def test_picker_kinds_match_previous_settings():
    """alone/parallel/sequential expressed with the model give the settings they always had."""
    from dataclasses import replace

    aspn, lysc = get_rule("Asp-N"), get_rule("Lys-C")
    prop = PTM.chemistry
    unmod = replace(prop, propionylation="none")
    assert S.Strategy("baseline").components(PTM) == (PTM,)
    assert S.Strategy("alone", aspn).components(PTM) == (
        replace(PTM, rules=(aspn,), name=None, chemistry=prop),)  # fmt: skip
    assert S.Strategy("parallel", lysc, "unmodified").components(PTM) == (
        PTM, replace(PTM, rules=(lysc,), name=None, chemistry=unmod))  # fmt: skip
    assert S.Strategy("sequential", aspn).components(PTM) == (
        replace(PTM, rules=(*PTM.rules, aspn), name=None, chemistry=prop),)  # fmt: skip
    # proteome recipes also go through the model and keep the baseline chemistry
    base = PRESETS["proteome"].baseline
    assert DigestRecipe("baseline", base).settings == (base,)
    assert DigestRecipe("parallel", base, aspn).settings == (base, base.with_rules(aspn))
    assert DigestRecipe("sequential", base, aspn).settings == (
        base.with_rules(*base.rules, aspn),)  # fmt: skip


def test_compare_against_saved_strategy(inputs):
    hs, ann, marks = inputs
    r = S.compare_workflows(hs, ann, marks, "P62805", "N-terminal tail", RIPUP_BASE, STANDARD,
                            against=RIPUP)  # fmt: skip
    assert r.ref_name == RIPUP.name
    s = r.all_goals_summary(STANDARD.name).set_index("goal")
    assert s.loc[S.GOALS["see"].title, "Compared-against → Custom"] == "5 → 4"
    assert s.loc[S.GOALS["see"].title, "lost_items"] == "K20"
    assert any(line == f"Compared against: {RIPUP.name}" for line in r.header_lines())
    assert any(line.startswith("Custom strategy: RIPUP") for line in r.header_lines())


def test_custom_strategies_rank_like_any_other(inputs):
    hs, ann, marks = inputs
    r = S.pick_strategies(hs, ann, marks, "P0C0S8", None, RIPUP_BASE, _candidates(),
                          custom=[RIPUP])  # fmt: skip
    labels = [o.label for o, _ in r.ranked("see")]
    assert RIPUP.name in labels
    card = r.card_for(RIPUP.name, "see")
    assert card.work == "2 sample preps, 2 enzymes, 2 LC-MS runs"
    assert "Work: 2 sample preps, 2 enzymes, 2 LC-MS runs" in card.text()
    table = r.advanced_table()
    row = table[table["code"] == RIPUP.name].iloc[0]
    assert row["kind"] == "custom" and row["variant"] == "custom"
    assert row["enzyme"] == "Arg-C Ultra (Literature), Recombinant chymotrypsin (Literature)"


# ----- separable-pairs correction and RIPUP regression ------------------------------------------
@pytest.fixture(scope="module")
def ripup_h4(inputs):
    hs, ann, marks = inputs
    return S.compare_workflows(hs, ann, marks, "P62805", (1, 25), RIPUP_BASE, RIPUP,
                               against=STANDARD)  # fmt: skip


def test_ripup_h4_tail(ripup_h4):
    r = ripup_h4
    assert r.scope.shift == 1  # Met removed, mature numbering
    s = r.get(RIPUP.name)
    seen_b, seen_a = r.goal_set(r.base, "see"), r.goal_set(r.items[s], "see")
    assert (len(seen_b), len(seen_a)) == (4, 5)
    assert r.names("see", seen_a - seen_b) == ["K20"]
    assert seen_b - seen_a == set()
    tog_b, tog_a = r.goal_set(r.base, "together"), r.goal_set(r.items[s], "together")
    assert r.names("together", tog_a - tog_b) == ["K12–K20", "K16–K20"]
    sep_b, sep_a = r.goal_set(r.base, "separate"), r.goal_set(r.items[s], "separate")
    assert r.names("separate", sep_a - sep_b) == ["K5–K12", "K5–K16", "K8–K12", "K8–K16"]


def test_separable_pairs_need_a_tangled_pair(ripup_h4):
    """K5-K20 and K8-K20 are separable with RIPUP but never shared a peptide in the standard
    digest (K20 is not seen there), so they are not counted."""
    r = ripup_h4
    s = r.get(RIPUP.name)
    raw = r.names("separate", r.items[s].separable)
    assert {"K5–K20", "K8–K20"} <= set(raw)
    counted = r.names("separate", r.goal_set(r.items[s], "separate"))
    assert "K5–K20" not in counted and "K8–K20" not in counted
    assert r.goal_total("separate") == len(r.base.together) == 6
    arcs = r.pair_states(S.Option([s], RIPUP.name, RIPUP.name, []))
    assert sorted(arcs.loc[arcs["category"] == "new separable", "pair"]) == [
        "K5–K12", "K5–K16", "K8–K12", "K8–K16"]  # fmt: skip


def test_ripup_h2a_whole_protein(inputs):
    hs, ann, marks = inputs
    r = S.compare_workflows(hs, ann, marks, "P0C0S8", None, RIPUP_BASE, RIPUP, against=STANDARD)
    s = r.get(RIPUP.name)
    seen_b, seen_a = r.goal_set(r.base, "see"), r.goal_set(r.items[s], "see")
    assert (len(seen_b), len(seen_a)) == (3, 13)
    assert r.names("see", seen_a - seen_b) == [
        "K13", "K15", "K74", "K75", "K99", "K118", "K119", "K125", "K127", "K129"]  # fmt: skip
    assert seen_b - seen_a == set()
    # the standard digest and the builder's standard preset are the same digest
    r2 = S.compare_workflows(hs, ann, marks, "P0C0S8", None, RIPUP_BASE, RIPUP)
    assert r2.base.key == r.base.key
