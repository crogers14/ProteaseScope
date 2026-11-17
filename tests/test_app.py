"""App layer: figures serialize for widgets (no NaN), and the app module imports."""

import json
import sys

import numpy as np
import pytest

from proteasescope import histones as H
from proteasescope import metrics as M
from proteasescope.presets import PRESETS, preset_recipes

from .conftest import ROOT

pytest.importorskip("plotly")
sys.path.insert(0, str(ROOT / "app"))
import figures as F  # noqa: E402


def _strict_json(fig):
    """shinywidgets serializes with allow_nan=False; NaN anywhere breaks the session."""
    return json.dumps(fig.to_plotly_json(), allow_nan=False,
                      default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))  # fmt: skip


@pytest.fixture(scope="module")
def histone_analyses():
    hs, ann, marks = H.load_histone_inputs(ROOT / "data" / "histones")
    cache: dict = {}
    an = [M.Analysis(r.run(hs, cache)) for r in preset_recipes(PRESETS["histone_ptm"])[:7]]
    return hs, ann, marks, an


def test_all_figures_are_strict_json(histone_analyses):
    hs, ann, marks, an = histone_analyses
    labels = [a.label for a in an]
    base = an[0]
    figs = [
        F.mark_grid(H.mark_visibility(an, marks, ann, entries=["P68431", "P62805"]), {}),
        F.isoform_matrix(
            H.discrimination_matrix(H.isoform_discrimination(an, ann), labels[0], "H2A"), "t"
        ),  # fmt: skip
        F.cooccurrence_matrix(
            H.lysine_cooccurrence(an, ann, entries=["P68431"])["matrix"],
            [4, 9, 14, 18, 23, 27, 36, 37, 56, 64, 79, 115, 122],
            "t",
        ),
        F.peptide_map(M.peptide_map(an, "P68431"), labels, "t"),
        F.peptide_map(M.peptide_map(an, "P62805", (1, 30)), labels, "t"),
        F.length_histogram(_hist_frame(base)),
        F.coverage_histogram(base.protein_stats["max_coverage_pct"].to_numpy()),
        F.coverage_distributions(
            {a.label: a.protein_stats["max_coverage_pct"].to_numpy() for a in an}
        ),  # fmt: skip
        F.empty("x"),
    ]
    for fig in figs:
        _strict_json(fig)


def _hist_frame(an):
    from proteasescope.results import LiveResults

    live = LiveResults(an.proteins, an.background, preset_recipes(PRESETS["histone_ptm"])[:1])
    return live.length_histogram(0)


def test_peptide_map_rows_and_lanes(histone_analyses):
    hs, ann, marks, an = histone_analyses
    pm = M.peptide_map(an[:2], "P68431")
    fig = F.peptide_map(pm, [a.label for a in an[:2]], "t")
    bars = [t for t in fig.data if t.type == "bar"]
    assert sum(len(b.y) for b in bars) == len(pm["peptides"])
    # lanes: overlapping peptides of one digest never share a lane
    for b in bars:
        spans = sorted(zip(b.y, b.base, np.asarray(b.base) + np.asarray(b.x), strict=True))
        for (y1, _s1, e1), (y2, s2, _e2) in zip(spans, spans[1:], strict=False):
            if y1 == y2:
                assert s2 >= e1


def test_app_imports():
    pytest.importorskip("shiny")
    pytest.importorskip("shinywidgets")
    import app  # noqa: F401

    assert hasattr(app, "app")


def test_strategy_picker_figures_are_strict_json(histone_analyses):
    from proteasescope import strategy as S
    from proteasescope.presets import PRESETS
    from proteasescope.rules import available_enzymes, get_rule

    hs, ann, marks, _ = histone_analyses
    rules = S.default_candidates([get_rule(n, "psi-ms") for n in available_enzymes("psi-ms")])
    base = PRESETS["histone_ptm"].baseline
    for acc, region, lys in [
        ("P68431", "N-terminal tail", None),
        ("P0C0S8", "C-terminal tail", None),
        ("P62805", None, [16, 20]),
        ("P62805", (30, 31), None),
    ]:
        res = S.pick_strategies(hs, ann, marks, acc, region, base, rules, lysines=lys)
        for goal in S.GOALS:
            cards = res.cards(goal)
            opt = cards[0].option if cards else None
            names = [res.scope.k_name(k) for k in res.scope.lysines]
            if S.GOALS[goal].unit == "lysine":
                fig = F.lysine_picture(res.lysine_states(goal, [c.option for c in cards]), "seen")
            else:
                fig = F.arc_view(res.pair_states(opt), names, opt.label if opt else None)
            _strict_json(fig)
            labels = ["Standard digest"] + ([opt.label] if opt else [])
            _strict_json(F.peptide_map(res.map_data(opt), labels, "t", "standard digest"))


def test_gained_lost_figure_labels_every_bar():
    import pandas as pd

    counts = pd.DataFrame([
        {"digest": "Trypsin | Asp-N", "peptides_gained": 1234567, "peptides_lost": 0,
         "K_sites_gained": 12, "K_sites_lost": 3, "STY_sites_gained": 0, "STY_sites_lost": 7,
         "median_coverage_pct": 82.25},
        {"digest": "Trypsin + Asp-N", "peptides_gained": 5, "peptides_lost": 2000,
         "K_sites_gained": 0, "K_sites_lost": 0, "STY_sites_gained": 4, "STY_sites_lost": 1,
         "median_coverage_pct": 79.9},
    ])  # fmt: skip
    base = {"peptides": 2844568, "K_sites": 100, "STY_sites": 200, "median_coverage_pct": 79.9}
    fig = F.diverging_changes(counts, base)
    _strict_json(fig)
    bars = {(t.xaxis, t.name): list(t.text) for t in fig.data}
    assert bars[("x", "Lost vs baseline")] == ["", "−2,000"]
    assert bars[("x", "Gained vs baseline")] == ["+1,234,567", "+5"]
    assert bars[("x2", "Lost vs baseline")] == ["−3", ""]
    assert bars[("x2", "Gained vs baseline")] == ["+12", "no change"]
    assert bars[("x3", "Gained vs baseline")] == ["", "+4"]
    cov = [t for t in fig.data if t.xaxis == "x4"][0]
    assert list(cov.text) == ["+2.3 pts", "no change"]
    pep_gain = [t for t in fig.data if t.xaxis == "x" and t.name == "Gained vs baseline"][0]
    assert "Net: +1,234,567 (+43.4% of baseline)" in pep_gain.hovertext[0]
    assert "Total: 2,844,568 → 4,079,135" in pep_gain.hovertext[0]
    assert "baseline 2,844,568" in fig.layout.annotations[0].text
    assert "baseline 79.9%" in fig.layout.annotations[3].text
    _strict_json(F.diverging_changes(counts.iloc[:0], base))


def test_builder_panel_ui():
    """One "Build your own strategy" button, one side panel with 3 digest rows and every
    control from CUSTOM_STRATEGY.md section 3."""
    import app

    html = str(app.app_ui)
    assert html.count(">Build your own strategy<") >= 1
    for i in range(3):
        assert f'id="bld_enz{i}"' in html and f'id="bld_chem{i}"' in html
        assert f'id="bld_rm{i}"' in html
    for ctrl in ("bld_name", "bld_start", "bld_add", "bld_combine", "bld_against", "bld_compare",
                 "bld_save", "bld_download", "bld_load"):
        assert f'id="{ctrl}"' in html, ctrl
    assert "Separate digests, results combined" in html and "Same tube, one after another" in html
    assert "Split the sample and run each digest on its own." in html
    assert "Standard digest (current preset)" in html
    assert list(app.CHEM_CHOICES.values()) == [
        "Propionylate before and after digestion", "Propionylate after digestion only",
        "No derivatization"]  # fmt: skip
    assert app.next_strategy_name({"My strategy 1", "My strategy 2"}) == "My strategy 3"
    assert app.next_strategy_name({}) == "My strategy 1"


def test_custom_comparison_figures_are_strict_json():
    from proteasescope import strategy as S
    from proteasescope import workflow as W

    hs, ann, marks = H.load_histone_inputs(ROOT / "data" / "histones")
    ripup = W.builder_preset("RIPUP (Turner et al. 2026)")
    res = S.compare_workflows(hs, ann, marks, "P62805", "N-terminal tail",
                              PRESETS["histone_ptm"].baseline, ripup,
                              against=W.builder_preset("Standard histone digest"))  # fmt: skip
    card = res.card_for(S.Strategy.custom(ripup), "see")
    assert card.headline == "**+0** nothing new for this goal"
    opt = card.option
    pm = res.map_data(opt, f"{opt.label} (combined)", split=True)
    assert pm["labels"] == ["Standard histone digest", "  Arg-C Ultra, no derivatization",
                            "  Recombinant chymotrypsin, no derivatization",
                            "RIPUP (Turner et al. 2026) (combined)"]  # fmt: skip
    # the combined row is the union of the two digest rows
    peps = pm["peptides"]
    parts = peps[peps["digest"].str.startswith("  ")][["start_mature", "end_mature"]]
    comb = peps[peps["digest"].str.endswith("(combined)")][["start_mature", "end_mature"]]
    assert set(map(tuple, parts.to_numpy())) == set(map(tuple, comb.to_numpy()))
    figs = [
        F.peptide_map(pm, pm["labels"], "Peptides", res.ref_name),
        F.arc_view(res.pair_states(opt), [res.scope.k_name(k) for k in res.scope.lysines],
                   opt.label, "pairs", res.ref_label),
        F.lysine_picture(res.lysine_states("see", [opt]), "seen", "lysines", res.ref_name),
    ]  # fmt: skip
    for fig in figs:
        _strict_json(fig)
    assert figs[1].layout.annotations[0].text == "Standard histone digest"
