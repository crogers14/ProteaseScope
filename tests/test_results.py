"""Results facade: precomputed and live results must give identical tables."""

import numpy as np
import pandas as pd
import pytest

from proteasescope import precomputed as P
from proteasescope.metrics import Analysis, Background, SharedPrefixIndex
from proteasescope.presets import PRESETS, preset_recipes, recipes_for
from proteasescope.results import LiveResults, PrecomputedResults
from proteasescope.rules import get_rule

from .conftest import ROOT

OUT = ROOT / "data" / "precomputed"


def _assert_same(pre, live, indices):
    for i in indices:
        pd.testing.assert_frame_equal(pre.entry_stats(i), live.entry_stats(i), check_dtype=False)
        pd.testing.assert_frame_equal(pre.entry_changes(i), live.entry_changes(i),
                                      check_dtype=False)  # fmt: skip
        assert pre.summary(i) == pytest.approx(live.summary(i))
        assert pre.comparison_counts(i) == live.comparison_counts(i)
        pd.testing.assert_frame_equal(pre.length_histogram(i), live.length_histogram(i),
                                      check_dtype=False)  # fmt: skip
        pd.testing.assert_frame_equal(pre.aa_detectability(i), live.aa_detectability(i),
                                      check_dtype=False)  # fmt: skip
        for key, df in pre.seen_vs_unseen(i).items():
            pd.testing.assert_frame_equal(df, live.seen_vs_unseen(i)[key], check_dtype=False)


@pytest.fixture(scope="module")
def histone_pair(histones, tmp_path_factory):
    preset = PRESETS["histone_ptm"]
    bg = Background(histones, index=SharedPrefixIndex.build(histones))
    out = tmp_path_factory.mktemp("pre")
    P.write_preset(P.build_preset(preset, histones, bg), preset, out)
    data = P.load_preset(out, preset.key, histones.accessions)
    recipes = preset_recipes(preset)
    return (PrecomputedResults(histones, bg, recipes, data), LiveResults(histones, bg, recipes))


def test_precomputed_equals_live_on_histones(histone_pair):
    pre, live = histone_pair
    assert pre.labels == live.labels and len(pre.labels) == 46
    _assert_same(pre, live, range(len(pre.recipes)))
    pd.testing.assert_frame_equal(pre.ranking_table(), live.ranking_table(), check_dtype=False)


def test_summary_table_and_coverage_pair(histone_pair):
    pre, _ = histone_pair
    t = pre.summary_table([5, 0, 7])
    assert t["digest"].tolist() == [pre.labels[0], pre.labels[5], pre.labels[7]]
    assert (t.loc[0, [c for c in t.columns if c.endswith("_diff")]] == 0).all()
    cp = pre.coverage_pair(5)
    assert len(cp) == len(pre.proteins)


def test_live_peptide_detail_matches_full_analysis(histone_pair, histones):
    pre, live = histone_pair
    full = live.analysis(6).peptide_table()
    key = ["accession", "start", "end"]
    sub = pre.peptide_table(6, ["P68431", "P62805"]).sort_values(key, ignore_index=True)
    want = full[full.accession.isin(["P68431", "P62805"])].sort_values(key, ignore_index=True)
    pd.testing.assert_frame_equal(sub, want, check_dtype=False)
    ans = pre.protein_analyses("P68431", [0, 6])
    assert [a.label for a in ans] == [pre.labels[0], pre.labels[6]]
    assert all(isinstance(a, Analysis) for a in ans)


def test_recipe_labels_match_digest_labels(histones):
    base = PRESETS["proteome"].baseline
    for r in recipes_for(base, [get_rule("Asp-N", "psi-ms")]):
        assert r.run(histones).label == r.label


def test_header_lists_selected_digests(histone_pair):
    pre, _ = histone_pair
    head = "\n".join(pre.header("Output 5", [0, 3]))
    assert "PREDICTED" in head and "Baseline:" in head and "Alternative 1:" in head


@pytest.mark.slow
@pytest.mark.skipif(not (OUT / "manifest.json").exists(), reason="run scripts/precompute.py")
def test_precomputed_equals_live_on_proteome():
    ref, idx, _ = P.load_reference(OUT)
    bg = Background(ref, index=idx)
    preset = PRESETS["middle_down"]
    recipes = preset_recipes(preset)
    pre = PrecomputedResults(ref, bg, recipes, P.load_preset(OUT, preset.key, ref.accessions))
    live = LiveResults(ref, bg, recipes)
    _assert_same(pre, live, [0, 2, 3])
    assert np.isfinite(pre.summary(0)["median_coverage_pct"])
