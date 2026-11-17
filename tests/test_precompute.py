"""Phase 4: shared-prefix uniqueness index, presets and precomputed data."""

import json
import random

import numpy as np
import pytest

from proteasescope import precomputed as P
from proteasescope.chemistry import Chemistry, Propionylation
from proteasescope.digest import DigestSettings, digest
from proteasescope.fasta import ProteinSet, parse_fasta
from proteasescope.metrics import (
    Analysis,
    Background,
    SharedPrefixIndex,
    compare,
    comparison_counts,
    peptide_mapping,
    ranking_table,
    uniqueness_background,
)
from proteasescope.presets import PRESETS, alternatives, preset_digests
from proteasescope.rules import available_enzymes, get_rule

from .conftest import ROOT

OUT = ROOT / "data" / "precomputed"


def _settings(name, mode="psi-ms", **kw):
    kw.setdefault("chemistry", Chemistry())
    return DigestSettings((get_rule(name, mode),), **kw)


# ----- index correctness ------------------------------------------------------------------
def _random_set(seed=3, n=40):
    """Random proteins sharing fragments (to create long shared prefixes), with I/L swaps."""
    rng = random.Random(seed)
    alpha = "ACDEFGHIKLMNPQRSTVWY"
    blocks = ["".join(rng.choice(alpha) for _ in range(rng.randint(5, 40))) for _ in range(15)]
    lines = []
    for i in range(n):
        parts = [rng.choice(blocks) if rng.random() < 0.5 else
                 "".join(rng.choice(alpha) for _ in range(rng.randint(3, 30)))
                 for _ in range(rng.randint(2, 8))]  # fmt: skip
        seq = "M" + "".join(parts)
        if i % 5 == 0:
            seq = seq.replace("L", "I")
        lines.append(f">sp|R{i}|R{i}_HUMAN Random OS=Homo sapiens GN=G{i // 3}\n{seq}")
    return parse_fasta("\n".join(lines))


def test_index_matches_brute_force():
    ps = _random_set()
    idx = SharedPrefixIndex.build(ps)
    folded = [s.replace("I", "L") for s in ps.sequences]
    genes = ps.meta["gene"].tolist()
    off = 0
    for e, s in enumerate(folded):
        others = [q for j, q in enumerate(folded) if j != e]
        other_genes = [q for j, q in enumerate(folded) if genes[j] != genes[e]]
        for i in range(len(s)):
            best = max((L for L in range(1, len(s) - i + 1)
                        if any(s[i:i + L] in q for q in others)), default=0)  # fmt: skip
            best_g = max((L for L in range(1, len(s) - i + 1)
                          if any(s[i:i + L] in q for q in other_genes)), default=0)  # fmt: skip
            assert idx.other_entry[off + i] == min(best, 255), (e, i)
            assert idx.other_gene[off + i] == min(best_g, 255), (e, i)
        off += len(s)


@pytest.mark.parametrize("name,kw", [
    ("Trypsin", dict(missed_cleavages=2)),
    ("Arg-C", dict(min_length=4, max_length=40, missed_cleavages=1,
                   chemistry=Chemistry(propionylation="before_and_after"))),
    ("Glu-C", dict(min_length=7, max_length=70)),
    ("Lys-N", dict(min_length=1, max_length=300, missed_cleavages=3)),
])  # fmt: skip
def test_index_flags_equal_substring_mapping(histones, name, kw):
    idx = SharedPrefixIndex.build(histones)
    d = digest(histones, _settings(name, **kw))
    up, ug = idx.unique_flags(d, d.identifiable)
    mp = peptide_mapping(d, histones)
    assert np.array_equal(up, mp["unique_protein"].to_numpy())
    assert np.array_equal(ug, mp["unique_gene"].to_numpy())


def test_index_save_load_roundtrip(histones, tmp_path):
    idx = SharedPrefixIndex.build(histones)
    idx.save(tmp_path / "i.npz")
    back = SharedPrefixIndex.load(tmp_path / "i.npz", histones)
    assert np.array_equal(back.other_entry, idx.other_entry)
    with pytest.raises(ValueError):
        SharedPrefixIndex.load(tmp_path / "i.npz", histones.subset(range(10)))


# ----- indexed background, including uploads not in the reference ---------------------------
def test_indexed_background_equals_live_mapping_with_uploads(histones):
    ref = histones.subset(range(40))
    ref.name = "ref"
    idx = SharedPrefixIndex.build(ref)
    # uploads: an exact reference entry, a novel protein containing a reference peptide,
    # and an edited copy of a reference entry under a new accession
    h3 = histones.sequence("P68431")
    up = parse_fasta(
        f">{ref.accessions[0]} dup\n{ref.sequences[0]}\n"
        ">sp|NOVEL1|N1_HUMAN Novel OS=Homo sapiens GN=NOV\nGGGGRKSTGGKAPRKQLATKAARGGGG\n"
        f">sp|EDIT1|E1_HUMAN Edited OS=Homo sapiens GN=ED\n{h3[:60]}W{h3[61:]}\n"
    )
    live = uniqueness_background(ref, up)
    fast = Background(live.proteins, index=idx)
    assert fast.novel_entries().accessions == ["NOVEL1", "EDIT1"]
    s = _settings("Trypsin", min_length=4, max_length=40, missed_cleavages=1,
                  chemistry=Chemistry(propionylation="before_and_after"))  # fmt: skip
    for query in (ref, up, live.proteins):
        d = digest(query, s)
        a = live.mapping(d)
        b = fast.mapping(d)
        assert np.array_equal(a["unique_protein"], b["unique_protein"])
        assert np.array_equal(a["unique_gene"], b["unique_gene"])
    # the novel protein makes a reference peptide non-unique
    d = digest(ref, s)
    peps = d.sequences(d.identifiable)
    flags = fast.mapping(d)["unique_protein"].to_numpy()
    assert "KSTGGKAPR" in peps
    assert not any(f for p, f in zip(peps, flags, strict=True) if p == "KSTGGKAPR")


def test_fast_comparison_counts_equal_compare(histones):
    bg = Background(histones)
    base = Analysis(digest(histones, _settings("Trypsin", missed_cleavages=1)), bg)
    for n in ("Asp-N", "Glu-C", "Chymotrypsin"):
        for s in (
            _settings(n),
            _settings("Trypsin").with_rules(get_rule("Trypsin", "psi-ms"), get_rule(n, "psi-ms")),
        ):
            alt = Analysis(digest(histones, s), bg)
            assert comparison_counts(base, alt) == compare(base, alt).counts


# ----- presets ------------------------------------------------------------------------------
def test_presets_match_claude_md():
    p = PRESETS
    assert p["proteome"].baseline.label == "Trypsin (PSI-MS)"
    assert (p["proteome"].baseline.min_length, p["proteome"].baseline.max_length,
            p["proteome"].baseline.missed_cleavages) == (7, 30, 1)  # fmt: skip
    h = p["histone_ptm"].baseline
    assert (h.min_length, h.max_length, h.missed_cleavages) == (4, 40, 1)
    assert h.chemistry.propionylation is Propionylation.BEFORE_AND_AFTER
    m = p["middle_down"].baseline
    assert (m.label, m.min_length, m.max_length, m.missed_cleavages) == ("Glu-C (PSI-MS)", 7, 70, 0)
    c = p["protein_cleaver"].baseline
    assert c.label == "Arg-C (ExPASy)" and p["protein_cleaver"].mode == "expasy"
    assert not c.chemistry.carbamidomethyl and not c.mz_filter and c.missed_cleavages == 0
    for key, n_enz in (("proteome", 16), ("protein_cleaver", 12)):
        assert len(alternatives(p[key])) == 3 * (n_enz - 1)
        assert len(available_enzymes(p[key].mode)) == n_enz


def test_precompute_roundtrip_matches_ranking_table(histones, tmp_path):
    preset = PRESETS["protein_cleaver"]
    idx = SharedPrefixIndex.build(histones)
    bg = Background(histones, index=idx)
    result = P.build_preset(preset, histones, bg)
    P.write_preset(result, preset, tmp_path)
    loaded = P.load_preset(tmp_path, preset.key, histones.accessions)
    ours = loaded.ranking_table()
    enzymes = [get_rule(n, "expasy") for n in available_enzymes("expasy")]
    ref = ranking_table(histones, preset.baseline, enzymes, bg)
    assert ours["digest"].tolist() == ref["digest"].tolist()
    for col in ref.columns:
        if col in ("digest", "kind", "enzyme", "flags"):
            assert ours[col].tolist() == ref[col].tolist(), col
        else:
            assert np.allclose(ours[col].astype(float), ref[col].astype(float)), col
    digests = preset_digests(preset, histones)
    for i in (0, 5):
        st = Analysis(digests[i][2], bg).protein_stats
        et = loaded.entry_table(i)
        for f in P.STAT_FIELDS:
            assert np.array_equal(et[f].to_numpy(), st[f].to_numpy()), f
    hist = loaded.info["digests"][0]["length_histogram"]
    assert sum(hist) == digests[0][2].n_total_peptides


# ----- shipped files ------------------------------------------------------------------------
@pytest.mark.skipif(not (OUT / "manifest.json").exists(), reason="run scripts/precompute.py")
def test_shipped_files_within_limits_and_verified():
    import hashlib

    manifest = json.loads((OUT / "manifest.json").read_text())
    assert set(manifest["presets"]) == set(PRESETS)
    for name, rec in manifest["files"].items():
        path = OUT / name
        assert path.stat().st_size == rec["bytes"] < 25 * 1024 * 1024, name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == rec["sha256"], name
    assert "PREDICTED" in manifest["notice"]


@pytest.mark.slow
@pytest.mark.skipif(not (OUT / "manifest.json").exists(), reason="run scripts/precompute.py")
def test_shipped_reference_and_preset_load():
    ref, idx, manifest = P.load_reference(OUT)
    assert len(ref) == manifest["reference"]["entries"]
    assert isinstance(ref, ProteinSet) and "UniProt release" in ref.name
    for key in PRESETS:
        pre = P.load_preset(OUT, key, ref.accessions)
        rt = pre.ranking_table()
        assert rt.loc[0, "kind"] == "baseline" and len(rt) == manifest["presets"][key]["n_digests"]
    # spot-check: proteome baseline per-entry stats equal a fresh digest
    pre = P.load_preset(OUT, "proteome", ref.accessions)
    an = Analysis(digest(ref, PRESETS["proteome"].baseline), Background(ref, index=idx))
    et = pre.entry_table(0)
    for f in P.STAT_FIELDS:
        assert np.array_equal(et[f].to_numpy(), an.protein_stats[f].to_numpy()), f
