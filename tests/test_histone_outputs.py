"""Histone outputs (CLAUDE.md section 8), checked against brute-force recomputation."""

import shutil
import time

import pandas as pd
import pytest

from proteasescope import histones as H
from proteasescope.chemistry import Chemistry
from proteasescope.digest import DigestSettings, digest, parallel
from proteasescope.fasta import parse_fasta
from proteasescope.metrics import Analysis, compare, peptide_map
from proteasescope.rules import get_rule

from .conftest import ROOT

HIST = ROOT / "data" / "histones"
PTM = DigestSettings(
    (get_rule("Trypsin", "psi-ms"),),
    missed_cleavages=1,
    min_length=4,
    max_length=40,
    chemistry=Chemistry(propionylation="before_and_after"),
    name="Histone PTM baseline",
)


@pytest.fixture(scope="module")
def inputs():
    return H.load_histone_inputs(HIST)


@pytest.fixture(scope="module")
def analyses(inputs):
    ps = inputs[0]
    base = digest(ps, PTM)

    def alt(name):
        return digest(ps, PTM.with_rules(get_rule(name, "psi-ms")))

    return [
        Analysis(base),
        Analysis(digest(ps, PTM.with_rules(*PTM.rules, get_rule("Glu-C", "psi-ms")))),
        Analysis(alt("Asp-N")),
        Analysis(parallel(base, alt("Chymotrypsin"))),
    ]


def _peps(an):
    """Brute force: {accession: [(start, end, sequence)]} UniProt numbering, identifiable."""
    out = {}
    for r in an.peptide_table().itertuples():
        out.setdefault(r.accession, []).append((r.start, r.end, r.peptide))
    return out


# ----- inputs ---------------------------------------------------------------------------
def test_inputs_require_reviewed_tables(tmp_path):
    for f in HIST.iterdir():
        if f.is_file():
            shutil.copy(f, tmp_path / f.name)
    ann = tmp_path / "histone_annotation.tsv"
    lines = ann.read_text().splitlines()
    lines[0] = "# status: DRAFT (test)"
    ann.write_text("\n".join(lines) + "\n")
    with pytest.raises(PermissionError):
        H.load_histone_inputs(tmp_path)


def test_inputs(inputs):
    ps, ann, marks = inputs
    assert len(ps) == len(ann) == 55
    assert set(marks.accession) <= set(ps.accessions)


# ----- 1. mark visibility -----------------------------------------------------------------
def test_mark_visibility_brute_force(inputs, analyses):
    ps, ann, marks = inputs
    vis = H.mark_visibility(analyses, marks, ann)
    defaults = set(H.default_entries(ann))
    assert set(vis.accession) <= defaults
    n_sites = marks[marks.accession.isin(defaults)][["accession", "position_uniprot"]]
    assert len(vis) == len(analyses) * len(n_sites.drop_duplicates())
    covered = {}
    for an in analyses:
        peps = _peps(an)
        sub = vis[vis.digest == an.label]
        for r in sub.itertuples():
            cov = [
                (s, e, p) for s, e, p in peps.get(r.accession, []) if s <= r.position_uniprot <= e
            ]
            assert r.covered == bool(cov)
            assert r.n_K == ";".join(str(p.count("K")) for _, _, p in cov)
            covered[(an.label, r.accession, r.position_uniprot)] = bool(cov)
    base = analyses[0].label
    for r in vis[vis.digest != base].itertuples():
        b = covered[(base, r.accession, r.position_uniprot)]
        expected = {(1, 1): "seen in both", (0, 1): "gained", (1, 0): "lost",
                    (0, 0): "unseen in both"}[(int(b), int(r.covered))]  # fmt: skip
        assert r.status == expected
    # parallel digests never lose a site
    assert not ((vis.digest == analyses[3].label) & (vis.status == "lost")).any()
    assert {"gained", "lost", "seen in both"} <= set(vis.status)


def test_h3k4_visible_on_tkqtar(inputs, analyses):
    ps, ann, marks = inputs
    vis = H.mark_visibility(analyses[:1], marks, ann, entries=["P68431"], modifications=["me3"])
    k4 = vis[vis.site == "K4"].iloc[0]
    assert k4.position_uniprot == 5 and k4.covered
    assert "TKQTAR(3-8)" in k4.covering_peptides.split("; ")
    assert "me3" in k4.modifications
    grid = H.mark_grid(vis)
    assert grid.loc[("H3.1", "K4"), PTM.label] == "covered"


# ----- 2. isoform discrimination ------------------------------------------------------------
def test_isoform_discrimination_brute_force(inputs, analyses):
    ps, ann, _ = inputs
    tab = H.isoform_discrimination(analyses, ann)
    folded = {a: s.replace("I", "L") for a, s in zip(ps.accessions, ps.sequences, strict=True)}
    for an in analyses:
        peps = _peps(an)
        for r in tab[tab.digest == an.label].itertuples():
            da = [p for _, _, p in peps.get(r.accession_a, [])
                  if p.replace("I", "L") not in folded[r.accession_b]]  # fmt: skip
            db = [p for _, _, p in peps.get(r.accession_b, [])
                  if p.replace("I", "L") not in folded[r.accession_a]]  # fmt: skip
            assert (r.n_discriminating_a, r.n_discriminating_b) == (len(da), len(db))
            expected = "resolved" if da and db else "one-way" if da or db else "unresolved"
            assert r.status == expected
    fam = ann.set_index("accession")["family"]
    assert (tab.accession_a.map(fam) == tab.accession_b.map(fam)).all()


def test_h31_vs_h32_resolved_by_c96(inputs, analyses):
    ps, ann, _ = inputs
    tab = H.isoform_discrimination(analyses[:1], ann)
    r = tab[(tab.entry_a == "H3.1") & (tab.entry_b == "H3.2")].iloc[0]
    assert r.status == "resolved"
    assert r.discriminating_peptides_a == "FQSSAVMALQEACEAYLVGLFEDTNLCAIHAKR(84-116)"
    assert r.discriminating_peptides_b == "FQSSAVMALQEASEAYLVGLFEDTNLCAIHAKR(84-116)"
    mat = H.discrimination_matrix(tab, PTM.label, "H3")
    assert mat.loc["H3.1", "H3.2"] == 1 and pd.isna(mat.loc["H3.1", "H3.1"])


def test_identical_and_one_way_statuses():
    ps = parse_fasta(
        ">a\nGGGGRAAAAAAR\n>b\nGGGGRAAAAAAR\n>c\nGGGGRAAAAAARCCCCCCR\n>d\nGGGGRALAAAAR\n>e\nGGGGRAIAAAAR\n"
    )
    ann = pd.DataFrame({
        "accession": list("abcde"), "display_name": list("ABCDE"), "family": "H3",
        "category": "canonical", "cooccurrence_default": False,
    })  # fmt: skip
    s = DigestSettings(
        (get_rule("Arg-C", "psi-ms"),),
        min_length=4,
        chemistry=Chemistry(carbamidomethyl=False, met_removal="never"),
    )
    tab = H.isoform_discrimination([Analysis(digest(ps, s))], ann).set_index(["entry_a", "entry_b"])
    assert tab.loc[("A", "B"), "status"] == "identical sequence"
    # C contains all of A's peptides; A contains none of C's CCCCCCR -> one-way
    assert tab.loc[("A", "C"), "status"] == "one-way"
    assert tab.loc[("A", "C"), "n_discriminating_b"] == 1
    # D (ALAAAAR) and E (AIAAAAR) are identical with I = L
    assert tab.loc[("D", "E"), "status"] == "identical sequence"
    assert tab.loc[("A", "D"), "status"] == "resolved"


# ----- 3. lysine co-occurrence -----------------------------------------------------------------
def _brute_sites_pairs(peps, acc):
    sites, pairs = set(), set()
    for s, _e, p in peps.get(acc, []):
        ks = [s + i for i, c in enumerate(p) if c == "K"]
        sites |= set(ks)
        pairs |= {(a, b) for i, a in enumerate(ks) for b in ks[i + 1 :]}
    return sites, pairs


def test_lysine_cooccurrence_brute_force(inputs, analyses):
    ps, ann, _ = inputs
    co = H.lysine_cooccurrence(analyses, ann)
    entries = H.cooccurrence_entries(ann)
    assert set(co["sites"].accession) == set(entries)
    assert "P84243" in entries  # H3.3 kept as explicit exception
    letters = dict(zip(co["letters"].digest, co["letters"].letter, strict=True))
    assert list(letters.values()) == ["A", "B", "C"]
    base_peps = _peps(analyses[0])
    seen = set()
    for acc in entries:
        bs, bp = _brute_sites_pairs(base_peps, acc)
        for an in analyses[1:]:
            a_sites, a_pairs = _brute_sites_pairs(_peps(an), acc)
            got = co["pairs"][(co["pairs"].accession == acc) & (co["pairs"].digest == an.label)]
            got = {(r.position_uniprot_i, r.position_uniprot_j): r.status for r in got.itertuples()}
            expected = {}
            for pr in bp | a_pairs:
                if pr in a_pairs:
                    expected[pr] = "kept" if pr in bp else "new"
                else:
                    expected[pr] = "split" if pr[0] in a_sites and pr[1] in a_sites else "lost"
            assert got == expected, (acc, an.label)
            seen |= set(expected.values())
            sites = co["sites"][(co["sites"].accession == acc) & (co["sites"].digest == an.label)]
            assert set(sites.loc[sites.covered, "position_uniprot"]) == a_sites
    assert {"kept", "new", "split", "lost"} <= seen
    m = co["matrix"]
    row = m[(m.display_name == "H3.1") & (m.K_i == 9) & (m.K_j == 14)].iloc[0]
    assert row.baseline_pair  # KSTGGKAPR (9-17)


def test_split_pair_definition():
    # Baseline (Arg-C) keeps K2 and K6 together; adding Glu-C cuts between them -> split.
    ps = parse_fasta(">x\nRAKAEAKAAR\n")
    chem = Chemistry(carbamidomethyl=False, met_removal="never", propionylation="before_and_after")
    s = DigestSettings((get_rule("Arg-C", "psi-ms"),), min_length=2, chemistry=chem)
    base = Analysis(digest(ps, s))
    alt = Analysis(digest(ps, s.with_rules(*s.rules, get_rule("Glu-C", "psi-ms"))))
    ann = pd.DataFrame({"accession": ["x"], "display_name": ["X"], "family": ["H3"],
                        "category": ["canonical"], "cooccurrence_default": [True]})  # fmt: skip
    co = H.lysine_cooccurrence([base, alt], ann)
    assert co["pairs"].status.tolist() == ["split"]
    assert co["matrix"].split_in.tolist() == ["A"]


# ----- 4. peptide maps ---------------------------------------------------------------------------
def test_tail_peptide_map(inputs, analyses):
    ps, ann, _ = inputs
    pm = H.tail_peptide_map(analyses, ann, "P68431")
    assert (pm["peptides"].start_mature <= 45).all()
    assert set(pm["lysines"].position_mature) == {4, 9, 14, 18, 23, 27, 36, 37}
    full = peptide_map(analyses, "P68431")
    cov = full["coverage"].set_index("digest")["coverage_pct"]
    for an in analyses:
        st = an.protein_stats.set_index("accession")
        assert cov[an.label] == pytest.approx(st.loc["P68431", "max_coverage_pct"])
    for an in analyses[1:]:
        res = compare(analyses[0], an).residues.set_index("accession")
        strips = full["strips"][full["strips"].digest == an.label]
        for status in ("gained", "lost"):
            n = int((strips[strips.status == status].end_mature
                     - strips[strips.status == status].start_mature + 1).sum())  # fmt: skip
            exp = int(res.loc["P68431", f"residues_{status}"]) if "P68431" in res.index else 0
            assert n == exp


# ----- performance ------------------------------------------------------------------------
def test_live_histone_digest_and_outputs_under_1s(inputs):
    ps, ann, marks = inputs
    t0 = time.perf_counter()
    base = digest(ps, PTM)
    an = [Analysis(base)] + [
        Analysis(digest(ps, PTM.with_rules(*PTM.rules, get_rule(n, "psi-ms"))))
        for n in ("Asp-N", "Glu-C", "Chymotrypsin")
    ]
    H.mark_visibility(an, marks, ann)
    H.isoform_discrimination(an, ann)
    H.lysine_cooccurrence(an, ann)
    H.tail_peptide_map(an, ann, "P68431")
    assert time.perf_counter() - t0 < 1.0
