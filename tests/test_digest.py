"""Unit tests for chemistry, FASTA parsing and the digestion engine (CLAUDE.md section 10)."""

import numpy as np
import pandas as pd
import pytest

from proteasescope.chemistry import (
    PROPIONYL_MONO,
    Chemistry,
    MetRemoval,
    Propionylation,
    met_removed,
    peptide_mass,
)
from proteasescope.digest import (
    DigestSettings,
    coverage,
    cut_sites,
    digest,
    parallel,
)
from proteasescope.fasta import ProteinSet, parse_fasta, parse_header
from proteasescope.rules import get_rule

NO_MODS = Chemistry(carbamidomethyl=False, met_removal="never")


def _settings(*names, mode="psi-ms", chem=NO_MODS, **kw):
    return DigestSettings(tuple(get_rule(n, mode) for n in names), chemistry=chem, **kw)


def _sites_by_protein(d):
    """Set of (entry, local cut position) from peptide boundaries of a 0-mc digest."""
    sel = d.missed == 0
    starts = set(zip(d.protein[sel].tolist(), d.start0[sel].tolist(), strict=True))
    mature_starts = set(enumerate(d.mature_offset.tolist()))
    return starts - mature_starts


# ----- FASTA ------------------------------------------------------------------------
def test_parse_uniprot_header_and_isoform():
    h = parse_header("sp|P68431|H31_HUMAN Histone H3.1 OS=Homo sapiens OX=9606 GN=H3C1 PE=1 SV=2")
    assert (h["accession"], h["gene"], h["protein_name"], h["pe"], h["sv"]) == (
        "P68431",
        "H3C1",
        "Histone H3.1",
        1,
        2,
    )
    assert not h["is_isoform"]
    iso = parse_header(
        "sp|Q71UI9-2|H2AV_HUMAN Isoform 2 of Histone H2A.V OS=Homo sapiens OX=9606 GN=H2AZ2"
    )
    assert iso["is_isoform"] and iso["canonical_accession"] == "Q71UI9" and iso["isoform"] == "2"
    assert iso["protein_name"] == "Isoform 2 of Histone H2A.V" and iso["pe"] is None
    other = parse_header("myprot some description")
    assert other["accession"] == "myprot" and other["gene"] is None


def test_parse_fasta_errors():
    with pytest.raises(ValueError, match="Duplicate"):
        parse_fasta(">a\nAAA\n>a\nCCC\n")
    with pytest.raises(ValueError, match="Invalid"):
        parse_fasta(">a\nAA1A\n")
    ps = parse_fasta(">a x\nmk\nAA*\n")
    assert ps.sequences == ["MKAA"]


def test_histone_fasta(histones):
    assert len(histones) == 54
    assert (histones.meta["pe"] == 5).sum() == 3


# ----- trypsin K-P --------------------------------------------------------------------
def test_trypsin_psims_does_not_cut_kp():
    tryp = get_rule("Trypsin", "psi-ms")
    assert tryp.cut_sites("AAKPAAKAARPAARA") == [7, 14]
    assert get_rule("Trypsin/P", "psi-ms").cut_sites("AAKPAAKAA") == [3, 7]


# ----- propionylation -------------------------------------------------------------------
def test_protein_level_propionylation_removes_all_k_sites(histones):
    prop = Chemistry(carbamidomethyl=False, propionylation="before_and_after", met_removal="never")
    for name in ["Trypsin", "Trypsin/P"]:
        d = digest(histones, _settings(name, chem=prop))
        sites = _sites_by_protein(d)
        assert sites and all(histones.sequences[i][p - 1] == "R" for i, p in sites)
    for name in ["Lys-C", "Lys-C/P", "Lys-N"]:
        assert digest(histones, _settings(name, chem=prop)).n_cut_sites == 0
        assert digest(histones, _settings(name)).n_cut_sites > 0
    # Trypsin with blocked K == PSI-MS Arg-C cut sites
    tryp = digest(histones, _settings("Trypsin", chem=prop))
    argc = digest(histones, _settings("Arg-C", chem=prop))
    assert _sites_by_protein(tryp) == _sites_by_protein(argc)
    # blocking never adds sites: ExPASy trypsin exception CR|K stays uncleaved
    exp_tryp = get_rule("Trypsin", "expasy")
    assert cut_sites("ACRKAA", exp_tryp, block_lysine=True).tolist() == []
    assert cut_sites("ACRKAA", exp_tryp, block_lysine=False).tolist() == [4]  # after K only


def test_propionyl_masses():
    pep = "KVLR"
    base = peptide_mass(pep, Chemistry(carbamidomethyl=False))
    both = peptide_mass(pep, Chemistry(carbamidomethyl=False, propionylation="before_and_after"))
    # one K + the peptide N-terminus
    assert both - base == pytest.approx(2 * PROPIONYL_MONO, abs=1e-9)
    assert peptide_mass("C", Chemistry(carbamidomethyl=True)) - peptide_mass("C") == pytest.approx(
        57.021464
    )


@pytest.mark.parametrize("old", ["protein-level", "protein + peptide N-term",
                                 "Propionylate before and after digestion"])  # fmt: skip
def test_old_propionylation_settings_migrate(old):
    assert Propionylation(old) is Propionylation.BEFORE_AND_AFTER
    assert Chemistry(propionylation=old).blocks_lysine
    assert Propionylation("none") is Propionylation.NONE
    assert Propionylation("No derivatization") is Propionylation.NONE
    assert [p.value for p in Propionylation] == ["before_and_after", "after_only", "none"]


def test_propionylate_after_digestion_only_does_not_block_k(histones):
    after = Chemistry(carbamidomethyl=False, propionylation="after_only", met_removal="never")
    none = Chemistry(carbamidomethyl=False, propionylation="none", met_removal="never")
    assert not after.blocks_lysine
    for name in ["Trypsin", "Lys-C", "Lys-N"]:
        a = digest(histones, _settings(name, chem=after))
        n = digest(histones, _settings(name, chem=none))
        assert _sites_by_protein(a) == _sites_by_protein(n)
        assert a.n_cut_sites > 0


def test_propionylate_after_digestion_only_labels_every_k_and_nterm(histones):
    after = Chemistry(carbamidomethyl=False, propionylation="after_only", met_removal="never")
    none = Chemistry(carbamidomethyl=False, propionylation="none", met_removal="never")
    for pep in ["KVLR", "GKGGKGLGKGGAKR", "AAR"]:
        nk = pep.count("K")
        assert peptide_mass(pep, after) - peptide_mass(pep, none) == pytest.approx(
            (nk + 1) * PROPIONYL_MONO, abs=1e-9)
        # no derivatization adds nothing
        assert peptide_mass(pep, none) == pytest.approx(peptide_mass(pep), abs=1e-12)
    # and in a digest table: Lys-C peptides end in a free K that is labeled after digestion
    a = digest(histones, _settings("Lys-C", chem=after))
    n = digest(histones, _settings("Lys-C", chem=none))
    assert len(a) == len(n) > 0
    n_k = np.array([s.count("K") for s in a.sequences()])
    ok = ~np.isnan(n.mono_mass)
    assert np.allclose((a.mono_mass - n.mono_mass)[ok], ((n_k + 1) * PROPIONYL_MONO)[ok],
                       atol=1e-6)  # fmt: skip
    assert (n_k > 0).mean() > 0.9  # all but C-terminal peptides end in K


# ----- Met removal ---------------------------------------------------------------------
@pytest.mark.parametrize(
    "seq,rule,expected",
    [
        ("MARTK", "MetAP rule", True),
        ("MCK", "MetAP rule", True),
        ("MGK", "MetAP rule", True),
        ("MPEP", "MetAP rule", True),
        ("MSGR", "MetAP rule", True),
        ("MTK", "MetAP rule", True),
        ("MVK", "MetAP rule", True),
        ("MKAA", "MetAP rule", False),
        ("MEAA", "MetAP rule", False),
        ("MLAA", "MetAP rule", False),
        ("M", "MetAP rule", False),
        ("SGRG", "MetAP rule", False),
        ("MKAA", "always", True),
        ("SKAA", "always", False),
        ("MAAA", "never", False),
    ],
)
def test_metap_rule(seq, rule, expected):
    assert met_removed(seq, MetRemoval(rule)) is expected


def test_met_removed_peptides_start_at_residue_2():
    ps = parse_fasta(">x\nMAKRGGGGGGGR\n")
    d = digest(ps, _settings("Arg-C", chem=Chemistry(met_removal="MetAP rule"), min_length=1))
    df = d.to_frame()
    assert df["peptide"].tolist() == ["AKR", "GGGGGGGR"]
    assert df["start"].tolist() == [2, 5] and df["start_mature"].tolist() == [1, 4]
    assert df["prev_aa"].tolist() == ["M", "R"] and df["next_aa"].tolist() == ["G", "-"]


# ----- numbering -------------------------------------------------------------------------
def test_h3_k4_mature_equals_uniprot_5(histones):
    histone_preset = Chemistry(propionylation="before_and_after", met_removal="MetAP rule")
    s = _settings("Trypsin", chem=histone_preset, min_length=4, max_length=40, missed_cleavages=1)
    df = digest(histones, s).to_frame()
    h3 = df[(df.accession == "P68431") & (df.missed_cleavages == 0)]
    row = h3[h3.peptide == "TKQTAR"].iloc[0]
    assert (row.start_mature, row.end_mature) == (3, 8)
    assert (row.start, row.end) == (4, 9)
    assert row.K_positions == "5"  # UniProt 5 = mature K4
    assert histones.sequence("P68431")[5 - 1] == "K"
    h4 = df[(df.accession == "P62805") & (df.peptide == "KVLR")].iloc[0]
    assert (h4.start_mature, h4.end_mature) == (20, 23)


# ----- sequential and parallel ------------------------------------------------------------
def test_sequential_is_union_of_cut_sites(histones):
    chem = Chemistry(carbamidomethyl=False, met_removal="MetAP rule")
    argc = digest(histones, _settings("Arg-C", mode="expasy", chem=chem))
    aspn = digest(histones, _settings("Asp-N", mode="expasy", chem=chem))
    seq = digest(histones, _settings("Arg-C", "Asp-N", mode="expasy", chem=chem))
    assert _sites_by_protein(seq) == _sites_by_protein(argc) | _sites_by_protein(aspn)
    assert seq.n_cut_sites == len(_sites_by_protein(argc) | _sites_by_protein(aspn))
    assert seq.settings[0].is_sequential
    # peptides are exactly the spans between consecutive union cut sites
    for i, s in enumerate(histones.sequences):
        off = seq.mature_offset[i]
        sites = sorted({p for j, p in _sites_by_protein(argc) | _sites_by_protein(aspn) if j == i})
        bounds = [off, *sites, len(s)]
        expected = list(zip(bounds[:-1], bounds[1:], strict=True))
        sel = (seq.protein == i) & (seq.missed == 0)
        got = list(zip(seq.start0[sel].tolist(), seq.end0[sel].tolist(), strict=True))
        assert got == expected


def test_missed_cleavages():
    ps = parse_fasta(">x\nAAKBBKCCKDD\n")
    d = digest(ps, _settings("Trypsin", min_length=1, missed_cleavages=2))
    got = sorted(zip(d.sequences(), d.missed.tolist(), strict=True))
    assert got == sorted(
        [
            ("AAK", 0),
            ("BBK", 0),
            ("CCK", 0),
            ("DD", 0),
            ("AAKBBK", 1),
            ("BBKCCK", 1),
            ("CCKDD", 1),
            ("AAKBBKCCK", 2),
            ("BBKCCKDD", 2),
        ]
    )
    assert d.n_total_peptides == 4


def test_parallel_is_union_of_peptides(histones):
    chem = Chemistry(carbamidomethyl=False, met_removal="MetAP rule")
    a = digest(histones, _settings("Arg-C", mode="expasy", chem=chem, missed_cleavages=1))
    b = digest(histones, _settings("Asp-N", mode="expasy", chem=chem, missed_cleavages=1))
    p = parallel(a, b)

    def peps(d, mask=None):
        mask = np.ones(len(d), bool) if mask is None else mask
        return set(
            zip(
                d.protein[mask].tolist(),
                d.start0[mask].tolist(),
                d.end0[mask].tolist(),
                strict=True,
            )
        )

    assert peps(p) == peps(a) | peps(b)
    assert len(p) == len(peps(a) | peps(b))
    assert peps(p, p.identifiable) == peps(a, a.identifiable) | peps(b, b.identifiable)
    both = peps(a) & peps(b)
    assert ((p.origin == 3).sum()) == len(both)
    # coverage of the parallel digest = union of covered residues
    ca, cb, cp = coverage(a), coverage(b), coverage(p)
    assert (cp.covered >= np.maximum(ca.covered, cb.covered)).all()
    assert (cp.covered <= ca.covered + cb.covered).all()


def test_masses_and_mz_filter():
    ps = ProteinSet(pd.DataFrame([parse_header("x")]), ["GGGGGGGRCAAAAAAAAR"])
    d = digest(
        ps,
        _settings(
            "Arg-C", chem=Chemistry(carbamidomethyl=True, met_removal="never"), mz_filter=True
        ),
    )
    df = d.to_frame(identifiable_only=False)
    assert df.mono_mass.tolist() == pytest.approx(
        [peptide_mass("GGGGGGGR"), peptide_mass("CAAAAAAAAR") + 57.021464], abs=1e-6
    )
    assert df.mz_2.iloc[0] == pytest.approx((df.mono_mass.iloc[0] + 2 * 1.007276467) / 2)
    # GGGGGGGR (~573 Da): 2+ = 278.6, 3+, 4+ smaller -> outside 350-1400; not identifiable
    assert d.identifiable.tolist() == [False, True]
