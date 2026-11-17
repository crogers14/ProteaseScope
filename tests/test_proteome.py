"""Outputs on the downloaded reviewed human proteome (skipped until fetch_data.py has run)."""

import random
import time

import pytest

from proteasescope.chemistry import Chemistry
from proteasescope.digest import DigestSettings, digest
from proteasescope.fasta import read_reference_fasta
from proteasescope.metrics import OUTPUT1_COLUMNS, Analysis, Background, uniqueness_background
from proteasescope.rules import get_rule

from .conftest import ROOT

REF = ROOT / "data" / "reference" / "uniprot_human_reviewed_isoforms.fasta"
pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not REF.exists(), reason="run scripts/fetch_data.py first"),
]

PROTEOME_PRESET = dict(missed_cleavages=1, min_length=7, max_length=30, chemistry=Chemistry())


@pytest.fixture(scope="module")
def ref():
    return read_reference_fasta(REF)


@pytest.fixture(scope="module")
def proteome_analysis(ref):
    s = DigestSettings((get_rule("Trypsin", "psi-ms"),), **PROTEOME_PRESET)
    t0 = time.perf_counter()
    an = Analysis(digest(ref, s), Background(ref))
    _ = an.summary, an.unique
    an.elapsed = time.perf_counter() - t0
    return an


def test_reference_provenance_in_name(ref):
    assert "UniProt release" in ref.name and "UP000005640" in ref.name
    assert ref.meta["is_isoform"].sum() > 10_000


def test_histone_fasta_matches_reference(ref, histones):
    seqs = dict(zip(ref.accessions, ref.sequences, strict=True))
    for acc, s in zip(histones.accessions, histones.sequences, strict=True):
        assert seqs.get(acc) == s, acc
    bg = uniqueness_background(ref, histones)
    assert len(bg.proteins) == len(ref)  # all histones already present


def test_proteome_preset_outputs(proteome_analysis):
    an = proteome_analysis
    s = an.summary
    n = len(an.proteins)
    print(f"\n{s} in {an.elapsed:.1f} s")
    assert an.elapsed < 60
    assert s["n_identifiable_proteins"] > 0.99 * n
    assert 0 < s["n_unique_peptides"] < s["n_identifiable_peptides"]
    df = an.peptide_table()
    assert list(df.columns) == OUTPUT1_COLUMNS and len(df) == s["n_identifiable_peptides"]
    assert df["mono_mass"].notna().mean() > 0.999  # only X-containing peptides lack a mass
    hist = an.length_histogram()
    assert hist["n_peptides"].sum() == s["n_total_peptides"]


def test_proteome_uniqueness_spot_check(ref, proteome_analysis):
    """200 random identifiable peptides: entry counts equal a plain substring search."""
    an = proteome_analysis
    df = an.peptide_table()
    counts = an.mapping["n_proteins"].to_numpy()
    folded = [s.replace("I", "L") for s in ref.sequences]
    accs = ref.accessions
    rng = random.Random(7)
    for i in rng.sample(range(len(df)), 200):
        pep = df.at[i, "peptide"].replace("I", "L")
        expected = len({a for a, s in zip(accs, folded, strict=True) if pep in s})
        assert counts[i] == expected, df.at[i, "peptide"]
