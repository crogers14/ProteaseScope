"""Peptide-to-protein mapping and uniqueness (substring default, digest-product option)."""

import numpy as np
import pandas as pd
import pytest

from proteasescope.chemistry import Chemistry
from proteasescope.digest import DigestSettings, digest
from proteasescope.fasta import parse_fasta
from proteasescope.metrics import peptide_mapping, peptide_matches
from proteasescope.rules import get_rule

NO_MODS = Chemistry(carbamidomethyl=False, met_removal="never")


def _trypsin(**kw):
    return DigestSettings((get_rule("Trypsin", "psi-ms"),), chemistry=NO_MODS, **kw)


def _table(d, background, method):
    return pd.concat([d.to_frame(), peptide_mapping(d, background, method=method)], axis=1)


def test_il_uniqueness():
    ps = parse_fasta(
        ">sp|A1|A_HUMAN Prot A OS=Homo sapiens GN=GA\nPEPTIDEKSAMPLEAAK\n"
        ">sp|A1-2|A_HUMAN Isoform 2 of Prot A OS=Homo sapiens GN=GA\nSAMPLEAAKGGGK\n"
        ">sp|B1|B_HUMAN Prot B OS=Homo sapiens GN=GB\nPEPTLDEKWWWWWWK\n"
        ">sp|C1|C_HUMAN Prot C OS=Homo sapiens GN=GC\nUNIQWEKMMMMMMMR\n"
    )
    d = digest(ps, _trypsin(min_length=6))
    for method in ("substring", "digest"):
        res = _table(d, d, method).set_index(["accession", "peptide"])
        # PEPTIDEK (A1) and PEPTLDEK (B1) differ only by I/L -> shared
        assert res.loc[("A1", "PEPTIDEK"), "n_proteins"] == 2
        assert not res.loc[("A1", "PEPTIDEK"), "unique_protein"]
        assert not res.loc[("B1", "PEPTLDEK"), "unique_gene"]
        # SAMPLEAAK in A1 and its isoform: not protein-unique, but gene-unique
        assert not res.loc[("A1", "SAMPLEAAK"), "unique_protein"]
        assert res.loc[("A1", "SAMPLEAAK"), "unique_gene"]
        assert res.loc[("C1", "UNIQWEK"), "unique_protein"]


def test_substring_and_digest_definitions_disagree():
    """SAMPLER is a tryptic product of A only. B contains SAMPIER (I/L-identical) after a G,
    so trypsin never releases it from B: digest mapping says unique, substring says shared."""
    ps = parse_fasta(
        ">sp|A1|A_HUMAN Prot A OS=Homo sapiens GN=GA\nGGGGKSAMPLERWWWWWWWK\n"
        ">sp|B1|B_HUMAN Prot B OS=Homo sapiens GN=GB\nHHHHGSAMPIERYYYYYYYK\n"
    )
    d = digest(ps, _trypsin())
    sub = _table(d, d, "substring").set_index(["accession", "peptide"])
    dig = _table(d, d, "digest").set_index(["accession", "peptide"])
    assert dig.loc[("A1", "SAMPLER"), "n_proteins"] == 1
    assert dig.loc[("A1", "SAMPLER"), "unique_protein"]
    assert sub.loc[("A1", "SAMPLER"), "n_proteins"] == 2
    assert not sub.loc[("A1", "SAMPLER"), "unique_protein"]
    assert not sub.loc[("A1", "SAMPLER"), "unique_gene"]
    # B's own tryptic product is unique under both definitions
    assert sub.loc[("B1", "HHHHGSAMPIER"), "unique_protein"]
    assert dig.loc[("B1", "HHHHGSAMPIER"), "unique_protein"]
    # default method is substring
    assert peptide_mapping(d, ps)["n_proteins"].tolist() == sub["n_proteins"].tolist()


def test_substring_mapping_ignores_cross_entry_junctions():
    # "KAAAAAAAAAAR" exists only across the a|b junction of the concatenated sequences;
    # its first 8 residues (the k-mer anchor) lie inside a, so the length check matters.
    ps = parse_fasta(">a\nGGGGGGGGGKAAAAAAA\n>b\nAAARGGGGGGG\n>c\nKAAAAAAAAAAR\n")
    d = digest(ps, DigestSettings((get_rule("Lys-N", "psi-ms"),), chemistry=NO_MODS))
    res = _table(d, ps, "substring").set_index(["accession", "peptide"])
    assert res.loc[("c", "KAAAAAAAAAAR"), "n_proteins"] == 1


@pytest.mark.parametrize("min_length,enzyme", [(7, "Arg-C"), (4, "Glu-C"), (4, "Lys-N")])
def test_substring_mapping_matches_brute_force(histones, min_length, enzyme):
    s = DigestSettings(
        (get_rule(enzyme, "expasy"),),
        missed_cleavages=1,
        min_length=min_length,
        max_length=40,
        chemistry=Chemistry(met_removal="MetAP rule"),
    )
    d = digest(histones, s)
    peps = d.sequences(d.identifiable)
    folded = [q.replace("I", "L") for q in histones.sequences]
    expected = {
        (i, j) for i, p in enumerate(peps) for j, q in enumerate(folded) if p.replace("I", "L") in q
    }
    pairs = peptide_matches(d, histones)
    assert set(zip(pairs["query"], pairs["protein"], strict=True)) == expected
    counts = peptide_mapping(d, histones)["n_proteins"].to_numpy()
    assert (counts >= 1).all()
    assert counts.tolist() == np.bincount([i for i, _ in expected], minlength=len(peps)).tolist()


def test_digest_mapping_matches_brute_force(histones):
    s = DigestSettings((get_rule("Arg-C", "expasy"),), missed_cleavages=1, chemistry=NO_MODS)
    d = digest(histones, s)
    all_products = {}
    for prot, pep in zip(d.protein.tolist(), d.sequences(), strict=True):
        all_products.setdefault(pep.replace("I", "L"), set()).add(prot)
    peps = d.sequences(d.identifiable)
    counts = peptide_mapping(d, d, method="digest")["n_proteins"].tolist()
    assert counts == [len(all_products[p.replace("I", "L")]) for p in peps]
    # every digest-product match is also a substring match
    sub = peptide_mapping(d, d)["n_proteins"].to_numpy()
    assert (sub >= np.array(counts)).all()


def test_mapping_argument_checks(histones):
    d = digest(histones, _trypsin())
    with pytest.raises(ValueError):
        peptide_mapping(d, d, method="exact")
    with pytest.raises(TypeError):
        peptide_mapping(d, histones, method="digest")
