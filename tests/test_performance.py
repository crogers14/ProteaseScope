"""Performance: trypsin, <=2 missed cleavages, whole reviewed human proteome + isoforms
in under 15 s natively (CLAUDE.md section 10).

Until Phase 2 downloads the real proteome, a synthetic proteome of the same size
(42,500 entries, ~24 M residues, human amino-acid composition) is used.
"""

import time

import numpy as np
import pandas as pd
import pytest

from proteasescope.chemistry import Chemistry
from proteasescope.digest import DigestSettings, coverage, digest
from proteasescope.fasta import ProteinSet, parse_header, read_fasta
from proteasescope.rules import get_rule

from .conftest import ROOT

REAL_PROTEOME = ROOT / "data" / "reference" / "uniprot_human_reviewed_isoforms.fasta"

# Approximate UniProtKB/Swiss-Prot human composition (%)
HUMAN_AA = {
    "A": 7.0,
    "R": 5.6,
    "N": 3.6,
    "D": 4.7,
    "C": 2.3,
    "Q": 4.8,
    "E": 7.1,
    "G": 6.6,
    "H": 2.6,
    "I": 4.3,
    "L": 10.0,
    "K": 5.7,
    "M": 2.1,
    "F": 3.7,
    "P": 6.3,
    "S": 8.3,
    "T": 5.4,
    "W": 1.2,
    "Y": 2.7,
    "V": 6.0,
}


def _synthetic_proteome(n=42_500, mean_len=565, seed=0):
    rng = np.random.default_rng(seed)
    letters = np.frombuffer("".join(HUMAN_AA).encode(), dtype=np.uint8)
    p = np.array(list(HUMAN_AA.values()))
    lengths = np.clip(rng.lognormal(np.log(mean_len) - 0.32, 0.8, n).astype(int), 30, 35_000)
    codes = rng.choice(letters, size=int(lengths.sum()), p=p / p.sum())
    text = codes.tobytes().decode()
    cuts = np.concatenate(([0], np.cumsum(lengths)))
    seqs = ["M" + text[a + 1 : b] for a, b in zip(cuts[:-1], cuts[1:], strict=True)]
    meta = pd.DataFrame([parse_header(f"sp|S{i}|SYN_HUMAN syn GN=G{i // 2}") for i in range(n)])
    return ProteinSet(meta, seqs)


@pytest.mark.slow
def test_proteome_trypsin_speed():
    proteins = read_fasta(REAL_PROTEOME) if REAL_PROTEOME.exists() else _synthetic_proteome()
    settings = DigestSettings(
        (get_rule("Trypsin", "psi-ms"),), missed_cleavages=2, chemistry=Chemistry()
    )
    t0 = time.perf_counter()
    d = digest(proteins, settings)
    _ = d.identifiable, d.mono_mass
    cov = coverage(d)
    elapsed = time.perf_counter() - t0
    print(
        f"\n{len(proteins)} entries, {int(proteins.lengths.sum()):,} residues, "
        f"{len(d):,} peptides ({int(d.identifiable.sum()):,} identifiable), "
        f"median coverage {cov.coverage_pct.median():.1f}% in {elapsed:.2f} s"
    )
    assert elapsed < 15


@pytest.mark.slow
def test_proteome_substring_uniqueness_speed():
    """Not a CLAUDE.md target: guards the substring mapping against regressions."""
    from proteasescope.metrics import peptide_mapping

    proteins = read_fasta(REAL_PROTEOME) if REAL_PROTEOME.exists() else _synthetic_proteome()
    settings = DigestSettings(
        (get_rule("Trypsin", "psi-ms"),), missed_cleavages=2, chemistry=Chemistry()
    )
    d = digest(proteins, settings)
    t0 = time.perf_counter()
    m = peptide_mapping(d, proteins)
    elapsed = time.perf_counter() - t0
    print(
        f"\nsubstring mapping of {len(m):,} identifiable peptides "
        f"(index build included): {elapsed:.2f} s"
    )
    assert (m["n_proteins"] >= 1).all()
    assert elapsed < 20
