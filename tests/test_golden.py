"""Golden tests against Protein Cleaver exports (CLAUDE.md section 10).

Settings: ExPASy mode, 7-30 aa, 0 missed cleavages, no Met removal, no fixed mods,
average mass. Peptide sets (accession, start, end, sequence) must match exactly and
average masses within 0.01 Da.
"""

import numpy as np
import pandas as pd
import pytest

from proteasescope.chemistry import Chemistry
from proteasescope.digest import DigestSettings, digest
from proteasescope.rules import get_rule

from .conftest import GOLDEN_DIR

GOLDEN = {
    # file, enzyme (ExPASy mode), expected identifiable-peptide count
    "Arg-C": ("All_Histone_ArgC_BaselineIdentifiable_peptides_2026-11-13.tsv", 307),
    "Lys-C": ("LysC_Identifiable_peptides_2026-11-13.tsv", 304),
    "Lys-N": ("LysN_Identifiable_peptides_2026-11-13.tsv", 301),
    "Asp-N": ("AspN_Identifiable_peptides_2026-11-13.tsv", 119),
    "Chymotrypsin high specificity": (
        "All_Histones_HSChemotrypsin_Identifiable_peptides_2026-11-13.tsv",
        156,
    ),
    "Chymotrypsin low specificity": ("Chemotrypsin_LS_Identifiable_peptides_2026-11-13.tsv", 345),
    "Proteinase K": ("ProtinaseK_Identifiable_peptides_2026-11-13.tsv", 107),
}

PC_CHEMISTRY = Chemistry(carbamidomethyl=False, propionylation="none", met_removal="never")


def _ours(histones, enzyme):
    settings = DigestSettings(
        (get_rule(enzyme, "expasy"),),
        missed_cleavages=0,
        min_length=7,
        max_length=30,
        chemistry=PC_CHEMISTRY,
    )
    return digest(histones, settings).to_frame()


def _theirs(fname):
    df = pd.read_csv(GOLDEN_DIR / fname, sep="\t")
    return df.rename(
        columns={
            "UniProtID": "accession",
            "Peptide": "peptide",
            "Start": "start",
            "End": "end",
            "Length": "length",
            "Peptide mass [Da]": "avg_mass",
        }
    )


KEY = ["accession", "start", "end", "peptide"]


@pytest.mark.parametrize("enzyme", list(GOLDEN))
def test_protein_cleaver_golden(histones, enzyme):
    fname, expected_n = GOLDEN[enzyme]
    theirs = _theirs(fname)
    ours = _ours(histones, enzyme)
    assert len(theirs) == expected_n
    assert len(ours) == expected_n

    ours_set = set(map(tuple, ours[KEY].to_numpy().tolist()))
    theirs_set = set(map(tuple, theirs[KEY].to_numpy().tolist()))
    assert ours_set - theirs_set == set(), "peptides predicted only by ProteaseScope"
    assert theirs_set - ours_set == set(), "peptides predicted only by Protein Cleaver"

    merged = ours.merge(theirs, on=KEY, suffixes=("", "_pc"), validate="one_to_one")
    assert len(merged) == expected_n
    assert (merged["length"] == merged["length_pc"]).all()
    diff = np.abs(merged["avg_mass"] - merged["avg_mass_pc"])
    assert diff.max() < 0.01, merged.loc[diff.idxmax(), KEY + ["avg_mass", "avg_mass_pc"]]
