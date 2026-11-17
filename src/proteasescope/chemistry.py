"""Masses, fixed modifications, propionylation and initiator-Met removal.

Units: all masses are neutral masses in daltons (Da); m/z values are for
[M + zH]^z+ ions using the proton mass. Numbering: UniProt (residue 1 = initiator
Met); "mature" numbering = UniProt - 1 when the initiator Met is removed.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

PROTON = 1.007276466812  # Da (CODATA 2018)
WATER_MONO = 18.0105646863  # H2O, monoisotopic
#: Average water mass of the ExPASy table below (matches Protein Cleaver output).
WATER_AVG = 18.01524

# Monoisotopic residue masses (Da), from elemental compositions
# (C 12, H 1.00782503207, N 14.0030740048, O 15.99491461956, S 31.97207100, Se 79.9165213).
RESIDUE_MONO = {
    "A": 71.037113805,
    "R": 156.101111050,
    "N": 114.042927470,
    "D": 115.026943065,
    "C": 103.009184505,
    "E": 129.042593135,
    "Q": 128.058577540,
    "G": 57.021463735,
    "H": 137.058911875,
    "I": 113.084064015,
    "L": 113.084064015,
    "K": 128.094963050,
    "M": 131.040484645,
    "F": 147.068413945,
    "P": 97.052763875,
    "S": 87.032028435,
    "T": 101.047678505,
    "W": 186.079312980,
    "Y": 163.063328575,
    "V": 99.068413945,
    "U": 150.953633405,
    "O": 237.147726925,
}

# Average residue masses (Da): ExPASy PeptideMass / Compute pI/Mw table. This table
# reproduces Protein Cleaver's "Peptide mass [Da]" exactly (see tests/test_golden.py).
RESIDUE_AVG = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
    "U": 150.0379,
    "O": 237.2982,
}
# B, Z, J, X (ambiguous/unknown) have no defined mass: peptides containing them get NaN.

CARBAMIDOMETHYL_MONO = 57.021464  # C2H3NO on Cys (Unimod 4)
CARBAMIDOMETHYL_AVG = 57.0513
PROPIONYL_MONO = 56.026215  # C3H4O (Unimod 58)
PROPIONYL_AVG = 56.0633

#: Second residues that allow initiator-Met removal by methionine aminopeptidase.
METAP_SECOND_RESIDUES = frozenset("ACGPSTV")


class Propionylation(StrEnum):
    """Derivatization of one digest: one plain-language choice (``CUSTOM_STRATEGY.md`` §2).

    It separates lysine blocking (which changes cleavage) from peptide labeling (which
    changes masses):

    * ``BEFORE_AND_AFTER`` ("Propionylate before and after digestion"): every K is blocked
      during digestion (no K-derived cleavage, for every rule) and carries +56.026215 Da;
      every peptide N-terminus is labeled after digestion (+56.026215 Da). Standard histone
      workflow.
    * ``AFTER_ONLY`` ("Propionylate after digestion only"): lysines are free during digestion
      (K cleavage unchanged); afterwards every free K and every peptide N-terminus is labeled
      (+56.026215 Da each).
    * ``NONE`` ("No derivatization"): no labels.

    Endogenous lysine marks are not modeled. Old settings migrate: ``protein-level`` and
    ``protein + peptide N-term`` -> ``BEFORE_AND_AFTER``; ``none`` -> ``NONE``. The display
    labels are accepted too.
    """

    BEFORE_AND_AFTER = "before_and_after"
    AFTER_ONLY = "after_only"
    NONE = "none"
    # Old names (aliases of BEFORE_AND_AFTER), kept so older code and settings still load.
    PROTEIN = "before_and_after"
    PROTEIN_PEPTIDE_NTERM = "before_and_after"

    @classmethod
    def _missing_(cls, value: object) -> Propionylation | None:
        if isinstance(value, str):
            key = value.strip().lower()
            if key in _PROPIONYLATION_ALIASES:
                return cls(_PROPIONYLATION_ALIASES[key])
        return None

    @property
    def label(self) -> str:
        """Plain-language choice shown to the user."""
        return PROPIONYLATION_LABELS[self.value]

    @property
    def blocks_lysine(self) -> bool:
        """True if lysines are blocked during digestion (K-derived cut sites removed)."""
        return self is Propionylation.BEFORE_AND_AFTER

    @property
    def labels_lysine(self) -> bool:
        """True if every K carries a propionyl in the measured peptide."""
        return self is not Propionylation.NONE

    @property
    def labels_nterm(self) -> bool:
        """True if every peptide N-terminus carries a propionyl."""
        return self is not Propionylation.NONE


#: Display label of each derivatization choice (by value).
PROPIONYLATION_LABELS = {
    "before_and_after": "Propionylate before and after digestion",
    "after_only": "Propionylate after digestion only",
    "none": "No derivatization",
}

#: Old setting names and display labels (lower case) -> value.
_PROPIONYLATION_ALIASES = {
    "protein-level": "before_and_after",
    "protein + peptide n-term": "before_and_after",
    "protein": "before_and_after",
    "propionylated": "before_and_after",
    "unmodified": "none",
    "after only": "after_only",
    "after-only": "after_only",
    **{label.lower(): value for value, label in PROPIONYLATION_LABELS.items()},
}


class MetRemoval(StrEnum):
    """Initiator-Met removal: never, always (when residue 1 is M), or by the MetAP rule."""

    NEVER = "never"
    ALWAYS = "always"
    METAP = "MetAP rule"


@dataclass(frozen=True)
class Chemistry:
    """Chemistry settings for a digest.

    Attributes:
        carbamidomethyl: Fixed carbamidomethylation of C (+57.021464 Da).
        propionylation: See :class:`Propionylation`.
        met_removal: See :class:`MetRemoval`.
    """

    carbamidomethyl: bool = True
    propionylation: Propionylation = Propionylation.NONE
    met_removal: MetRemoval = MetRemoval.METAP

    def __post_init__(self) -> None:
        object.__setattr__(self, "propionylation", Propionylation(self.propionylation))
        object.__setattr__(self, "met_removal", MetRemoval(self.met_removal))

    @property
    def blocks_lysine(self) -> bool:
        """True if lysines are blocked during digestion (no K-derived cleavage)."""
        return self.propionylation.blocks_lysine

    def describe(self) -> dict[str, str]:
        """Settings as display strings."""
        return {
            "Fixed modifications": "Carbamidomethyl (C)" if self.carbamidomethyl else "none",
            "Propionylation": self.propionylation.label,
            "Initiator Met removal": str(self.met_removal),
        }


def met_removed(sequence: str, rule: MetRemoval | str = MetRemoval.METAP) -> bool:
    """Whether the initiator Met (UniProt position 1) is removed.

    ``always`` removes residue 1 only if it is M; ``MetAP rule`` removes it when residue 2
    is A, C, G, P, S, T or V.
    """
    rule = MetRemoval(rule)
    if rule is MetRemoval.NEVER or not sequence.startswith("M") or len(sequence) < 2:
        return False
    if rule is MetRemoval.ALWAYS:
        return True
    return sequence[1] in METAP_SECOND_RESIDUES


def residue_mass_table(chem: Chemistry, average: bool = False) -> np.ndarray:
    """Residue mass lookup indexed by ASCII code (length 256), including fixed mods.

    Includes carbamidomethyl C if enabled and +propionyl on K whenever lysines are
    propionylated (before and after digestion, or after digestion only). Codes without a
    defined mass are NaN.
    """
    base = RESIDUE_AVG if average else RESIDUE_MONO
    table = np.full(256, np.nan)
    for aa, m in base.items():
        table[ord(aa)] = m
    if chem.carbamidomethyl:
        table[ord("C")] += CARBAMIDOMETHYL_AVG if average else CARBAMIDOMETHYL_MONO
    if chem.propionylation.labels_lysine:
        table[ord("K")] += PROPIONYL_AVG if average else PROPIONYL_MONO
    return table


def nterm_mass(chem: Chemistry, average: bool = False) -> float:
    """Extra mass on every peptide N-terminus (Da): propionyl if peptides are propionylated
    after digestion, else 0."""
    if chem.propionylation.labels_nterm:
        return PROPIONYL_AVG if average else PROPIONYL_MONO
    return 0.0


def mz(neutral_mass: np.ndarray | float, charge: int) -> np.ndarray | float:
    """m/z of [M + zH]^z+ for a neutral monoisotopic mass (Da)."""
    return (np.asarray(neutral_mass) + charge * PROTON) / charge


def peptide_mass(sequence: str, chem: Chemistry | None = None, average: bool = False) -> float:
    """Neutral mass (Da) of a single peptide string under ``chem`` (default: no mods)."""
    chem = chem or Chemistry(carbamidomethyl=False)
    table = residue_mass_table(chem, average)
    codes = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    water = WATER_AVG if average else WATER_MONO
    return float(table[codes].sum() + water + nterm_mass(chem, average))
