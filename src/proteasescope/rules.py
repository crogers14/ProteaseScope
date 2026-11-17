"""Cleavage rules, their sources, and rule-source (compatibility) modes.

Every rule is stored as data. Regexes are written in **zero-width** form (only
lookbehind / lookahead assertions), so a match at string index ``i`` means
"cut between residue ``i - 1`` and residue ``i``" (0-based), i.e. between
UniProt positions ``i`` and ``i + 1``. Matches at the sequence ends are ignored
by the digestion engine.

Two rule sources are supported:

* ``"psi-ms"`` (default): HUPO-PSI PSI-MS controlled vocabulary cleavage-agent
  regexes (``psi-ms.obo``, version recorded in :data:`PSI_MS_VERSION`).
* ``"expasy"``: ExPASy PeptideCutter rules, as used by Protein Cleaver and the
  R package ``cleaver``. Rules vendored from ``pyteomics.parser.expasy_rules``
  (which writes them in residue-consuming form) were rewritten in zero-width
  form; ``tests/test_rules.py`` checks identical cut positions.

Enzymes that exist only in ExPASy (chymotrypsin high/low specificity,
proteinase K) are offered in PSI-MS mode too, but keep their ExPASy source label
so the source printed in every export is always the rule actually used.

Enzymes in neither source (ProAlanase, and the Promega products Arg-C Ultra and
recombinant chymotrypsin used in Turner et al. 2026) are defined from the literature
(``source="Literature"``, ``mode="literature"``) and offered in both modes.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import cached_property

#: PSI-MS CV release the PSI-MS regexes were checked against.
PSI_MS_VERSION = {
    "data_version": "4.2.17",
    "obo_date": "2026-11-11",
    "git_commit": "4d2a2346f7dd20d2da5eadfe93c0a897efe47a6b",
    "url": "https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo",
}

#: pyteomics version the vendored rules were verified against (test-time only).
PYTEOMICS_VERIFIED_VERSION = "5.0.2"

MODES = ("psi-ms", "expasy")

CITE_PSI_MS = (
    "HUPO-PSI PSI-MS controlled vocabulary (psi-ms.obo v{v}); Mayer G et al. (2013) "
    "Database 2013:bat009.".format(v=PSI_MS_VERSION["data_version"])
)
CITE_EXPASY = (
    "ExPASy PeptideCutter cleavage rules; Gasteiger E et al. (2005) The Proteomics "
    "Protocols Handbook, Humana Press, 571-607. As implemented in Protein Cleaver and "
    "the R/Bioconductor package cleaver."
)
CITE_LYSN = "Taouatas N et al. (2008) Nat Methods 5:405-407 (Lys-N cleaves N-terminal to K)."
CITE_PROALANASE = (
    "[to verify] Samodova D et al. (2020) ProAlanase is an effective alternative to trypsin "
    "for proteomics applications and disulfide bond mapping. Mol Cell Proteomics 19:2139-2157. "
    "Promega ProAlanase, Mass Spec Grade, product information."
)

CITE_TURNER_2026 = (
    "Turner et al. (2026) Anal Chem 98:18694, doi:10.1021/acs.analchem.6c01147 (search "
    "settings; citation as provided by the user, to verify)."
)
CITE_ARGC_ULTRA = (
    "Promega Arg-C Ultra, MS grade, product information. " + CITE_TURNER_2026 + " Also "
    "Ryzhaya et al. [to verify: full citation not yet recorded]."
)
CITE_RCHYMO = (
    "Promega rChymoselect (recombinant chymotrypsin; Early Access at time of writing), product "
    "information. " + CITE_TURNER_2026
)
#: Label of the evidence notes from Turner et al. 2026.
TURNER_EVIDENCE_LABEL = "Observed in Turner et al. 2026 (HEK293T histones)"

STANDARD_RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


@dataclass(frozen=True)
class CleavageRule:
    """One cleavage rule.

    Attributes:
        name: Enzyme name shown to the user (e.g. ``"Arg-C"``).
        regex: Zero-width regex; a match at index ``i`` = cut before 0-based residue ``i``.
        source: ``"PSI-MS"``, ``"ExPASy"``, ``"Protein Cleaver"`` or ``"Custom"``.
        mode: Rule-source mode the rule belongs to (``"psi-ms"``, ``"expasy"``, ``"custom"``).
        citation: Literature / database citation.
        accession: PSI-MS CV accession of the cleavage agent, if any.
        notes: Free-text notes (differences between sources, caveats).
        flags: Short warning labels shown with results (e.g. ``"low confidence"``).
        ui_note: Practical note for the enzyme's info tooltip in the app (bench conditions,
            suggested settings); not part of the rule.
        suggested_missed_cleavages: Suggested maximum missed cleavages for this enzyme (shown
            in its tooltip and used by builder presets); None if no suggestion. The sidebar
            setting still applies unless a preset overrides it.
        evidence: One-line observation from the literature for the tooltip, already labeled
            with its source (e.g. "Observed in Turner et al. 2026 (HEK293T histones): ...").
    """

    name: str
    regex: str
    source: str
    mode: str
    citation: str
    accession: str | None = None
    notes: str = ""
    flags: tuple[str, ...] = field(default_factory=tuple)
    ui_note: str = ""
    suggested_missed_cleavages: int | None = None
    evidence: str = ""

    @cached_property
    def pattern(self) -> re.Pattern[str]:
        """Compiled :attr:`regex`."""
        return re.compile(self.regex)

    @property
    def label(self) -> str:
        """Display label including the source, e.g. ``"Arg-C (ExPASy)"``."""
        return f"{self.name} ({self.source})"

    def cut_sites(self, sequence: str) -> list[int]:
        """Return interior cut positions in ``sequence``.

        Positions are 0-based string indices ``i`` with ``0 < i < len(sequence)``;
        a cut at ``i`` separates residue ``i - 1`` from residue ``i``.
        """
        n = len(sequence)
        return [p for p in (m.start() for m in self.pattern.finditer(sequence)) if 0 < p < n]


def _rule(name, regex, source, mode, citation, accession=None, notes="", flags=(), ui_note="",
          suggested_missed_cleavages=None, evidence=""):  # fmt: skip
    return CleavageRule(name, regex, source, mode, citation, accession, notes, tuple(flags),
                        ui_note, suggested_missed_cleavages, evidence)  # fmt: skip


# ---------------------------------------------------------------------------
# PSI-MS rules (default). Regexes copied from psi-ms.obo (with "\!" unescaped).
# ---------------------------------------------------------------------------
_PSI = [
    _rule(
        "Trypsin",
        r"(?<=[KR])(?!P)",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS,
        "MS:1001251",
        "After K or R, not before P.",
    ),
    _rule(
        "Trypsin/P",
        r"(?<=[KR])",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS,
        "MS:1001313",
        "After K or R, including before P.",
    ),
    _rule(
        "Arg-C",
        r"(?<=R)(?!P)",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS,
        "MS:1001303",
        "After R, not before P (ExPASy Arg-C also cuts before P).",
    ),
    _rule(
        "Lys-C",
        r"(?<=K)(?!P)",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS,
        "MS:1001309",
        "After K, not before P (ExPASy Lys-C also cuts before P).",
    ),
    _rule(
        "Lys-C/P",
        r"(?<=K)",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS,
        "MS:1001310",
        "After K, including before P.",
    ),
    _rule(
        "Lys-N",
        r"(?=K)",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS + " " + CITE_LYSN,
        "MS:1003093",
        "Before K. psi-ms.obo v4.2.17 term MS:1003093 states 'cleaves proteins on the amino "
        "side of lysine residues' and annotates its has_regexp as '(?=K)', but the linked "
        "regex ID is MS:1001335, the Lys-C regex '(?<=K)(?!P)' (a CV inconsistency). "
        "We use '(?=K)', matching the definition and Taouatas et al. 2008.",
    ),
    _rule(
        "Asp-N",
        r"(?=[BD])",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS,
        "MS:1001304",
        "Before D or B (ExPASy Asp-N: before D only).",
    ),
    _rule(
        "Asp-N ambic",
        r"(?=[DE])",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS,
        "MS:1001305",
        "Before D or E (Asp-N in ammonium bicarbonate).",
    ),
    _rule(
        "Glu-C",
        r"(?<=[^E]E)",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS,
        "MS:1001917",
        "PSI-MS 'glutamyl endopeptidase' (synonym Glu-C): after E, but not after EE "
        "(nor after an E at the protein N-terminus). ExPASy: after every E.",
    ),
    _rule(
        "Chymotrypsin",
        r"(?<=[FYWL])(?!P)",
        "PSI-MS",
        "psi-ms",
        CITE_PSI_MS,
        "MS:1001306",
        "After F, Y, W or L, not before P. ExPASy splits chymotrypsin into high and low "
        "specificity rules.",
    ),
]

# ---------------------------------------------------------------------------
# ExPASy PeptideCutter rules (compatibility mode), zero-width form.
# pyteomics key in brackets.
# ---------------------------------------------------------------------------
_TRYPSIN_EXPASY_MAIN = r"(?<=[KR])(?=[^P])|(?<=WK)(?=P)|(?<=MR)(?=P)"
_TRYPSIN_EXPASY_EXCEPTIONS = r"(?<=[CD]K)(?=D)|(?<=CK)(?=[HY])|(?<=CR)(?=K)|(?<=RR)(?=[HR])"

_EXPASY = [
    _rule(
        "Trypsin",
        rf"(?!{_TRYPSIN_EXPASY_EXCEPTIONS})(?:{_TRYPSIN_EXPASY_MAIN})",
        "ExPASy",
        "expasy",
        CITE_EXPASY,
        None,
        "ExPASy trypsin [trypsin] with its exceptions [trypsin_exception]: after K/R not "
        "before P; also WK|P and MR|P; no cut at CK|D, DK|D, CK|H, CK|Y, CR|K, RR|H, RR|R.",
    ),
    _rule(
        "Arg-C",
        r"(?<=R)",
        "ExPASy",
        "expasy",
        CITE_EXPASY,
        None,
        "After R, including before P [arg-c].",
    ),
    _rule(
        "Lys-C",
        r"(?<=K)",
        "ExPASy",
        "expasy",
        CITE_EXPASY,
        None,
        "After K, including before P [lysc].",
    ),
    _rule(
        "Lys-N",
        r"(?=K)",
        "Protein Cleaver",
        "expasy",
        CITE_LYSN,
        None,
        "Before K. Not in ExPASy PeptideCutter; Protein Cleaver defines Lys-N as before K.",
    ),
    _rule("Asp-N", r"(?=D)", "ExPASy", "expasy", CITE_EXPASY, None, "Before D [asp-n]."),
    _rule(
        "Glu-C",
        r"(?<=E)",
        "ExPASy",
        "expasy",
        CITE_EXPASY,
        None,
        "After every E [glutamyl endopeptidase].",
    ),
    _rule(
        "Chymotrypsin high specificity",
        r"(?<=[FY])(?=[^P])|(?<=W)(?=[^MP])",
        "ExPASy",
        "expasy",
        CITE_EXPASY,
        None,
        "After F or Y not before P; after W not before M or P [chymotrypsin high specificity].",
    ),
    _rule(
        "Chymotrypsin low specificity",
        r"(?<=[FLY])(?=[^P])|(?<=W)(?=[^MP])|(?<=M)(?=[^PY])|(?<=H)(?=[^DMPW])",
        "ExPASy",
        "expasy",
        CITE_EXPASY,
        None,
        "After F, L or Y not before P; W not before M/P; M not before P/Y; H not before "
        "D/M/P/W [chymotrypsin low specificity].",
    ),
    _rule(
        "Proteinase K",
        r"(?<=[AEFILTVWY])",
        "ExPASy",
        "expasy",
        CITE_EXPASY,
        None,
        "After A, E, F, I, L, T, V, W or Y [proteinase k]. Broad specificity: predictions "
        "are unreliable. Not used from PSI-MS: psi-ms.obo v4.2.17 'peptidase K' "
        "(MS:1003981) gives '(?=[AEFILTVWY])', i.e. cleavage *before* these residues, "
        "which contradicts the P1 specificity in ExPASy/BRENDA.",
        flags=("nonspecific", "low confidence"),
    ),
]

# ---------------------------------------------------------------------------
# Literature-defined rules: in neither psi-ms.obo nor ExPASy; offered in both modes.
# ---------------------------------------------------------------------------
_LITERATURE = [
    _rule(
        "ProAlanase",
        r"(?<=[PA])",
        "Literature",
        "literature",
        CITE_PROALANASE,
        None,
        "After P or A, with no restriction on the following residue (it also cuts before P). "
        "Behavior at P-P bonds is unverified: the rule cuts between every P-P pair, which may "
        "over-predict cleavage in proline-rich stretches. Not in psi-ms.obo v4.2.17 (no "
        "ProAlanase cleavage-agent term; checked 2026-11-15) nor in ExPASy PeptideCutter. "
        "Does not cut at K, so it works on propionylated protein.",
        ui_note="Acidic digestion (around pH 1.5–2); more missed cleavages expected than "
        "trypsin, so consider ≤1–2 missed cleavages.",
    ),
    _rule(
        "Arg-C Ultra",
        r"(?<=R)(?!P)",
        "Literature",
        "literature",
        CITE_ARGC_ULTRA,
        None,
        "After R, not before P; same cut sites as the PSI-MS Arg-C rule, kept under the "
        "product name. Rule from the search settings of Turner et al. 2026; whether 'not "
        "before P' applies to Arg-C Ultra is ambiguous in the paper's text (assumed yes; to "
        "confirm). The paper reports about 84% of peptides fully cleaved on histones. Does not "
        "cut at K, so no derivatization is needed for Arg-C-like peptides.",
        ui_note="Suggested: ≤1 missed cleavage.",
        suggested_missed_cleavages=1,
        evidence=f"{TURNER_EVIDENCE_LABEL}: clean, efficient cleavage; peptides like "
        "propionylated trypsin.",
    ),
    _rule(
        "Recombinant chymotrypsin",
        r"(?<=[FYLM])(?!P)",
        "Literature",
        "literature",
        CITE_RCHYMO,
        None,
        "After F, Y, L or M, not before P; never after W or H (unlike the PSI-MS and ExPASy "
        "chymotrypsin rules). Rule from the search settings of Turner et al. 2026. The paper "
        "reports about 20-25% of peptides fully cleaved, so allow more missed cleavages. Does "
        "not cut at K.",
        ui_note="Suggested: ≤3 missed cleavages.",
        suggested_missed_cleavages=3,
        evidence=f"{TURNER_EVIDENCE_LABEL}: very good H2A, H2A.Z and H1 coverage; poor H3 "
        "coverage (10–20%) despite higher predicted coverage.",
    ),
]

#: Enzymes available only from ExPASy; offered in PSI-MS mode with their ExPASy source.
EXPASY_ONLY = ("Chymotrypsin high specificity", "Chymotrypsin low specificity", "Proteinase K")

#: Literature-defined enzymes; offered in every mode with their ``"Literature"`` source.
LITERATURE_ONLY = tuple(r.name for r in _LITERATURE)

RULES: dict[tuple[str, str], CleavageRule] = {
    (r.mode, r.name): r for r in _PSI + _EXPASY + _LITERATURE
}


def available_enzymes(mode: str = "psi-ms") -> list[str]:
    """Enzyme names selectable in a rule-source mode, in display order."""
    _check_mode(mode)
    names = [name for (m, name) in RULES if m == mode]
    if mode == "psi-ms":
        names += list(EXPASY_ONLY)
    return names + list(LITERATURE_ONLY)


def get_rule(name: str, mode: str = "psi-ms") -> CleavageRule:
    """Return the rule for ``name`` in ``mode``.

    In PSI-MS mode the ExPASy-only enzymes (:data:`EXPASY_ONLY`) return their ExPASy rule,
    whose ``source`` stays ``"ExPASy"``; in every mode the literature-defined enzymes
    (:data:`LITERATURE_ONLY`) return their ``"Literature"`` rule. Raises ``KeyError`` for
    unknown names.
    """
    _check_mode(mode)
    if (mode, name) in RULES:
        return RULES[(mode, name)]
    if mode == "psi-ms" and name in EXPASY_ONLY:
        return RULES[("expasy", name)]
    if name in LITERATURE_ONLY:
        return RULES[("literature", name)]
    raise KeyError(f"No rule named {name!r} in mode {mode!r}. Available: {available_enzymes(mode)}")


def custom_rule(
    after: str = "", before: str = "", except_before_p: bool = False, name: str = "Custom"
) -> CleavageRule:
    """Build a custom rule from the simple form ("cut after", "cut before", "except before P").

    Args:
        after: Residues after which to cut (C-terminal side), e.g. ``"KR"``.
        before: Residues before which to cut (N-terminal side), e.g. ``"D"``.
        except_before_p: If True, never cut when the next residue is P.
        name: Display name.

    Returns:
        A :class:`CleavageRule` whose ``regex`` is the compiled zero-width regex shown to
        the user.
    """
    after = _clean_residues(after)
    before = _clean_residues(before)
    if not after and not before:
        raise ValueError("A custom enzyme needs at least one 'cut after' or 'cut before' residue.")
    clauses = []
    if after:
        clauses.append(f"(?<=[{after}])")
    if before:
        clauses.append(f"(?=[{before}])")
    regex = "|".join(clauses)
    if except_before_p:
        regex = f"(?:{regex})(?!P)" if len(clauses) > 1 else f"{regex}(?!P)"
    desc = []
    if after:
        desc.append(f"after {after}")
    if before:
        desc.append(f"before {before}")
    note = "Custom rule: cut " + " or ".join(desc) + (", not before P" if except_before_p else "")
    return CleavageRule(name, regex, "Custom", "custom", "User-defined rule.", None, note)


def _clean_residues(residues: str) -> str:
    out = "".join(sorted(set(residues.upper().replace(" ", "").replace(",", ""))))
    bad = set(out) - set(STANDARD_RESIDUES)
    if bad:
        raise ValueError(f"Unknown residue code(s) {sorted(bad)}; use {STANDARD_RESIDUES}.")
    return out


def _check_mode(mode: str) -> None:
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
