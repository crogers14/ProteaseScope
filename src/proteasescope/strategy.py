"""Histone protease strategy picker (Phase 6, ``STRATEGY_PICKER_2.md``).

The user picks a histone entry, where to look (a region, or specific lysines) and a goal.
The picker compares a small set of protease strategies with the standard digest (normally
the Histone PTM preset: protein-level propionylation + trypsin) and answers with up to three
recommendations, each described in plain sentences. Nothing is optimized: every strategy is
scored, strategies that change nothing are dropped, identical ones are merged, and the rest
are sorted.

Strategies, per candidate enzyme:

* ``alone``: the enzyme instead of the standard digest (1 prep, 1 enzyme step, 1 run).
* ``parallel``: the standard digest plus a separate digest with the enzyme; union of
  identifiable peptides (2 preps, 2 enzyme steps, 2 runs).
* ``sequential``: the standard digest, then the enzyme in the same tube; union of cut sites
  (1 prep, 2 enzyme steps, 1 run).

* ``custom``: a strategy from the builder (:class:`~proteasescope.workflow.Workflow`: 1-3
  digests, each with its own chemistry, run separately or in one tube).

Every strategy, the standard digest included, is expressed with the strategy model of
:mod:`proteasescope.workflow` (:meth:`Strategy.workflow_for`), so there is one code path.
A strategy can also be compared against a saved custom strategy instead of the standard
digest (the "compared-against" digest, :attr:`StrategyResults.reference`).

Chemistry: ``sequential`` is always "propionylate before and after digestion" (same tube,
lysines blocked for both enzymes). In ``sequential`` the second enzyme acts before the
peptide-level N-terminal propionylation, so every peptide N-terminus, including new ones,
carries the peptide-level propionyl; a single digest with the union of cut sites computes
exactly this. ``alone`` and ``parallel`` are propionylated before and after, plus an
``unmodified`` variant (separate prep, no derivatization) for enzymes that cleave at K; for
other enzymes the unmodified variant gives the same peptides.

Definitions (all restricted to the scope: a region in mature numbering, or selected lysines):

* **Peptide in scope:** an identifiable peptide that overlaps the region by >=1 residue
  (with selected lysines: one that carries >=1 selected lysine).
* **Lysine seen:** a lysine in scope on >=1 identifiable peptide.
* **Measurable alone (<= N):** on >=1 identifiable peptide carrying <= N lysines in total
  (lysines outside the scope count toward N).
* **Pair read together:** two lysines in scope on one identifiable peptide.
* **Pair newly separable** (``CUSTOM_STRATEGY.md`` §8): the two lysines share >=1
  identifiable peptide in the compared-against digest ("tangled"), and with the strategy
  each lysine also lies on >=1 identifiable peptide without the other. Pairs not tangled in
  the compared-against digest are never counted, so the goal's items are the tangled pairs
  that are separable (``Items.separable`` holds the raw per-digest separable pairs).
* **Known-mark lysine:** a lysine with >=1 entry in ``known_marks.tsv`` (optionally of one
  type); each position counts once.
* **Variant told apart:** another entry B of the same family is told apart if the chosen entry
  has >=1 identifiable peptide in scope whose sequence (I = L) is not contained anywhere in
  B's full sequence. Entries with an identical (I = L) sequence are excluded.

For every goal the items are sets, so Net = Strategy - Standard = Gained - Lost.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from itertools import combinations
from math import comb

import pandas as pd

from . import workflow as W
from .chemistry import Propionylation
from .digest import DigestSettings, digest
from .fasta import ProteinSet
from .rules import CleavageRule
from .workflow import Workflow, cleaves_at_lysine

# ----- where to look --------------------------------------------------------------------------
#: N-terminal tail per family (mature numbering, inclusive).
N_TAILS = {"H2A": (1, 20), "H2B": (1, 35), "H3": (1, 40), "H4": (1, 25)}
#: C-terminal tail per family, offered for entries of canonical length only.
C_TAILS = {"H2A": (110, 129)}
#: An entry gets the C-terminal tail preset if its mature length is within this many
#: residues of the tail's end (excludes e.g. macroH2A, whose residue 110-129 is not a tail).
C_TAIL_LENGTH_TOLERANCE = 10
N_TAIL, C_TAIL, WHOLE = "N-terminal tail", "C-terminal tail", "Whole protein"

#: Rules left out of the default candidate list.
EXCLUDED_BY_DEFAULT = frozenset({"Proteinase K"})


def region_presets(family: str, mature_length: int) -> dict[str, tuple[int, int] | None]:
    """Region presets offered for an entry (``None`` = whole protein). H1 and other entries
    without a tail preset get only the whole protein (plus a custom range in the UI)."""
    out: dict[str, tuple[int, int] | None] = {}
    if family in N_TAILS:
        out[N_TAIL] = N_TAILS[family]
    if family in C_TAILS and abs(mature_length - C_TAILS[family][1]) <= C_TAIL_LENGTH_TOLERANCE:
        out[C_TAIL] = C_TAILS[family]
    out[WHOLE] = None
    return out


# ----- strategies ----------------------------------------------------------------------------
KINDS = ("alone", "parallel", "sequential")
#: Every kind, including the standard digest and builder strategies.
ALL_KINDS = ("baseline", *KINDS, "custom")
VARIANTS = ("propionylated", "unmodified")

#: (sample preps, enzyme steps, LC-MS runs) per strategy kind.
COSTS = {"baseline": (1, 1, 1), "alone": (1, 1, 1), "parallel": (2, 2, 2),
         "sequential": (1, 2, 1)}  # fmt: skip

MIDDLE_DOWN_LENGTH = 40

#: Warning badges: key -> (badge text, tooltip).
WARNINGS = {
    "unmodified": ("Unmodified lysines only",
                   "This enzyme cuts at lysine and is used without propionylation. Acetylated "
                   "or methylated lysines block cleavage, so modified forms are missed."),
    "middle_down": ("Long peptides (middle-down)",
                    f"At least one peptide used is longer than {MIDDLE_DOWN_LENGTH} residues; "
                    "it needs middle-down methods."),
    "custom": ("Custom rule", "Uses a cleavage rule you defined."),
    "low_confidence": ("Low confidence",
                       "Nonspecific cleavage rule; its predictions are unreliable."),
}  # fmt: skip

_ENZYME_NAMES = {
    "Trypsin": "trypsin",
    "Trypsin/P": "trypsin/P",
    "Chymotrypsin": "chymotrypsin",
    "Chymotrypsin high specificity": "chymotrypsin (high specificity)",
    "Chymotrypsin low specificity": "chymotrypsin (low specificity)",
    "Asp-N ambic": "Asp-N (ambic)",
    "Proteinase K": "proteinase K",
    "Recombinant chymotrypsin": "recombinant chymotrypsin",
    "Custom": "the custom enzyme",
}


def enzyme_name(rule: CleavageRule) -> str:
    """Enzyme name as written in sentences ("trypsin", "Asp-N", "chymotrypsin (high
    specificity)")."""
    return _ENZYME_NAMES.get(rule.name, rule.name)


@dataclass(frozen=True)
class Strategy:
    """One way to digest the sample, relative to the standard digest.

    ``kind`` is ``baseline``, ``alone``, ``parallel``, ``sequential`` or ``custom``;
    ``variant`` is ``propionylated`` or ``unmodified`` (only ``alone``/``parallel`` can be
    unmodified). A ``custom`` strategy carries its :class:`~proteasescope.workflow.Workflow`
    and no rule.
    """

    kind: str
    rule: CleavageRule | None = None
    variant: str = "propionylated"
    workflow: Workflow | None = None

    def __post_init__(self) -> None:
        if self.kind not in ALL_KINDS:
            raise ValueError(f"unknown strategy kind {self.kind!r}")
        if self.variant not in VARIANTS:
            raise ValueError(f"unknown variant {self.variant!r}")
        if self.kind == "sequential" and self.variant != "propionylated":
            raise ValueError("sequential strategies are always propionylated (same tube)")
        if (self.kind in ("baseline", "custom")) != (self.rule is None):
            raise ValueError("only the baseline and custom strategies have no rule")
        if (self.kind == "custom") != (self.workflow is not None):
            raise ValueError("custom strategies (and only they) need a workflow")

    @classmethod
    def custom(cls, workflow: Workflow) -> Strategy:
        """A builder strategy."""
        return cls("custom", workflow=workflow)

    @property
    def code(self) -> str:
        """Compact id, e.g. ``parallel + Asp-N``, ``Lys-C alone (unmodified)`` or a custom
        strategy's name."""
        if self.kind == "baseline":
            return "standard digest"
        if self.kind == "custom":
            return self.workflow.name
        name = self.rule.name
        text = {"alone": f"{name} alone", "parallel": f"parallel + {name}",
                "sequential": f"sequential + {name}"}[self.kind]  # fmt: skip
        return text + (" (unmodified)" if self.variant == "unmodified" else "")

    @property
    def chemistry(self) -> Propionylation:
        """Chemistry of the added digest: before and after (propionylated) or none."""
        if self.variant == "unmodified":
            return Propionylation.NONE
        return Propionylation.BEFORE_AND_AFTER

    def workflow_for(self, baseline: DigestSettings) -> Workflow:
        """This strategy in the strategy model."""
        if self.kind == "baseline":
            return W.standard(baseline)
        if self.kind == "custom":
            return self.workflow
        if self.kind == "alone":
            return W.alone(self.rule, self.chemistry, self.code)
        if self.kind == "parallel":
            return W.parallel_with(baseline, self.rule, self.chemistry, self.code)
        return W.sequential_after(baseline, self.rule, Propionylation.BEFORE_AND_AFTER,
                                  self.code)  # fmt: skip

    @property
    def costs(self) -> tuple[int, int, int]:
        """(sample preps, enzyme steps, LC-MS runs)."""
        if self.kind == "custom":
            return self.workflow.costs
        return COSTS[self.kind]

    @property
    def enzymes(self) -> tuple[CleavageRule, ...]:
        """Enzymes this strategy adds (custom: every enzyme of the workflow)."""
        if self.kind == "custom":
            return self.workflow.enzymes
        return () if self.rule is None else (self.rule,)

    def components(self, baseline: DigestSettings) -> tuple[DigestSettings, ...]:
        """Settings of the digest(s) whose identifiable peptides make up this strategy.

        Several settings = parallel digests (union of identifiable peptides, each
        identifiable under its own settings).
        """
        return self.workflow_for(baseline).components(baseline)


def propionylated(chem):
    """The chemistry with lysines blocked ("propionylate before and after digestion")."""
    if chem.blocks_lysine:
        return chem
    return replace(chem, propionylation=Propionylation.BEFORE_AND_AFTER)


def unmodified(chem):
    """The chemistry without derivatization (separate prep on unmodified protein)."""
    return replace(chem, propionylation=Propionylation.NONE)


def enumerate_strategies(
    baseline: DigestSettings,
    rules: Sequence[CleavageRule],
    kinds: Sequence[str] = KINDS,
    all_unmodified: bool = False,
) -> list[Strategy]:
    """Every strategy for every candidate rule, in a fixed order (rule, kind, variant).

    Unmodified variants of ``alone``/``parallel`` are generated only for rules that cleave at
    K (``all_unmodified=True`` generates them for every rule). Propionylated strategies whose
    rule is already in the standard digest repeat it and are skipped.
    """
    out = []
    for rule in rules:
        in_base = rule in baseline.rules
        k_dep = all_unmodified or cleaves_at_lysine(rule)
        for kind in kinds:
            variants = ("propionylated", "unmodified") if kind != "sequential" and k_dep \
                else ("propionylated",)  # fmt: skip
            for variant in variants:
                if in_base and variant == "propionylated":
                    continue
                out.append(Strategy(kind, rule, variant))
    return out


def default_candidates(rules: Sequence[CleavageRule]) -> list[CleavageRule]:
    """Candidate list without Proteinase K."""
    return [r for r in rules if r.name not in EXCLUDED_BY_DEFAULT]


# ----- goals ------------------------------------------------------------------------------------
#: Known-mark filters: name -> modification abbreviations (``known_marks.tsv`` mod_abbrev).
MARK_FILTERS: dict[str, frozenset[str] | None] = {
    "any": None,
    "acetylation": frozenset({"ac", "acme"}),
    "methylation": frozenset({"me", "me1", "me2", "me3", "acme"}),
    "ubiquitination": frozenset({"ub"}),
}


@dataclass(frozen=True)
class Goal:
    """A goal as shown to the user, with the sentence templates that describe it.

    Templates use ``{n}`` (a count), ``{s}`` (plural suffix of the counted noun), ``{list}``
    and, for messages, ``{total}``, ``{b}``, ``{where}``.
    """

    key: str
    title: str
    short: str  # column name in the advanced table
    unit: str  # "lysine" | "pair" | "variant"
    gain: str
    loss: str
    complete: str
    status: str
    nothing: str


GOALS: dict[str, Goal] = {g.key: g for g in (
    Goal("see", "See lysines the standard digest misses", "Lysines seen", "lysine",
         "+{n} lysine{s} seen ({list})",
         "{n} lysine{s} no longer seen ({list})",
         "sees {all} lysine{ts}", "sees {b} of {total} lysine{ts}",
         "there are no lysines"),
    Goal("alone", "Measure lysines one at a time", "Lysines measurable alone", "lysine",
         "+{n} lysine{s} measurable on {their} own ({list})",
         "{n} lysine{s} no longer measurable on {their} own ({list})",
         "can measure {all} lysine{ts} on {their_t} own ({nk})",
         "can measure {b} of {total} lysine{ts} on their own ({nk})",
         "there are no lysines"),
    Goal("separate", "Separate lysines that share a peptide", "Separable lysine pairs", "pair",
         "+{n} lysine pair{s} that share a peptide can now be measured separately ({list})",
         "{n} lysine pair{s} that share a peptide can no longer be measured separately "
         "({list})",
         "can also measure {all} lysine pair{ts} that share a peptide separately",
         "can also measure {b} of the {total} lysine pair{ts} that share a peptide separately",
         "no two lysines share an identifiable peptide in the {ref}"),
    Goal("together", "Read lysines together (combinations / crosstalk)",
         "Lysine pairs on one peptide", "pair",
         "+{n} new lysine combination{s} on one peptide ({list})",
         "{n} lysine combination{s} no longer on one peptide ({list})",
         "reads {all} lysine pair{ts} together", "reads {b} of {total} lysine pair{ts} together",
         "there are fewer than two lysines"),
    Goal("marks", "Cover known modification sites", "Known modification sites seen", "lysine",
         "+{n} known {mark}site{s} seen ({list})",
         "{n} known {mark}site{s} no longer seen ({list})",
         "sees {all} known {mark}site{ts}", "sees {b} of {total} known {mark}site{ts}",
         "there are no known {mark}sites"),
    Goal("variants", "Tell variants and isoforms apart", "Variants told apart", "variant",
         "+{n} more variant{s} told apart ({list})",
         "{n} variant{s} no longer told apart ({list})",
         "tells {entry} apart from {all} other {family} entr{ies}",
         "tells {entry} apart from {b} of {total} other {family} entr{ies}",
         "no other {family} entry has a different sequence"),
)}  # fmt: skip


def plural(n: int) -> str:
    """Plural suffix for ``n`` items."""
    return "" if n == 1 else "s"


def short_list(items: Sequence[str], limit: int = 4) -> str:
    """``K27, K28 and K30``; lists longer than ``limit`` become ``K27, K28, K30 and 2 more``."""
    items = list(items)
    if len(items) > limit:
        return f"{', '.join(items[:3])} and {len(items) - 3} more"
    if len(items) <= 1:
        return "".join(items)
    return f"{', '.join(items[:-1])} and {items[-1]}"


# ----- scope and per-strategy items -----------------------------------------------------------
@dataclass
class Scope:
    """One histone entry and where to look in it.

    Positions are 0-based indices into the UniProt sequence; ``shift`` converts to mature
    numbering (mature = index + 1 - shift). With ``selected`` lysines the scope is those
    lysines only; otherwise the region ``r0:r1``.
    """

    accession: str
    display_name: str
    family: str
    sequence: str
    shift: int
    region: tuple[int, int]  # mature numbering, inclusive, clipped to the entry
    region_name: str
    r0: int
    r1: int
    lysines: list[int]
    selected: bool
    marks: dict[int, list[str]]
    others: list[tuple[str, str, str]]  # (accession, display name, I/L-folded sequence)
    identical: list[str]  # display names of entries with an identical (I = L) sequence

    def k_name(self, i: int) -> str:
        """``K<n>`` in mature numbering for 0-based index ``i``."""
        return f"K{i + 1 - self.shift}"

    def in_scope(self, start0: int, end0: int, ks: Sequence[int]) -> bool:
        """Whether a peptide (0-based half-open, carrying in-scope lysines ``ks``) is in
        scope."""
        if self.selected:
            return bool(ks)
        return start0 < self.r1 and end0 > self.r0

    def mark_sites(self, mark_filter: str = "any") -> set[int]:
        """Lysines in scope with >=1 known mark of the filtered type."""
        allowed = MARK_FILTERS[mark_filter]
        return {k for k, mods in self.marks.items()
                if allowed is None or allowed.intersection(mods)}  # fmt: skip

    @property
    def where(self) -> str:
        """``H2B type 1-C/E/F/G/I 1-35`` or ``H4 K16 and K20``."""
        if self.selected:
            return f"{self.display_name} {short_list([self.k_name(k) for k in self.lysines])}"
        a, b = self.region
        return f"{self.display_name} {a}–{b}"

    @property
    def description(self) -> str:
        """Scope for export headers."""
        if self.selected:
            ks = ", ".join(self.k_name(k) for k in self.lysines)
            return f"{self.display_name} ({self.accession}), selected lysines {ks}"
        a, b = self.region
        return (f"{self.display_name} ({self.accession}), residues {a}-{b} mature numbering"
                f" ({self.region_name})")  # fmt: skip

    @property
    def loses_label(self) -> str:
        """Lead-in of the "what you give up" line."""
        if self.selected:
            return "Loses (of the lysines you selected)"
        return "Loses in this protein" if self.region_name == WHOLE else "Loses in this region"


def make_scope(
    proteins: ProteinSet,
    annotation: pd.DataFrame,
    marks: pd.DataFrame,
    accession: str,
    region: tuple[int, int] | str | None,
    shift: int,
    lysines: Sequence[int] | None = None,
) -> Scope:
    """Build the :class:`Scope` for ``accession``.

    Args:
        region: (start, end) in mature numbering, inclusive; a key of
            :func:`region_presets`; or None for the whole protein. The end is clipped to the
            entry. Ignored when ``lysines`` is given.
        shift: 1 if the initiator Met is removed under the standard chemistry, else 0.
        lysines: Specific lysines (mature numbering); the scope is then those lysines only.
    """
    idx = proteins.index_of(accession)
    seq = proteins.sequences[idx]
    ann = annotation.set_index("accession")
    if accession not in ann.index:
        raise KeyError(f"{accession} is not in the histone annotation")
    family = ann.loc[accession, "family"]
    name = ann.loc[accession, "display_name"]
    mature_len = len(seq) - shift
    if lysines:
        ks = sorted({k - 1 + shift for k in lysines})
        bad = [k + 1 - shift for k in ks if not (shift <= k < len(seq) and seq[k] == "K")]
        if bad:
            raise ValueError(f"not lysines of {accession} (mature numbering): {bad}")
        region_name, a, b = "selected lysines", ks[0] + 1 - shift, ks[-1] + 1 - shift
    else:
        if isinstance(region, str):
            presets = region_presets(family, mature_len)
            if region not in presets:
                raise KeyError(f"no {region!r} preset for {accession}")
            region_name, region = region, presets[region]
        else:
            region_name = "custom" if region is not None else WHOLE
        a, b = (1, mature_len) if region is None else region
        if not 1 <= a <= b:
            raise ValueError(f"invalid region {a}-{b}")
        if a > mature_len:
            raise ValueError(f"region starts after the end of {accession} ({mature_len} aa)")
        b = min(b, mature_len)
    r0, r1 = a - 1 + shift, b + shift
    lys = ks if lysines else [i for i in range(r0, r1) if seq[i] == "K"]
    m = marks[(marks["accession"] == accession) & (marks["residue"] == "K")]
    by_pos: dict[int, list[str]] = {}
    for pos, abbrev in zip(m["position_uniprot"].astype(int), m["mod_abbrev"], strict=True):
        if pos - 1 in lys and abbrev not in by_pos.setdefault(pos - 1, []):
            by_pos[pos - 1].append(abbrev)
    folded = seq.replace("I", "L")
    others, identical = [], []
    accs = set(proteins.accessions)
    for acc, nm in zip(ann.index[ann["family"] == family],
                       ann.loc[ann["family"] == family, "display_name"], strict=True):  # fmt: skip
        if acc == accession or acc not in accs:
            continue
        f = proteins.sequences[proteins.index_of(acc)].replace("I", "L")
        if f == folded:
            identical.append(nm)
        else:
            others.append((acc, nm, f))
    return Scope(accession, name, family, seq, shift, (a, b), region_name, r0, r1, lys,
                 bool(lysines), by_pos, others, identical)  # fmt: skip


@dataclass
class Items:
    """What one strategy sees in the scope."""

    peptides: list[tuple[int, int, int]]  # (start0, end0, component bit mask), in scope
    key: frozenset[tuple[int, int]] = frozenset()  # in-scope peptides (identity)
    seen: set[int] = field(default_factory=set)
    min_k: dict[int, int] = field(default_factory=dict)  # fewest K on a peptide carrying k
    together: set[tuple[int, int]] = field(default_factory=set)
    separable: set[tuple[int, int]] = field(default_factory=set)
    variants: set[str] = field(default_factory=set)
    discriminating: dict[str, list[str]] = field(default_factory=dict)


def _items(peptides: dict[tuple[int, int], int], scope: Scope) -> Items:
    """Score a set of identifiable peptides ``{(start0, end0): component mask}``."""
    seq, lys = scope.sequence, set(scope.lysines)
    out = Items([])
    carriers: dict[int, set[int]] = {}
    for (s, e), mask in sorted(peptides.items()):
        ks_all = [i for i in range(s, e) if seq[i] == "K"]
        ks = [i for i in ks_all if i in lys]
        if not scope.in_scope(s, e, ks):
            continue
        pid = len(out.peptides)
        out.peptides.append((s, e, mask))
        for k in ks:
            out.seen.add(k)
            out.min_k[k] = min(out.min_k.get(k, len(ks_all)), len(ks_all))
            carriers.setdefault(k, set()).add(pid)
        out.together.update(combinations(ks, 2))
        pep = seq[s:e].replace("I", "L")
        for acc, _, fold in scope.others:
            if pep not in fold:
                out.variants.add(acc)
                out.discriminating.setdefault(acc, []).append(
                    f"{seq[s:e]} ({s + 1 - scope.shift}-{e - scope.shift})"
                )
    for i, j in combinations(scope.lysines, 2):
        ci, cj = carriers.get(i, set()), carriers.get(j, set())
        if ci - cj and cj - ci:
            out.separable.add((i, j))
    out.key = frozenset((s, e) for s, e, _ in out.peptides)
    return out


# ----- options, cards, messages -------------------------------------------------------------
@dataclass
class Option:
    """One recommendation: strategies of one kind that give identical peptides in scope
    (merged), described in plain words."""

    members: list[Strategy]
    title: str  # with **bold** enzyme names
    short: str
    warnings: list[str]  # keys of WARNINGS shared by every member

    @property
    def kind(self) -> str:
        return self.members[0].kind

    @property
    def label(self) -> str:
        """``title`` without markup (the sort name)."""
        return self.title.replace("**", "")

    @property
    def costs(self) -> tuple[int, int, int]:
        return self.members[0].costs


def extra_work(costs: tuple[int, int, int]) -> str:
    """Extra bench work relative to the standard digest (1 prep, 1 enzyme step, 1 run)."""
    preps, steps, runs = (c - b for c, b in zip(costs, COSTS["baseline"], strict=True))
    parts = []
    if preps:
        parts.append(f"+{preps} sample prep{plural(preps)}")
    elif steps:
        parts.append(f"+{steps} enzyme step{plural(steps)}")
    if runs:
        parts.append(f"+{runs} LC-MS run{plural(runs)}")
    elif parts:
        parts.append("no extra run")
    return ", ".join(parts) if parts else "same work"


def _variant_note(variants: set[str]) -> str:
    if variants == {"unmodified"}:
        return " (no propionylation)"
    if variants == {"propionylated", "unmodified"}:
        return " (with or without propionylation)"
    return ""


def option_title(members: Sequence[Strategy], baseline: DigestSettings) -> tuple[str, str]:
    """(title with **bold** enzymes, short name) for strategies of one kind merged into one
    recommendation."""
    kind = members[0].kind
    if kind == "custom":
        names = [m.workflow.name for m in members]
        bold = [f"**{x}**" for x in names]
        title = bold[0] if len(bold) == 1 else f"{', '.join(bold[:-1])} or {bold[-1]}"
        return title, " / ".join(names)
    by_rule: dict[CleavageRule, set[str]] = {}
    for m in members:
        by_rule.setdefault(m.rule, set()).add(m.variant)
    notes = {frozenset(v) for v in by_rule.values()}
    shared = _variant_note(set(next(iter(notes)))) if len(notes) == 1 else ""
    names, shorts = [], []
    for rule, variants in by_rule.items():
        note = "" if len(notes) == 1 else _variant_note(variants)
        names.append(f"**{enzyme_name(rule)}**{note}")
        shorts.append(enzyme_name(rule) + note.replace("no propionylation", "no prop.")
                      .replace("with or without propionylation", "±prop."))  # fmt: skip
    enz = names[0] if len(names) == 1 else f"{', '.join(names[:-1])} or {names[-1]}"
    enz_short = " / ".join(shorts)
    base = " + ".join(enzyme_name(r) for r in baseline.rules)
    short_shared = shared.replace("no propionylation", "no prop.").replace(
        "with or without propionylation", "±prop."
    )
    if kind == "alone":
        return f"Use {enz} instead{shared}", f"{enz_short} instead{short_shared}"
    if kind == "parallel":
        return f"Add a separate {enz} digest{shared}", f"+ separate {enz_short}{short_shared}"
    return (f"Follow {base} with {enz} in the same tube",
            f"{base} → {enz_short}")  # fmt: skip


@dataclass
class Card:
    """One answer card (plain sentences; ``**`` marks bold)."""

    rank: int
    option: Option
    title: str
    short: str
    gain_n: int
    headline: str  # "**+4** lysine pairs can now be measured separately (K5–K12, ...)"
    headline_items: list[str]  # full list (hover)
    goal_loss: str  # sentence when the goal also loses items, else ""
    gives_up: str  # always present
    gives_up_items: list[str]
    extra_work: str
    warnings: list[tuple[str, str]]  # (badge, tooltip)
    net: int
    hover: dict[str, list[str]] = field(default_factory=dict)  # item -> details
    work: str = ""  # absolute work, e.g. "2 sample preps, 2 enzymes, 2 LC-MS runs"

    def text(self) -> str:
        """The card as plain text."""
        lines = [f"{self.rank}. {self.title}", f"   {self.headline}"]
        if self.goal_loss:
            lines.append(f"   {self.goal_loss}")
        lines.append(f"   What you give up: {self.gives_up}")
        if self.option.kind == "custom":
            lines.append(f"   Work: {self.work}")
        else:
            lines.append(f"   Extra work: {self.extra_work}")
        if self.warnings:
            lines.append("   Warnings: " + "; ".join(f"[{b}]" for b, _ in self.warnings))
        return "\n".join(lines)


@dataclass
class StrategyResults:
    """Scores of every strategy for one scope. Build with :func:`pick_strategies`."""

    scope: Scope
    baseline: DigestSettings
    strategies: list[Strategy]
    items: dict[Strategy, Items]
    warnings: dict[Strategy, list[str]]
    #: The compared-against digest: the standard digest, or a saved custom strategy.
    reference: Strategy = field(default_factory=lambda: Strategy("baseline"))

    @property
    def base(self) -> Items:
        """Items seen by the compared-against digest (normally the standard digest)."""
        return self.items[self.reference]

    @property
    def ref_name(self) -> str:
        """``standard digest`` or the compared-against custom strategy's name."""
        return self.reference.code

    @property
    def ref_label(self) -> str:
        """Row label of the compared-against digest: ``Standard digest`` or its name."""
        return "Standard digest" if self.reference.kind == "baseline" else self.reference.code

    def get(self, strategy: Strategy | str) -> Strategy:
        """Strategy by object, code (``parallel + Asp-N``) or plain title."""
        if isinstance(strategy, Strategy):
            return strategy
        if strategy in ("standard digest", "baseline"):
            return Strategy("baseline")
        if strategy == self.reference.code:
            return self.reference
        for s in self.strategies:
            if strategy in (s.code, option_title([s], self.baseline)[0].replace("**", "")):
                return s
        raise KeyError(strategy)

    # -- goal values ----------------------------------------------------------------------
    def goal_set(self, it: Items, goal: str, n: int = 1, marks: str = "any") -> set:
        """The items a digest achieves for ``goal`` (N for "alone", mark filter for
        "marks")."""
        if goal == "see":
            return set(it.seen)
        if goal == "alone":
            return {k for k, c in it.min_k.items() if c <= n}
        if goal == "separate":
            # newly separable pairs only count if tangled in the compared-against digest
            return set(it.separable) & self.base.together
        if goal == "together":
            return set(it.together)
        if goal == "marks":
            return it.seen & self.scope.mark_sites(marks)
        if goal == "variants":
            return set(it.variants)
        raise KeyError(goal)

    def goal_total(self, goal: str, marks: str = "any") -> int:
        """Most items the goal can have in this scope."""
        nk = len(self.scope.lysines)
        return {"see": nk, "alone": nk, "separate": len(self.base.together),
                "together": comb(nk, 2),
                "marks": len(self.scope.mark_sites(marks)),
                "variants": len(self.scope.others)}[goal]  # fmt: skip

    def names(self, goal: str, items) -> list[str]:
        """Item names in reading order (``K5``, ``K5–K12``, variant names)."""
        sc = self.scope
        if GOALS[goal].unit == "pair":
            return [f"{sc.k_name(i)}–{sc.k_name(j)}" for i, j in sorted(items)]
        if goal == "variants":
            nm = {acc: n for acc, n, _ in sc.others}
            order = {acc: i for i, (acc, _, _) in enumerate(sc.others)}
            return [nm[a] for a in sorted(items, key=order.get)]
        return [sc.k_name(k) for k in sorted(items)]

    def _fill(self, template: str, goal: str, n: int = 0, names=(), *, total: int = 0,
              b: int = 0, nk: int = 1, marks: str = "any") -> str:  # fmt: skip
        return template.format(
            n=n,
            s=plural(n),
            list=short_list(list(names)),
            total=total,
            b=b,
            ts=plural(total),
            all=f"all {total}" if total != 1 else "the one",
            their_t="their" if total != 1 else "its",
            ies="ies" if total != 1 else "y",
            their="its" if n == 1 else "their",
            entry=self.scope.display_name,
            family=self.scope.family,
            where=self.scope.where,
            nk=f"peptides with at most {nk} lysine{plural(nk)}",
            mark="modification " if marks == "any" else f"{marks} ",
            ref=self.ref_name,
        )

    # -- options and ranking ----------------------------------------------------------------
    def options(self) -> list[Option]:
        """Strategies that differ from the standard digest in scope, merged when one kind
        gives identical peptides in scope (enumeration order)."""
        groups: dict[tuple[str, frozenset], list[Strategy]] = {}
        for s in self.strategies:
            key = self.items[s].key
            if key == self.base.key:
                continue
            groups.setdefault((s.kind, key), []).append(s)
        out = []
        for members in groups.values():
            title, short = option_title(members, self.baseline)
            shared = [w for w in self.warnings[members[0]]
                      if all(w in self.warnings[m] for m in members)]  # fmt: skip
            out.append(Option(members, title, short, shared))
        return out

    def ranked(self, goal: str, n: int = 1, marks: str = "any",
               positive_only: bool = True) -> list[tuple[Option, dict]]:  # fmt: skip
        """Options sorted by the ranking rules, with their goal numbers.

        Sort: goal net (descending), lysines lost (ascending), LC-MS runs, sample preps,
        name. With ``positive_only`` (the cards), only options with goal net > 0.
        """
        rows = []
        for opt in self.options():
            nums = self.numbers(opt.members[0], goal, n, marks)
            if positive_only and nums["net"] <= 0:
                continue
            rows.append((opt, nums))
        rows.sort(key=lambda r: (-r[1]["net"], len(r[1]["lysines_lost"]), r[0].costs[2],
                                 r[0].costs[0], r[0].label))  # fmt: skip
        return rows

    def numbers(self, strategy: Strategy, goal: str, n: int = 1, marks: str = "any") -> dict:
        """Goal numbers of one strategy versus the compared-against digest: ``standard``,
        ``strategy`` (counts), ``gained``, ``lost`` (sets), ``net``, ``lysines_lost``."""
        b = self.goal_set(self.base, goal, n, marks)
        it = self.items[strategy]
        a = self.goal_set(it, goal, n, marks)
        return {"standard": len(b), "strategy": len(a), "gained": a - b, "lost": b - a,
                "net": len(a) - len(b), "lysines_lost": self.base.seen - it.seen}  # fmt: skip

    def _card(self, rank: int, opt: Option, nums: dict, goal: str, n: int,
              marks: str) -> Card:  # fmt: skip
        g = GOALS[goal]
        gained = self.names(goal, nums["gained"])
        lost = self.names(goal, nums["lost"])
        lys_lost = self.names("see", nums["lysines_lost"])
        headline = self._fill(g.gain, goal, len(gained), gained, nk=n, marks=marks)
        headline = headline.replace(f"+{len(gained)}", f"**+{len(gained)}**", 1)
        if not gained:  # only for card_for (ranked cards always gain)
            headline = "**+0** nothing new for this goal"
        goal_loss = ("" if not lost or goal == "see"
                     else "but " + self._fill(g.loss, goal, len(lost), lost, nk=n,
                                              marks=marks))  # fmt: skip
        if self.scope.selected:
            gives = f"{self.scope.loses_label}: {short_list(lys_lost) or 'nothing'}"
        elif lys_lost:
            gives = f"{self.scope.loses_label}: {short_list(lys_lost)}"
        else:
            gives = "Nothing lost in this " + (
                "protein" if self.scope.region_name == WHOLE else "region"
            )
        hover: dict[str, list[str]] = {}
        if goal == "variants":
            it = self.items[opt.members[0]]
            nm = {acc: x for acc, x, _ in self.scope.others}
            hover = {nm[a]: it.discriminating[a] for a in nums["gained"]}
        elif goal == "marks":
            hover = {self.scope.k_name(k): self.scope.marks[k] for k in nums["gained"]}
        work = opt.members[0].workflow_for(self.baseline).work
        return Card(rank, opt, opt.title, opt.short, len(gained), headline, gained, goal_loss,
                    gives, lys_lost, extra_work(opt.costs), [WARNINGS[w] for w in opt.warnings],
                    nums["net"], hover, work)  # fmt: skip

    def cards(self, goal: str, n: int = 1, marks: str = "any", k: int = 3) -> list[Card]:
        """Up to ``k`` answer cards for ``goal`` (none if nothing improves it)."""
        top = self.ranked(goal, n, marks)[:k]
        return [self._card(rank, opt, nums, goal, n, marks)
                for rank, (opt, nums) in enumerate(top, start=1)]  # fmt: skip

    def card_for(self, strategy: Strategy | str, goal: str, n: int = 1,
                 marks: str = "any") -> Card:  # fmt: skip
        """The answer card of one strategy (e.g. a custom one), whatever its net change."""
        s = self.get(strategy)
        title, short = option_title([s], self.baseline)
        opt = Option([s], title, short, list(self.warnings.get(s, [])))
        return self._card(1, opt, self.numbers(s, goal, n, marks), goal, n, marks)

    def all_goals_summary(self, strategy: Strategy | str, n: int = 1,
                          marks: str = "any") -> pd.DataFrame:  # fmt: skip
        """One row per goal for one strategy versus the compared-against digest.

        Columns: ``goal``, ``Compared-against → Custom`` (e.g. ``4 → 5``), ``Gained``,
        ``Lost``, ``Net``, the numbers behind them (``compared_against``, ``strategy``,
        ``total``) and the items gained and lost.
        """
        s = self.get(strategy)
        rows = []
        for key, g in GOALS.items():
            nums = self.numbers(s, key, n, marks)
            rows.append({
                "goal": g.title,
                "Compared-against → Custom": f"{nums['standard']} → {nums['strategy']}",
                "Gained": len(nums["gained"]),
                "Lost": len(nums["lost"]),
                "Net": nums["net"],
                "compared_against": nums["standard"],
                "strategy": nums["strategy"],
                "total": self.goal_total(key, marks),
                "gained_items": ", ".join(self.names(key, nums["gained"])),
                "lost_items": ", ".join(self.names(key, nums["lost"])),
            })  # fmt: skip
        return pd.DataFrame(rows)

    def status(self, goal: str, n: int = 1, marks: str = "any") -> str:
        """``best``, ``already complete``, ``no improvement`` or ``nothing to measure``."""
        total = self.goal_total(goal, marks)
        if total == 0:
            return "nothing to measure"
        if len(self.goal_set(self.base, goal, n, marks)) == total:
            return "already complete"
        return "best" if self.ranked(goal, n, marks) else "no improvement"

    def message(self, goal: str, n: int = 1, marks: str = "any") -> str | None:
        """Plain message shown instead of cards when nothing improves the goal (else None)."""
        g = GOALS[goal]
        st = self.status(goal, n, marks)
        if st == "best":
            return None
        total = self.goal_total(goal, marks)
        b = len(self.goal_set(self.base, goal, n, marks))
        kw = {"total": total, "b": b, "nk": n, "marks": marks}
        where = self.scope.where
        if st == "nothing to measure":
            reason = self._fill(g.nothing, goal, **kw)
            tail = "Try another goal or region."
            return f"Nothing to measure for this goal in {where}: {reason}. {tail}"
        if st == "already complete":
            return (f"The {self.ref_name} already **{self._fill(g.complete, goal, **kw)}** in "
                    f"{where} with these settings, so no strategy can improve this. Try "
                    "another goal, region or settings.")  # fmt: skip
        return (f"No strategy improves this in {where} with these settings. The "
                f"{self.ref_name} {self._fill(g.status, goal, **kw)}, and every alternative "
                "gives the same or less. Try another goal, region or settings.")  # fmt: skip

    # -- advanced table ---------------------------------------------------------------------
    def advanced_table(self, n: int = 1, marks: str = "any") -> pd.DataFrame:
        """Every strategy (enumeration order) with every goal's Standard, Strategy, Gained,
        Lost and Net, costs, warnings, and how it appears among the recommendations."""
        shown = {}
        for opt in self.options():
            for m in opt.members:
                shown[m] = opt.label
        rows = []
        for s in self.strategies:
            it = self.items[s]
            preps, steps, runs = s.costs
            row = {"strategy": option_title([s], self.baseline)[0].replace("**", ""),
                   "code": s.code, "kind": s.kind,
                   "enzyme": ", ".join(r.label for r in s.enzymes),
                   "variant": "custom" if s.kind == "custom" else s.variant,
                   "same_as_standard": it.key == self.base.key,
                   "listed_as": shown.get(s, "")}  # fmt: skip
            for key, g in GOALS.items():
                b = self.goal_set(self.base, key, n, marks)
                a = self.goal_set(it, key, n, marks)
                row |= {f"{g.short}: standard": len(b), f"{g.short}: strategy": len(a),
                        f"{g.short}: gained": len(a - b), f"{g.short}: lost": len(b - a),
                        f"{g.short}: net": len(a) - len(b)}  # fmt: skip
            row |= {"lysines lost": len(self.base.seen - it.seen), "sample preps": preps,
                    "enzyme steps": steps, "LC-MS runs": runs,
                    "warnings": "; ".join(WARNINGS[w][0] for w in self.warnings[s])}  # fmt: skip
            rows.append(row)
        return pd.DataFrame(rows)

    # -- data for the visuals ------------------------------------------------------------------
    def _items_of(self, option: Option | None) -> Items:
        return self.base if option is None else self.items[option.members[0]]

    def lysine_states(self, goal: str, options: Sequence[Option], n: int = 1,
                      marks: str = "any") -> pd.DataFrame:  # fmt: skip
        """The lysine picture for a lysine-level goal ("see", "alone", "marks").

        One row per digest (``Standard digest`` first, then ``options``) and lysine in
        scope. ``state``: ``standard`` (achieved, as with the standard digest), ``gained``,
        ``lost``, ``neither``, or ``not counted`` (no known mark of the filtered type, for
        the "marks" goal).
        """
        if GOALS[goal].unit != "lysine":
            raise ValueError(f"{goal!r} is not a lysine-level goal")
        sc = self.scope
        counted = sc.mark_sites(marks) if goal == "marks" else set(sc.lysines)
        b = self.goal_set(self.base, goal, n, marks)
        rows = []
        for label, opt in [(self.ref_label, None)] + [(o.label, o) for o in options]:
            a = self.goal_set(self._items_of(opt), goal, n, marks)
            for k in sc.lysines:
                if k not in counted:
                    state = "not counted"
                elif opt is None:
                    state = "standard" if k in b else "neither"
                else:
                    state = {(1, 1): "standard", (0, 1): "gained", (1, 0): "lost",
                             (0, 0): "neither"}[(int(k in b), int(k in a))]  # fmt: skip
                rows.append({"row": label, "lysine": sc.k_name(k),
                             "position_mature": k + 1 - sc.shift, "state": state,
                             "marks": ", ".join(sc.marks.get(k, []))})  # fmt: skip
        return pd.DataFrame(rows)

    def pair_states(self, option: Option | None) -> pd.DataFrame:
        """Arcs for the pair goals: one row per lysine pair with a category.

        ``panel`` is ``standard`` (pairs the standard digest reads together) or ``strategy``
        (relative to the standard digest): ``together`` (read together by both), ``new
        together``, ``lost together``, ``new separable``, ``lost separable``. Positions are
        indices into the scope's lysines.
        """
        sc = self.scope
        pos = {k: i for i, k in enumerate(sc.lysines)}
        rows = []

        def add(panel, cat, pairs):
            for i, j in sorted(pairs):
                rows.append({"panel": panel, "category": cat, "i": pos[i], "j": pos[j],
                             "pair": f"{sc.k_name(i)}–{sc.k_name(j)}"})  # fmt: skip

        b = self.base
        add("standard", "together", b.together)
        if option is not None:
            a = self._items_of(option)
            add("strategy", "together", a.together & b.together)
            add("strategy", "new together", a.together - b.together)
            add("strategy", "lost together", b.together - a.together)
            b_sep = self.goal_set(b, "separate")
            a_sep = self.goal_set(a, "separate")
            add("strategy", "new separable", a_sep - b_sep)
            add("strategy", "lost separable", b_sep - a_sep)
        return pd.DataFrame(rows, columns=["panel", "category", "i", "j", "pair"])

    def variant_states(self, option: Option | None) -> pd.DataFrame:
        """Family entries with ``status``: ``told apart`` (also by the standard digest),
        ``newly told apart``, ``no longer told apart``, ``not told apart`` or ``identical
        sequence``; ``peptides`` lists the discriminating peptides (strategy, else
        standard)."""
        b = self.base
        a = self._items_of(option)
        rows = []
        for acc, name, _ in self.scope.others:
            in_b, in_a = acc in b.variants, acc in a.variants
            status = {(1, 1): "told apart", (0, 1): "newly told apart",
                      (1, 0): "no longer told apart", (0, 0): "not told apart"}[
                (int(in_b), int(in_a))]  # fmt: skip
            peps = a.discriminating.get(acc) or b.discriminating.get(acc) or []
            rows.append({"entry": name, "accession": acc, "status": status,
                         "peptides": "; ".join(peps)})  # fmt: skip
        for name in self.scope.identical:
            rows.append({"entry": name, "accession": "", "status": "identical sequence",
                         "peptides": ""})  # fmt: skip
        return pd.DataFrame(rows, columns=["entry", "accession", "status", "peptides"])

    def map_window(self) -> tuple[int, int]:
        """Residue window (mature numbering) of the scope's peptide map: the region, or for
        selected lysines the span of every in-scope peptide."""
        if not self.scope.selected:
            return self.scope.region
        sh = self.scope.shift
        spans = [(s + 1 - sh, e - sh) for it in self.items.values() for s, e, _ in it.peptides]
        a, b = self.scope.region
        return (min([a] + [x for x, _ in spans]), max([b] + [y for _, y in spans]))

    def map_data(self, option: Option | None, label: str = "",
                 split: bool = False) -> dict[str, pd.DataFrame | list[str]]:  # fmt: skip
        """Peptide map data (same layout as :func:`proteasescope.metrics.peptide_map`) for
        the compared-against digest and ``option``: in-scope peptides, lysines in scope,
        residues gained/lost inside the window (combined row only), and coverage of the
        window. ``labels`` lists the rows in order.

        With ``split``, a strategy of several separate digests also gets one row per digest
        (``  <enzyme>, <chemistry>``) before its combined row.
        """
        sc, (w0, w1) = self.scope, self.map_window()
        seq = sc.sequence
        rows = [(self.ref_label, self.base)]
        if option is not None:
            it = self._items_of(option)
            wf = option.members[0].workflow_for(self.baseline)
            if split and wf.combine == "separate" and len(wf.digests) > 1:
                for bit, step in enumerate(wf.digests):
                    part = [(s, e, m) for s, e, m in it.peptides if m >> bit & 1]
                    rows.append((f"  {step.label}", Items(part)))
            rows.append((label or option.label, it))
        peps, cov, covered = [], [], {}
        for lab, it in rows:
            res = set()
            for s, e, _ in it.peptides:
                ks = sum(1 for i in range(s, e) if seq[i] == "K")
                peps.append({"digest": lab, "peptide": seq[s:e], "start_mature": s + 1 - sc.shift,
                             "end_mature": e - sc.shift, "length": e - s, "n_K": ks})  # fmt: skip
                res.update(r for r in range(s + 1 - sc.shift, e + 1 - sc.shift) if w0 <= r <= w1)
            covered[lab] = res
            cov.append({"digest": lab, "coverage_window_pct": 100 * len(res) / (w1 - w0 + 1)})
        strips = []
        if option is not None:
            b_res, (al, a_res) = covered[rows[0][0]], list(covered.items())[-1]
            for status, resid in (("gained", a_res - b_res), ("lost", b_res - a_res)):
                for start, end in _runs(sorted(resid)):
                    strips.append({"digest": al, "status": status, "start_mature": start,
                                   "end_mature": end})  # fmt: skip
        lys = pd.DataFrame({"position_mature": [k + 1 - sc.shift for k in sc.lysines]})
        return {
            "peptides": pd.DataFrame(peps, columns=["digest", "peptide", "start_mature",
                                                    "end_mature", "length", "n_K"]),
            "lysines": lys,
            "strips": pd.DataFrame(strips, columns=["digest", "status", "start_mature",
                                                    "end_mature"]),
            "coverage": pd.DataFrame(cov),
            "labels": [lab for lab, _ in rows],
        }  # fmt: skip

    def answer_text(self, goal: str, n: int = 1, marks: str = "any") -> str:
        """The answer (cards or message) as plain text, with the settings."""
        g = GOALS[goal]
        head = ["ProteaseScope strategy picker - PREDICTED (in silico)",
                standard_line(self.baseline), f"Histone: {self.scope.description}",
                f"Goal: {g.title}"]  # fmt: skip
        if goal == "alone":
            head.append(f"Measurable on its own: peptides with at most {n} lysine(s)")
        if goal == "marks":
            head.append(f"Known modification sites: {marks}")
        customs = [s for s in [self.reference, *self.strategies] if s.kind == "custom"]
        head += [f"Custom strategy: {s.workflow.description}" for s in dict.fromkeys(customs)]
        if self.reference.kind == "custom":
            head.append(f"Compared against: {self.reference.workflow.name}")
        cards = self.cards(goal, n, marks)
        body = [c.text().replace("**", "") for c in cards] or [
            (self.message(goal, n, marks) or "").replace("**", "")
        ]
        return "\n".join(head) + "\n\n" + "\n\n".join(body) + "\n"

    # -- peptides and headers -----------------------------------------------------------------
    def peptides(self, strategy: Strategy | str) -> pd.DataFrame:
        """In-scope identifiable peptides of a strategy (mature numbering), with the K they
        carry, for the peptide map."""
        s = self.get(strategy)
        seq, sc = self.scope.sequence, self.scope
        comps = [c.label for c in s.components(self.baseline)]
        base = self.base.key
        rows = []
        for a, b, mask in self.items[s].peptides:
            ks = [i for i in range(a, b) if seq[i] == "K"]
            rows.append({
                "strategy": s.code,
                "start_mature": a + 1 - sc.shift,
                "end_mature": b - sc.shift,
                "peptide": seq[a:b],
                "length": b - a,
                "n_K": len(ks),
                "K_positions": ",".join(sc.k_name(i) for i in ks),
                "from": " | ".join(c for i, c in enumerate(comps) if mask >> i & 1),
                "in_standard": (a, b) in base,
            })  # fmt: skip
        return pd.DataFrame(rows)

    def header_lines(self, title: str = "Strategy picker", n: int = 1,
                     marks: str = "any") -> list[str]:  # fmt: skip
        """Settings header (without ``#``) for exports."""
        from . import report

        lines = report.header_lines(
            [(self.baseline,)],
            proteins_name=self.scope.description,
            background_description=(
                f"variants: the other {self.scope.family} entries of the histone annotation; "
                "a peptide tells entry B apart if B's full sequence does not contain it "
                "(I = L)"),
            title=title,
        )  # fmt: skip
        cands = list(dict.fromkeys(s.rule.label for s in self.strategies if s.rule is not None))
        lines += [
            standard_line(self.baseline),
            f"Measurable alone: peptides with at most {n} lysine(s) (all K on the peptide "
            f"count). Known-mark filter: {marks}.",
            "Strategies: alone = enzyme instead of the standard digest; parallel = standard + "
            "separate enzyme digest (union of identifiable peptides); sequential = standard "
            "then enzyme in the same tube (union of cut sites, always propionylated). "
            "Unmodified = separate prep without propionylation (lysine-cleaving enzymes).",
            *(["Candidate enzymes: " + ", ".join(cands)] if cands else []),
            "Separable pairs: only pairs that share an identifiable peptide in the "
            f"compared-against digest ({self.ref_name}) are counted.",
            "Net = Strategy - Standard = Gained - Lost.",
        ]  # fmt: skip
        customs = [s for s in [self.reference, *self.strategies] if s.kind == "custom"]
        for s in dict.fromkeys(customs):
            lines.append(f"Custom strategy: {s.workflow.description}")
        if self.reference.kind == "custom":
            lines.append(f"Compared against: {self.reference.workflow.name}")
        return lines


def _runs(values: Sequence[int]) -> list[tuple[int, int]]:
    """Consecutive runs of sorted integers as (first, last)."""
    out: list[tuple[int, int]] = []
    for v in values:
        if out and v == out[-1][1] + 1:
            out[-1] = (out[-1][0], v)
        else:
            out.append((v, v))
    return out


def standard_line(baseline: DigestSettings) -> str:
    """``Compared with: propionylation + trypsin · 4–40 aa · ≤1 missed cleavage``."""
    enz = " + ".join(enzyme_name(r) for r in baseline.rules)
    parts = [("propionylation + " if baseline.chemistry.blocks_lysine else "") + enz,
             f"{baseline.min_length}–{baseline.max_length} aa"]  # fmt: skip
    mc = baseline.missed_cleavages
    parts.append("no missed cleavages" if mc == 0 else f"≤{mc} missed cleavage{plural(mc)}")
    if baseline.mz_filter:
        parts.append(f"m/z {baseline.mz_min:g}–{baseline.mz_max:g}")
    return "Compared with: " + " · ".join(parts)


def pick_strategies(
    proteins: ProteinSet,
    annotation: pd.DataFrame,
    marks: pd.DataFrame,
    accession: str,
    region: tuple[int, int] | str | None,
    baseline: DigestSettings,
    rules: Sequence[CleavageRule],
    lysines: Sequence[int] | None = None,
    kinds: Sequence[str] = KINDS,
    all_unmodified: bool = False,
    custom: Sequence[Workflow] = (),
    against: Workflow | None = None,
) -> StrategyResults:
    """Enumerate and score every strategy for one histone entry and scope.

    Args:
        proteins: The histone set (all family members are needed for variants).
        annotation: Reviewed ``histone_annotation.tsv`` (family, display names).
        marks: Reviewed ``known_marks.tsv``.
        accession: The entry to study.
        region: (start, end) mature numbering inclusive, a :func:`region_presets` key, or
            None for the whole protein.
        baseline: Standard digest settings (normally the Histone PTM preset with the
            sidebar's length, missed-cleavage, rule-source and m/z settings).
        rules: Candidate enzymes (see :func:`default_candidates`).
        lysines: Specific lysines (mature numbering) instead of a region.
        custom: Saved builder strategies, scored like any other strategy.
        against: Compare against this saved strategy instead of the standard digest.
    """
    for w in [*custom, *([against] if against else [])]:
        w.validate()
    idx = proteins.index_of(accession)
    sub = proteins.subset([idx])
    shift = int(digest(sub, baseline).mature_offset[0])
    scope = make_scope(proteins, annotation, marks, accession, region, shift, lysines)
    strategies = enumerate_strategies(baseline, rules, kinds, all_unmodified)
    strategies += [Strategy.custom(w) for w in custom]
    reference = Strategy("baseline") if against is None else Strategy.custom(against)

    cache_: dict[DigestSettings, list[tuple[int, int]]] = {}

    def identifiable(s: DigestSettings) -> list[tuple[int, int]]:
        if s not in cache_:
            d = digest(sub, s)
            m = d.identifiable
            cache_[s] = list(zip(d.start0[m].tolist(), d.end0[m].tolist(), strict=True))
        return cache_[s]

    items: dict[Strategy, Items] = {}
    warnings: dict[Strategy, list[str]] = {}
    for st in dict.fromkeys([Strategy("baseline"), reference, *strategies]):
        peps: dict[tuple[int, int], int] = {}
        for bit, comp in enumerate(st.components(baseline)):
            for p in identifiable(comp):
                peps[p] = peps.get(p, 0) | (1 << bit)
        it = _items(peps, scope)
        items[st] = it
        if st.kind == "baseline":
            continue
        w = []
        added = st.components(baseline)
        if st.kind == "parallel":
            added = added[1:]
        if any(not c.chemistry.blocks_lysine and any(cleaves_at_lysine(r) for r in c.rules)
               for c in added):  # fmt: skip
            w.append("unmodified")
        if any(e - s > MIDDLE_DOWN_LENGTH for s, e, _ in it.peptides):
            w.append("middle_down")
        if any(r.mode == "custom" for r in st.enzymes):
            w.append("custom")
        if any({"nonspecific", "low confidence"} & set(r.flags) for r in st.enzymes):
            w.append("low_confidence")
        warnings[st] = w
    return StrategyResults(scope, baseline, strategies, items, warnings, reference)


def compare_workflows(
    proteins: ProteinSet,
    annotation: pd.DataFrame,
    marks: pd.DataFrame,
    accession: str,
    region: tuple[int, int] | str | None,
    baseline: DigestSettings,
    workflow: Workflow,
    against: Workflow | None = None,
    lysines: Sequence[int] | None = None,
) -> StrategyResults:
    """Score one custom strategy against the standard digest (or a saved strategy).

    ``baseline`` holds the shared settings (length window, missed cleavages, m/z filter,
    fixed mods, Met removal) and defines the standard digest. Use
    :meth:`StrategyResults.card_for` and :meth:`StrategyResults.all_goals_summary` with
    ``workflow.name``.
    """
    return pick_strategies(proteins, annotation, marks, accession, region, baseline, [],
                           lysines=lysines, custom=[workflow], against=against)  # fmt: skip


# ----- compare across all histones ------------------------------------------------------------
def grid_rows(proteins: ProteinSet, annotation: pd.DataFrame,
              baseline: DigestSettings) -> list[tuple[str, str]]:  # fmt: skip
    """(accession, region preset) rows of the "compare across all histones" grid: every
    default canonical histone with its tail presets (N-terminal; C-terminal where offered)."""
    from .chemistry import met_removed
    from .histones import cooccurrence_entries

    ann = annotation.set_index("accession")
    rows = []
    for acc in cooccurrence_entries(annotation):
        seq = proteins.sequences[proteins.index_of(acc)]
        shift = int(met_removed(seq, baseline.chemistry.met_removal))
        for name in region_presets(ann.loc[acc, "family"], len(seq) - shift):
            if name != WHOLE:
                rows.append((acc, name))
    return rows


def compare_all_histones(
    proteins: ProteinSet,
    annotation: pd.DataFrame,
    marks: pd.DataFrame,
    baseline: DigestSettings,
    rules: Sequence[CleavageRule],
    goals: Sequence[str] | None = None,
    n: int = 1,
    mark_filter: str = "any",
) -> pd.DataFrame:
    """One cell per canonical histone region and goal: the best strategy's short name and
    headline gain, or "already complete" / "no improvement" / "nothing to measure".

    Returns:
        Columns ``goal, accession, histone, region_name, region, status, best, gain, cell``.
    """
    goals = list(GOALS) if goals is None else list(goals)
    out = []
    for acc, region in grid_rows(proteins, annotation, baseline):
        res = pick_strategies(proteins, annotation, marks, acc, region, baseline, rules)
        for goal in goals:
            st = res.status(goal, n, mark_filter)
            best, gain = "", ""
            if st == "best":
                card = res.cards(goal, n, mark_filter, k=1)[0]
                best = card.short
                gain = card.headline.replace("**", "").split(" (")[0]
            out.append({
                "goal": goal,
                "accession": acc,
                "histone": res.scope.display_name,
                "region_name": region,
                "region": "{}–{}".format(*res.scope.region),
                "status": st,
                "best": best,
                "gain": gain,
                "cell": f"{best}: {gain}" if best else st,
            })  # fmt: skip
    return pd.DataFrame(out)
