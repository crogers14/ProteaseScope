"""Strategy model: digests, the chemistry of each, and how they are combined
(``CUSTOM_STRATEGY.md`` §7).

A strategy (here :class:`Workflow`, to keep it apart from the picker's enumerated
:class:`~proteasescope.strategy.Strategy`) is::

    {"name": "...",
     "digests": [{"enzyme": "...", "chemistry": "before_and_after | after_only | none"}],
     "combine": "separate | same_tube"}

* **separate:** each digest is run on its own (split sample) with its own chemistry; the
  strategy's identifiable peptides are the union of each digest's identifiable peptides.
* **same tube:** one chemistry for all digests (the first digest's); the cut sites of every
  enzyme are united (K sites removed if lysines are blocked) and the sample is digested once.

Every other setting (length window, missed cleavages, m/z filter, fixed modifications, Met
removal) comes from the shared digest settings (the sidebar or the preset). The standard
digest and the picker's ``alone``, ``parallel`` and ``sequential`` strategies are all
expressed with this model (:func:`standard`, :func:`alone`, :func:`parallel_with`,
:func:`sequential_after`), so there is one code path.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from functools import cache

from .chemistry import Propionylation
from .digest import DigestSettings, cut_sites
from .rules import CleavageRule, get_rule

COMBINES = ("separate", "same_tube")
COMBINE_LABELS = {
    "separate": "Separate digests, results combined",
    "same_tube": "Same tube, one after another",
}
COMBINE_TOOLTIPS = {
    "separate": "Split the sample and run each digest on its own.",
    "same_tube": "All enzymes cut the same sample; one LC-MS run.",
}
MAX_DIGESTS = 3

#: Guard-rail notes (``CUSTOM_STRATEGY.md`` §3).
NOTE_SAME_TUBE = "In one tube, all digests share the same chemistry."
NOTE_K_BLOCKED = "Lysines are blocked, so this enzyme will barely cut at K."
NOTE_NO_DERIVATIZATION = "Short, very basic peptides are often lost without derivatization."
NOTE_DUPLICATE = "This digest is already in the strategy"


@cache
def cleaves_at_lysine(rule: CleavageRule) -> bool:
    """True if any cut site of ``rule`` depends on K (it disappears when K is blocked).

    Tested on every residue pair around K (``xKy``), so the answer does not depend on the
    entry being digested.
    """
    aas = "ACDEFGHIKLMNPQRSTVWY"
    probe = "".join(f"G{a}K{b}G" for a in aas for b in aas)
    return len(cut_sites(probe, rule, block_lysine=False)) != len(
        cut_sites(probe, rule, block_lysine=True)
    )


@dataclass(frozen=True)
class DigestStep:
    """One digest of a strategy: its enzyme(s) and its chemistry.

    A builder row has exactly one enzyme. Several enzymes in one step mean one tube (used to
    express a standard digest whose settings hold several rules).
    """

    enzymes: tuple[CleavageRule, ...]
    chemistry: Propionylation = Propionylation.BEFORE_AND_AFTER

    def __post_init__(self) -> None:
        if isinstance(self.enzymes, CleavageRule):
            object.__setattr__(self, "enzymes", (self.enzymes,))
        object.__setattr__(self, "enzymes", tuple(self.enzymes))
        if not self.enzymes:
            raise ValueError("a digest needs at least one enzyme")
        object.__setattr__(self, "chemistry", Propionylation(self.chemistry))

    @property
    def enzyme(self) -> CleavageRule:
        """The (first) enzyme."""
        return self.enzymes[0]

    @property
    def label(self) -> str:
        """``Arg-C Ultra, no derivatization``."""
        enz = " + ".join(r.name for r in self.enzymes)
        return f"{enz}, {self.chemistry.label[0].lower()}{self.chemistry.label[1:]}"

    def settings(self, shared: DigestSettings) -> DigestSettings:
        """Digest settings for this step alone: ``shared`` with this step's enzymes and
        chemistry (label cleared). A step equal to the shared digest returns ``shared``
        itself."""
        if self.enzymes == shared.rules and self.chemistry == shared.chemistry.propionylation:
            return shared
        chem = replace(shared.chemistry, propionylation=self.chemistry)
        return replace(shared, rules=self.enzymes, name=None, chemistry=chem)

    def to_dict(self) -> dict:
        """JSON-ready dict (one ``enzyme``; ``enzymes`` when there are several)."""
        enz = [_rule_dict(r) for r in self.enzymes]
        out: dict = {"enzyme": enz[0]["enzyme"]} if len(enz) == 1 else {"enzymes": enz}
        if len(enz) == 1:
            out |= {k: v for k, v in enz[0].items() if k != "enzyme"}
        out["chemistry"] = self.chemistry.value
        return out


def _rule_dict(rule: CleavageRule) -> dict:
    out = {"enzyme": rule.name, "source": rule.mode}
    if rule.mode == "custom":
        out |= {"regex": rule.regex, "notes": rule.notes}
    return out


def _rule_from(d: dict, mode: str) -> CleavageRule:
    src = d.get("source") or mode
    if src == "custom":
        return CleavageRule(d["enzyme"], d["regex"], "Custom", "custom", "User-defined rule.",
                            None, d.get("notes", "Custom rule"))  # fmt: skip
    if src == "literature":
        src = mode
    return get_rule(d["enzyme"], src)


@dataclass(frozen=True)
class Workflow:
    """A strategy: 1-3 digests (builder; more for internal use) and how they are combined.

    For ``same_tube`` the chemistry of every digest is set to the first digest's when the
    strategy is created (``synced`` records that it changed something).
    """

    name: str
    digests: tuple[DigestStep, ...]
    combine: str = "separate"
    synced: bool = field(default=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "digests", tuple(self.digests))
        if not self.digests:
            raise ValueError("a strategy needs at least one digest")
        if self.combine not in COMBINES:
            raise ValueError(f"combine must be one of {COMBINES}, not {self.combine!r}")
        if self.combine == "same_tube":
            chem = self.digests[0].chemistry
            if any(d.chemistry != chem for d in self.digests):
                object.__setattr__(self, "digests",
                                   tuple(replace(d, chemistry=chem) for d in self.digests))
                object.__setattr__(self, "synced", True)  # fmt: skip

    # -- what is digested -----------------------------------------------------------------
    @property
    def enzymes(self) -> tuple[CleavageRule, ...]:
        """Every enzyme, in order."""
        return tuple(r for d in self.digests for r in d.enzymes)

    def components(self, shared: DigestSettings) -> tuple[DigestSettings, ...]:
        """Settings of the digest(s) whose identifiable peptides make up the strategy.

        ``separate``: one settings object per digest (union of identifiable peptides, each
        identifiable under its own chemistry). ``same_tube``: one settings object holding
        every enzyme (union of cut sites) with the shared chemistry.
        """
        if self.combine == "same_tube" or len(self.digests) == 1:
            step = DigestStep(self.enzymes, self.digests[0].chemistry)
            return (step.settings(shared),)
        return tuple(d.settings(shared) for d in self.digests)

    @property
    def costs(self) -> tuple[int, int, int]:
        """(sample preps, enzyme steps, LC-MS runs): separate = N, N, N; same tube = 1, N,
        1 (N = number of enzymes)."""
        n = len(self.enzymes)
        if self.combine == "separate":
            return (len(self.digests), n, len(self.digests))
        return (1, n, 1)

    @property
    def work(self) -> str:
        """``2 sample preps, 2 enzymes, 2 LC-MS runs``."""
        preps, steps, runs = self.costs
        return (f"{preps} sample prep{'s' * (preps != 1)}, {steps} enzyme{'s' * (steps != 1)}, "
                f"{runs} LC-MS run{'s' * (runs != 1)}")  # fmt: skip

    # -- guard rails ----------------------------------------------------------------------
    def duplicates(self) -> list[int]:
        """Indices of digests identical to an earlier one (same enzymes and chemistry)."""
        seen, out = set(), []
        for i, d in enumerate(self.digests):
            key = (d.enzymes, d.chemistry)
            if key in seen:
                out.append(i)
            seen.add(key)
        return out

    def notes(self) -> list[tuple[int | None, str, bool]]:
        """Plain-language guard-rail notes as (digest index or None, note, blocks Compare).

        Only identical digests block Compare.
        """
        out: list[tuple[int | None, str, bool]] = []
        if self.synced:
            out.append((None, NOTE_SAME_TUBE, False))
        dup = set(self.duplicates())
        for i, d in enumerate(self.digests):
            if d.chemistry.blocks_lysine and any(cleaves_at_lysine(r) for r in d.enzymes):
                out.append((i, NOTE_K_BLOCKED, False))
            if d.chemistry is Propionylation.NONE:
                out.append((i, NOTE_NO_DERIVATIZATION, False))
            if i in dup:
                out.append((i, NOTE_DUPLICATE, True))
        return out

    @property
    def can_compare(self) -> bool:
        """False if a guard rail blocks Compare (identical digests)."""
        return not self.duplicates()

    def validate(self) -> None:
        """Raise ``ValueError`` if Compare is blocked."""
        if not self.can_compare:
            raise ValueError(f"{self.name}: {NOTE_DUPLICATE.lower()}")

    # -- text and JSON --------------------------------------------------------------------
    @property
    def description(self) -> str:
        """For exports and Methods: ``RIPUP: Arg-C Ultra, no derivatization; recombinant
        chymotrypsin, no derivatization (separate digests, results combined)``."""
        steps = "; ".join(d.label for d in self.digests)
        if len(self.digests) == 1:
            return f"{self.name}: {steps}"
        return f"{self.name}: {steps} ({COMBINE_LABELS[self.combine].lower()})"

    def to_dict(self) -> dict:
        """JSON-ready dict in the documented format."""
        return {"name": self.name, "digests": [d.to_dict() for d in self.digests],
                "combine": self.combine}  # fmt: skip

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, ensure_ascii=False)

    @classmethod
    def from_dict(cls, data: dict, mode: str = "psi-ms") -> Workflow:
        """Load from a dict. Enzymes are looked up in their recorded ``source`` mode, else
        in ``mode``; old chemistry names migrate (see
        :class:`~proteasescope.chemistry.Propionylation`)."""
        steps = []
        for d in data["digests"]:
            rules = ([_rule_from(x, mode) for x in d["enzymes"]] if "enzymes" in d
                     else [_rule_from(d, mode)])  # fmt: skip
            steps.append(DigestStep(tuple(rules), Propionylation(d.get("chemistry", "none"))))
        return cls(data.get("name", "My strategy 1"), tuple(steps),
                   data.get("combine", "separate"))  # fmt: skip

    @classmethod
    def from_json(cls, text: str, mode: str = "psi-ms") -> Workflow:
        return cls.from_dict(json.loads(text), mode)


# ----- the picker's strategies, expressed with the model ------------------------------------
def standard(shared: DigestSettings, name: str = "Standard digest") -> Workflow:
    """The standard digest: the shared settings' rules (one tube) and chemistry."""
    return Workflow(name, (DigestStep(shared.rules, shared.chemistry.propionylation),))


def alone(rule: CleavageRule, chemistry: Propionylation | str, name: str = "") -> Workflow:
    """``rule`` instead of the standard digest."""
    return Workflow(name or f"{rule.name} alone", (DigestStep((rule,), chemistry),))


def parallel_with(shared: DigestSettings, rule: CleavageRule, chemistry: Propionylation | str,
                  name: str = "") -> Workflow:  # fmt: skip
    """The standard digest plus a separate digest with ``rule`` (union of peptides)."""
    steps = (*standard(shared).digests, DigestStep((rule,), chemistry))
    return Workflow(name or f"parallel + {rule.name}", steps, "separate")


def sequential_after(shared: DigestSettings, rule: CleavageRule,
                     chemistry: Propionylation | str | None = None,
                     name: str = "") -> Workflow:  # fmt: skip
    """The standard digest, then ``rule`` in the same tube (union of cut sites); chemistry
    defaults to the standard digest's."""
    chem = Propionylation(chemistry or shared.chemistry.propionylation)
    steps = (DigestStep(shared.rules, chem), DigestStep((rule,), chem))
    return Workflow(name or f"sequential + {rule.name}", steps, "same_tube")


# ----- builder presets ("Start from...", CUSTOM_STRATEGY.md §5) ------------------------------
@dataclass(frozen=True)
class BuilderPreset:
    """A "Start from..." entry: digests (enzyme name, chemistry), combination, citation."""

    name: str
    digests: tuple[tuple[str, str], ...]
    combine: str = "separate"
    citation: str = ""

    def workflow(self, mode: str = "psi-ms", name: str | None = None) -> Workflow:
        """The preset as a :class:`Workflow` (enzymes looked up in ``mode``)."""
        steps = tuple(DigestStep((get_rule(e, mode),), c) for e, c in self.digests)
        return Workflow(name or self.name, steps, self.combine)


TURNER_2026 = "Turner et al. 2026, Anal Chem 98:18694, doi:10.1021/acs.analchem.6c01147"

BUILDER_PRESETS: dict[str, BuilderPreset] = {p.name: p for p in (
    BuilderPreset("Standard histone digest", (("Trypsin", "before_and_after"),),
                  citation="Propionylation before and after trypsin digestion (standard "
                  "bottom-up histone workflow)."),
    BuilderPreset("RIPUP (Turner et al. 2026)",
                  (("Arg-C Ultra", "none"), ("Recombinant chymotrypsin", "none")), "separate",
                  citation=TURNER_2026),
    BuilderPreset("Arg-C Ultra + propionylation (Turner et al. 2026)",
                  (("Arg-C Ultra", "after_only"),), citation=TURNER_2026),
)}  # fmt: skip


def builder_preset(name: str, mode: str = "psi-ms") -> Workflow:
    """A "Start from..." preset as a :class:`Workflow`."""
    return BUILDER_PRESETS[name].workflow(mode)


def builder_enzymes(available: Sequence[str]) -> list[str]:
    """Enzyme dropdown order for histones: Arg-C Ultra and recombinant chymotrypsin first."""
    top = [e for e in ("Arg-C Ultra", "Recombinant chymotrypsin") if e in available]
    return top + [e for e in available if e not in top]
