"""Presets (CLAUDE.md section 5) and the digests compared under each preset.

Each preset defines a baseline digest. Its alternatives are every other enzyme of the
preset's rule-source mode, run alone, in parallel with the baseline, and sequentially after
the baseline (union of cut sites), all with the baseline's length, missed-cleavage, m/z and
chemistry settings.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import workflow as W
from .chemistry import Chemistry
from .digest import Digest, DigestSettings, digest, parallel
from .fasta import ProteinSet
from .rules import CleavageRule, available_enzymes, get_rule
from .workflow import Workflow


@dataclass(frozen=True)
class Preset:
    """A named baseline plus the rule-source mode used for its alternatives."""

    key: str
    name: str
    mode: str
    baseline: DigestSettings
    notes: str = ""


PRESETS: dict[str, Preset] = {
    p.key: p
    for p in (
        Preset(
            "proteome",
            "Proteome",
            "psi-ms",
            DigestSettings(
                (get_rule("Trypsin", "psi-ms"),),
                missed_cleavages=1,
                min_length=7,
                max_length=30,
                chemistry=Chemistry(),
            ),  # fmt: skip
            "Trypsin (PSI-MS), 7-30 aa, <=1 missed cleavage.",
        ),
        Preset(
            "histone_ptm",
            "Histone PTM",
            "psi-ms",
            DigestSettings(
                (get_rule("Trypsin", "psi-ms"),),
                missed_cleavages=1,
                min_length=4,
                max_length=40,
                chemistry=Chemistry(propionylation="before_and_after"),
            ),  # fmt: skip
            "Trypsin with protein-level propionylation (Arg-C-like), 4-40 aa, <=1 missed "
            "cleavage, 'propionylate before and after digestion' (formerly 'protein + peptide "
            "N-term').",
        ),
        Preset(
            "middle_down",
            "Middle-down",
            "psi-ms",
            DigestSettings(
                (get_rule("Glu-C", "psi-ms"),),
                missed_cleavages=0,
                min_length=7,
                max_length=70,
                chemistry=Chemistry(),
            ),  # fmt: skip
            "Glu-C (PSI-MS), 7-70 aa, 0 missed cleavages.",
        ),
        Preset(
            "protein_cleaver",
            "Protein Cleaver compatible",
            "expasy",
            DigestSettings(
                (get_rule("Arg-C", "expasy"),),
                missed_cleavages=0,
                min_length=7,
                max_length=30,
                mz_filter=False,
                chemistry=Chemistry(carbamidomethyl=False, met_removal="never"),
            ),  # fmt: skip
            "Arg-C (ExPASy), ExPASy rules, 7-30 aa, 0 missed cleavages, m/z filter off, no "
            "fixed mods. Initiator Met is not removed, as in Protein Cleaver (the preset "
            "table in CLAUDE.md does not specify Met removal).",
        ),
    )
}


@dataclass(frozen=True)
class DigestRecipe:
    """How to make one digest relative to a baseline.

    ``kind`` is ``baseline``, ``alone`` (the rule with the baseline's other settings),
    ``parallel`` (baseline digest and the ``alone`` digest, union of peptides) or
    ``sequential`` (baseline rules + the rule, union of cut sites).
    """

    kind: str
    baseline: DigestSettings
    rule: CleavageRule | None = None

    @property
    def workflow(self) -> Workflow:
        """The recipe in the strategy model (every digest with the baseline's chemistry)."""
        b = self.baseline
        chem = b.chemistry.propionylation
        if self.kind == "baseline":
            return W.standard(b)
        if self.kind == "alone":
            return W.alone(self.rule, chem)
        if self.kind == "parallel":
            return W.parallel_with(b, self.rule, chem)
        if self.kind == "sequential":
            return W.sequential_after(b, self.rule, chem)
        raise ValueError(f"unknown recipe kind {self.kind!r}")

    @property
    def settings(self) -> tuple[DigestSettings, ...]:
        """Settings of the digest(s) involved (two for a parallel digest)."""
        return self.workflow.components(self.baseline)

    @property
    def label(self) -> str:
        """Digest label, identical to :attr:`Digest.label` of the digest it makes."""
        return " | ".join(s.label for s in self.settings)

    @property
    def enzyme(self) -> str:
        """Label of the added rule ("" for the baseline)."""
        return self.rule.label if self.rule is not None else ""

    def run(self, proteins: ProteinSet, cache: dict | None = None) -> Digest:
        """Digest ``proteins``; ``cache`` (settings -> Digest) reuses single digests."""
        cache = {} if cache is None else cache

        def single(s: DigestSettings) -> Digest:
            if s not in cache:
                cache[s] = digest(proteins, s)
            return cache[s]

        ds = [single(s) for s in self.settings]
        return ds[0] if len(ds) == 1 else parallel(*ds)


def recipes_for(baseline: DigestSettings, rules: list[CleavageRule],
                kinds: tuple[str, ...] = ("alone", "parallel", "sequential"),
                ) -> list[DigestRecipe]:  # fmt: skip
    """Baseline recipe followed by ``kinds`` for each rule not already in the baseline."""
    out = [DigestRecipe("baseline", baseline)]
    for rule in rules:
        if rule in baseline.rules:
            continue
        out += [DigestRecipe(k, baseline, rule) for k in kinds]
    return out


def preset_recipes(preset: Preset) -> list[DigestRecipe]:
    """All digests of a preset: baseline, then every other enzyme of its mode alone, parallel
    and sequential."""
    rules = [get_rule(n, preset.mode) for n in available_enzymes(preset.mode)]
    return recipes_for(preset.baseline, rules)


def alternatives(preset: Preset) -> list[tuple[str, str]]:
    """(kind, enzyme name) for every alternative digest of ``preset``, in order."""
    return [(r.kind, r.rule.name) for r in preset_recipes(preset)[1:]]


def preset_digests(preset: Preset, proteins: ProteinSet) -> list[tuple[str, str, Digest]]:
    """Baseline and all alternative digests: list of (kind, enzyme label, Digest)."""
    cache: dict = {}
    return [(r.kind, r.enzyme, r.run(proteins, cache)) for r in preset_recipes(preset)]
