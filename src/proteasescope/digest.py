"""Digestion engine: single, sequential and parallel digests.

Internal representation
-----------------------
All sequences of a :class:`~proteasescope.fasta.ProteinSet` are laid end to end in
one ASCII byte array. Peptides are stored as integer arrays only: entry index,
global start (inclusive) and global end (exclusive) into that array, plus the number
of missed cleavages. Sequences are sliced on demand.

Numbering conventions
---------------------
* ``start``/``end`` in exported tables: UniProt numbering, 1-based, inclusive.
* ``start_mature``/``end_mature``: UniProt - 1 when the initiator Met is removed,
  otherwise equal to UniProt.

Definitions (CLAUDE.md section 6)
---------------------------------
* Single digest: cut sites of one rule, after chemistry (K-derived sites removed when
  lysines are blocked).
* Sequential digest: union of the cut sites of several rules (one :class:`DigestSettings`
  with several rules).
* Parallel digests: union of peptide sets (:func:`parallel`).
* Identifiable: passes the length window and, if enabled, the m/z filter.
* Total peptides: all cleavage products with 0 missed cleavages, before filters.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from functools import cached_property

import numpy as np
import pandas as pd

from .chemistry import (
    WATER_AVG,
    WATER_MONO,
    Chemistry,
    met_removed,
    mz,
    nterm_mass,
    residue_mass_table,
)
from .fasta import ProteinSet
from .rules import CleavageRule


def bincount(values: np.ndarray, minlength: int = 0) -> np.ndarray:
    """``np.bincount`` that also works in 32-bit WebAssembly (Pyodide), where numpy refuses
    to cast int64 input to its int32 index type."""
    return np.bincount(np.asarray(values).astype(np.intp, copy=False), minlength=minlength)


#: Placeholder written over K to find K-derived cut sites (J appears in no rule).
_BLOCKED_K = "J"
_K_TO_J = str.maketrans("K", _BLOCKED_K)


@dataclass(frozen=True)
class DigestSettings:
    """Settings for one digest (one or more rules; several rules = sequential digest).

    Attributes:
        rules: Cleavage rules; with more than one rule the cut sites are united
            (sequential digest).
        missed_cleavages: Maximum missed cleavages per peptide.
        min_length, max_length: Identifiable length window (residues, inclusive).
        mz_filter: If True, identifiable peptides need >=1 charge state in ``charges``
            with monoisotopic m/z inside [``mz_min``, ``mz_max``].
        chemistry: Fixed mods, propionylation and Met removal.
        name: Optional label; defaults to the rule labels joined by " + ".
    """

    rules: tuple[CleavageRule, ...]
    missed_cleavages: int = 0
    min_length: int = 7
    max_length: int = 30
    mz_filter: bool = False
    mz_min: float = 350.0
    mz_max: float = 1400.0
    charges: tuple[int, ...] = (2, 3, 4)
    chemistry: Chemistry = field(default_factory=Chemistry)
    name: str | None = None

    def __post_init__(self) -> None:
        rules = (self.rules,) if isinstance(self.rules, CleavageRule) else tuple(self.rules)
        if not rules:
            raise ValueError("At least one cleavage rule is required")
        object.__setattr__(self, "rules", rules)
        if self.missed_cleavages < 0:
            raise ValueError("missed_cleavages must be >= 0")
        if not 1 <= self.min_length <= self.max_length:
            raise ValueError("need 1 <= min_length <= max_length")

    @property
    def label(self) -> str:
        """Digest label used in the ``digest`` column."""
        if self.name:
            return self.name
        return " + ".join(r.label for r in self.rules)

    @property
    def is_sequential(self) -> bool:
        """True when several rules are combined (union of cut sites)."""
        return len(self.rules) > 1

    @cached_property
    def pattern(self) -> re.Pattern[str]:
        """One compiled regex whose matches are the union of all rules' cut sites."""
        if len(self.rules) == 1:
            return self.rules[0].pattern
        return re.compile("|".join(f"(?:{r.regex})" for r in self.rules))

    def with_rules(self, *rules: CleavageRule, name: str | None = None) -> DigestSettings:
        """Copy with different rules (and label)."""
        return replace(self, rules=tuple(rules), name=name)


class _Concat:
    """Concatenated sequences of a ProteinSet plus per-entry offsets (cached per set)."""

    def __init__(self, proteins: ProteinSet):
        self.text = "".join(proteins.sequences)
        self.codes = np.frombuffer(self.text.encode("ascii"), dtype=np.uint8)
        lengths = proteins.lengths
        self.offsets = np.zeros(len(lengths) + 1, dtype=np.int64)
        np.cumsum(lengths, out=self.offsets[1:])


_CONCAT_CACHE: dict[int, tuple[ProteinSet, _Concat]] = {}


def _concat(proteins: ProteinSet) -> _Concat:
    hit = _CONCAT_CACHE.get(id(proteins))
    if hit is None or hit[0] is not proteins:
        hit = (proteins, _Concat(proteins))
        _CONCAT_CACHE[id(proteins)] = hit
        if len(_CONCAT_CACHE) > 8:
            _CONCAT_CACHE.pop(next(iter(_CONCAT_CACHE)))
    return hit[1]


def mature_offsets(proteins: ProteinSet, chem: Chemistry) -> np.ndarray:
    """Per entry: 1 if the initiator Met is removed under ``chem``, else 0 (int64)."""
    return np.fromiter(
        (met_removed(s, chem.met_removal) for s in proteins.sequences),
        dtype=np.int64,
        count=len(proteins),
    )


def cut_sites(
    sequence: str,
    settings: DigestSettings | CleavageRule | Sequence[CleavageRule],
    block_lysine: bool | None = None,
) -> np.ndarray:
    """Interior cut positions for one sequence (0-based index of the first residue after
    each cut; ``0 < i < len(sequence)``).

    Args:
        sequence: The sequence as digested (already Met-trimmed if applicable).
        settings: Digest settings, a rule, or a list of rules (union of cut sites).
        block_lysine: Remove K-derived sites. Defaults to the settings' chemistry; False
            when only rules are given.

    A cut site is K-derived if it disappears when every K is replaced by an unreactive
    residue: blocking lysines can only remove sites, never create them.
    """
    if not isinstance(settings, DigestSettings):
        rules = (settings,) if isinstance(settings, CleavageRule) else tuple(settings)
        settings = DigestSettings(rules, chemistry=Chemistry(carbamidomethyl=False))
        block_lysine = bool(block_lysine)
    elif block_lysine is None:
        block_lysine = settings.chemistry.blocks_lysine
    n = len(sequence)
    pat = settings.pattern
    sites = [p for p in (m.start() for m in pat.finditer(sequence)) if 0 < p < n]
    if block_lysine and "K" in sequence:
        blocked = sequence.translate(_K_TO_J)
        keep = {m.start() for m in pat.finditer(blocked)}
        sites = [p for p in sites if p in keep]
    return np.asarray(sites, dtype=np.int64)


@dataclass
class Digest:
    """Result of digesting a :class:`ProteinSet`.

    Peptide arrays (one element per peptide, all cleavage products up to the missed-cleavage
    limit, before length/m/z filters):

    * ``protein``: entry index into ``proteins``.
    * ``gstart`` / ``gend``: global start (inclusive) / end (exclusive) into the concatenated
      sequence array; local 0-based start = ``gstart - offsets[protein]``.
    * ``missed``: missed cleavages.
    * ``origin``: bit mask of contributing digests (bit 0 for a single digest; see
      :func:`parallel`).

    ``mature_offset`` holds, per entry, 1 if the initiator Met was removed, else 0.
    """

    proteins: ProteinSet
    settings: tuple[DigestSettings, ...]
    protein: np.ndarray
    gstart: np.ndarray
    gend: np.ndarray
    missed: np.ndarray
    origin: np.ndarray
    mature_offset: np.ndarray
    n_cut_sites: int = 0

    # ----- basic derived arrays -------------------------------------------------
    @property
    def label(self) -> str:
        """Digest label (parallel digests joined by " | ")."""
        return " | ".join(s.label for s in self.settings)

    @property
    def chemistry(self) -> Chemistry:
        """Chemistry shared by all component digests."""
        return self.settings[0].chemistry

    @property
    def offsets(self) -> np.ndarray:
        """Global offset of each entry (length n_entries + 1)."""
        return _concat(self.proteins).offsets

    def __len__(self) -> int:
        return len(self.protein)

    @cached_property
    def length(self) -> np.ndarray:
        """Peptide lengths (residues)."""
        return self.gend - self.gstart

    @cached_property
    def start0(self) -> np.ndarray:
        """Local 0-based start (inclusive) within the entry sequence."""
        return (self.gstart - self.offsets[self.protein]).astype(self.gstart.dtype)

    @cached_property
    def end0(self) -> np.ndarray:
        """Local 0-based end (exclusive) within the entry sequence."""
        return (self.gend - self.offsets[self.protein]).astype(self.gend.dtype)

    @cached_property
    def _masses(self) -> tuple[np.ndarray, np.ndarray]:
        cat = _concat(self.proteins)
        chem = self.chemistry
        out = []
        for average, water in ((False, WATER_MONO), (True, WATER_AVG)):
            res = residue_mass_table(chem, average)[cat.codes]
            unknown = np.isnan(res)
            cum = np.concatenate(([0.0], np.cumsum(np.where(unknown, 0.0, res))))
            ucum = np.concatenate(([0], np.cumsum(unknown)))
            m = cum[self.gend] - cum[self.gstart] + water + nterm_mass(chem, average)
            m[(ucum[self.gend] - ucum[self.gstart]) > 0] = np.nan
            out.append(m)
        return out[0], out[1]

    @property
    def mono_mass(self) -> np.ndarray:
        """Monoisotopic neutral mass (Da), incl. fixed mods; NaN if a residue is ambiguous."""
        return self._masses[0]

    @property
    def avg_mass(self) -> np.ndarray:
        """Average neutral mass (Da), ExPASy average residue table, incl. fixed mods."""
        return self._masses[1]

    def mz(self, charge: int) -> np.ndarray:
        """Monoisotopic m/z at ``charge``."""
        return mz(self.mono_mass, charge)

    def identifiable_mask(self, settings: DigestSettings | None = None) -> np.ndarray:
        """Boolean mask of identifiable peptides.

        For a parallel digest a peptide is identifiable if it is identifiable under the
        settings of any digest it came from.
        """
        if settings is not None:
            return self._identifiable_under(settings)
        mask = np.zeros(len(self), dtype=bool)
        for i, s in enumerate(self.settings):
            from_i = (self.origin >> i) & 1 == 1
            mask |= from_i & self._identifiable_under(s)
        return mask

    def _identifiable_under(self, s: DigestSettings) -> np.ndarray:
        ok = (self.length >= s.min_length) & (self.length <= s.max_length)
        if s.mz_filter:
            in_window = np.zeros(len(self), dtype=bool)
            for z in s.charges:
                v = self.mz(z)
                in_window |= (v >= s.mz_min) & (v <= s.mz_max)
            ok &= in_window
        return ok

    @cached_property
    def identifiable(self) -> np.ndarray:
        """Cached :meth:`identifiable_mask` for this digest's own settings."""
        return self.identifiable_mask()

    @property
    def n_total_peptides(self) -> int:
        """All cleavage products with 0 missed cleavages, before filters."""
        return int((self.missed == 0).sum())

    # ----- selection and export ----------------------------------------------------
    def select(self, mask: np.ndarray) -> Digest:
        """New digest restricted to peptides where ``mask`` is True."""
        return Digest(
            self.proteins,
            self.settings,
            self.protein[mask],
            self.gstart[mask],
            self.gend[mask],
            self.missed[mask],
            self.origin[mask],
            self.mature_offset,
            self.n_cut_sites,
        )

    def sequences(self, mask: np.ndarray | None = None) -> list[str]:
        """Peptide sequences (optionally for ``mask`` only), sliced on demand."""
        text = _concat(self.proteins).text
        gs, ge = (self.gstart, self.gend) if mask is None else (self.gstart[mask], self.gend[mask])
        return [text[a:b] for a, b in zip(gs.tolist(), ge.tolist(), strict=True)]

    def to_frame(
        self, identifiable_only: bool = True, digest_label: str | None = None
    ) -> pd.DataFrame:
        """Peptide table (Output 1 columns available from the core digest).

        Positions use UniProt numbering (1-based, inclusive); ``*_mature`` subtract 1 when the
        initiator Met was removed. ``prev_aa``/``next_aa`` are the flanking residues in the
        UniProt sequence (``-`` at the sequence ends). ``K_positions`` lists lysines in
        UniProt numbering, ``;``-separated.
        """
        mask = self.identifiable if identifiable_only else np.ones(len(self), dtype=bool)
        d = self.select(mask)
        meta = self.proteins.meta
        text = _concat(self.proteins).text
        peps = d.sequences()
        gs, ge = d.gstart.tolist(), d.gend.tolist()
        first = self.offsets[d.protein].tolist()
        last = self.offsets[d.protein + 1].tolist()
        prev_aa = [text[a - 1] if a > f else "-" for a, f in zip(gs, first, strict=True)]
        next_aa = [text[b] if b < e else "-" for b, e in zip(ge, last, strict=True)]
        start = d.start0 + 1
        shift = self.mature_offset[d.protein]
        k_pos = [
            ";".join(str(s + i) for i, c in enumerate(p) if c == "K")
            for p, s in zip(peps, start.tolist(), strict=True)
        ]
        df = pd.DataFrame(
            {
                "accession": meta["accession"].to_numpy()[d.protein],
                "gene": meta["gene"].to_numpy()[d.protein],
                "protein_name": meta["protein_name"].to_numpy()[d.protein],
                "digest": digest_label or self.label,
                "peptide": peps,
                "prev_aa": prev_aa,
                "next_aa": next_aa,
                "start": start,
                "end": d.end0,
                "start_mature": start - shift,
                "end_mature": d.end0 - shift,
                "length": d.length,
                "missed_cleavages": d.missed,
                "mono_mass": d.mono_mass,
                "avg_mass": d.avg_mass,
                "mz_2": d.mz(2),
                "mz_3": d.mz(3),
                "mz_4": d.mz(4),
                "n_K": [p.count("K") for p in peps],
                "K_positions": k_pos,
            }
        )
        return df


def digest(proteins: ProteinSet, settings: DigestSettings) -> Digest:
    """Digest every entry of ``proteins`` with ``settings``.

    Returns all cleavage products with up to ``settings.missed_cleavages`` missed cleavages;
    use :attr:`Digest.identifiable` for the filtered set. When the initiator Met is removed
    the rules are applied to the mature sequence (residue 2 onward).
    """
    cat = _concat(proteins)
    offsets = cat.offsets
    n = len(proteins)
    shift = mature_offsets(proteins, settings.chemistry)
    starts = offsets[:-1] + shift
    ends = offsets[1:]

    chunks: list[np.ndarray] = []
    n_sites = np.zeros(n, dtype=np.int64)
    block = settings.chemistry.blocks_lysine
    for i, seq in enumerate(proteins.sequences):
        s = shift[i]
        sites = cut_sites(seq[s:] if s else seq, settings, block_lysine=block)
        if len(sites):
            chunks.append(sites + starts[i])
            n_sites[i] = len(sites)
    cuts = np.concatenate(chunks) if chunks else np.zeros(0, dtype=np.int64)
    cprot = np.repeat(np.arange(n, dtype=np.int64), n_sites.astype(np.intp))

    # Boundaries per entry: mature start, cut sites, end -> sorted by (entry, position).
    idx = np.arange(n, dtype=np.int64)
    pos = np.concatenate([starts, cuts, ends])
    prot = np.concatenate([idx, cprot, idx])
    order = np.lexsort((pos, prot))
    pos, prot = pos[order], prot[order]

    pieces = []
    for mc in range(settings.missed_cleavages + 1):
        d = mc + 1
        if d >= len(pos):
            break
        same = prot[:-d] == prot[d:]
        gs, ge = pos[:-d][same], pos[d:][same]
        keep = ge > gs
        pieces.append(
            (
                prot[:-d][same][keep],
                gs[keep],
                ge[keep],
                np.full(int(keep.sum()), mc, dtype=np.int16),
            )
        )
    if pieces:
        p, gs, ge, mcs = (np.concatenate(x) for x in zip(*pieces, strict=True))
    else:
        p = gs = ge = np.zeros(0, dtype=np.int64)
        mcs = np.zeros(0, dtype=np.int16)
    order = np.lexsort((ge, gs, p))
    it = index_dtype(int(offsets[-1]))
    return Digest(
        proteins,
        (settings,),
        p[order].astype(it),
        gs[order].astype(it),
        ge[order].astype(it),
        mcs[order],
        np.ones(len(order), dtype=np.int32),
        shift,
        int(len(cuts)),
    )


def index_dtype(n: int) -> type:
    """int32 when positions up to ``n`` fit (halves memory; native index type in Pyodide),
    else int64."""
    return np.int32 if n < np.iinfo(np.int32).max else np.int64


def parallel(*digests: Digest) -> Digest:
    """Parallel digests: union of peptide sets over the same proteins.

    Peptides are identical when entry, start and end match. ``origin`` bit ``i`` marks
    peptides produced by ``digests[i]`` (bits of already-parallel inputs are shifted).
    A peptide produced by several digests keeps its smallest missed-cleavage count.
    All inputs must share the same ProteinSet object and chemistry.
    """
    if not digests:
        raise ValueError("need at least one digest")
    base = digests[0]
    for d in digests[1:]:
        if d.proteins is not base.proteins:
            raise ValueError("parallel digests must use the same ProteinSet object")
        if d.chemistry != base.chemistry:
            raise ValueError("parallel digests must share chemistry settings")
    settings: list[DigestSettings] = []
    parts = []
    for d in digests:
        parts.append((d.protein, d.gstart, d.gend, d.missed,
                      (d.origin << len(settings)).astype(np.int32)))  # fmt: skip
        settings.extend(d.settings)
    p, gs, ge, mcs, org = (np.concatenate(x) for x in zip(*parts, strict=True))
    order = np.lexsort((ge, gs, p))
    p, gs, ge, mcs, org = p[order], gs[order], ge[order], mcs[order], org[order]
    new = np.ones(len(p), dtype=bool)
    new[1:] = (gs[1:] != gs[:-1]) | (ge[1:] != ge[:-1])  # gstart is global -> entry implied
    group = np.cumsum(new) - 1
    merged_origin = np.zeros(int(new.sum()), dtype=np.int32)
    np.bitwise_or.at(merged_origin, group, org)
    merged_missed = np.full(len(merged_origin), np.iinfo(np.int16).max, dtype=np.int16)
    np.minimum.at(merged_missed, group, mcs.astype(np.int16))
    return Digest(
        base.proteins,
        tuple(settings),
        p[new],
        gs[new],
        ge[new],
        merged_missed,
        merged_origin,
        base.mature_offset,
        sum(d.n_cut_sites for d in digests),
    )


def coverage(d: Digest, mask: np.ndarray | None = None) -> pd.DataFrame:
    """Per-entry residue coverage by peptides (default: identifiable peptides).

    Uses a difference array over the concatenated sequences (vectorized across entries).
    The denominator is the mature length (sequence length minus a removed initiator Met).

    Returns:
        DataFrame indexed like ``d.proteins.meta`` with ``accession``, ``length``,
        ``mature_length``, ``covered`` (residues) and ``coverage_pct``.
    """
    covered = covered_residues(d, mask)
    offs = d.offsets
    per = (
        np.add.reduceat(covered.astype(np.int64), offs[:-1].astype(np.intp))
        if len(covered)
        else np.zeros(0)
    )
    lengths = np.diff(offs)
    per = np.where(lengths > 0, per, 0)
    mature_len = lengths - d.mature_offset
    with np.errstate(invalid="ignore", divide="ignore"):
        pct = np.where(mature_len > 0, 100.0 * per / mature_len, 0.0)
    return pd.DataFrame(
        {
            "accession": d.proteins.meta["accession"].to_numpy(),
            "length": lengths,
            "mature_length": mature_len,
            "covered": per,
            "coverage_pct": pct,
        }
    )


def covered_residues(d: Digest, mask: np.ndarray | None = None) -> np.ndarray:
    """Boolean array over the concatenated sequences: residue lies on a selected peptide."""
    mask = d.identifiable if mask is None else mask
    total = int(d.offsets[-1])
    diff = np.zeros(total + 1, dtype=np.int32)
    np.add.at(diff, d.gstart[mask], 1)
    np.add.at(diff, d.gend[mask], -1)
    return np.cumsum(diff[:-1], dtype=np.int32) > 0
