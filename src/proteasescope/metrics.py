"""Peptide-to-protein mapping, uniqueness, and outputs 1-5 (CLAUDE.md section 7).

Mapping definitions
-------------------
* ``"substring"`` (default): a peptide maps to every entry of the reference database whose
  full UniProt sequence contains it, with I and L treated as identical, whether or not that
  entry's digest would produce the peptide. This definition is also used for histone
  isoform discrimination.
* ``"digest"`` (option): a peptide maps to an entry only if the I/L-folded peptide is one of
  that entry's cleavage products under the background digest's settings.

Entries are counted by distinct accession. Genes are counted by distinct GN= value;
an entry without a gene name counts as its own gene.

Substring search uses an exact index of the 8 residues starting at every background position
(packed into a uint64, so there are no hash collisions; residues past an entry's end are
coded 0). A peptide of <= 8 residues is a range lookup in that index; longer peptides anchor
on their first 8 residues and every candidate is verified residue by residue with numpy.
Results are exact. Full UniProt sequences are searched (a removed initiator Met is kept).
"""

from __future__ import annotations

import string
from collections.abc import Sequence
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

from .digest import (
    Digest,
    DigestSettings,
    _concat,
    bincount,
    covered_residues,
    digest,
    parallel,
)
from .fasta import ProteinSet
from .rules import CleavageRule

MAPPING_METHODS = ("substring", "digest")

#: Residues per index key (5 bits each, packed into a uint64).
KMER = 8
#: Maximum candidate hits verified per batch (bounds memory use).
_CANDIDATE_CHUNK = 4_000_000

# Letter -> small code (A=1 ... Z=26); I folded onto L. Non-letters -> 0.
_FOLD = np.zeros(256, dtype=np.uint8)
for _i, _ch in enumerate(string.ascii_uppercase):
    _FOLD[ord(_ch)] = _i + 1
_FOLD[ord("I")] = _FOLD[ord("L")]


def il_key(sequences: Sequence[str]) -> list[str]:
    """Sequences with I replaced by L (I/L are isobaric and treated as identical)."""
    return [s.replace("I", "L") for s in sequences]


class _KmerIndex:
    """Exact 8-mer keys of every background position, sorted.

    The key of position ``p`` packs residues ``p .. p+7`` (5 bits each, first residue in the
    high bits). Residues past the end of ``p``'s entry are coded 0, so keys never span two
    entries, and an L-residue prefix (L < 8) matches a contiguous range of keys.
    """

    def __init__(self, proteins: ProteinSet):
        cat = _concat(proteins)
        n = len(cat.codes)
        self.codes = _FOLD[cat.codes]
        self.offsets = cat.offsets
        lengths = np.diff(cat.offsets)
        self.prot_of = np.repeat(np.arange(len(lengths), dtype=np.int32), lengths.astype(np.intp))
        pos = np.arange(n, dtype=np.int64)
        remaining = (self.offsets[1:][self.prot_of] - pos).astype(np.int32)
        padded = np.concatenate([self.codes, np.zeros(KMER, dtype=np.uint8)])
        keys = np.zeros(n, dtype=np.uint64)
        for j in range(KMER):
            c = np.where(remaining > j, padded[j : n + j], 0).astype(np.uint64)
            keys = (keys << np.uint64(5)) | c
        del pos, remaining
        order = np.argsort(keys)
        self.keys = keys[order]
        self.pos = order.astype(np.int64 if n >= 2**31 else np.int32)


_INDEX_CACHE: dict[int, tuple[ProteinSet, _KmerIndex]] = {}


def _kmer_index(proteins: ProteinSet) -> _KmerIndex:
    hit = _INDEX_CACHE.get(id(proteins))
    if hit is None or hit[0] is not proteins:
        hit = (proteins, _KmerIndex(proteins))
        _INDEX_CACHE[id(proteins)] = hit
        while len(_INDEX_CACHE) > 2:
            _INDEX_CACHE.pop(next(iter(_INDEX_CACHE)))
    return hit[1]


def _substring_matches(
    query: Digest, background: ProteinSet, mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    qcodes = _FOLD[_concat(query.proteins).codes]
    qstart = query.gstart[mask]
    qlen = query.length[mask]
    if not len(qstart):
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    idx = _kmer_index(background)
    # Key of the first min(L, 8) query residues, left-aligned like the index keys.
    k = np.minimum(qlen, KMER)
    qkey = np.zeros(len(qstart), dtype=np.uint64)
    for j in range(KMER):
        c = np.where(k > j, qcodes[np.minimum(qstart + j, len(qcodes) - 1)], 0)
        qkey = (qkey << np.uint64(5)) | c.astype(np.uint64)
    span = np.left_shift(np.uint64(1), (5 * (KMER - k)).astype(np.uint64))
    lo = np.searchsorted(idx.keys, qkey, side="left")
    hi = np.searchsorted(idx.keys, qkey + span, side="left")
    counts = hi - lo

    out_q, out_p = [], []
    bounds = np.searchsorted(np.cumsum(counts), np.arange(0, counts.sum(), _CANDIDATE_CHUNK))
    bounds = np.unique(np.concatenate(([0], bounds, [len(counts)])))
    for a, b in zip(bounds[:-1], bounds[1:], strict=True):
        c = counts[a:b]
        if not c.sum():
            continue
        qi = np.repeat(np.arange(a, b), c)
        first = np.repeat(lo[a:b] - (np.cumsum(c) - c), c) + np.arange(int(c.sum()))
        p = idx.pos[first].astype(np.int64)
        length = qlen[qi]
        prot = idx.prot_of[p]
        # Peptides of <= 8 residues are exact matches already; longer ones must fit in the
        # entry and match residue by residue beyond the anchor.
        ok = p + length <= idx.offsets[prot + 1]
        for j in range(KMER, int(length.max())):
            sel = np.flatnonzero(ok & (length > j))
            if not len(sel):
                break
            ok[sel] = idx.codes[p[sel] + j] == qcodes[qstart[qi[sel]] + j]
        out_q.append(qi[ok])
        out_p.append(prot[ok].astype(np.int64))
    if not out_q:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(out_q), np.concatenate(out_p)


def _digest_matches(
    query: Digest, background: Digest, mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    q = pd.DataFrame({"key": il_key(query.sequences(mask))})
    q["qi"] = np.arange(len(q))
    bsel = np.isin(background.length, np.unique(query.length[mask]))
    bg = pd.DataFrame({"key": il_key(background.sequences(bsel)), "prot": background.protein[bsel]})
    pairs = q.merge(bg[bg["key"].isin(set(q["key"]))], on="key")
    return pairs["qi"].to_numpy(np.int64), pairs["prot"].to_numpy(np.int64)


def peptide_matches(
    query: Digest,
    background: Digest | ProteinSet,
    mask: np.ndarray | None = None,
    method: str = "substring",
) -> pd.DataFrame:
    """All (query peptide, background entry) pairs where the peptide maps to the entry.

    Args:
        query: Digest whose peptides are mapped.
        background: Reference database: a ProteinSet (``"substring"``) or a Digest
            (``"substring"`` uses its proteins; ``"digest"`` uses its cleavage products).
        mask: Query peptides to map (default: identifiable).
        method: ``"substring"`` (default) or ``"digest"``; see the module docstring.

    Returns:
        DataFrame with ``query`` (position among the selected query peptides) and
        ``protein`` (background entry index), one row per distinct pair.
    """
    if method not in MAPPING_METHODS:
        raise ValueError(f"method must be one of {MAPPING_METHODS}, got {method!r}")
    mask = query.identifiable if mask is None else mask
    if method == "substring":
        proteins = background.proteins if isinstance(background, Digest) else background
        qi, prot = _substring_matches(query, proteins, mask)
    else:
        if not isinstance(background, Digest):
            raise TypeError("method='digest' needs a Digest of the background")
        qi, prot = _digest_matches(query, background, mask)
    proteins = background.proteins if isinstance(background, Digest) else background
    qi, prot = _unique_pairs(qi, prot, len(proteins))
    return pd.DataFrame({"query": qi, "protein": prot})


def _unique_pairs(a: np.ndarray, b: np.ndarray, nb: int) -> tuple[np.ndarray, np.ndarray]:
    """Distinct (a, b) pairs sorted by a then b; ``b`` must lie in ``[0, nb)``."""
    combined = np.sort(a.astype(np.int64) * nb + b)  # np.unique's hash path is slow here
    if len(combined):
        combined = combined[np.concatenate(([True], combined[1:] != combined[:-1]))]
    return combined // nb, combined % nb


def peptide_mapping(
    query: Digest,
    background: Digest | ProteinSet,
    mask: np.ndarray | None = None,
    method: str = "substring",
) -> pd.DataFrame:
    """Count background entries and genes each query peptide maps to.

    See :func:`peptide_matches` for arguments. The background should include the query's
    own entries (CLAUDE.md: reference proteome + isoforms plus the uploaded FASTA).

    Returns:
        DataFrame (one row per selected query peptide, in order) with ``n_proteins``,
        ``n_genes``, ``unique_protein`` (maps to exactly one entry) and ``unique_gene``
        (maps to exactly one gene).
    """
    mask = query.identifiable if mask is None else mask
    pairs = peptide_matches(query, background, mask, method)
    proteins = background.proteins if isinstance(background, Digest) else background
    meta = proteins.meta
    gene = meta["gene"].where(meta["gene"].notna(), "acc:" + meta["accession"].astype(str))
    gene_id = pd.factorize(gene)[0]
    acc_id = pd.factorize(meta["accession"])[0]
    n = int(mask.sum())
    # entries are counted by distinct accession
    aq, _ = _unique_pairs(
        pairs["query"].to_numpy(), acc_id[pairs["protein"].to_numpy()], int(acc_id.max()) + 1
    )
    n_prot = bincount(aq, minlength=n)
    gq, _ = _unique_pairs(
        pairs["query"].to_numpy(), gene_id[pairs["protein"].to_numpy()], int(gene_id.max()) + 1
    )
    n_gene = bincount(gq, minlength=n)
    return pd.DataFrame(
        {
            "n_proteins": n_prot,
            "n_genes": n_gene,
            "unique_protein": n_prot == 1,
            "unique_gene": n_gene == 1,
        }
    )


# =========================================================================================
# Uniqueness background
# =========================================================================================
@dataclass(eq=False)
class Background:
    """The uniqueness background (reference database) and mapping method.

    Attributes:
        proteins: Reference entries (CLAUDE.md default: reviewed human proteome + isoforms,
            plus the uploaded FASTA; see :func:`uniqueness_background`).
        method: ``"substring"`` (default) or ``"digest"``.
        index: Optional precomputed :class:`SharedPrefixIndex` of the reference part of
            ``proteins``. With ``method="substring"`` uniqueness is then looked up instead of
            searched; entries not in the index (uploads) are handled by live mapping. With an
            index, ``n_proteins``/``n_genes`` are 1 for unique peptides and <NA> otherwise.
    """

    proteins: ProteinSet
    method: str = "substring"
    index: SharedPrefixIndex | None = None
    _digests: dict = field(default_factory=dict, repr=False)
    _novel: ProteinSet | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if self.method not in MAPPING_METHODS:
            raise ValueError(f"method must be one of {MAPPING_METHODS}")

    @property
    def description(self) -> str:
        """Background and method, as printed in export headers."""
        how = {
            "substring": "a peptide maps to every entry whose sequence contains it (I = L)",
            "digest": "a peptide maps to every entry whose cleavage products under the same "
            "settings include it (I = L)",
        }[self.method]
        return f"{self.proteins.name or f'{len(self.proteins)} entries'}; {how}"

    def _background_digest(self, settings: DigestSettings) -> Digest:
        if settings not in self._digests:
            self._digests[settings] = digest(self.proteins, settings)
        return self._digests[settings]

    def mapping(self, d: Digest, mask: np.ndarray | None = None) -> pd.DataFrame:
        """:func:`peptide_mapping` of ``d``'s peptides against this background."""
        if self.method == "substring" and self.index is not None:
            return self._indexed_mapping(d, d.identifiable if mask is None else mask)
        if self.method == "substring":
            return peptide_mapping(d, self.proteins, mask)
        bgs = [self._background_digest(s) for s in d.settings]
        bg = bgs[0] if len(bgs) == 1 else parallel(*bgs)
        return peptide_mapping(d, bg, mask, method="digest")

    def novel_entries(self) -> ProteinSet:
        """Background entries absent from the index (accession + sequence), e.g. uploads."""
        if self._novel is None:
            known = set(zip(self.index.proteins.accessions, self.index.proteins.sequences,
                            strict=True))  # fmt: skip
            keep = [i for i, k in enumerate(zip(self.proteins.accessions,
                                                self.proteins.sequences, strict=True))
                    if k not in known]  # fmt: skip
            self._novel = self.proteins.subset(keep)
        return self._novel

    def _indexed_mapping(self, d: Digest, mask: np.ndarray) -> pd.DataFrame:
        idx = self.index
        sel = np.flatnonzero(mask)
        q_off = idx.offsets_for(d.proteins, strict=False)
        entry = d.protein[sel]
        in_ref = q_off[entry] >= 0
        up = np.zeros(len(sel), dtype=bool)
        ug = np.zeros(len(sel), dtype=bool)
        gs = q_off[entry[in_ref]] + d.start0[sel][in_ref]
        length = d.length[sel][in_ref]
        up[in_ref] = idx.other_entry[gs] < length
        ug[in_ref] = idx.other_gene[gs] < length
        novel = self.novel_entries()
        if len(novel) and in_ref.any():
            # reference peptides must also be absent from uploaded entries of another
            # accession (protein) / gene (gene)
            m_ref = np.zeros(len(d), dtype=bool)
            m_ref[sel[in_ref]] = True
            pairs = peptide_matches(d, novel, m_ref)
            ref_pos = np.flatnonzero(in_ref)[pairs["query"].to_numpy()]
            q_meta = d.proteins.meta
            q_acc = q_meta["accession"].to_numpy()[entry[ref_pos]]
            q_gene = _gene_keys(q_meta)[entry[ref_pos]]
            n_acc = novel.meta["accession"].to_numpy()[pairs["protein"].to_numpy()]
            n_gene = _gene_keys(novel.meta)[pairs["protein"].to_numpy()]
            up[ref_pos[n_acc != q_acc]] = False
            ug[ref_pos[n_gene != q_gene]] = False
        if (~in_ref).any():
            # peptides of entries outside the index: live mapping against the whole background
            m_new = np.zeros(len(d), dtype=bool)
            m_new[sel[~in_ref]] = True
            mp = peptide_mapping(d, self.proteins, m_new)
            up[~in_ref] = mp["unique_protein"].to_numpy()
            ug[~in_ref] = mp["unique_gene"].to_numpy()

        def ones_or_na(flags: np.ndarray) -> pd.Series:
            return pd.Series(np.ones(len(flags), dtype=np.int64), dtype="Int64").where(flags)

        return pd.DataFrame(
            {
                "n_proteins": ones_or_na(up),
                "n_genes": ones_or_na(ug),
                "unique_protein": up,
                "unique_gene": ug,
            }
        )


def _gene_keys(meta: pd.DataFrame) -> np.ndarray:
    """Gene name per entry, or ``acc:<accession>`` when the entry has none."""
    return (
        meta["gene"].where(meta["gene"].notna(), "acc:" + meta["accession"].astype(str))
    ).to_numpy()


def uniqueness_background(
    reference: ProteinSet, *uploaded: ProteinSet, method: str = "substring"
) -> Background:
    """Reference database plus uploaded FASTA(s); exact duplicates (same accession and
    sequence) are dropped."""
    return Background(ProteinSet.concat([reference, *uploaded]), method)


# =========================================================================================
# Per-digest analysis: outputs 1-4
# =========================================================================================
OUTPUT1_COLUMNS = [
    "accession", "gene", "protein_name", "digest", "peptide", "prev_aa", "next_aa", "start",
    "end", "start_mature", "end_mature", "length", "missed_cleavages", "mono_mass", "avg_mass",
    "mz_2", "mz_3", "mz_4", "n_K", "K_positions", "unique_protein", "unique_gene",
]  # fmt: skip
OUTPUT3_COLUMNS = [
    "accession", "gene", "protein_name", "length", "n_identifiable_peptides",
    "n_unique_peptides", "n_total_peptides", "max_coverage_pct", "unique_coverage_pct",
    "confident",
]  # fmt: skip
SUMMARY_METRICS = [
    "n_identifiable_proteins", "n_identifiable_peptides", "n_unique_peptides",
    "n_confident_proteins", "median_coverage_pct", "overall_coverage_pct",
]  # fmt: skip
#: Residues highlighted in amino-acid detectability and tracked as sites.
HIGHLIGHT_RESIDUES = "KSTYC"
SITE_RESIDUES = "KSTY"
#: Keys of :attr:`Comparison.counts`, in order.
COMPARISON_COUNTS = [
    "rescued_proteins", "lost_proteins", "proteins_gaining_unique", "residues_gained",
    "residues_lost", "peptides_gained", "peptides_lost",
    *(f"{r}_sites_{s}" for r in SITE_RESIDUES for s in ("gained", "lost")),
]  # fmt: skip


class Analysis:
    """Outputs 1-4 for one digest (single, sequential or parallel).

    Peptide counts count peptides by position (entry, start, end), as Protein Cleaver does:
    with missed cleavages, overlapping peptides count separately. Coverage is over the mature
    sequence (initiator Met excluded when removed).

    Args:
        d: The digest.
        background: Uniqueness background; defaults to the digested entries alone.
    """

    def __init__(self, d: Digest, background: Background | None = None):
        self.digest = d
        self.background = background or Background(d.proteins)

    @property
    def proteins(self) -> ProteinSet:
        """The digested entries."""
        return self.digest.proteins

    @property
    def label(self) -> str:
        """Digest label."""
        return self.digest.label

    # ----- peptide-level arrays (length = all cleavage products) -----------------------------
    @cached_property
    def identifiable(self) -> np.ndarray:
        """Mask of identifiable peptides."""
        return self.digest.identifiable

    @cached_property
    def mapping(self) -> pd.DataFrame:
        """Mapping counts for the identifiable peptides, in order."""
        return self.background.mapping(self.digest, self.identifiable)

    def _expand(self, column: str) -> np.ndarray:
        out = np.zeros(len(self.digest), dtype=bool)
        out[self.identifiable] = self.mapping[column].to_numpy()
        return out

    @cached_property
    def unique(self) -> np.ndarray:
        """Mask of identifiable peptides that map to exactly one entry."""
        return self._expand("unique_protein")

    @cached_property
    def unique_gene(self) -> np.ndarray:
        """Mask of identifiable peptides that map to exactly one gene."""
        return self._expand("unique_gene")

    @cached_property
    def covered(self) -> np.ndarray:
        """Residues (concatenated sequences) on >=1 identifiable peptide."""
        return covered_residues(self.digest, self.identifiable)

    @cached_property
    def covered_unique(self) -> np.ndarray:
        """Residues on >=1 unique identifiable peptide."""
        return covered_residues(self.digest, self.unique)

    # ----- Output 1 and 2 -----------------------------------------------------------------
    def peptide_table(self, unique_only: bool = False) -> pd.DataFrame:
        """Output 1 (identifiable peptides) or, with ``unique_only``, Output 2.

        For a parallel digest, ``digest`` names the digest(s) that produce each peptide.
        """
        df = self.digest.to_frame(identifiable_only=True)
        df["unique_protein"] = self.mapping["unique_protein"].to_numpy()
        df["unique_gene"] = self.mapping["unique_gene"].to_numpy()
        if len(self.digest.settings) > 1:
            origin = self.digest.origin[self.identifiable]
            labels = [s.label for s in self.digest.settings]
            names = {
                o: "; ".join(lab for i, lab in enumerate(labels) if (o >> i) & 1)
                for o in np.unique(origin).tolist()
            }
            df["digest"] = [names[o] for o in origin.tolist()]
        df = df[OUTPUT1_COLUMNS]
        return df[df["unique_protein"]].reset_index(drop=True) if unique_only else df

    def unique_peptide_table(self) -> pd.DataFrame:
        """Output 2: identifiable peptides with ``unique_protein == True``."""
        return self.peptide_table(unique_only=True)

    # ----- Output 3 -----------------------------------------------------------------------
    @cached_property
    def protein_stats(self) -> pd.DataFrame:
        """Per-entry statistics for every entry (identifiable or not)."""
        d = self.digest
        n = len(self.proteins)
        offs = d.offsets
        lengths = np.diff(offs)
        mature = lengths - d.mature_offset
        n_ident = bincount(d.protein[self.identifiable], minlength=n)
        n_uniq = bincount(d.protein[self.unique], minlength=n)
        n_total = bincount(d.protein[d.missed == 0], minlength=n)
        cov = _per_entry_sum(self.covered, offs)
        ucov = _per_entry_sum(self.covered_unique, offs)
        with np.errstate(invalid="ignore", divide="ignore"):
            pct = np.where(mature > 0, 100.0 * cov / mature, 0.0)
            upct = np.where(mature > 0, 100.0 * ucov / mature, 0.0)
        meta = self.proteins.meta
        return pd.DataFrame(
            {
                "accession": meta["accession"].to_numpy(),
                "gene": meta["gene"].to_numpy(),
                "protein_name": meta["protein_name"].to_numpy(),
                "length": lengths,
                "mature_length": mature,
                "n_identifiable_peptides": n_ident,
                "n_unique_peptides": n_uniq,
                "n_total_peptides": n_total,
                "covered_residues": cov,
                "max_coverage_pct": pct,
                "unique_covered_residues": ucov,
                "unique_coverage_pct": upct,
                "identifiable": n_ident >= 1,
                "confident": n_uniq >= 2,
            }
        )

    def protein_table(self, identifiable_only: bool = True) -> pd.DataFrame:
        """Output 3: identifiable proteins (``identifiable_only=False``: every entry)."""
        st = self.protein_stats
        if identifiable_only:
            st = st[st["identifiable"]]
        return st[OUTPUT3_COLUMNS].reset_index(drop=True)

    # ----- Summary -------------------------------------------------------------------------
    @cached_property
    def summary(self) -> dict[str, float]:
        """Headline metrics (see :data:`SUMMARY_METRICS`) plus ``n_total_peptides``.

        ``median_coverage_pct`` is the median over all entries (non-identifiable entries count
        as 0 %); ``overall_coverage_pct`` = covered residues / mature residues, all entries.
        """
        st = self.protein_stats
        mature = int(st["mature_length"].sum())
        return {
            "n_identifiable_proteins": int(st["identifiable"].sum()),
            "n_identifiable_peptides": int(self.identifiable.sum()),
            "n_unique_peptides": int(self.unique.sum()),
            "n_confident_proteins": int(st["confident"].sum()),
            "median_coverage_pct": float(st["max_coverage_pct"].median()),
            "overall_coverage_pct": 100.0 * float(st["covered_residues"].sum()) / mature
            if mature
            else 0.0,
            "n_total_peptides": int(self.digest.n_total_peptides),
        }

    # ----- Output 4 -----------------------------------------------------------------------
    def length_histogram(self) -> pd.DataFrame:
        """Peptide length frequency of all cleavage products (0 missed cleavages).

        ``in_window`` marks lengths inside the identifiable length window (any component
        digest's window for a parallel digest).
        """
        d = self.digest
        lengths = d.length[d.missed == 0]
        counts = bincount(lengths) if len(lengths) else np.zeros(1, dtype=np.int64)
        length = np.arange(len(counts))
        in_window = np.zeros(len(counts), dtype=bool)
        for s in d.settings:
            in_window |= (length >= s.min_length) & (length <= s.max_length)
        df = pd.DataFrame({"length": length, "n_peptides": counts, "in_window": in_window})
        return df[df["length"] >= 1].reset_index(drop=True)

    def coverage_summary(self) -> dict[str, float]:
        """Proteome-wide coverage: overall % of mature residues on identifiable peptides,
        and the median / mean of per-entry coverage."""
        st = self.protein_stats
        return {
            "overall_coverage_pct": self.summary["overall_coverage_pct"],
            "median_coverage_pct": float(st["max_coverage_pct"].median()),
            "mean_coverage_pct": float(st["max_coverage_pct"].mean()),
            "n_entries": len(st),
        }

    def aa_detectability(self) -> pd.DataFrame:
        """For each residue type: occurrences in mature sequences, occurrences on identifiable
        peptides, and the fraction. ``highlight`` marks K, S, T, Y and C."""
        cat = _concat(self.proteins)
        in_mature = np.ones(len(cat.codes), dtype=bool)
        in_mature[cat.offsets[:-1][self.digest.mature_offset == 1]] = False
        total = bincount(cat.codes[in_mature], minlength=256)
        cov = bincount(cat.codes[in_mature & self.covered], minlength=256)
        present = np.flatnonzero(total)
        res = [chr(c) for c in present]
        return pd.DataFrame(
            {
                "residue": res,
                "occurrences": total[present],
                "on_identifiable_peptides": cov[present],
                "fraction_detectable": cov[present] / total[present],
                "highlight": [r in HIGHLIGHT_RESIDUES for r in res],
            }
        )

    def digest_statistics(self) -> dict[str, object]:
        """Output 4: ``length_histogram``, ``coverage_summary``, ``coverage_per_protein``
        and ``aa_detectability``."""
        return {
            "length_histogram": self.length_histogram(),
            "coverage_summary": self.coverage_summary(),
            "coverage_per_protein": self.protein_stats[
                ["accession", "gene", "mature_length", "covered_residues", "max_coverage_pct"]
            ],
            "aa_detectability": self.aa_detectability(),
        }


def _per_entry_sum(mask: np.ndarray, offsets: np.ndarray) -> np.ndarray:
    csum = np.concatenate(([0], np.cumsum(mask, dtype=np.int64)))
    return csum[offsets[1:]] - csum[offsets[:-1]]


def analyze(d: Digest, background: Background | None = None) -> Analysis:
    """Convenience constructor for :class:`Analysis`."""
    return Analysis(d, background)


# =========================================================================================
# Output 5: comparison against the baseline
# =========================================================================================
@dataclass
class Comparison:
    """Seen vs unseen between a baseline and an alternative digest (Output 5).

    Attributes:
        baseline, alternative: The two analyses (same ProteinSet object).
        summary: Metric table: ``metric``, ``baseline``, ``alternative``, ``difference``.
        rescued_proteins: Not identifiable in the baseline, identifiable in the alternative.
        lost_proteins: Identifiable in the baseline, not in the alternative.
        proteins_gaining_unique: Entries with more unique peptides in the alternative.
        residues: Per entry with changes: residues gained / lost and their UniProt ranges.
        sites: K, S, T and Y residues gained or lost (covered by identifiable peptides).
        peptides_gained, peptides_lost: Identifiable peptides (by position) only in one digest.
        per_protein: Coverage of every entry in both digests (for the scatter plot).
        counts: Headline counts of the lists above.
    """

    baseline: Analysis
    alternative: Analysis
    summary: pd.DataFrame
    rescued_proteins: pd.DataFrame
    lost_proteins: pd.DataFrame
    proteins_gaining_unique: pd.DataFrame
    residues: pd.DataFrame
    sites: pd.DataFrame
    peptides_gained: int
    peptides_lost: int
    per_protein: pd.DataFrame
    counts: dict[str, int]


def compare(baseline: Analysis, alternative: Analysis) -> Comparison:
    """Compare ``alternative`` against ``baseline`` (both on the same ProteinSet object)."""
    if baseline.proteins is not alternative.proteins:
        raise ValueError("baseline and alternative must digest the same ProteinSet object")
    b, a = baseline.protein_stats, alternative.protein_stats
    meta = b[["accession", "gene", "protein_name"]]

    summary = pd.DataFrame(
        {
            "metric": SUMMARY_METRICS,
            "baseline": [baseline.summary[m] for m in SUMMARY_METRICS],
            "alternative": [alternative.summary[m] for m in SUMMARY_METRICS],
        }
    )
    summary["difference"] = summary["alternative"] - summary["baseline"]

    rescued = ~b["identifiable"] & a["identifiable"]
    lost = b["identifiable"] & ~a["identifiable"]
    gaining = a["n_unique_peptides"] > b["n_unique_peptides"]

    def _ptable(sel: pd.Series) -> pd.DataFrame:
        out = meta[sel].copy()
        out["n_identifiable_baseline"] = b.loc[sel, "n_identifiable_peptides"]
        out["n_identifiable_alternative"] = a.loc[sel, "n_identifiable_peptides"]
        out["coverage_pct_baseline"] = b.loc[sel, "max_coverage_pct"]
        out["coverage_pct_alternative"] = a.loc[sel, "max_coverage_pct"]
        return out.reset_index(drop=True)

    gain_tab = meta[gaining].copy()
    gain_tab["n_unique_baseline"] = b.loc[gaining, "n_unique_peptides"]
    gain_tab["n_unique_alternative"] = a.loc[gaining, "n_unique_peptides"]
    gain_tab["gain"] = gain_tab["n_unique_alternative"] - gain_tab["n_unique_baseline"]
    gain_tab = gain_tab.reset_index(drop=True)

    residues, sites = _residue_changes(baseline, alternative)

    bset = _peptide_keys(baseline)
    aset = _peptide_keys(alternative)
    per_protein = meta.copy()
    per_protein["coverage_pct_baseline"] = b["max_coverage_pct"]
    per_protein["coverage_pct_alternative"] = a["max_coverage_pct"]
    per_protein["identifiable_baseline"] = b["identifiable"]
    per_protein["identifiable_alternative"] = a["identifiable"]

    status = sites["status"] if len(sites) else pd.Series(dtype=object)
    counts = {
        "rescued_proteins": int(rescued.sum()),
        "lost_proteins": int(lost.sum()),
        "proteins_gaining_unique": int(gaining.sum()),
        "residues_gained": int(residues["residues_gained"].sum()) if len(residues) else 0,
        "residues_lost": int(residues["residues_lost"].sum()) if len(residues) else 0,
        "peptides_gained": int((~np.isin(aset, bset)).sum()),
        "peptides_lost": int((~np.isin(bset, aset)).sum()),
    }
    for r in SITE_RESIDUES:
        for st in ("gained", "lost"):
            counts[f"{r}_sites_{st}"] = (
                int(((sites["residue"] == r) & (status == st)).sum()) if len(sites) else 0
            )
    return Comparison(
        baseline,
        alternative,
        summary,
        _ptable(rescued),
        _ptable(lost),
        gain_tab,
        residues,
        sites,
        counts["peptides_gained"],
        counts["peptides_lost"],
        per_protein,
        counts,
    )


def _peptide_keys(an: Analysis) -> np.ndarray:
    """Unique int64 keys of identifiable peptide spans (global start, end)."""
    d = an.digest
    m = an.identifiable
    span = int(d.offsets[-1]) + 1
    return np.unique(d.gstart[m].astype(np.int64) * span + d.gend[m])


def comparison_counts(baseline: Analysis, alternative: Analysis) -> dict[str, int]:
    """:attr:`Comparison.counts` computed without building the lists (fast, vectorized)."""
    per = comparison_per_protein(baseline, alternative)
    b, a = baseline.protein_stats, alternative.protein_stats
    bk, ak = _peptide_keys(baseline), _peptide_keys(alternative)
    out = {
        "rescued_proteins": int((~b["identifiable"] & a["identifiable"]).sum()),
        "lost_proteins": int((b["identifiable"] & ~a["identifiable"]).sum()),
        "proteins_gaining_unique": int((a["n_unique_peptides"] > b["n_unique_peptides"]).sum()),
        "peptides_gained": int((~np.isin(ak, bk)).sum()),
        "peptides_lost": int((~np.isin(bk, ak)).sum()),
    }
    out.update({c: int(per[c].sum()) for c in per.columns if c != "accession"})
    return {k: out[k] for k in COMPARISON_COUNTS}


def comparison_per_protein(baseline: Analysis, alternative: Analysis) -> pd.DataFrame:
    """Per entry: residues and K/S/T/Y sites gained and lost (identifiable-peptide coverage)."""
    if baseline.proteins is not alternative.proteins:
        raise ValueError("baseline and alternative must digest the same ProteinSet object")
    cat = _concat(baseline.proteins)
    offs = cat.offsets
    gained = alternative.covered & ~baseline.covered
    lost = baseline.covered & ~alternative.covered
    cols = {
        "accession": baseline.proteins.meta["accession"].to_numpy(),
        "residues_gained": _per_entry_sum(gained, offs),
        "residues_lost": _per_entry_sum(lost, offs),
    }
    for r in SITE_RESIDUES:
        is_r = cat.codes == ord(r)
        cols[f"{r}_sites_gained"] = _per_entry_sum(gained & is_r, offs)
        cols[f"{r}_sites_lost"] = _per_entry_sum(lost & is_r, offs)
    return pd.DataFrame(cols)


def _ranges(positions: np.ndarray) -> str:
    """'3-8;12;20-25' from sorted 1-based positions."""
    if not len(positions):
        return ""
    breaks = np.flatnonzero(np.diff(positions) != 1)
    starts = np.concatenate(([positions[0]], positions[breaks + 1]))
    ends = np.concatenate((positions[breaks], [positions[-1]]))
    return ";".join(str(s) if s == e else f"{s}-{e}" for s, e in zip(starts, ends, strict=True))


def _residue_changes(baseline: Analysis, alt: Analysis) -> tuple[pd.DataFrame, pd.DataFrame]:
    cat = _concat(baseline.proteins)
    offs = cat.offsets
    gained = alt.covered & ~baseline.covered
    lost = baseline.covered & ~alt.covered
    changed = gained | lost
    meta = baseline.proteins.meta
    prot_of = np.searchsorted(offs, np.flatnonzero(changed), side="right") - 1
    shift = baseline.digest.mature_offset

    rows = []
    for i in np.unique(prot_of).tolist():
        lo, hi = offs[i], offs[i + 1]
        g = np.flatnonzero(gained[lo:hi]) + 1
        lo_ = np.flatnonzero(lost[lo:hi]) + 1
        rows.append(
            {
                "accession": meta.at[i, "accession"],
                "gene": meta.at[i, "gene"],
                "residues_gained": len(g),
                "residues_lost": len(lo_),
                "gained_ranges": _ranges(g),
                "lost_ranges": _ranges(lo_),
            }
        )
    residues = pd.DataFrame(
        rows,
        columns=[
            "accession", "gene", "residues_gained", "residues_lost", "gained_ranges",
            "lost_ranges",
        ],
    )  # fmt: skip

    is_site = np.isin(cat.codes, np.frombuffer(SITE_RESIDUES.encode(), dtype=np.uint8))
    site_rows = []
    for status, mask in (("gained", gained), ("lost", lost)):
        pos = np.flatnonzero(mask & is_site)
        prot = np.searchsorted(offs, pos, side="right") - 1
        uni = pos - offs[prot] + 1
        site_rows.append(
            pd.DataFrame(
                {
                    "accession": meta["accession"].to_numpy()[prot],
                    "gene": meta["gene"].to_numpy()[prot],
                    "residue": [chr(c) for c in cat.codes[pos]],
                    "position": uni,
                    "position_mature": uni - shift[prot],
                    "status": status,
                }
            )
        )
    sites = pd.concat(site_rows, ignore_index=True)
    return residues, sites.sort_values(["accession", "position"], ignore_index=True)


def summary_table(baseline: Analysis, alternatives: Sequence[Analysis]) -> pd.DataFrame:
    """Output 5 summary: one row per digest (baseline first) with each metric and its
    difference from the baseline (``<metric>_diff``)."""
    rows = []
    for an in [baseline, *alternatives]:
        row = {"digest": an.label}
        for m in SUMMARY_METRICS:
            row[m] = an.summary[m]
            row[f"{m}_diff"] = an.summary[m] - baseline.summary[m]
        rows.append(row)
    return pd.DataFrame(rows)


def ranking_table(
    proteins: ProteinSet,
    baseline: DigestSettings,
    enzymes: Sequence[CleavageRule],
    background: Background | None = None,
) -> pd.DataFrame:
    """All-pairs ranking table (Output 5).

    Rows: the baseline, then for every enzyme whose rule differs from the baseline rules:
    the enzyme alone, the enzyme in parallel with the baseline, and the sequential digest
    baseline + enzyme (union of cut sites). Alternative digests use the baseline's length,
    missed-cleavage, m/z and chemistry settings. Sort the returned DataFrame by any column.
    """
    background = background or Background(proteins)
    base = Analysis(digest(proteins, baseline), background)
    rows = [_ranking_row("baseline", None, base, base)]
    for rule in enzymes:
        if rule in baseline.rules:
            continue
        alone_d = digest(proteins, baseline.with_rules(rule))
        variants = {
            "alone": alone_d,
            "parallel": parallel(base.digest, alone_d),
            "sequential": digest(proteins, baseline.with_rules(*baseline.rules, rule)),
        }
        for kind, d in variants.items():
            rows.append(_ranking_row(kind, rule, Analysis(d, background), base))
    return pd.DataFrame(rows)


def _ranking_row(kind: str, rule: CleavageRule | None, an: Analysis, base: Analysis) -> dict:
    row: dict[str, object] = {
        "digest": an.label,
        "kind": kind,
        "enzyme": rule.label if rule else "",
        "flags": ", ".join(rule.flags) if rule else "",
    }
    for m in SUMMARY_METRICS:
        row[m] = an.summary[m]
        row[f"{m}_diff"] = an.summary[m] - base.summary[m]
    row.update(dict.fromkeys(COMPARISON_COUNTS, 0) if an is base else comparison_counts(base, an))
    return row


# =========================================================================================
# Peptide map (protein viewer, CLAUDE.md section 7)
# =========================================================================================
def peptide_map(
    analyses: Sequence[Analysis],
    accession: str,
    window: tuple[int, int] | None = None,
) -> dict[str, pd.DataFrame]:
    """Data for the peptide map of one protein: one row per digest (baseline first).

    Args:
        analyses: Analyses on the same ProteinSet object; ``analyses[0]`` is the baseline.
        accession: Entry to show.
        window: Optional (first, last) residue window in **mature** numbering; peptides and
            strips are limited to it and ``coverage_window_pct`` is computed over it.

    Returns:
        ``peptides`` (digest, peptide, start/end UniProt and mature, length, n_K),
        ``lysines`` (position_uniprot, position_mature),
        ``strips`` (digest, status gained/lost vs the baseline, start_mature, end_mature),
        ``coverage`` (digest, coverage_pct over the mature sequence, coverage_window_pct).
    """
    base = analyses[0]
    proteins = base.proteins
    for an in analyses[1:]:
        if an.proteins is not proteins:
            raise ValueError("all analyses must use the same ProteinSet object")
    idx = proteins.index_of(accession)
    seq = proteins.sequences[idx]
    shift = int(base.digest.mature_offset[idx])
    lo, hi = (1, len(seq) - shift) if window is None else window
    in_win = np.zeros(len(seq), dtype=bool)  # by 0-based UniProt index
    in_win[max(lo + shift - 1, 0) : min(hi + shift, len(seq))] = True
    off = int(_concat(proteins).offsets[idx])

    pep_rows, strip_rows, cov_rows = [], [], []
    base_cov = base.covered[off : off + len(seq)]
    for an in analyses:
        d = an.digest
        sel = an.identifiable & (d.protein == idx)
        for s, e in zip(d.start0[sel].tolist(), d.end0[sel].tolist(), strict=True):
            if not in_win[s:e].any():
                continue
            p = seq[s:e]
            pep_rows.append(
                {
                    "digest": an.label,
                    "peptide": p,
                    "start": s + 1,
                    "end": e,
                    "start_mature": s + 1 - shift,
                    "end_mature": e - shift,
                    "length": e - s,
                    "n_K": p.count("K"),
                }
            )
        cov = an.covered[off : off + len(seq)]
        if an is not base:
            for status, mask in (("gained", cov & ~base_cov), ("lost", base_cov & ~cov)):
                pos = np.flatnonzero(mask & in_win) + 1 - shift
                for r in _ranges(pos).split(";") if len(pos) else []:
                    a, _, b = r.partition("-")
                    strip_rows.append(
                        {
                            "digest": an.label,
                            "status": status,
                            "start_mature": int(a),
                            "end_mature": int(b or a),
                        }
                    )
        mature = len(seq) - shift
        cov_rows.append(
            {
                "digest": an.label,
                "coverage_pct": 100.0 * cov[shift:].sum() / mature if mature else 0.0,
                "coverage_window_pct": 100.0 * (cov & in_win).sum() / in_win.sum()
                if in_win.any()
                else 0.0,
            }
        )
    lys = [i + 1 for i, c in enumerate(seq) if c == "K" and in_win[i]]
    return {
        "peptides": pd.DataFrame(pep_rows),
        "lysines": pd.DataFrame(
            {"position_uniprot": lys, "position_mature": [p - shift for p in lys]}
        ),
        "strips": pd.DataFrame(
            strip_rows, columns=["digest", "status", "start_mature", "end_mature"]
        ),
        "coverage": pd.DataFrame(cov_rows),
    }


# =========================================================================================
# Precomputed uniqueness index (shared-prefix lengths)
# =========================================================================================
#: Largest shared-prefix length stored (uint8); peptides longer than this are treated as
#: shared whenever the stored value is 255.
SHARED_PREFIX_CAP = 255


def _dense_rank(key: np.ndarray) -> np.ndarray:
    """1-based dense ranks of ``key`` (int32)."""
    order = np.argsort(key)
    s = key[order]
    new = np.empty(len(s), dtype=bool)
    new[:1] = True
    np.not_equal(s[1:], s[:-1], out=new[1:])
    rank = np.empty(len(key), dtype=np.int32)
    rank[order] = np.cumsum(new, dtype=np.int32)
    return rank


def _segmented_cummin(values: np.ndarray, segment: np.ndarray) -> np.ndarray:
    """Running minimum restarting at every new ``segment`` id (ids non-decreasing)."""
    big = int(values.max(initial=0)) + 1
    shifted = values.astype(np.int64) - segment.astype(np.int64) * big
    return np.minimum.accumulate(shifted) + segment.astype(np.int64) * big


@dataclass(eq=False)
class SharedPrefixIndex:
    """For every residue of a reference ProteinSet, the length of the longest prefix starting
    there that also occurs (I = L) in another entry (``other_entry``) or in an entry of another
    gene (``other_gene``), capped at :data:`SHARED_PREFIX_CAP`.

    A peptide of length L starting at residue s of reference entry e is protein-unique iff
    ``other_entry[s] < L`` and gene-unique iff ``other_gene[s] < L``. This is exactly the
    substring mapping definition, for any enzyme or settings.

    Attributes:
        proteins: The reference entries the index was built on.
        other_entry, other_gene: uint8 arrays over the concatenated reference sequences.
    """

    proteins: ProteinSet
    other_entry: np.ndarray
    other_gene: np.ndarray

    @classmethod
    def build(cls, proteins: ProteinSet) -> SharedPrefixIndex:
        """Build by prefix doubling over the entries joined with unique sentinels."""
        cat = _concat(proteins)
        n_ent = len(proteins)
        lengths = np.diff(cat.offsets)
        n_res = int(cat.offsets[-1])
        n = n_res + n_ent
        # text = entry residues followed by a sentinel unique to the entry
        res_text_pos = np.arange(n_res, dtype=np.int64) + np.repeat(
            np.arange(n_ent, dtype=np.int64), lengths.astype(np.intp)
        )
        text = np.empty(n, dtype=np.int32)
        sentinel_pos = cat.offsets[1:] + np.arange(n_ent)
        text[sentinel_pos] = 27 + np.arange(n_ent, dtype=np.int32)
        text[res_text_pos] = _FOLD[cat.codes]
        ranks = {1: text}
        k = 1
        while k < SHARED_PREFIX_CAP + 1:
            r = ranks[k]
            nxt = np.zeros(n, dtype=np.int64)
            nxt[: n - k] = r[k:]
            ranks[2 * k] = _dense_rank(r.astype(np.int64) * (int(r.max()) + 1) + nxt)
            k *= 2
        sa = np.argsort(ranks[k], kind="stable")
        del ranks[k]
        # adjacent LCP (capped) by binary lifting over the power-of-two ranks
        a, b = sa[:-1], sa[1:]
        lcp_adj = np.zeros(n - 1, dtype=np.int32)
        for step in sorted(ranks, reverse=True):
            r = ranks[step]
            ia, ib = a + lcp_adj, b + lcp_adj
            ok = (ia < n) & (ib < n)
            eq = np.zeros(n - 1, dtype=bool)
            eq[ok] = r[ia[ok]] == r[ib[ok]]
            lcp_adj += np.where(eq, step, 0).astype(np.int32)
        del ranks
        lcp_adj = np.minimum(lcp_adj, SHARED_PREFIX_CAP)
        owner = np.empty(n, dtype=np.int64)
        owner[res_text_pos] = np.repeat(np.arange(n_ent), lengths.astype(np.intp))
        owner[sentinel_pos] = np.arange(n_ent)
        meta = proteins.meta
        gene = meta["gene"].where(meta["gene"].notna(), "acc:" + meta["accession"].astype(str))
        gene_id = pd.factorize(gene)[0]
        acc_id = pd.factorize(meta["accession"])[0]
        out = []
        for group in (acc_id[owner], gene_id[owner]):
            g = group[sa]
            new_run = np.empty(n, dtype=bool)
            new_run[0] = True
            np.not_equal(g[1:], g[:-1], out=new_run[1:])
            seg = np.cumsum(new_run) - 1
            # left: min LCP back to the nearest suffix of another group (0 for the first run)
            left_vals = np.concatenate(([0], lcp_adj))
            left_vals[seg == 0] = 0
            left = _segmented_cummin(left_vals, seg)
            # right: min LCP forward to the nearest suffix of another group (0 for last run)
            right_vals = np.concatenate((lcp_adj, [0]))
            right_vals[seg == seg[-1]] = 0
            rseg = (seg[-1] - seg)[::-1]
            right = _segmented_cummin(right_vals[::-1], rseg)[::-1]
            m_sa = np.maximum(left, right)
            m_text = np.empty(n, dtype=np.uint8)
            m_text[sa] = m_sa.astype(np.uint8)
            out.append(m_text[res_text_pos])
        return cls(proteins, out[0], out[1])

    # ----- lookup -------------------------------------------------------------------------
    def unique_flags(self, d: Digest, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Protein- and gene-uniqueness of ``d``'s selected peptides.

        ``d`` must digest the indexed ProteinSet itself, or a set whose entries all occur in
        it with identical accession and sequence (e.g. the histone FASTA in the proteome).
        Raises ``KeyError`` for entries absent from the reference.
        """
        if d.proteins is self.proteins:
            gstart = d.gstart[mask]
        else:
            ref_off = self.offsets_for(d.proteins)
            gstart = ref_off[d.protein[mask]] + d.start0[mask]
        length = d.length[mask]
        return self.other_entry[gstart] < length, self.other_gene[gstart] < length

    def offsets_for(self, query: ProteinSet, strict: bool = True) -> np.ndarray:
        """Reference offset of each query entry (same accession and sequence). Missing
        entries raise ``KeyError`` (``strict``) or get -1."""
        if query is self.proteins:
            return _concat(self.proteins).offsets[:-1].copy()
        ref_seq = dict(zip(self.proteins.accessions, self.proteins.sequences, strict=True))
        ref_idx = {a: i for i, a in enumerate(self.proteins.accessions)}
        offs = _concat(self.proteins).offsets
        out = np.full(len(query), -1, dtype=np.int64)
        for i, (acc, seq) in enumerate(zip(query.accessions, query.sequences, strict=True)):
            if ref_seq.get(acc) == seq:
                out[i] = offs[ref_idx[acc]]
            elif strict:
                raise KeyError(f"{acc} is not in the indexed reference with this sequence")
        return out

    def save(self, path: str | Path) -> None:
        """Save both arrays (compressed npz) with the accession order for validation."""
        np.savez_compressed(
            path,
            other_entry=self.other_entry,
            other_gene=self.other_gene,
            accessions=np.array(self.proteins.accessions),
            lengths=self.proteins.lengths,
        )

    @classmethod
    def load(cls, path: str | Path, proteins: ProteinSet) -> SharedPrefixIndex:
        """Load an index saved by :meth:`save` for ``proteins`` (same entries, same order)."""
        z = np.load(path)
        if list(z["accessions"]) != proteins.accessions or not np.array_equal(
            z["lengths"], proteins.lengths
        ):
            raise ValueError("index does not match the given reference entries")
        return cls(proteins, z["other_entry"], z["other_gene"])
