"""One interface over precomputed and live results, used by the app.

``Results`` exposes per-digest summaries, per-entry statistics and changes versus the
baseline, Output 4 data, the ranking table, the "seen vs unseen" lists, and live
peptide-level detail for selected entries. :class:`PrecomputedResults` reads the
precomputed arrays (whole proteome, preset settings); :class:`LiveResults` digests in place
(uploaded FASTA or changed settings). Both return identical tables for the same input; the
tests check this.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd

from .digest import Digest, mature_offsets
from .fasta import ProteinSet
from .metrics import (
    COMPARISON_COUNTS,
    HIGHLIGHT_RESIDUES,
    SUMMARY_METRICS,
    Analysis,
    Background,
    comparison_counts,
    comparison_per_protein,
)
from .precomputed import CHANGE_FIELDS, HIST_MAX, PrecomputedPreset
from .presets import DigestRecipe
from .report import header_lines

ENTRY_STAT_COLUMNS = [
    "accession", "gene", "protein_name", "length", "mature_length", "n_identifiable_peptides",
    "n_unique_peptides", "n_total_peptides", "covered_residues", "max_coverage_pct",
    "unique_covered_residues", "unique_coverage_pct", "identifiable", "confident",
]  # fmt: skip
Progress = Callable[[int, int, str], None]


class Results:
    """Shared behaviour; subclasses provide the per-digest primitives."""

    mode: str
    proteins: ProteinSet
    background: Background
    recipes: list[DigestRecipe]

    # ----- primitives (subclasses) ------------------------------------------------------
    def summary(self, i: int) -> dict[str, float]:  # pragma: no cover - abstract
        raise NotImplementedError

    def comparison_counts(self, i: int) -> dict[str, int]:  # pragma: no cover
        raise NotImplementedError

    def entry_stats(self, i: int) -> pd.DataFrame:  # pragma: no cover
        raise NotImplementedError

    def entry_changes(self, i: int) -> pd.DataFrame:  # pragma: no cover
        raise NotImplementedError

    def _hist_counts(self, i: int) -> np.ndarray:  # pragma: no cover
        raise NotImplementedError

    def _aa_counts(self, i: int) -> pd.DataFrame:  # pragma: no cover
        raise NotImplementedError

    # ----- shared --------------------------------------------------------------------------
    @property
    def labels(self) -> list[str]:
        """Digest labels, baseline first."""
        return [r.label for r in self.recipes]

    def header(self, title: str, indices: Sequence[int] = (0,)) -> list[str]:
        """Export header for a table computed from the given digests."""
        return header_lines(
            [self.recipes[i].settings for i in indices],
            str(self.proteins.name or f"{len(self.proteins)} entries"),
            self.background.description,
            title,
        )

    def length_histogram(self, i: int) -> pd.DataFrame:
        """Lengths of all cleavage products (0 missed cleavages); last bin = > HIST_MAX aa."""
        counts = self._hist_counts(i)
        lengths = np.arange(1, HIST_MAX + 2)
        in_window = np.zeros(len(lengths), dtype=bool)
        for s in self.recipes[i].settings:
            in_window |= (lengths >= s.min_length) & (lengths <= s.max_length)
        return pd.DataFrame({
            "length": [str(x) for x in lengths[:-1]] + [f">{HIST_MAX}"],
            "length_value": lengths,
            "n_peptides": counts,
            "in_window": in_window,
        })  # fmt: skip

    def aa_detectability(self, i: int) -> pd.DataFrame:
        """Output 4 amino-acid detectability."""
        df = self._aa_counts(i).copy()
        df["fraction_detectable"] = df["on_identifiable_peptides"] / df["occurrences"]
        df["highlight"] = df["residue"].isin(list(HIGHLIGHT_RESIDUES))
        return df

    def summary_table(self, indices: Sequence[int]) -> pd.DataFrame:
        """Output 5 summary: metrics and differences from the baseline, baseline first."""
        base = self.summary(0)
        rows = []
        for i in [0, *[j for j in indices if j != 0]]:
            s = self.summary(i)
            row: dict[str, object] = {"digest": self.recipes[i].label,
                                      "kind": self.recipes[i].kind}  # fmt: skip
            for m in SUMMARY_METRICS:
                row[m] = s[m]
                row[f"{m}_diff"] = s[m] - base[m]
            rows.append(row)
        return pd.DataFrame(rows)

    def ranking_table(self, progress: Progress | None = None) -> pd.DataFrame:
        """All-pairs ranking table over every recipe."""
        base = self.summary(0)
        rows = []
        for i, r in enumerate(self.recipes):
            if progress:
                progress(i, len(self.recipes), r.label)
            s = self.summary(i)
            row: dict[str, object] = {
                "digest": r.label,
                "kind": r.kind,
                "enzyme": r.enzyme,
                "flags": ", ".join(r.rule.flags) if r.rule is not None else "",
            }
            for m in SUMMARY_METRICS:
                row[m] = s[m]
                row[f"{m}_diff"] = s[m] - base[m]
            row.update(self.comparison_counts(i))
            rows.append(row)
        return pd.DataFrame(rows)

    def seen_vs_unseen(self, i: int) -> dict[str, pd.DataFrame]:
        """Output 5 lists for digest ``i`` versus the baseline (per-entry counts)."""
        b, a = self.entry_stats(0), self.entry_stats(i)
        meta = b[["accession", "gene", "protein_name"]]

        def ptable(sel: pd.Series) -> pd.DataFrame:
            out = meta[sel].copy()
            out["n_identifiable_baseline"] = b.loc[sel, "n_identifiable_peptides"]
            out["n_identifiable_alternative"] = a.loc[sel, "n_identifiable_peptides"]
            out["coverage_pct_baseline"] = b.loc[sel, "max_coverage_pct"]
            out["coverage_pct_alternative"] = a.loc[sel, "max_coverage_pct"]
            return out.reset_index(drop=True)

        gaining = a["n_unique_peptides"] > b["n_unique_peptides"]
        gain = meta[gaining].copy()
        gain["n_unique_baseline"] = b.loc[gaining, "n_unique_peptides"]
        gain["n_unique_alternative"] = a.loc[gaining, "n_unique_peptides"]
        gain["gain"] = gain["n_unique_alternative"] - gain["n_unique_baseline"]
        ch = self.entry_changes(i)
        changed = ch[CHANGE_FIELDS].sum(axis=1) > 0
        changes = pd.concat([meta[changed], ch.loc[changed, CHANGE_FIELDS]], axis=1)
        return {
            "rescued_proteins": ptable(~b["identifiable"] & a["identifiable"]),
            "lost_proteins": ptable(b["identifiable"] & ~a["identifiable"]),
            "proteins_gaining_unique": gain.reset_index(drop=True),
            "residue_and_site_changes": changes.reset_index(drop=True),
        }

    def coverage_pair(self, i: int) -> pd.DataFrame:
        """Per-entry baseline vs alternative coverage (scatter data)."""
        b, a = self.entry_stats(0), self.entry_stats(i)
        return pd.DataFrame({
            "accession": b["accession"], "gene": b["gene"], "protein_name": b["protein_name"],
            "coverage_pct_baseline": b["max_coverage_pct"],
            "coverage_pct_alternative": a["max_coverage_pct"],
        })  # fmt: skip

    # ----- live peptide-level detail -------------------------------------------------------
    def protein_analyses(self, accession: str, indices: Sequence[int]) -> list[Analysis]:
        """Analyses of one entry (digested alone) for the given recipes, in order."""
        sub = self.proteins.subset([self.proteins.index_of(accession)])
        sub.name = self.proteins.name
        cache: dict = {}
        return [Analysis(self.recipes[i].run(sub, cache), self.background) for i in indices]

    def peptide_table(self, i: int, accessions: Sequence[str]) -> pd.DataFrame:
        """Output 1 for the given entries under digest ``i`` (computed live)."""
        idx = [self.proteins.index_of(a) for a in accessions]
        sub = self.proteins.subset(idx)
        sub.name = self.proteins.name
        return Analysis(self.recipes[i].run(sub), self.background).peptide_table()


def _stats_frame(proteins: ProteinSet, mature: np.ndarray, n_ident, n_uniq, n_total, cov,
                 ucov) -> pd.DataFrame:  # fmt: skip
    lengths = proteins.lengths
    with np.errstate(invalid="ignore", divide="ignore"):
        pct = np.where(mature > 0, 100.0 * cov / np.maximum(mature, 1), 0.0)
        upct = np.where(mature > 0, 100.0 * ucov / np.maximum(mature, 1), 0.0)
    meta = proteins.meta
    return pd.DataFrame({
        "accession": meta["accession"].to_numpy(),
        "gene": meta["gene"].to_numpy(),
        "protein_name": meta["protein_name"].to_numpy(),
        "length": lengths,
        "mature_length": mature,
        "n_identifiable_peptides": np.asarray(n_ident, dtype=np.int64),
        "n_unique_peptides": np.asarray(n_uniq, dtype=np.int64),
        "n_total_peptides": np.asarray(n_total, dtype=np.int64),
        "covered_residues": np.asarray(cov, dtype=np.int64),
        "max_coverage_pct": pct,
        "unique_covered_residues": np.asarray(ucov, dtype=np.int64),
        "unique_coverage_pct": upct,
        "identifiable": np.asarray(n_ident) >= 1,
        "confident": np.asarray(n_uniq) >= 2,
    })[ENTRY_STAT_COLUMNS]  # fmt: skip


class PrecomputedResults(Results):
    """Results read from a precomputed preset (whole reference, preset settings)."""

    mode = "precomputed"

    def __init__(self, reference: ProteinSet, background: Background,
                 recipes: list[DigestRecipe], data: PrecomputedPreset):  # fmt: skip
        labels = [d["label"] for d in data.info["digests"]]
        if labels != [r.label for r in recipes]:
            raise ValueError("precomputed digests do not match the preset recipes")
        self.proteins = reference
        self.background = background
        self.recipes = recipes
        self.data = data
        self._mature = reference.lengths - mature_offsets(reference, recipes[0].baseline.chemistry)
        self._stats: dict[int, pd.DataFrame] = {}

    def summary(self, i: int) -> dict[str, float]:
        return dict(self.data.info["digests"][i]["summary"])

    def comparison_counts(self, i: int) -> dict[str, int]:
        return {k: int(self.data.info["digests"][i]["comparison"][k]) for k in COMPARISON_COUNTS}

    def entry_stats(self, i: int) -> pd.DataFrame:
        if i not in self._stats:
            a = self.data.arrays
            self._stats[i] = _stats_frame(
                self.proteins, self._mature, a["n_identifiable_peptides"][i],
                a["n_unique_peptides"][i], a["n_total_peptides"][i], a["covered_residues"][i],
                a["unique_covered_residues"][i],
            )  # fmt: skip
        return self._stats[i]

    def entry_changes(self, i: int) -> pd.DataFrame:
        df = pd.DataFrame({f: self.data.arrays[f][i].astype(np.int64) for f in CHANGE_FIELDS})
        df.insert(0, "accession", self.proteins.accessions)
        return df

    def _hist_counts(self, i: int) -> np.ndarray:
        return np.asarray(self.data.info["digests"][i]["length_histogram"], dtype=np.int64)

    def _aa_counts(self, i: int) -> pd.DataFrame:
        aa = self.data.info["digests"][i]["aa_detectability"]
        return pd.DataFrame(aa)


class LiveResults(Results):
    """Results computed in place by digesting ``proteins`` with every recipe (lazily).

    Each recipe is reduced to a compact record (summary, per-entry statistics and changes,
    comparison counts, histogram, amino-acid counts). For large inputs (``lean``; default
    when ``proteins`` has more than :data:`LEAN_RESIDUES` residues) only the baseline digest
    stays in memory and every other digest is released once summarized: a whole-proteome
    digest takes several hundred MB, more than the browser (Pyodide) heap can hold for many
    digests at once.
    """

    mode = "live"

    def __init__(self, proteins: ProteinSet, background: Background,
                 recipes: list[DigestRecipe], lean: bool | None = None):  # fmt: skip
        self.proteins = proteins
        self.background = background
        self.recipes = recipes
        self.lean = bool(proteins.lengths.sum() > LEAN_RESIDUES) if lean is None else lean
        self._singles: dict = {}  # DigestSettings -> Digest (baseline + reusable singles)
        self._an: dict[int, Analysis] = {}
        self._records: dict[int, dict] = {}

    def _single(self, s) -> Digest:
        from .digest import digest

        if s not in self._singles:
            if self.lean:  # keep the baseline only
                base = self.recipes[0].baseline
                self._singles = {k: v for k, v in self._singles.items() if k == base}
            self._singles[s] = digest(self.proteins, s)
        return self._singles[s]

    def analysis(self, i: int) -> Analysis:
        """The analysis of recipe ``i`` (cached unless ``lean`` and ``i`` is not 0)."""
        if i in self._an:
            return self._an[i]
        from .digest import parallel

        ds = [self._single(s) for s in self.recipes[i].settings]
        an = Analysis(ds[0] if len(ds) == 1 else parallel(*ds), self.background)
        if not self.lean or i == 0:
            self._an[i] = an
        return an

    def record(self, i: int) -> dict:
        """Compact results of recipe ``i`` (computed once)."""
        if i in self._records:
            return self._records[i]
        an = self.analysis(i)
        base = self.analysis(0)
        st = an.protein_stats
        d = an.digest
        from .digest import bincount

        lengths = d.length[d.missed == 0]
        aa = an.aa_detectability()
        if i == 0:
            changes = pd.DataFrame({f: np.zeros(len(self.proteins), dtype=np.int64)
                                    for f in CHANGE_FIELDS})  # fmt: skip
            changes.insert(0, "accession", self.proteins.accessions)
            counts = dict.fromkeys(COMPARISON_COUNTS, 0)
        else:
            changes = comparison_per_protein(base, an)
            counts = comparison_counts(base, an)
        rec = {
            "summary": dict(an.summary),
            "stats": _stats_frame(
                self.proteins, st["mature_length"].to_numpy(), st["n_identifiable_peptides"],
                st["n_unique_peptides"], st["n_total_peptides"], st["covered_residues"],
                st["unique_covered_residues"],
            ),
            "changes": changes,
            "counts": counts,
            "hist": bincount(np.minimum(lengths, HIST_MAX + 1), minlength=HIST_MAX + 2)[1:],
            "aa": aa[["residue", "occurrences", "on_identifiable_peptides"]],
        }  # fmt: skip
        self._records[i] = rec
        return rec

    def summary(self, i: int) -> dict[str, float]:
        return dict(self.record(i)["summary"])

    def comparison_counts(self, i: int) -> dict[str, int]:
        return dict(self.record(i)["counts"])

    def entry_stats(self, i: int) -> pd.DataFrame:
        return self.record(i)["stats"]

    def entry_changes(self, i: int) -> pd.DataFrame:
        return self.record(i)["changes"]

    def _hist_counts(self, i: int) -> np.ndarray:
        return self.record(i)["hist"]

    def _aa_counts(self, i: int) -> pd.DataFrame:
        return self.record(i)["aa"]


#: Inputs larger than this (residues) use lean live results by default.
LEAN_RESIDUES = 2_000_000
