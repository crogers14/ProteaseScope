"""Precomputed app data: build, write and load.

Layout of ``data/precomputed/``:

* ``manifest.json``: software and data versions, presets, digests, file sizes, sha256, timings.
* ``reference.fasta.gz``: the reference proteome (gzip FASTA).
* ``uniqueness_index.npz``: :class:`~proteasescope.metrics.SharedPrefixIndex` arrays.
* ``preset_<key>.npz``: per-digest x per-entry arrays (uint16), see :data:`ENTRY_FIELDS`.
* ``preset_<key>.json``: per-digest tables: ranking/summary rows, length histograms,
  amino-acid detectability.

Per-entry arrays are protein-level summaries; peptide-level detail is computed live for the
protein being viewed (CLAUDE.md section 2).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .digest import bincount
from .fasta import ProteinSet, parse_fasta
from .metrics import (
    COMPARISON_COUNTS,
    SUMMARY_METRICS,
    Analysis,
    Background,
    SharedPrefixIndex,
    comparison_counts,
    comparison_per_protein,
)
from .presets import Preset, preset_digests
from .report import digest_lines
from .rules import PSI_MS_VERSION

#: Per-entry statistics of each digest.
STAT_FIELDS = [
    "n_identifiable_peptides", "n_unique_peptides", "n_total_peptides", "covered_residues",
    "unique_covered_residues",
]  # fmt: skip
#: Per-entry changes versus the preset baseline.
CHANGE_FIELDS = [
    "residues_gained", "residues_lost", "K_sites_gained", "K_sites_lost", "S_sites_gained",
    "S_sites_lost", "T_sites_gained", "T_sites_lost", "Y_sites_gained", "Y_sites_lost",
]  # fmt: skip
ENTRY_FIELDS = STAT_FIELDS + CHANGE_FIELDS
#: Length histogram bins: 1..HIST_MAX aa, plus one overflow bin (> HIST_MAX).
HIST_MAX = 100


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_preset(preset: Preset, reference: ProteinSet, background: Background) -> dict:
    """Compute every digest of ``preset`` on ``reference``.

    Returns a dict with ``arrays`` (field -> uint16/uint32 array, shape digests x entries),
    ``digests`` (one record per digest: label, kind, enzyme, settings lines, summary metrics,
    comparison counts, length histogram, amino-acid detectability) and ``seconds``.
    """
    t0 = time.perf_counter()
    digests = preset_digests(preset, reference)
    base = Analysis(digests[0][2], background)
    n_d, n_e = len(digests), len(reference)
    arrays = {f: np.zeros((n_d, n_e), dtype=np.int64) for f in ENTRY_FIELDS}
    records = []
    for i, (kind, enzyme, d) in enumerate(digests):
        an = base if i == 0 else Analysis(d, background)
        st = an.protein_stats
        for f in STAT_FIELDS:
            arrays[f][i] = st[f].to_numpy()
        if i:
            per = comparison_per_protein(base, an)
            for f in CHANGE_FIELDS:
                arrays[f][i] = per[f].to_numpy()
            counts = comparison_counts(base, an)
        else:
            counts = dict.fromkeys(COMPARISON_COUNTS, 0)
        lengths = d.length[d.missed == 0]
        hist = bincount(np.minimum(lengths, HIST_MAX + 1), minlength=HIST_MAX + 2)[1:]
        aa = an.aa_detectability()
        records.append(
            {
                "index": i,
                "label": an.label,
                "kind": kind,
                "enzyme": enzyme,
                "flags": list(
                    dict.fromkeys(f for s in d.settings for r in s.rules for f in r.flags)
                ),
                "settings": [line for s in d.settings for line in digest_lines(s)],
                "summary": an.summary,
                "comparison": counts,
                "length_histogram": hist.tolist(),
                "aa_detectability": {
                    "residue": aa["residue"].tolist(),
                    "occurrences": aa["occurrences"].tolist(),
                    "on_identifiable_peptides": aa["on_identifiable_peptides"].tolist(),
                },
            }
        )
    out_arrays = {}
    for f, a in arrays.items():
        dtype = np.uint16 if a.max(initial=0) < 2**16 else np.uint32
        out_arrays[f] = a.astype(dtype)
    return {"arrays": out_arrays, "digests": records, "seconds": time.perf_counter() - t0}


def write_preset(result: dict, preset: Preset, outdir: Path) -> dict:
    """Write ``preset_<key>.npz`` and ``preset_<key>.json``; return their manifest records."""
    npz = outdir / f"preset_{preset.key}.npz"
    js = outdir / f"preset_{preset.key}.json"
    np.savez_compressed(npz, **result["arrays"])
    js.write_text(
        json.dumps(
            {
                "preset": preset.key,
                "name": preset.name,
                "mode": preset.mode,
                "notes": preset.notes,
                "hist_max": HIST_MAX,
                "digests": result["digests"],
            }
        )
    )
    return {p.name: {"bytes": p.stat().st_size, "sha256": _sha256(p)} for p in (npz, js)}


def write_reference(reference: ProteinSet, fasta_path: Path, index: SharedPrefixIndex,
                    outdir: Path) -> dict:  # fmt: skip
    """Write ``reference.fasta.gz`` (from the original FASTA bytes) and the index."""
    ref = outdir / "reference.fasta.gz"
    ref.write_bytes(gzip.compress(Path(fasta_path).read_bytes(), compresslevel=9, mtime=0))
    idx = outdir / "uniqueness_index.npz"
    index.save(idx)
    return {p.name: {"bytes": p.stat().st_size, "sha256": _sha256(p)} for p in (ref, idx)}


def write_manifest(outdir: Path, files: dict, reference_meta: dict, presets: dict,
                   timings: dict) -> Path:  # fmt: skip
    """Write ``manifest.json``."""
    manifest = {
        "software": f"ProteaseScope {__version__}",
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "notice": "PREDICTED (in silico) results; not observed data.",
        "psi_ms": PSI_MS_VERSION,
        "reference": reference_meta,
        "uniqueness": "substring mapping (I = L) via a precomputed shared-prefix index",
        "presets": presets,
        "files": files,
        "timings_seconds": timings,
        "entry_fields": ENTRY_FIELDS,
    }
    path = outdir / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


# ----- loading --------------------------------------------------------------------------
@dataclass
class PrecomputedPreset:
    """A loaded preset: per-digest tables and per-entry arrays."""

    key: str
    info: dict
    arrays: dict[str, np.ndarray]
    accessions: list[str]

    @property
    def digests(self) -> pd.DataFrame:
        """One row per digest: index, label, kind, enzyme, flags."""
        return pd.DataFrame(
            [{k: r[k] for k in ("index", "label", "kind", "enzyme")} | {
                "flags": ", ".join(r["flags"])} for r in self.info["digests"]]
        )  # fmt: skip

    def ranking_table(self) -> pd.DataFrame:
        """All-pairs ranking table (same columns as :func:`metrics.ranking_table`)."""
        rows = []
        base = self.info["digests"][0]["summary"]
        for r in self.info["digests"]:
            row = {"digest": r["label"], "kind": r["kind"], "enzyme": r["enzyme"],
                   "flags": ", ".join(r["flags"])}  # fmt: skip
            for m in SUMMARY_METRICS:
                row[m] = r["summary"][m]
                row[f"{m}_diff"] = r["summary"][m] - base[m]
            row.update(r["comparison"])
            rows.append(row)
        return pd.DataFrame(rows)

    def entry_table(self, digest_index: int) -> pd.DataFrame:
        """Per-entry arrays of one digest as a DataFrame (accession + ENTRY_FIELDS)."""
        df = pd.DataFrame({f: self.arrays[f][digest_index] for f in self.arrays})
        df.insert(0, "accession", self.accessions)
        return df


def load_reference_bytes(manifest: bytes, fasta_gz: bytes, index_npz: bytes
                         ) -> tuple[ProteinSet, SharedPrefixIndex, dict]:  # fmt: skip
    """Load the shipped reference from raw file contents (e.g. fetched over HTTP)."""
    import io

    man = json.loads(manifest)
    ref = parse_fasta(gzip.decompress(fasta_gz).decode())
    r = man["reference"]
    ref.name = (
        f"{r.get('description')}, UniProt release {r.get('release')} ({r.get('release_date')};"
        f" downloaded {r.get('download_date')}), {len(ref):,} entries"
    )
    return ref, SharedPrefixIndex.load(io.BytesIO(index_npz), ref), man


def load_reference(outdir: Path) -> tuple[ProteinSet, SharedPrefixIndex, dict]:
    """Load the shipped reference, its uniqueness index and the manifest."""
    outdir = Path(outdir)
    return load_reference_bytes(
        (outdir / "manifest.json").read_bytes(),
        (outdir / "reference.fasta.gz").read_bytes(),
        (outdir / "uniqueness_index.npz").read_bytes(),
    )


def load_preset_bytes(key: str, json_bytes: bytes, npz_bytes: bytes, accessions: list[str]
                      ) -> PrecomputedPreset:  # fmt: skip
    """Load a preset from raw file contents."""
    import io

    info = json.loads(json_bytes)
    z = np.load(io.BytesIO(npz_bytes))
    arrays = {f: z[f] for f in z.files}
    if any(a.shape[1] != len(accessions) for a in arrays.values()):
        raise ValueError("preset arrays do not match the reference entries")
    return PrecomputedPreset(key, info, arrays, accessions)


def load_preset(outdir: Path, key: str, accessions: list[str]) -> PrecomputedPreset:
    """Load ``preset_<key>`` files; ``accessions`` = reference entry order."""
    outdir = Path(outdir)
    return load_preset_bytes(
        key,
        (outdir / f"preset_{key}.json").read_bytes(),
        (outdir / f"preset_{key}.npz").read_bytes(),
        accessions,
    )
