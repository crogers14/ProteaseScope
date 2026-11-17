"""FASTA reading and UniProt header parsing (including isoforms).

Sequences are stored as uppercase one-letter strings indexed 0-based in memory;
all user-facing positions use UniProt numbering (residue 1 = initiator Met).
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

_UNIPROT_HEADER = re.compile(r"^(?P<db>sp|tr)\|(?P<acc>[^|\s]+)\|(?P<entry>\S+)\s*(?P<desc>.*)$")
_KEYS = ("OS", "OX", "GN", "PE", "SV")
_KEY_SPLIT = re.compile(r"\s(OS|OX|GN|PE|SV)=")

META_COLUMNS = [
    "accession",
    "entry_name",
    "protein_name",
    "gene",
    "organism",
    "organism_id",
    "pe",
    "sv",
    "db",
    "canonical_accession",
    "isoform",
    "is_isoform",
    "header",
]


def parse_header(header: str) -> dict:
    """Parse one FASTA header (without the leading ``>``).

    UniProt headers (``sp|ACC|ENTRY Name OS=.. OX=.. GN=.. PE=.. SV=..``, isoforms as
    ``ACC-2``) are split into fields; other headers use the first word as accession and
    the rest as protein name.

    Returns:
        Dict with the keys in :data:`META_COLUMNS`. Missing fields are ``None``.
    """
    header = header.strip()
    out = dict.fromkeys(META_COLUMNS)
    out["header"] = header
    m = _UNIPROT_HEADER.match(header)
    if m:
        out["db"] = m["db"]
        out["accession"] = m["acc"]
        out["entry_name"] = m["entry"]
        parts = _KEY_SPLIT.split(" " + m["desc"])
        out["protein_name"] = parts[0].strip() or None
        fields = dict(zip(parts[1::2], (p.strip() for p in parts[2::2]), strict=True))
        out["organism"] = fields.get("OS")
        out["organism_id"] = fields.get("OX")
        out["gene"] = fields.get("GN")
        out["pe"] = int(fields["PE"]) if fields.get("PE", "").isdigit() else None
        out["sv"] = int(fields["SV"]) if fields.get("SV", "").isdigit() else None
    else:
        first, _, rest = header.partition(" ")
        out["accession"] = first
        out["protein_name"] = rest.strip() or None
    acc = out["accession"]
    base, dash, iso = acc.rpartition("-")
    if dash and base and iso.isdigit():
        out["canonical_accession"], out["isoform"], out["is_isoform"] = base, iso, True
    else:
        out["canonical_accession"], out["isoform"], out["is_isoform"] = acc, None, False
    return out


@dataclass
class ProteinSet:
    """A set of protein entries.

    Attributes:
        meta: One row per entry, columns :data:`META_COLUMNS` (row order = entry index).
        sequences: Uppercase sequences, same order as ``meta``.
        name: Description of the source (file name, database release), printed in exports.
    """

    meta: pd.DataFrame
    sequences: list[str]
    name: str = ""

    def __post_init__(self) -> None:
        if len(self.meta) != len(self.sequences):
            raise ValueError("meta and sequences differ in length")
        self.meta = self.meta.reset_index(drop=True)

    def __len__(self) -> int:
        return len(self.sequences)

    @property
    def accessions(self) -> list[str]:
        """Entry accessions in order."""
        return self.meta["accession"].tolist()

    @property
    def lengths(self) -> np.ndarray:
        """Sequence lengths (residues), int64."""
        return np.fromiter((len(s) for s in self.sequences), dtype=np.int64, count=len(self))

    def index_of(self, accession: str) -> int:
        """Entry index of ``accession`` (first match). Raises ``KeyError`` if absent."""
        hits = np.flatnonzero(self.meta["accession"].to_numpy() == accession)
        if not len(hits):
            raise KeyError(accession)
        return int(hits[0])

    def sequence(self, accession: str) -> str:
        """Sequence of ``accession``."""
        return self.sequences[self.index_of(accession)]

    def subset(self, indices: Iterable[int]) -> ProteinSet:
        """New set with the given entry indices, in that order."""
        idx = list(indices)
        return ProteinSet(self.meta.iloc[idx].copy(), [self.sequences[i] for i in idx], self.name)

    @staticmethod
    def concat(sets: Iterable[ProteinSet], drop_exact_duplicates: bool = True) -> ProteinSet:
        """Concatenate sets.

        By default an entry is dropped when an earlier entry has the same accession *and* the
        same sequence. Entries sharing an accession but differing in sequence are kept.
        """
        sets = list(sets)
        meta = pd.concat([s.meta for s in sets], ignore_index=True)
        seqs = [q for s in sets for q in s.sequences]
        keep = np.ones(len(seqs), dtype=bool)
        if drop_exact_duplicates:
            key = pd.Series(list(zip(meta["accession"], seqs, strict=True)))
            keep = ~key.duplicated().to_numpy()
            meta = meta[keep]
            seqs = [q for q, k in zip(seqs, keep, strict=True) if k]
        names, start = [], 0
        for i, s in enumerate(sets):
            label = s.name or f"{len(s)} entries"
            added = int(keep[start : start + len(s)].sum())
            if i and added < len(s):
                label += f" ({added} not already present)"
            names.append(label)
            start += len(s)
        return ProteinSet(meta, seqs, " + ".join(names))


def parse_fasta(text: str) -> ProteinSet:
    """Parse FASTA text into a :class:`ProteinSet`.

    Whitespace and a trailing ``*`` are removed and residues are uppercased. Raises
    ``ValueError`` on duplicate accessions, non-letter residues or entries without
    sequence lines.
    """
    headers: list[str] = []
    chunks: list[list[str]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            headers.append(line[1:])
            chunks.append([])
        elif not headers:
            raise ValueError("FASTA text does not start with a '>' header line")
        else:
            chunks[-1].append(line)
    seqs = []
    for h, c in zip(headers, chunks, strict=True):
        s = "".join(c).replace(" ", "").upper().rstrip("*")
        if not s:
            raise ValueError(f"Entry has no sequence: >{h}")
        if not s.isalpha() or not s.isascii():
            bad = sorted({ch for ch in s if not (ch.isascii() and ch.isalpha())})
            raise ValueError(f"Invalid residue characters {bad} in >{h}")
        seqs.append(s)
    meta = pd.DataFrame([parse_header(h) for h in headers], columns=META_COLUMNS)
    dup = meta["accession"][meta["accession"].duplicated()]
    if len(dup):
        raise ValueError(f"Duplicate accessions in FASTA: {sorted(set(dup))[:10]}")
    return ProteinSet(meta, seqs)


def read_fasta(path: str | Path, name: str | None = None) -> ProteinSet:
    """Read a FASTA file (plain text) into a :class:`ProteinSet` named after the file."""
    path = Path(path)
    ps = parse_fasta(path.read_text())
    ps.name = name or f"{path.name} ({len(ps)} entries)"
    return ps


def read_reference_fasta(path: str | Path, metadata: str | Path | None = None) -> ProteinSet:
    """Read a downloaded reference FASTA and name it with its provenance.

    ``metadata`` defaults to ``metadata.json`` next to the file (written by
    ``scripts/fetch_data.py``); its record supplies the description and UniProt release.
    """
    import json

    path = Path(path)
    meta_path = Path(metadata) if metadata else path.with_name("metadata.json")
    ps = read_fasta(path)
    if meta_path.exists():
        rec = json.loads(meta_path.read_text())["files"].get(path.name)
        if rec:
            ps.name = (
                f"{rec.get('description', path.name)}, UniProt release {rec.get('release')} "
                f"({rec.get('release_date')}; downloaded {rec.get('download_date')}), "
                f"{len(ps):,} entries"
            )
    return ps
