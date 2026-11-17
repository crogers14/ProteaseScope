"""Build data/histones/known_marks.tsv from the downloaded UniProt features.

Keeps lysine features (MOD_RES and CROSSLNK) for the entries in the histone FASTA, adds
mature numbering (MetAP rule), and writes the table with a DRAFT status header for review.

Usage:
    uv run python scripts/build_known_marks.py
"""

from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from proteasescope.histones import (  # noqa: E402
    known_marks_from_features,
    load_histone_set,
    parse_uniprot_features,
    transfer_marks,
)

FEATURES = ROOT / "data" / "reference" / "histone_features.tsv"
METADATA = ROOT / "data" / "reference" / "metadata.json"
HISTONES = ROOT / "data" / "histones"
OUT = ROOT / "data" / "histones" / "known_marks.tsv"


def main() -> None:
    rec = json.loads(METADATA.read_text())["files"][FEATURES.name]
    proteins = load_histone_set(HISTONES)  # user FASTA + histone_isoforms.fasta
    marks = known_marks_from_features(parse_uniprot_features(FEATURES.read_text()), proteins)
    marks = transfer_marks(marks, proteins)
    header = [
        "status: DRAFT (generated; awaiting user review)",
        "Known lysine marks of the histone FASTA entries, from UniProt features "
        "ft_mod_res (modified residues) and ft_crosslnk (ubiquitin/SUMO/UFM1 isopeptides).",
        f"UniProt release {rec['release']} ({rec['release_date']}), downloaded "
        f"{rec['download_date']}; source file data/reference/{FEATURES.name}.",
        f"Generated {dt.date.today().isoformat()} by scripts/build_known_marks.py.",
        "position_uniprot: initiator Met = 1. position_mature: UniProt - 1 when the initiator "
        "Met is removed (MetAP rule).",
        "Splice isoforms (histone_isoforms.fasta) inherit the marks of their UniProt entry "
        "where the site lies in identical sequence (evidence_type 'transferred from ...').",
        "One row per (site, modification). evidence = UniProt evidence string; evidence_type = "
        "strongest class (experimental > large-scale proteomics > by similarity > curator "
        "inference).",
    ]
    with OUT.open("w", newline="") as fh:
        for line in header:
            fh.write(f"# {line}\n")
        marks.to_csv(fh, sep="\t", index=False)
    n_sites = marks[["accession", "position_uniprot"]].drop_duplicates().shape[0]
    print(
        f"{OUT.relative_to(ROOT)}: {len(marks)} rows, {n_sites} lysine sites, "
        f"{marks.accession.nunique()} entries"
    )


if __name__ == "__main__":
    main()
