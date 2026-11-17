"""Download UniProt reference data and record provenance in data/reference/metadata.json.

Downloads:
  1. Reviewed human proteome (UP000005640) with isoforms, FASTA.
  2. UniProt features (modified residues, cross-links) for the accessions in the histone FASTA.

Usage:
    uv run python scripts/fetch_data.py [--force]
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from proteasescope.fasta import parse_fasta, read_fasta  # noqa: E402

REF = ROOT / "data" / "reference"
METADATA = REF / "metadata.json"
HISTONE_FASTA = ROOT / "data" / "histones" / "human_histones.fasta"

PROTEOME_URL = (
    "https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true"
    "&query=(proteome:UP000005640)+AND+(reviewed:true)"
)
PROTEOME_FILE = REF / "uniprot_human_reviewed_isoforms.fasta"

FEATURE_FIELDS = "accession,gene_names,protein_name,length,ft_mod_res,ft_crosslnk"
FEATURES_FILE = REF / "histone_features.tsv"

USER_AGENT = "ProteaseScope-fetch/0.1 (python-urllib)"


def _download(url: str) -> tuple[bytes, dict[str, str]]:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=600) as resp:
        return resp.read(), {k.lower(): v for k, v in resp.headers.items()}


def _record(path: Path, url: str, headers: dict[str, str], data: bytes, **extra) -> dict:
    return {
        "source_url": url,
        "release": headers.get("x-uniprot-release"),
        "release_date": headers.get("x-uniprot-release-date"),
        "download_date": dt.date.today().isoformat(),
        "sha256": hashlib.sha256(data).hexdigest(),
        "bytes": len(data),
        **extra,
    }


def fetch_proteome() -> dict:
    data, headers = _download(PROTEOME_URL)
    proteins = parse_fasta(data.decode("ascii"))  # validates before writing
    PROTEOME_FILE.write_bytes(data)
    n_iso = int(proteins.meta["is_isoform"].sum())
    return _record(
        PROTEOME_FILE,
        PROTEOME_URL,
        headers,
        data,
        description="UniProtKB/Swiss-Prot reviewed human proteome UP000005640 incl. isoforms",
        entries=len(proteins),
        canonical_entries=len(proteins) - n_iso,
        isoform_entries=n_iso,
        residues=int(proteins.lengths.sum()),
    )


def fetch_histone_features() -> dict:
    accessions = read_fasta(HISTONE_FASTA).accessions
    query = " OR ".join(f"(accession:{a})" for a in accessions)
    url = (
        "https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields="
        + FEATURE_FIELDS
        + "&query="
        + urllib.parse.quote(query)
    )
    data, headers = _download(url)
    lines = data.decode("utf-8").splitlines()
    got = {line.split("\t", 1)[0] for line in lines[1:] if line}
    missing = sorted(set(accessions) - got)
    if missing:
        raise RuntimeError(f"UniProt returned no features for {missing}")
    FEATURES_FILE.write_bytes(data)
    return _record(
        FEATURES_FILE,
        url,
        headers,
        data,
        description="UniProt features (ft_mod_res, ft_crosslnk) for the histone FASTA accessions",
        fields=FEATURE_FIELDS,
        entries=len(lines) - 1,
    )


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--force", action="store_true", help="re-download existing files")
    args = ap.parse_args(argv)
    REF.mkdir(parents=True, exist_ok=True)
    meta = json.loads(METADATA.read_text()) if METADATA.exists() else {"files": {}}
    for path, fetch in ((PROTEOME_FILE, fetch_proteome), (FEATURES_FILE, fetch_histone_features)):
        if path.exists() and path.name in meta["files"] and not args.force:
            print(f"{path.name}: present, skipping (use --force to re-download)")
            continue
        print(f"{path.name}: downloading ...", flush=True)
        meta["files"][path.name] = fetch()
        METADATA.write_text(json.dumps(meta, indent=2) + "\n")
        rec = meta["files"][path.name]
        print(f"  release {rec['release']} ({rec['release_date']}), {rec['bytes']:,} bytes")


if __name__ == "__main__":
    main()
