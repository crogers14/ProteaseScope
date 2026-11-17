"""Export Phase 2 outputs (1-5) as TSV files with settings headers.

Writes:
  results/phase2/histones/  Histone PTM preset on the histone FASTA, all outputs, plus
                            comparisons against every alternative digest (alone, parallel,
                            sequential) and the all-pairs ranking table.
  results/phase2/proteome/  Proteome preset on the reviewed human proteome: summary, Output 3,
                            Output 4 (Output 1 has ~2.8 M rows; use --proteome-peptides).

Usage:
    uv run python scripts/export_results.py [--proteome-peptides] [--proteome-ranking]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from proteasescope import metrics, report  # noqa: E402
from proteasescope.chemistry import Chemistry  # noqa: E402
from proteasescope.digest import DigestSettings, digest, parallel  # noqa: E402
from proteasescope.fasta import read_fasta, read_reference_fasta  # noqa: E402
from proteasescope.rules import available_enzymes, get_rule  # noqa: E402

OUT = ROOT / "results" / "phase2"
REF = ROOT / "data" / "reference" / "uniprot_human_reviewed_isoforms.fasta"
HISTONES = ROOT / "data" / "histones" / "human_histones.fasta"

HISTONE_PTM = DigestSettings(
    (get_rule("Trypsin", "psi-ms"),),
    missed_cleavages=1,
    min_length=4,
    max_length=40,
    chemistry=Chemistry(propionylation="before_and_after"),
    name="Histone PTM baseline: Trypsin (PSI-MS), propionylated K",
)
PROTEOME = DigestSettings(
    (get_rule("Trypsin", "psi-ms"),),
    missed_cleavages=1,
    min_length=7,
    max_length=30,
    chemistry=Chemistry(),
    name="Proteome baseline: Trypsin (PSI-MS)",
)


def _write(df: pd.DataFrame, path: Path, analyses, title: str) -> None:
    report.write_tsv(df, path, report.settings_header(analyses, title))
    print(f"  {path.relative_to(ROOT)}  ({len(df):,} rows)")


def _per_digest(an: metrics.Analysis, folder: Path, prefix: str, peptides: bool = True) -> None:
    if peptides:
        _write(
            an.peptide_table(),
            folder / f"{prefix}output1_identifiable_peptides.tsv",
            an,
            "Output 1: Identifiable peptides",
        )
        _write(
            an.unique_peptide_table(),
            folder / f"{prefix}output2_unique_peptides.tsv",
            an,
            "Output 2: Unique peptides",
        )
    _write(
        an.protein_table(),
        folder / f"{prefix}output3_identifiable_proteins.tsv",
        an,
        "Output 3: Identifiable proteins",
    )
    st = an.digest_statistics()
    _write(
        st["length_histogram"],
        folder / f"{prefix}output4_length_histogram.tsv",
        an,
        "Output 4: Cleaved peptide frequency (all products, 0 missed cleavages)",
    )
    _write(
        st["coverage_per_protein"],
        folder / f"{prefix}output4_coverage_per_protein.tsv",
        an,
        "Output 4: Cleaved peptide coverage per protein",
    )
    _write(
        pd.DataFrame([st["coverage_summary"]]),
        folder / f"{prefix}output4_coverage_summary.tsv",
        an,
        "Output 4: Cleaved peptide coverage, proteome-wide",
    )
    _write(
        st["aa_detectability"],
        folder / f"{prefix}output4_aa_detectability.tsv",
        an,
        "Output 4: Amino acid detectability",
    )


def export_histones(ref, background) -> None:
    folder = OUT / "histones"
    folder.mkdir(parents=True, exist_ok=True)
    h = read_fasta(HISTONES)
    print("Histones (Histone PTM preset):")
    base = metrics.Analysis(digest(h, HISTONE_PTM), background)
    _per_digest(base, folder, "baseline_")

    alts = []
    for name in available_enzymes("psi-ms"):
        rule = get_rule(name, "psi-ms")
        if rule in HISTONE_PTM.rules:
            continue
        alone = digest(h, HISTONE_PTM.with_rules(rule))
        alts += [
            ("alone", rule, metrics.Analysis(alone, background)),
            ("parallel", rule, metrics.Analysis(parallel(base.digest, alone), background)),
            (
                "sequential",
                rule,
                metrics.Analysis(
                    digest(h, HISTONE_PTM.with_rules(*HISTONE_PTM.rules, rule)), background
                ),
            ),
        ]
    analyses = [base] + [a for _, _, a in alts]

    lists: dict[str, list[pd.DataFrame]] = {k: [] for k in (
        "rescued_proteins", "lost_proteins", "proteins_gaining_unique", "residues", "sites",
        "per_protein")}  # fmt: skip
    for kind, _rule, an in alts:
        c = metrics.compare(base, an)
        for key in lists:
            df = getattr(c, key).copy()
            df.insert(0, "kind", kind)
            df.insert(0, "alternative", an.label)
            lists[key].append(df)
    summary = metrics.summary_table(base, [a for _, _, a in alts])
    summary.insert(1, "kind", ["baseline"] + [k for k, _, _ in alts])
    _write(summary, folder / "output5_summary.tsv", analyses, "Output 5: Summary table")
    for key, frames in lists.items():
        _write(
            pd.concat(frames, ignore_index=True),
            folder / f"output5_{key}.tsv",
            analyses,
            f"Output 5: {key.replace('_', ' ')}",
        )
    rank = metrics.ranking_table(
        h, HISTONE_PTM, [get_rule(n, "psi-ms") for n in available_enzymes("psi-ms")], background
    )
    _write(rank, folder / "output5_ranking_table.tsv", analyses, "Output 5: All-pairs ranking")
    (folder / "methods.txt").write_text(report.methods_text(base) + "\n")


def export_proteome(ref, background, peptides: bool, ranking: bool) -> None:
    folder = OUT / "proteome"
    folder.mkdir(parents=True, exist_ok=True)
    print("Proteome (Proteome preset):")
    base = metrics.Analysis(digest(ref, PROTEOME), background)
    _write(pd.DataFrame([base.summary]), folder / "baseline_summary.tsv", base, "Summary metrics")
    _per_digest(base, folder, "baseline_", peptides=peptides)
    if ranking:
        rank = metrics.ranking_table(
            ref,
            PROTEOME,
            [get_rule(n, "psi-ms") for n in available_enzymes("psi-ms")],
            background,
        )
        _write(rank, folder / "output5_ranking_table.tsv", base, "Output 5: All-pairs ranking")
    (folder / "methods.txt").write_text(report.methods_text(base) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--proteome-peptides",
        action="store_true",
        help="also write proteome Outputs 1-2 (~2.8 M rows, several hundred MB)",
    )
    ap.add_argument(
        "--proteome-ranking",
        action="store_true",
        help="also compute the proteome all-pairs ranking table (several minutes)",
    )
    args = ap.parse_args()
    ref = read_reference_fasta(REF)
    h = read_fasta(HISTONES)
    export_histones(ref, metrics.uniqueness_background(ref, h))
    export_proteome(ref, metrics.Background(ref), args.proteome_peptides, args.proteome_ranking)


if __name__ == "__main__":
    main()
