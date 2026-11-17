"""Export Phase 3 histone outputs as TSV files with settings headers.

Baseline: Histone PTM preset (Trypsin PSI-MS, propionylation before and after digestion,
4-40 aa, <=1 missed cleavage). Alternatives: every other PSI-MS-mode enzyme alone, in
parallel with the baseline, and sequential after it (36 digests, letters A..AJ).

Writes results/phase3/histones/. Usage:
    uv run python scripts/export_histone_results.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from proteasescope import histones as H  # noqa: E402
from proteasescope import report  # noqa: E402
from proteasescope.chemistry import Chemistry  # noqa: E402
from proteasescope.digest import DigestSettings, digest, parallel  # noqa: E402
from proteasescope.metrics import Analysis  # noqa: E402
from proteasescope.rules import available_enzymes, get_rule  # noqa: E402

OUT = ROOT / "results" / "phase3" / "histones"
PTM = DigestSettings(
    (get_rule("Trypsin", "psi-ms"),),
    missed_cleavages=1,
    min_length=4,
    max_length=40,
    chemistry=Chemistry(propionylation="before_and_after"),
    name="Histone PTM baseline: Trypsin (PSI-MS), propionylated K",
)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ps, ann, marks = H.load_histone_inputs(ROOT / "data" / "histones")
    base = Analysis(digest(ps, PTM))
    analyses, kinds = [base], [("baseline", "")]
    for name in available_enzymes("psi-ms"):
        rule = get_rule(name, "psi-ms")
        if rule in PTM.rules:
            continue
        alone = digest(ps, PTM.with_rules(rule))
        for kind, d in (
            ("alone", alone),
            ("parallel", parallel(base.digest, alone)),
            ("sequential", digest(ps, PTM.with_rules(*PTM.rules, rule))),
        ):
            analyses.append(Analysis(d))
            kinds.append((kind, rule.label))
    header = report.settings_header(analyses, "Histone outputs (Phase 3)")

    def write(df: pd.DataFrame, name: str, title: str) -> None:
        path = OUT / name
        report.write_tsv(df, path, [header[0] + f" - {title}", *header[1:]])
        print(f"  {path.relative_to(ROOT)}  ({len(df):,} rows)")

    co = H.lysine_cooccurrence(analyses, ann)
    legend = pd.DataFrame(
        {
            "letter": [""] + co["letters"]["letter"].tolist(),
            "digest": [a.label for a in analyses],
            "kind": [k for k, _ in kinds],
            "enzyme": [e for _, e in kinds],
        }
    )
    write(legend, "digests.tsv", "Digest legend (letters used in co-occurrence matrix)")

    vis = H.mark_visibility(analyses, marks, ann)
    write(vis, "mark_visibility.tsv", "Known-mark visibility (canonical entries)")
    write(H.mark_grid(vis).reset_index(), "mark_grid.tsv", "Known marks x digests grid")

    iso = H.isoform_discrimination(analyses, ann)
    write(iso, "isoform_discrimination.tsv", "Isoform discrimination (canonical entries)")
    for fam in iso[iso.digest == base.label]["family"].unique():
        mat = H.discrimination_matrix(iso, base.label, fam)
        write(mat.reset_index(names="entry"), f"isoform_matrix_baseline_{fam}.tsv",
              f"Isoform discrimination matrix, baseline, {fam}")  # fmt: skip

    write(co["sites"], "cooccurrence_sites.tsv", "Lysine sites covered")
    write(co["pairs"], "cooccurrence_pairs.tsv", "Lysine pairs vs baseline")
    write(co["matrix"], "cooccurrence_matrix.tsv", "Lysine co-occurrence matrix data")

    maps = {"peptides": [], "strips": [], "coverage": []}
    for acc in H.cooccurrence_entries(ann):
        pm = H.tail_peptide_map(analyses, ann, acc)
        for key in maps:
            df = pm[key].copy()
            df.insert(0, "accession", acc)
            df.insert(1, "display_name", ann.set_index("accession").loc[acc, "display_name"])
            maps[key].append(df)
    for key, frames in maps.items():
        write(pd.concat(frames, ignore_index=True), f"tail_map_{key}.tsv",
              f"Histone tail peptide maps: {key}")  # fmt: skip
    (OUT / "methods.txt").write_text(report.methods_text(analyses[:1]) + "\n")


if __name__ == "__main__":
    main()
