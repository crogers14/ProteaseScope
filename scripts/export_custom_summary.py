"""Export "All goals" summaries of a custom strategy versus the standard digest.

Default: RIPUP (Turner et al. 2026) vs the standard histone digest for canonical H2A type 1,
H2B type 1-C, H3.1 and H4, at the N-terminal tail and the whole protein, under the
CUSTOM_STRATEGY.md section 9 regression settings (PSI-MS rules, 7-30 aa, 0 missed cleavages,
Met removed by the MetAP rule, mature numbering). ``--preset histone_ptm`` uses the Histone
PTM preset's shared settings (4-40 aa, <=1 missed cleavage) instead.

    uv run python scripts/export_custom_summary.py [--preset histone_ptm] [--out results/phase7]
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import pandas as pd

from proteasescope import histones as H
from proteasescope import strategy as S
from proteasescope import workflow as W
from proteasescope.chemistry import Chemistry
from proteasescope.digest import DigestSettings
from proteasescope.presets import PRESETS
from proteasescope.rules import get_rule

ROOT = Path(__file__).resolve().parents[1]
ENTRIES = [("P0C0S8", "H2A type 1"), ("P62807", "H2B type 1-C"), ("P68431", "H3.1"),
           ("P62805", "H4")]  # fmt: skip
REGIONS = [S.N_TAIL, S.WHOLE]
REGRESSION = DigestSettings(
    (get_rule("Trypsin", "psi-ms"),), missed_cleavages=0, min_length=7, max_length=30,
    chemistry=Chemistry(propionylation="before_and_after"),
)  # fmt: skip


def summaries(shared: DigestSettings, custom: W.Workflow, against: W.Workflow):
    hs, ann, marks = H.load_histone_inputs(ROOT / "data" / "histones")
    out = []
    for acc, _ in ENTRIES:
        for region in REGIONS:
            r = S.compare_workflows(hs, ann, marks, acc, region, shared, custom, against=against)
            s = r.all_goals_summary(custom.name)
            s.insert(0, "region", f"{region} ({r.scope.region[0]}–{r.scope.region[1]})")
            s.insert(0, "histone", f"{r.scope.display_name} ({acc})")
            out.append((r, s))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", default=None, choices=[None, "histone_ptm"])
    ap.add_argument("--out", default=str(ROOT / "results" / "phase7"))
    args = ap.parse_args()
    shared = REGRESSION if args.preset is None else replace(
        PRESETS[args.preset].baseline, chemistry=REGRESSION.chemistry)  # fmt: skip
    custom = W.builder_preset("RIPUP (Turner et al. 2026)")
    std = W.builder_preset("Standard histone digest")
    res = summaries(shared, custom, std)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = args.preset or "regression"
    header = res[0][0].header_lines(title=f"{custom.name} vs {std.name}: all goals")
    scopes = "; ".join(dict.fromkeys(r.scope.description for r, _ in res))
    header = [f"Histone scopes: {scopes}" if line.startswith("Proteins digested") else line
              for line in header if not line.startswith("Uniqueness background")]  # fmt: skip
    table = pd.concat([s for _, s in res], ignore_index=True)
    tsv = out / f"ripup_vs_standard_{tag}.tsv"
    with tsv.open("w") as fh:
        fh.writelines(f"# {line}\n" for line in header)
        table.to_csv(fh, sep="\t", index=False)
    md = [f"# {custom.name} vs {std.name}: all goals (PREDICTED, in silico)", ""]
    md += [f"- {line}" for line in header[1:]] + [""]
    for _, s in res:
        md += [f"## {s['histone'].iloc[0]}, {s['region'].iloc[0]}", ""]
        md += ["| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |",
               "| --- | --- | --- | --- | --- | --- | --- |"]  # fmt: skip
        for _, row in s.iterrows():
            md.append(f"| {row['goal']} | {row['Compared-against → Custom']} | {row['Gained']} | "
                      f"{row['Lost']} | {row['Net']:+d} | {row['gained_items'] or '—'} | "
                      f"{row['lost_items'] or '—'} |")  # fmt: skip
        md.append("")
    (out / f"ripup_vs_standard_{tag}.md").write_text("\n".join(md))
    print(f"wrote {tsv} and {tsv.with_suffix('.md')}")


if __name__ == "__main__":
    main()
