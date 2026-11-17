"""Strategy picker examples (STRATEGY_PICKER_2.md step 1).

For H4 1-25, H3.1 1-40 and H2B type 1-C/E/F/G/I 1-35 under the Histone PTM preset, writes
the top-3 answer cards (or the no-improvement message) for every goal as plain text, plus
the advanced table and the "compare across all histones" grid as TSV with the settings
header, into results/phase6/.

Usage:
    uv run python scripts/export_strategy_example.py
"""

from __future__ import annotations

import time
from pathlib import Path

from proteasescope import histones as H
from proteasescope import strategy as S
from proteasescope.presets import PRESETS
from proteasescope.report import write_tsv
from proteasescope.rules import available_enzymes, get_rule

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase6"
EXAMPLES = [("P62805", "H4"), ("P68431", "H3.1"), ("P62807", "H2B")]


def main() -> None:
    hs, ann, marks = H.load_histone_inputs(ROOT / "data" / "histones")
    rules = S.default_candidates([get_rule(n, "psi-ms") for n in available_enzymes("psi-ms")])
    base = PRESETS["histone_ptm"].baseline
    OUT.mkdir(parents=True, exist_ok=True)
    text = [S.standard_line(base), ""]
    for acc, short in EXAMPLES:
        t = time.perf_counter()
        res = S.pick_strategies(hs, ann, marks, acc, S.N_TAIL, base, rules)
        blocks = []
        for key, goal in S.GOALS.items():
            cards = res.cards(key)
            body = [c.text() for c in cards] or [res.message(key)]
            blocks.append(f"### {goal.title}\n" + "\n".join(body))
        elapsed = time.perf_counter() - t
        text += [f"## {res.scope.where}  ({elapsed:.2f} s for all goals)", *blocks, ""]
        write_tsv(res.advanced_table(), OUT / f"{short}_advanced_table.tsv",
                  res.header_lines("Strategy picker: all strategies (advanced)"))  # fmt: skip
    (OUT / "cards.txt").write_text("\n\n".join(text).replace("\n\n\n", "\n\n"))
    grid = S.compare_all_histones(hs, ann, marks, base, rules)
    write_tsv(grid, OUT / "compare_all_histones.tsv",
              [S.standard_line(base), "Strategy picker: compare across all histones"])  # fmt: skip
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
