"""Open-time of precomputed views inside Pyodide (CLAUDE.md target: < 2 s after first load)."""

import sys
import time

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].endswith(".py") else "/proj"
sys.path.insert(0, f"{ROOT}/src")

from proteasescope import precomputed as P  # noqa: E402
from proteasescope.digest import digest  # noqa: E402
from proteasescope.metrics import Analysis, Background  # noqa: E402
from proteasescope.presets import PRESETS  # noqa: E402

OUT = f"{ROOT}/data/precomputed"


def timed(label, fn):
    t = time.perf_counter()
    out = fn()
    print(f"{label}: {time.perf_counter() - t:.2f} s", flush=True)
    return out


ref, idx, manifest = timed("first load: reference + uniqueness index + manifest",
                           lambda: P.load_reference(OUT))
for key in PRESETS:
    pre = timed(f"open preset {key}: load arrays + ranking table",
                lambda k=key: (lambda p: (p, p.ranking_table()))(P.load_preset(OUT, k, ref.accessions))[0])
    timed(f"  {key}: per-protein table of one alternative", lambda p=pre: p.entry_table(5))

bg = Background(ref, index=idx)
one = ref.subset([ref.index_of("P68431")])
for key in PRESETS:
    timed(f"protein viewer (H3.1), live peptides + uniqueness, {key}",
          lambda k=key: Analysis(digest(one, PRESETS[k].baseline), Background(ref, index=idx)).peptide_table())
titin = ref.subset([ref.index_of("Q8WZ42")])
timed("protein viewer (titin), live peptides + uniqueness, proteome",
      lambda: Analysis(digest(titin, PRESETS["proteome"].baseline), bg).peptide_table())
