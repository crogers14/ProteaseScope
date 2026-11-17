"""Live whole-proteome digests in Pyodide: lean LiveResults memory/time (Histone PTM, 6 aa)."""
import gc
import sys
import time
from dataclasses import replace

sys.path.insert(0, "/proj/src")
from proteasescope import precomputed as P  # noqa: E402
from proteasescope.metrics import Background  # noqa: E402
from proteasescope.presets import PRESETS, recipes_for  # noqa: E402
from proteasescope.results import LiveResults  # noqa: E402
from proteasescope.rules import get_rule  # noqa: E402

t = time.perf_counter()
ref, idx, _ = P.load_reference("/proj/data/precomputed")
print(f"load reference + index: {time.perf_counter() - t:.1f} s", flush=True)
base = replace(PRESETS["histone_ptm"].baseline, min_length=6)
rules = [get_rule(n, "psi-ms") for n in ("Asp-N", "Glu-C", "Chymotrypsin")]
res = LiveResults(ref, Background(ref, index=idx), recipes_for(base, rules, ("parallel", "sequential")))
print("lean:", res.lean, "recipes:", len(res.recipes), flush=True)
t0 = time.perf_counter()
for i, r in enumerate(res.recipes):
    t = time.perf_counter()
    res.record(i)
    gc.collect()
    import js  # noqa: E402

    heap = js.eval("globalThis.__pyodideHeapMB ? globalThis.__pyodideHeapMB() : -1")
    print(f"  {i} {r.label}: {time.perf_counter() - t:.1f} s, wasm heap {heap} MB", flush=True)
print(f"total {time.perf_counter() - t0:.1f} s; summary(0) {res.summary(0)['n_identifiable_peptides']:,}")
