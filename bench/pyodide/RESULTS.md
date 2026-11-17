# Pyodide benchmark results (Phases 4–5, 2026-11-14)

The Phase 4 tables below used Pyodide 314.0.7. The app ships on **Pyodide 0.27.7**, bundled with Shinylive 0.10.16 (Python 3.12, numpy 2.0.2, pandas 2.2.3); its results are in the Phase 5 section at the end. `package.json` now pins 0.27.7.

Pyodide 314.0.7 (Python 3.14) under Node 25.9 on an Apple M3 Pro, so expect somewhat slower times in a browser. Native timings are from the same machine (CPython 3.13). Reference: UniProt 2026_04 reviewed human proteome with isoforms, 42,545 entries and 24.5 M residues.

Run with `npm install`, then `node run.mjs <bench_*.py> ../..` from this folder.

## Live uniqueness options (`bench_live.py`, `bench_index.py`)

| Operation | Native | Pyodide |
|---|---|---|
| Histone digest ×4 + all histone outputs | 0.05 s | 0.17 s |
| Read reference FASTA | 0.21 s | 0.43 s |
| Proteome digest, trypsin ≤1 missed cleavage | 0.96 s | 1.7 s |
| **A.** Build k-mer index + map all 2.84 M peptides | 6.0 s | 16.5 s |
| **C.** `str.find` scan: one typical protein / titin | 0.03 s / 21 s | 0.04 s / 26 s |
| **B.** Load shared-prefix index (6.5 MB) | — | 0.10 s |
| **B.** Lookup, all 2.84 M peptides | 0.03 s | 0.14 s |
| **B.** Lookup, histone set / titin | — | 0.02 s |

**Chosen: B**, the precomputed shared-prefix index. It gives exactly the substring-mapping result for any enzyme or settings; this was verified on 5.6 M proteome peptides and against brute force. Entries uploaded by the user that are not in the reference fall back to live mapping.

## Precomputed views (`bench_precomputed.py`)

| Operation | Pyodide |
|---|---|
| First load: reference + index + manifest (after download) | 0.66 s |
| Open a preset (arrays + ranking table) | 0.09–0.14 s |
| Per-protein table of one digest | < 0.01 s |
| Protein viewer, live peptides + uniqueness (H3.1; titin) | 0.03–0.05 s |

Download time for the ~80 MB of shipped data is not included.

## wasm32 portability

In Pyodide `np.intp` is int32. `np.repeat` with int64 counts and `ufunc.reduceat` with int64 indices raise `TypeError`, so the core casts both to `np.intp`. `bincount`, fancy indexing, `ufunc.at`, `searchsorted` and `lexsort` accept int64.

## Phase 5: the app's runtime (Pyodide 0.27.7) and the real browser

| Operation | Pyodide 0.27.7 (Node) | Browser (Shinylive, local server) |
|---|---|---|
| Histone digests ×4 + all histone outputs | 0.11 s | — |
| Proteome digest, trypsin ≤1 missed cleavage | 1.7 s | — |
| App ready with precomputed Proteome view | — | about 4 s after page load |
| Switch preset (precomputed) | — | about 1 s |
| Upload a 2-entry FASTA; live digests, uniqueness against the proteome | — | 5.0 s |
| **Live whole proteome**, Histone PTM baseline at 6 aa + 3 enzymes × 2 ways (7 digests) | 37.8 s, heap peak 1.84 GB | 39 s |

**Memory.** Before Phase 5's changes the live whole-proteome run failed in the browser with `_ArrayMemoryError` (all 7 digests were kept in memory). The fixes:
- `LiveResults` is lean for inputs over 2 M residues: it keeps the baseline digest plus compact per-protein records, and releases every other digest.
- Digest positions are stored as int32.
- Coverage counters are int32.
- Cut sites are collected as numpy chunks rather than Python lists.

Together these cut the heap peak from 2.43 GB to 1.84 GB.

**numpy 2.0.2 in wasm32** additionally rejects int64 input to `np.bincount` and `np.take`, not only `np.repeat`/`reduceat`. The core uses `digest.bincount()`, which casts to `np.intp`.

**Test suite inside Pyodide 0.27.7** (`run_tests.py`): 98 passed. The pyteomics module is skipped (test-only dependency), as is the subprocess-based import test.
