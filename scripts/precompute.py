"""Build the precomputed app data in data/precomputed/ (CLAUDE.md section 11, Phase 4).

For every preset: the baseline and every enzyme alone / in parallel with the baseline /
sequential after it, on the whole reviewed human proteome with isoforms. Uniqueness uses the
precomputed shared-prefix index (substring mapping, I = L).

Usage:
    uv run python scripts/precompute.py [--presets proteome histone_ptm ...]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from proteasescope import precomputed as P  # noqa: E402
from proteasescope.fasta import read_reference_fasta  # noqa: E402
from proteasescope.metrics import Background, SharedPrefixIndex  # noqa: E402
from proteasescope.presets import PRESETS  # noqa: E402

REF = ROOT / "data" / "reference" / "uniprot_human_reviewed_isoforms.fasta"
OUT = ROOT / "data" / "precomputed"
SIZE_LIMIT = 25 * 1024 * 1024


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--presets", nargs="*", default=list(PRESETS), choices=list(PRESETS))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    timings: dict[str, float] = {}

    t = time.perf_counter()
    ref = read_reference_fasta(REF)
    rec = json.loads((REF.parent / "metadata.json").read_text())["files"][REF.name]
    print(f"reference: {len(ref):,} entries ({time.perf_counter() - t:.1f} s)", flush=True)

    t = time.perf_counter()
    index = SharedPrefixIndex.build(ref)
    timings["uniqueness_index_build"] = round(time.perf_counter() - t, 1)
    print(f"uniqueness index: {timings['uniqueness_index_build']} s", flush=True)
    files = P.write_reference(ref, REF, index, OUT)
    background = Background(ref, index=index)

    presets = {}
    for key in args.presets:
        preset = PRESETS[key]
        result = P.build_preset(preset, ref, background)
        files |= P.write_preset(result, preset, OUT)
        timings[f"preset_{key}"] = round(result["seconds"], 1)
        presets[key] = {
            "name": preset.name,
            "mode": preset.mode,
            "notes": preset.notes,
            "n_digests": len(result["digests"]),
            "baseline": result["digests"][0]["settings"],
        }
        print(f"preset {key}: {len(result['digests'])} digests in {result['seconds']:.0f} s",
              flush=True)  # fmt: skip
    ref_meta = {k: rec[k] for k in ("description", "source_url", "release", "release_date",
                                    "download_date", "sha256", "entries")}  # fmt: skip
    P.write_manifest(OUT, files, ref_meta, presets, timings)
    print("\nfiles:")
    for name, f in sorted(files.items()):
        flag = "  OVER 25 MB LIMIT" if f["bytes"] > SIZE_LIMIT else ""
        print(f"  {name:32s} {f['bytes'] / 1e6:7.2f} MB{flag}")
    print(f"  total {sum(f['bytes'] for f in files.values()) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
