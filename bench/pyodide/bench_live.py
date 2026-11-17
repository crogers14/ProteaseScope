"""Live-mode timings inside Pyodide: histone outputs, proteome digest, live uniqueness.

Run natively with `python bench_live.py /path/to/project` or in Pyodide via run.mjs.
"""

import sys
import time

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].endswith(".py") else "/proj"
sys.path.insert(0, f"{ROOT}/src")

import numpy as np  # noqa: E402

from proteasescope import histones as H  # noqa: E402
from proteasescope import metrics  # noqa: E402
from proteasescope.chemistry import Chemistry  # noqa: E402
from proteasescope.digest import DigestSettings, digest  # noqa: E402
from proteasescope.fasta import read_fasta  # noqa: E402
from proteasescope.rules import get_rule  # noqa: E402


def timed(label, fn):
    t = time.perf_counter()
    out = fn()
    print(f"{label}: {time.perf_counter() - t:.2f} s", flush=True)
    return out


PTM = DigestSettings((get_rule("Trypsin", "psi-ms"),), missed_cleavages=1, min_length=4,
                     max_length=40, chemistry=Chemistry(propionylation="before_and_after"))
PROTEOME = DigestSettings((get_rule("Trypsin", "psi-ms"),), missed_cleavages=1)

ps, ann, marks = timed("histone inputs", lambda: H.load_histone_inputs(f"{ROOT}/data/histones"))


def histone_all():
    an = [metrics.Analysis(digest(ps, PTM))] + [
        metrics.Analysis(digest(ps, PTM.with_rules(*PTM.rules, get_rule(n, "psi-ms"))))
        for n in ("Asp-N", "Glu-C", "Chymotrypsin")
    ]
    H.mark_visibility(an, marks, ann)
    H.isoform_discrimination(an, ann)
    H.lysine_cooccurrence(an, ann)
    H.tail_peptide_map(an, ann, "P68431")


timed("histone digest x4 + all histone outputs (histone-only background)", histone_all)

ref = timed("read reference FASTA (42,545 entries)",
            lambda: read_fasta(f"{ROOT}/data/reference/uniprot_human_reviewed_isoforms.fasta"))
d = timed("proteome digest, trypsin <=1 missed cleavage", lambda: digest(ref, PROTEOME))
an = metrics.Analysis(d, metrics.Background(ref))
timed("identifiable + coverage per protein (no uniqueness)",
      lambda: (an.identifiable, metrics.covered_residues(d, an.identifiable)))

# Option C: one protein's peptides, plain str.find over the folded proteome (no index)
folded = "\n".join(s.replace("I", "L") for s in ref.sequences)
idx = ref.index_of("P68431")
sel = an.identifiable & (d.protein == idx)
peps = [p.replace("I", "L") for p in d.sequences(sel)]


def scan():
    counts = []
    for p in peps:
        n, i = 0, folded.find(p)
        while i != -1 and n < 2:
            n += 1
            i = folded.find(p, i + 1)
        counts.append(n)
    return counts


timed(f"option C: str.find scan, {len(peps)} peptides of one protein", scan)
big = ref.index_of("Q8WZ42")  # titin
selb = an.identifiable & (d.protein == big)
peps = [p.replace("I", "L") for p in d.sequences(selb)]
timed(f"option C: str.find scan, {len(peps)} peptides of titin", scan)

# Option A: k-mer index + substring mapping for the whole proteome
try:
    timed("option A: build k-mer index (whole proteome)", lambda: metrics._kmer_index(ref))
    timed("option A: map all identifiable peptides", lambda: an.mapping)
    print("option A unique fraction:", float(np.mean(an.mapping["unique_protein"])))
except MemoryError as e:  # pragma: no cover
    print("option A: MemoryError", e)
