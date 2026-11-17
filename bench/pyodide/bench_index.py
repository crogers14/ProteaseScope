"""Option B in Pyodide: load gzipped reference + shared-prefix index, look up uniqueness."""
import gzip
import sys
import time

ROOT = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].endswith(".py") else "/proj"
sys.path.insert(0, f"{ROOT}/src")

from proteasescope import metrics  # noqa: E402
from proteasescope.chemistry import Chemistry  # noqa: E402
from proteasescope.digest import DigestSettings, digest  # noqa: E402
from proteasescope.fasta import parse_fasta, read_fasta  # noqa: E402
from proteasescope.rules import get_rule  # noqa: E402


def timed(label, fn):
    t = time.perf_counter()
    out = fn()
    print(f"{label}: {time.perf_counter() - t:.2f} s", flush=True)
    return out


raw = open(f"{ROOT}/data/precomputed/reference.fasta.gz", "rb").read()
text = timed("gunzip reference FASTA", lambda: gzip.decompress(raw).decode("ascii"))
ref = timed("parse reference FASTA", lambda: parse_fasta(text))
idx = timed("load shared-prefix index (npz)",
            lambda: metrics.SharedPrefixIndex.load(f"{ROOT}/data/precomputed/uniqueness_index.npz", ref))
d = timed("proteome digest trypsin <=1 mc",
          lambda: digest(ref, DigestSettings((get_rule("Trypsin", "psi-ms"),), missed_cleavages=1,
                                             chemistry=Chemistry())))
up, ug = timed("uniqueness lookup, all identifiable peptides", lambda: idx.unique_flags(d, d.identifiable))
print("unique fraction", float(up.mean()))
h = read_fasta(f"{ROOT}/data/histones/human_histones.fasta")
hd = digest(h, DigestSettings((get_rule("Trypsin", "psi-ms"),), missed_cleavages=1, min_length=4,
                              max_length=40, chemistry=Chemistry(propionylation="before_and_after")))
timed("histone peptides vs proteome (lookup incl. accession matching)",
      lambda: idx.unique_flags(hd, hd.identifiable))
one = ref.subset([ref.index_of("Q8WZ42")])
timed("titin alone: digest + lookup", lambda: idx.unique_flags(*(lambda x: (x, x.identifiable))(
    digest(one, DigestSettings((get_rule("Trypsin", "psi-ms"),), missed_cleavages=1)))))
