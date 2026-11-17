"""Regression checks against the earlier prototype (CLAUDE.md section 10).

ExPASy rules, 7-30 aa, 0 missed cleavages, mature numbering (MetAP rule), coverage over
the mature sequence, rounded to whole percent.
"""

import pytest

from proteasescope.chemistry import Chemistry
from proteasescope.digest import DigestSettings, coverage, digest
from proteasescope.rules import get_rule

CHEM = Chemistry(carbamidomethyl=False, met_removal="MetAP rule")


def _cov(histones, *enzymes):
    s = DigestSettings(
        tuple(get_rule(e, "expasy") for e in enzymes),
        missed_cleavages=0,
        min_length=7,
        max_length=30,
        chemistry=CHEM,
    )
    return coverage(digest(histones, s)).set_index("accession")["coverage_pct"]


@pytest.mark.parametrize(
    "accession,argc,argc_aspn",
    [
        ("P62807", 61, 84),  # H2B type 1-C/E/F/G/I
        ("P68431", 53, 61),  # H3.1
    ],
)
def test_prototype_coverage(histones, accession, argc, argc_aspn):
    assert round(_cov(histones, "Arg-C")[accession]) == argc
    assert round(_cov(histones, "Arg-C", "Asp-N")[accession]) == argc_aspn
