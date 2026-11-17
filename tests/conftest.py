"""Shared fixtures and paths."""

from pathlib import Path

import pytest

from proteasescope.fasta import read_fasta

ROOT = Path(__file__).resolve().parents[1]
HISTONE_FASTA = ROOT / "data" / "histones" / "human_histones.fasta"
GOLDEN_DIR = ROOT / "tests" / "data" / "protein_cleaver"
PSI_MS_OBO = ROOT / "data" / "reference" / "psi-ms.obo"


@pytest.fixture(scope="session")
def histones():
    """The 54-entry human histone FASTA."""
    return read_fasta(HISTONE_FASTA)
