"""Histone annotation and known-mark tables (Phase 3, part 1)."""

import pandas as pd
import pytest

from proteasescope.histones import (
    CATEGORIES,
    KNOWN_MARK_COLUMNS,
    cooccurrence_entries,
    default_entries,
    known_marks_from_features,
    load_annotation,
    load_histone_set,
    load_known_marks,
    parse_uniprot_features,
    read_status,
    transfer_marks,
    validate_annotation,
)

from .conftest import ROOT

HIST = ROOT / "data" / "histones"
ANNOTATION = HIST / "histone_annotation.tsv"
MARKS = HIST / "known_marks.tsv"
FEATURES = ROOT / "data" / "reference" / "histone_features.tsv"
needs_features = pytest.mark.skipif(not FEATURES.exists(), reason="run fetch_data.py")


@pytest.fixture(scope="module")
def hset():
    """Histone set: user FASTA (54) + histone_isoforms.fasta (macroH2A1.2)."""
    return load_histone_set(HIST)


@pytest.fixture(scope="module")
def marks(hset):
    m = known_marks_from_features(parse_uniprot_features(FEATURES.read_text()), hset)
    return transfer_marks(m, hset)


def _has(marks, acc, mature, abbrev):
    sel = (marks.accession == acc) & (marks.position_mature == mature)
    return abbrev in set(marks.loc[sel, "mod_abbrev"])


@needs_features
def test_known_landmark_marks(marks):
    assert _has(marks, "P68431", 4, "me3")  # H3.1 K4me3
    assert _has(marks, "P68431", 27, "me3")  # H3K27me3
    assert _has(marks, "P68431", 36, "me3")  # H3K36me3
    assert _has(marks, "P84243", 9, "ac")  # H3.3 K9ac
    assert _has(marks, "P62805", 16, "ac")  # H4K16ac
    assert _has(marks, "P62805", 20, "me1")  # H4K20me1
    assert _has(marks, "P0C0S8", 119, "ub")  # H2AK119ub
    assert _has(marks, "P62807", 120, "ub")  # H2BK120ub
    h3k4 = marks[(marks.accession == "P68431") & (marks.position_mature == 4)].iloc[0]
    assert h3k4.position_uniprot == 5


@needs_features
def test_known_marks_are_lysines_with_consistent_numbering(hset, marks):
    seqs = dict(zip(hset.accessions, hset.sequences, strict=True))
    assert list(marks.columns) == KNOWN_MARK_COLUMNS
    for row in marks.itertuples():
        assert seqs[row.accession][row.position_uniprot - 1] == "K" == row.residue
        shift = row.position_uniprot - row.position_mature
        assert shift in (0, 1)
    # H2B.W (Q7Z2G1, starts MA...): Met removed under the MetAP rule
    assert (
        marks.loc[marks.accession == "Q7Z2G1", "position_uniprot"]
        - marks.loc[marks.accession == "Q7Z2G1", "position_mature"]
        == 1
    ).all()
    assert (marks["mod_abbrev"] != "").all()
    assert set(marks["evidence_type"]) <= {
        "experimental",
        "large-scale proteomics",
        "by similarity",
        "curator inference",
        "transferred from O75367 (identical sequence)",
    }


def test_histone_set_has_both_macroh2a1_splice_forms(hset):
    assert len(hset) == 55
    s11, s12 = hset.sequence("O75367"), hset.sequence("O75367-1")
    assert (len(s11), len(s12)) == (369, 372)
    assert s11[197:226] == "QVVQADIASIDSDAVVHPTNTDFYIGGEV"  # mH2A1.1 exon (UniProt 198-226)
    assert s12[197:229] == "NLIHSEISNLAGFEVEAIINPTNADIDLKDDL"  # mH2A1.2 exon
    assert s11[:197] == s12[:197] and s11[226:] == s12[229:]


@needs_features
def test_marks_transferred_to_macroh2a12_by_identity(hset, marks):
    m11 = marks[marks.accession == "O75367"]
    m12 = marks[marks.accession == "O75367-1"]
    assert len(m12) == len(m11) > 0
    shift = (m12.position_uniprot.to_numpy() - m11.position_uniprot.to_numpy()).tolist()
    # identical before the swapped exon, +3 after it (32-aa exon replaces 29 aa)
    for p11, d in zip(m11.position_uniprot, shift, strict=True):
        assert d == (0 if p11 < 198 else 3)
    assert (m12.evidence_type == "transferred from O75367 (identical sequence)").all()


@needs_features
def test_known_marks_file_is_up_to_date(marks):
    on_disk = load_known_marks(MARKS, allow_draft=True)
    pd.testing.assert_frame_equal(on_disk, marks, check_dtype=False)


def test_feature_mismatch_is_rejected(histones):
    text = (
        "Entry\tModified residue\tCross-link\n"
        'P68431\tMOD_RES 3; /note="N6-acetyllysine"; /evidence="ECO:0000269"\t\n'
    )
    with pytest.raises(ValueError, match="do not match"):
        known_marks_from_features(parse_uniprot_features(text), histones)


def test_drafts_are_refused_until_reviewed(tmp_path):
    for path, loader in ((ANNOTATION, load_annotation), (MARKS, load_known_marks)):
        if read_status(path).upper().startswith("DRAFT"):
            with pytest.raises(PermissionError, match="DRAFT"):
                loader(path)
        reviewed = tmp_path / path.name
        reviewed.write_text(
            path.read_text().replace(
                "# status: DRAFT", "# status: REVIEWED by tester on 2026-11-14 (was DRAFT", 1
            )
        )
        assert len(loader(reviewed)) > 0


def test_annotation_draft_is_consistent(hset):
    ann = load_annotation(ANNOTATION, allow_draft=True)
    assert "is_canonical" not in ann.columns
    assert validate_annotation(ann, hset) == []
    assert set(ann["category"]) == set(CATEGORIES)
    names = dict(zip(hset.accessions, hset.meta.protein_name, strict=True))
    for row in ann.itertuples():
        if row.family in ("H2A", "H2B", "H3", "H4"):
            assert row.family in names[row.accession] or row.accession == "P49450"
    pe5 = set(hset.meta.loc[hset.meta.pe == 5, "accession"])
    assert (ann.loc[ann.accession.isin(pe5), "category"] == "putative").all()


def test_user_category_decisions():
    ann = load_annotation(ANNOTATION, allow_draft=True).set_index("gene")
    assert ann.loc["H2AC1", "category"] == "tissue_specific"
    assert ann.loc["H2AC25", "category"] == "canonical_divergent"  # differs from H2A type 1
    assert "legacy alias" in ann.loc["H2AC25", "notes"]
    assert ann.loc["H2BC26", "category"] == "canonical_divergent"
    assert ann.loc["H2BC12L", "category"] == "replacement_variant"
    assert "Arg-C" in ann.loc["H2BC12L", "notes"]
    assert ann.loc["H3-7", "category"] == "putative"
    by_acc = ann.reset_index().set_index("accession")
    assert by_acc.loc["O75367", "display_name"] == "macroH2A1.1"
    assert by_acc.loc["O75367-1", "display_name"] == "macroH2A1.2"


def test_h2bc12l_k28r_adds_arg_c_site(hset):
    from proteasescope.rules import get_rule

    argc = get_rule("Arg-C", "psi-ms")
    assert 28 in argc.cut_sites(hset.sequence("P57053"))  # cut after UniProt R28
    assert 28 not in argc.cut_sites(hset.sequence("O60814"))  # K28 in H2B type 1-K


def test_default_views_use_canonical_only():
    ann = load_annotation(ANNOTATION, allow_draft=True)
    cat = dict(zip(ann.accession, ann.category, strict=True))
    defaults = default_entries(ann)
    assert defaults and all(cat[a] == "canonical" for a in defaults)
    assert "P84243" not in defaults  # H3.3 is a replacement variant
    # co-occurrence: canonical representatives plus H3.3 by explicit exception
    assert sorted(cooccurrence_entries(ann)) == sorted(
        ["P0C0S8", "P62807", "P68431", "P84243", "P62805"]
    )
    assert cat["P84243"] == "replacement_variant"
    bad = ann.copy()
    bad.loc[bad.accession == "Q6NXT2", "cooccurrence_default"] = True  # H3.5, not an exception
    assert any("must be canonical" in p for p in validate_annotation(bad, load_histone_set(HIST)))


def test_validate_annotation_reports_problems():
    ps = load_histone_set(HIST)
    ann = load_annotation(ANNOTATION, allow_draft=True)
    bad = ann.iloc[1:].copy()
    bad.loc[bad.index[0], "family"] = "H5"
    problems = validate_annotation(bad, ps)
    assert any("without annotation" in p for p in problems)
    assert any("unknown families" in p for p in problems)
