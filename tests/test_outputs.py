"""Outputs 1-5 and export headers, tested on the histone FASTA (CLAUDE.md sections 6-7).

Each output is checked against an independent brute-force computation in plain Python.
"""

import io
import json
from collections import Counter

import pandas as pd
import pytest

from proteasescope import report
from proteasescope.chemistry import Chemistry
from proteasescope.digest import DigestSettings, digest, parallel
from proteasescope.fasta import parse_fasta
from proteasescope.metrics import (
    OUTPUT1_COLUMNS,
    OUTPUT3_COLUMNS,
    SUMMARY_METRICS,
    Analysis,
    Background,
    compare,
    ranking_table,
    summary_table,
    uniqueness_background,
)
from proteasescope.rules import available_enzymes, get_rule

from .conftest import ROOT

CHEM = Chemistry(carbamidomethyl=False, met_removal="MetAP rule")


def _settings(*names, mode="expasy", **kw):
    kw.setdefault("chemistry", CHEM)
    return DigestSettings(tuple(get_rule(n, mode) for n in names), **kw)


@pytest.fixture(scope="module")
def bg(histones):
    return Background(histones)


@pytest.fixture(scope="module")
def argc(histones, bg):
    return Analysis(digest(histones, _settings("Arg-C", missed_cleavages=1)), bg)


@pytest.fixture(scope="module")
def argc_aspn(histones, bg):
    return Analysis(digest(histones, _settings("Arg-C", "Asp-N", missed_cleavages=1)), bg)


def _covered_positions(an):
    """Brute force: {(accession, UniProt position)} on identifiable peptides."""
    out = set()
    for row in an.peptide_table().itertuples():
        out.update((row.accession, p) for p in range(row.start, row.end + 1))
    return out


# ----- Output 1 and 2 ---------------------------------------------------------------------
def test_output1_columns_and_uniqueness(histones, argc):
    df = argc.peptide_table()
    assert list(df.columns) == OUTPUT1_COLUMNS
    assert len(df) == argc.summary["n_identifiable_peptides"]
    folded = {
        a: s.replace("I", "L") for a, s in zip(histones.accessions, histones.sequences, strict=True)
    }
    genes = dict(zip(histones.accessions, histones.meta["gene"], strict=True))
    for row in df.itertuples():
        hits = [a for a, s in folded.items() if row.peptide.replace("I", "L") in s]
        assert row.unique_protein == (len(hits) == 1)
        assert row.unique_gene == (len({genes[a] for a in hits}) == 1)
        assert histones.sequence(row.accession)[row.start - 1 : row.end] == row.peptide


def test_output2_is_unique_subset(argc):
    out2 = argc.unique_peptide_table()
    out1 = argc.peptide_table()
    assert list(out2.columns) == OUTPUT1_COLUMNS
    assert out2["unique_protein"].all()
    assert len(out2) == out1["unique_protein"].sum() == argc.summary["n_unique_peptides"]


def test_parallel_digest_column_names_components(histones, bg):
    a = digest(histones, _settings("Arg-C"))
    b = digest(histones, _settings("Asp-N"))
    df = Analysis(parallel(a, b), bg).peptide_table()
    labels = set(df["digest"])
    assert labels <= {"Arg-C (ExPASy)", "Asp-N (ExPASy)", "Arg-C (ExPASy); Asp-N (ExPASy)"}
    assert "Arg-C (ExPASy)" in labels and "Asp-N (ExPASy)" in labels


# ----- Output 3 ----------------------------------------------------------------------------
def test_output3_matches_brute_force(histones, argc):
    out3 = argc.protein_table(identifiable_only=False)
    assert list(out3.columns) == OUTPUT3_COLUMNS
    pep = argc.peptide_table()
    allp = argc.digest.to_frame(identifiable_only=False)
    n_ident = Counter(pep["accession"])
    n_uniq = Counter(pep.loc[pep["unique_protein"], "accession"])
    n_total = Counter(allp.loc[allp["missed_cleavages"] == 0, "accession"])
    covered = _covered_positions(argc)
    ucovered = set()
    for row in pep[pep["unique_protein"]].itertuples():
        ucovered.update((row.accession, p) for p in range(row.start, row.end + 1))
    for row in out3.itertuples():
        seq = histones.sequence(row.accession)
        mature = len(seq) - (1 if seq[1] in "ACGPSTV" and seq[0] == "M" else 0)
        assert row.n_identifiable_peptides == n_ident[row.accession]
        assert row.n_unique_peptides == n_uniq[row.accession]
        assert row.n_total_peptides == n_total[row.accession]
        cov = sum(1 for a, _ in covered if a == row.accession)
        ucov = sum(1 for a, _ in ucovered if a == row.accession)
        assert row.max_coverage_pct == pytest.approx(100 * cov / mature)
        assert row.unique_coverage_pct == pytest.approx(100 * ucov / mature)
        assert row.confident == (n_uniq[row.accession] >= 2)
    ident = argc.protein_table()
    assert (ident["n_identifiable_peptides"] >= 1).all()
    assert len(ident) == argc.summary["n_identifiable_proteins"]


# ----- Output 4 ----------------------------------------------------------------------------
def test_output4_statistics(histones, argc):
    stats = argc.digest_statistics()
    hist = stats["length_histogram"]
    assert hist["n_peptides"].sum() == argc.summary["n_total_peptides"]
    assert (hist["in_window"] == hist["length"].between(7, 30)).all()

    covered = _covered_positions(argc)
    total, cov = Counter(), Counter()
    for acc, seq in zip(histones.accessions, histones.sequences, strict=True):
        start = 2 if seq[0] == "M" and seq[1] in "ACGPSTV" else 1
        for pos in range(start, len(seq) + 1):
            total[seq[pos - 1]] += 1
            cov[seq[pos - 1]] += (acc, pos) in covered
    aa = stats["aa_detectability"].set_index("residue")
    for r in total:
        assert aa.loc[r, "occurrences"] == total[r]
        assert aa.loc[r, "fraction_detectable"] == pytest.approx(cov[r] / total[r])
    assert set(aa.index[aa["highlight"]]) == set("KSTYC") & set(total)
    s = stats["coverage_summary"]
    assert s["overall_coverage_pct"] == pytest.approx(100 * len(covered) / sum(total.values()))
    assert len(stats["coverage_per_protein"]) == len(histones)


# ----- Output 5 ----------------------------------------------------------------------------
def test_compare_with_itself_is_empty(argc):
    c = compare(argc, argc)
    assert (c.summary["difference"] == 0).all()
    assert all(v == 0 for v in c.counts.values())
    assert c.rescued_proteins.empty and c.lost_proteins.empty and c.sites.empty


def test_compare_matches_brute_force(argc, argc_aspn):
    c = compare(argc, argc_aspn)
    base_cov, alt_cov = _covered_positions(argc), _covered_positions(argc_aspn)
    gained, lost = alt_cov - base_cov, base_cov - alt_cov
    assert c.counts["residues_gained"] == len(gained)
    assert c.counts["residues_lost"] == len(lost)
    res = c.residues.set_index("accession")
    for acc in {a for a, _ in gained | lost}:
        assert res.loc[acc, "residues_gained"] == sum(1 for a, _ in gained if a == acc)
        assert res.loc[acc, "residues_lost"] == sum(1 for a, _ in lost if a == acc)
    for r in "KSTY":
        seq = {a: s for a, s in zip(argc.proteins.accessions, argc.proteins.sequences, strict=True)}
        g = {(a, p) for a, p in gained if seq[a][p - 1] == r}
        lo = {(a, p) for a, p in lost if seq[a][p - 1] == r}
        assert c.counts[f"{r}_sites_gained"] == len(g)
        assert c.counts[f"{r}_sites_lost"] == len(lo)
        got = c.sites[c.sites["residue"] == r]
        assert (
            set(
                zip(
                    got.loc[got.status == "gained", "accession"],
                    got.loc[got.status == "gained", "position"],
                    strict=True,
                )
            )
            == g
        )
    shift = {a: int(s[0] == "M" and s[1] in "ACGPSTV") for a, s in seq.items()}
    expected_mature = c.sites["position"] - c.sites["accession"].map(shift)
    assert (c.sites["position_mature"] == expected_mature).all()
    assert set(shift.values()) == {0, 1}  # both cases occur in the histone set

    b_pep = set(map(tuple, argc.peptide_table()[["accession", "start", "end"]].to_numpy()))
    a_pep = set(map(tuple, argc_aspn.peptide_table()[["accession", "start", "end"]].to_numpy()))
    assert c.peptides_gained == len(a_pep - b_pep) and c.peptides_lost == len(b_pep - a_pep)
    s = c.summary.set_index("metric")
    for m in SUMMARY_METRICS:
        assert s.loc[m, "difference"] == pytest.approx(argc_aspn.summary[m] - argc.summary[m])


def test_rescued_lost_and_gaining_unique(histones, bg):
    # Asp-N alone leaves some histones without identifiable peptides; Arg-C covers all.
    aspn = Analysis(digest(histones, _settings("Asp-N")), bg)
    argc0 = Analysis(digest(histones, _settings("Arg-C")), bg)
    c = compare(aspn, argc0)
    b_id = set(aspn.protein_table()["accession"])
    a_id = set(argc0.protein_table()["accession"])
    assert set(c.rescued_proteins["accession"]) == a_id - b_id
    assert set(c.lost_proteins["accession"]) == b_id - a_id
    assert len(a_id - b_id) > 0
    bu = aspn.protein_stats.set_index("accession")["n_unique_peptides"]
    au = argc0.protein_stats.set_index("accession")["n_unique_peptides"]
    assert set(c.proteins_gaining_unique["accession"]) == set(au.index[au > bu])
    assert (c.proteins_gaining_unique["gain"] > 0).all()


def test_summary_table(argc, argc_aspn):
    t = summary_table(argc, [argc_aspn])
    assert t["digest"].tolist() == [argc.label, argc_aspn.label]
    for m in SUMMARY_METRICS:
        assert t.loc[0, f"{m}_diff"] == 0
        assert t.loc[1, m] == argc_aspn.summary[m]


def test_ranking_table(histones, bg):
    base = _settings("Arg-C")
    enzymes = [get_rule(n, "expasy") for n in available_enzymes("expasy")]
    t = ranking_table(histones, base, enzymes, bg)
    n_other = len([e for e in enzymes if e not in base.rules])
    assert len(t) == 1 + 3 * n_other
    assert t.loc[0, "kind"] == "baseline"
    assert (t.filter(like="_diff").iloc[0] == 0).all()
    aspn = t[(t.enzyme == "Asp-N (ExPASy)")].set_index("kind")
    alone = Analysis(digest(histones, _settings("Asp-N")), bg)
    seq = Analysis(digest(histones, _settings("Arg-C", "Asp-N")), bg)
    par = Analysis(parallel(digest(histones, base), alone.digest), bg)
    for kind, an in (("alone", alone), ("sequential", seq), ("parallel", par)):
        for m in SUMMARY_METRICS:
            assert aspn.loc[kind, m] == pytest.approx(an.summary[m])
    assert "low confidence" in t.loc[t.enzyme == "Proteinase K (ExPASy)", "flags"].iloc[0]
    # sortable by any metric
    assert t.sort_values("K_sites_gained", ascending=False).iloc[0]["K_sites_gained"] >= 0


# ----- background ---------------------------------------------------------------------------
def test_uniqueness_background_concat():
    ref = parse_fasta(">a\nAAAAKCCCCK\n>b\nDDDDK\n")
    ref.name = "ref"
    up = parse_fasta(">a\nAAAAKCCCCK\n>c\nEEEEK\n>b\nDDDDKX\n")
    up.name = "upload"
    bg = uniqueness_background(ref, up)
    assert bg.proteins.accessions == ["a", "b", "c", "b"]  # same acc, different seq kept
    assert bg.proteins.name == "ref + upload (2 not already present)"


def test_digest_mapping_background(histones):
    s = _settings("Arg-C")
    d = digest(histones, s)
    sub = Analysis(d, Background(histones))
    dig = Analysis(d, Background(histones, method="digest"))
    assert dig.summary["n_unique_peptides"] >= sub.summary["n_unique_peptides"]
    assert "cleavage products" in dig.background.description


def test_parallel_keeps_smallest_missed_cleavages():
    # AAAARAAAAK has 1 missed cleavage under trypsin but 0 under Lys-C.
    ps = parse_fasta(">x\nAAAAKAAAARAAAAK\n")
    kw = dict(missed_cleavages=1, min_length=1, chemistry=CHEM)
    tryp = digest(ps, DigestSettings((get_rule("Trypsin", "psi-ms"),), **kw))
    lysc = digest(ps, DigestSettings((get_rule("Lys-C", "psi-ms"),), **kw))
    p = parallel(tryp, lysc)  # trypsin first: "first digest wins" would give 1
    got = dict(zip(p.sequences(), p.missed.tolist(), strict=True))
    assert got["AAAARAAAAK"] == 0
    assert got["AAAAKAAAAR"] == 1
    assert p.n_total_peptides == 4  # AAAAK, AAAAR, AAAAK (trypsin) + AAAARAAAAK (Lys-C)


# ----- export header and methods ----------------------------------------------------------
def test_settings_header_and_tsv(argc, argc_aspn, tmp_path):
    header = report.settings_header([argc, argc_aspn], "Output 5: Comparison")
    text = "\n".join(header)
    assert report.PREDICTION_NOTICE in text
    assert "Baseline:" in text and "Alternative 1:" in text
    for r in ("(?<=R)", "(?=D)"):
        assert r in text
    assert "Peptide length: 7-30 aa" in text and "Max missed cleavages: 1" in text
    assert "Initiator Met removal: MetAP rule" in text
    assert "maps to every entry whose sequence contains it" in text
    path = report.write_tsv(argc.peptide_table(), tmp_path / "out1.tsv", header)
    raw = path.read_text()
    assert raw.startswith("# ProteaseScope")
    back = pd.read_csv(io.StringIO(raw), sep="\t", comment="#")
    assert list(back.columns) == OUTPUT1_COLUMNS and len(back) == len(argc.peptide_table())
    methods = report.methods_text([argc, argc_aspn])
    assert "predictions" in methods and "(?<=R)" in methods and "MetAP rule" in methods


# ----- provenance ---------------------------------------------------------------------------
META = ROOT / "data" / "reference" / "metadata.json"


@pytest.mark.skipif(not META.exists(), reason="no downloaded reference data")
def test_downloaded_files_have_provenance():
    import hashlib

    files = json.loads(META.read_text())["files"]
    for name, rec in files.items():
        path = META.parent / name
        assert rec["source_url"].startswith("https://") and rec["download_date"]
        assert rec.get("release"), name
        if not path.exists():
            # Large raw FASTA downloads are not committed (.gitignore); fetch_data.py re-creates
            # them. Their provenance record above is still required.
            assert path.suffix == ".fasta", name
            continue
        if "sha256" in rec:
            assert hashlib.sha256(path.read_bytes()).hexdigest() == rec["sha256"], name


@pytest.mark.skipif(__import__("sys").platform == "emscripten", reason="no subprocesses in wasm")
def test_core_has_no_ui_imports():
    """The core must not import Shiny/plotly (CLAUDE.md section 2); checked in a fresh
    interpreter because pytest plugins may load shiny into this one."""
    import subprocess
    import sys

    code = (
        "import sys, importlib, pkgutil, proteasescope\n"
        "for m in pkgutil.iter_modules(proteasescope.__path__):\n"
        "    importlib.import_module('proteasescope.' + m.name)\n"
        "bad = [m for m in ('shiny', 'shinywidgets', 'plotly') if m in sys.modules]\n"
        "print(','.join(bad))\n"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                         check=True)  # fmt: skip
    assert out.stdout.strip() == ""
