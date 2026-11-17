"""Histone module: annotation and known-mark tables (Phase 3, part 1).

Positions: ``position_uniprot`` uses UniProt numbering (initiator Met = 1);
``position_mature`` = UniProt - 1 when the initiator Met is removed (MetAP rule by default),
e.g. H3.1 K4 = UniProt 5.

Both curated tables carry a ``# status:`` header line. Files with ``status: DRAFT`` are refused
by the loaders until the user has reviewed them (CLAUDE.md section 8).

The histone set is the user-provided FASTA plus ``histone_isoforms.fasta`` (splice forms added
by user decision, e.g. macroH2A1.2 = O75367-1). Default histone views use entries whose
``category`` is ``canonical`` only; the one exception is H3.3 as a lysine co-occurrence default.
"""

from __future__ import annotations

import csv
import difflib
import io
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .chemistry import MetRemoval, met_removed
from .fasta import ProteinSet, read_fasta

ANNOTATION_COLUMNS = [
    "accession", "display_name", "family", "variant", "category", "notes",
    "gene", "pe", "cooccurrence_default",
]  # fmt: skip
#: Annotation categories (user decision 2026-11-14).
CATEGORIES = (
    "canonical",  # replication-dependent, reference sequence of its family
    "canonical_divergent",  # replication-dependent, sequence diverges from the reference
    "tissue_specific",  # expressed mainly in one tissue (e.g. testis)
    "replacement_variant",  # replication-independent variant
    "putative",  # uncertain existence or pseudogene product
)
#: Non-canonical entries allowed as lysine co-occurrence defaults (user decision 2026-11-14:
#: keep H3.3, as in CLAUDE.md section 8). All other default views stay canonical-only.
COOCCURRENCE_EXCEPTIONS = frozenset({"P84243"})
KNOWN_MARK_COLUMNS = [
    "accession", "position_uniprot", "position_mature", "residue", "modification",
    "evidence", "mod_abbrev", "feature", "evidence_type", "uniprot_note",
]  # fmt: skip
FAMILIES = ("H1", "H2A", "H2B", "H3", "H4")

#: Short names for lysine modifications (UniProt controlled vocabulary -> abbreviation).
MOD_ABBREV = {
    "N6-acetyllysine": "ac",
    "N6-methyllysine": "me1",
    "N6,N6-dimethyllysine": "me2",
    "N6,N6,N6-trimethyllysine": "me3",
    "N6-methylated lysine": "me",
    "N6-acetyl-N6-methyllysine": "acme",
    "N6-crotonyllysine": "cr",
    "N6-(2-hydroxyisobutyryl)lysine": "hib",
    "N6-(beta-hydroxybutyryl)lysine": "bhb",
    "N6-lactoyllysine": "la",
    "N6-glutaryllysine": "glu",
    "N6-methacryllysine": "mac",
    "N6-isonicotinyllysine": "inic",
    "N6-benzoyllysine": "bz",
    "N6-succinyllysine": "succ",
    "N6-butyryllysine": "bu",
    "N6-propionyllysine": "pr",
    "N6-malonyllysine": "ma",
    "Allysine": "allysine",
    "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin)": "ub",
    "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)": "sumo2",
    "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in UFM1)": "ufm1",
}

#: UniProt evidence codes, strongest first.
EVIDENCE_TYPES = {
    "ECO:0000269": "experimental",
    "ECO:0007744": "large-scale proteomics",
    "ECO:0000250": "by similarity",
    "ECO:0000305": "curator inference",
}

_FEATURE = re.compile(r'(MOD_RES|CROSSLNK) (\d+); /note="([^"]*)"(?:; /evidence="([^"]*)")?')


def parse_uniprot_features(text: str) -> pd.DataFrame:
    """Parse a UniProt TSV with ``Modified residue`` and ``Cross-link`` columns.

    Returns one row per feature: ``accession``, ``position_uniprot`` (int), ``feature``
    (MOD_RES/CROSSLNK), ``modification`` (first part of the note), ``uniprot_note``,
    ``evidence`` (UniProt evidence string) and ``evidence_type`` (strongest evidence class).
    Features with uncertain or ranged positions are skipped.
    """
    rows = []
    for rec in csv.DictReader(io.StringIO(text), delimiter="\t"):
        for col in ("Modified residue", "Cross-link"):
            for m in _FEATURE.finditer(rec.get(col) or ""):
                kind, pos, note, ev = m.groups()
                ev = ev or ""
                rows.append(
                    {
                        "accession": rec["Entry"],
                        "position_uniprot": int(pos),
                        "feature": kind,
                        "modification": note.split(";")[0].strip(),
                        "uniprot_note": note,
                        "evidence": ev,
                        "evidence_type": _evidence_type(ev),
                    }
                )
    return pd.DataFrame(rows)


def _evidence_type(evidence: str) -> str:
    for code, label in EVIDENCE_TYPES.items():
        if code in evidence:
            return label
    return "unspecified"


def known_marks_from_features(
    features: pd.DataFrame,
    proteins: ProteinSet,
    residues: str = "K",
    met_removal: MetRemoval | str = MetRemoval.METAP,
) -> pd.DataFrame:
    """Known marks on ``residues`` (default lysine), with UniProt and mature numbering.

    Raises ``ValueError`` if a feature position does not hold the expected residue in the
    sequence (i.e. FASTA and feature table come from different UniProt versions).
    """
    seqs = dict(zip(proteins.accessions, proteins.sequences, strict=True))
    df = features[features["accession"].isin(seqs)].copy()
    df["residue"] = [seqs[a][p - 1] for a, p in zip(df.accession, df.position_uniprot, strict=True)]
    bad = df[(df["feature"] == "CROSSLNK") & (df["residue"] != "K")]
    bad = pd.concat([bad, df[df["modification"].str.contains("lysine") & (df["residue"] != "K")]])
    if len(bad):
        raise ValueError(f"Feature positions do not match the sequences:\n{bad.head()}")
    df = df[df["residue"].isin(list(residues))]
    shift = {a: int(met_removed(s, met_removal)) for a, s in seqs.items()}
    df["position_mature"] = df["position_uniprot"] - df["accession"].map(shift)
    df["mod_abbrev"] = df["modification"].map(MOD_ABBREV).fillna("")
    order = {a: i for i, a in enumerate(proteins.accessions)}
    df = df.sort_values(
        ["accession", "position_uniprot", "modification"],
        key=lambda s: s.map(order) if s.name == "accession" else s,
    )
    return df[KNOWN_MARK_COLUMNS].drop_duplicates().reset_index(drop=True)


# ----- curated-table I/O with review status ---------------------------------------------
def read_status(path: str | Path) -> str:
    """The value of the ``# status:`` header line ("" if absent)."""
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#"):
            break
        if line.lower().startswith("# status:"):
            return line.split(":", 1)[1].strip()
    return ""


def _check_reviewed(path: Path, allow_draft: bool) -> None:
    status = read_status(path)
    if status.upper().startswith("DRAFT") and not allow_draft:
        raise PermissionError(
            f"{path.name} is a DRAFT awaiting user review. Edit its '# status:' line to "
            "'REVIEWED by <name> on <date>' once checked."
        )


def load_annotation(path: str | Path, allow_draft: bool = False) -> pd.DataFrame:
    """Load ``histone_annotation.tsv`` (refuses DRAFT files unless ``allow_draft``)."""
    path = Path(path)
    _check_reviewed(path, allow_draft)
    df = pd.read_csv(path, sep="\t", comment="#", dtype={"pe": "Int64"})
    df["cooccurrence_default"] = (
        df["cooccurrence_default"].astype(str).str.lower().map({"true": True, "false": False})
    )
    return df


def load_known_marks(path: str | Path, allow_draft: bool = False) -> pd.DataFrame:
    """Load ``known_marks.tsv`` (refuses DRAFT files unless ``allow_draft``)."""
    path = Path(path)
    _check_reviewed(path, allow_draft)
    return pd.read_csv(
        path,
        sep="\t",
        comment="#",
        keep_default_na=False,
        dtype={"position_uniprot": int, "position_mature": int},
    )


def validate_annotation(annotation: pd.DataFrame, proteins: ProteinSet) -> list[str]:
    """Problems found in an annotation table (empty list = consistent)."""
    problems = []
    accs = set(proteins.accessions)
    ann = set(annotation["accession"])
    if missing := sorted(accs - ann):
        problems.append(f"entries without annotation: {missing}")
    if extra := sorted(ann - accs):
        problems.append(f"annotated accessions not in FASTA: {extra}")
    if annotation["accession"].duplicated().any():
        problems.append("duplicate accessions")
    if bad := sorted(set(annotation["family"]) - set(FAMILIES)):
        problems.append(f"unknown families: {bad}")
    if annotation["display_name"].duplicated().any():
        dup = annotation.loc[annotation["display_name"].duplicated(), "display_name"].tolist()
        problems.append(f"duplicate display names: {dup}")
    if bad := sorted(set(annotation["category"]) - set(CATEGORIES)):
        problems.append(f"unknown categories: {bad} (allowed: {list(CATEGORIES)})")
    if annotation["cooccurrence_default"].isna().any():
        problems.append("cooccurrence_default must be true/false")
    co = annotation[annotation["cooccurrence_default"] == True]  # noqa: E712
    off = co[(co["category"] != "canonical") & ~co["accession"].isin(COOCCURRENCE_EXCEPTIONS)]
    if len(off):
        problems.append(
            "co-occurrence defaults must be canonical or a listed exception "
            f"(COOCCURRENCE_EXCEPTIONS): {off['accession'].tolist()}"
        )
    return problems


def default_entries(annotation: pd.DataFrame) -> list[str]:
    """Accessions shown by default in histone views: ``category == "canonical"`` only."""
    return annotation.loc[annotation["category"] == "canonical", "accession"].tolist()


def cooccurrence_entries(annotation: pd.DataFrame) -> list[str]:
    """Default entries for lysine co-occurrence: one canonical representative per core family,
    plus the explicit exceptions in :data:`COOCCURRENCE_EXCEPTIONS` (H3.3)."""
    ok = (annotation["category"] == "canonical") | annotation["accession"].isin(
        COOCCURRENCE_EXCEPTIONS
    )
    return annotation.loc[ok & annotation["cooccurrence_default"], "accession"].tolist()


def load_histone_set(directory: str | Path) -> ProteinSet:
    """The histone set: ``human_histones.fasta`` plus ``histone_isoforms.fasta`` if present."""
    directory = Path(directory)
    sets = [read_fasta(directory / "human_histones.fasta")]
    iso = directory / "histone_isoforms.fasta"
    if iso.exists():
        sets.append(read_fasta(iso))
    return ProteinSet.concat(sets) if len(sets) > 1 else sets[0]


def transfer_marks(
    marks: pd.DataFrame, proteins: ProteinSet, met_removal: MetRemoval | str = MetRemoval.METAP
) -> pd.DataFrame:
    """Copy marks of a canonical UniProt entry onto its splice isoforms in ``proteins``.

    A mark is transferred only when its residue lies in a block of identical sequence (aligned
    with :class:`difflib.SequenceMatcher`), so it keeps its sequence context. Transferred rows
    get ``evidence_type = "transferred from <accession> (identical sequence)"``.
    """
    seqs = dict(zip(proteins.accessions, proteins.sequences, strict=True))
    out = [marks]
    for acc, seq in seqs.items():
        base = acc.rsplit("-", 1)[0]
        if "-" not in acc or base not in seqs or acc in set(marks["accession"]):
            continue
        src = marks[marks["accession"] == base]
        sm = difflib.SequenceMatcher(None, seqs[base], seq, autojunk=False)
        pos_map = {}
        for a, b, n in sm.get_matching_blocks():
            pos_map.update({a + k + 1: b + k + 1 for k in range(n)})
        moved = src[src["position_uniprot"].isin(pos_map)].copy()
        moved["accession"] = acc
        moved["position_uniprot"] = moved["position_uniprot"].map(pos_map)
        moved["position_mature"] = moved["position_uniprot"] - int(met_removed(seq, met_removal))
        moved["evidence_type"] = f"transferred from {base} (identical sequence)"
        assert all(seq[p - 1] == "K" for p in moved["position_uniprot"])
        out.append(moved)
    return pd.concat(out, ignore_index=True)


# =========================================================================================
# Histone outputs (CLAUDE.md section 8). ``analyses[0]`` is the baseline in every function.
# =========================================================================================
#: Display windows (mature numbering) for histone tail peptide maps. View defaults only;
#: they do not affect any computation.
TAIL_WINDOWS = {"H1": (1, 40), "H2A": (1, 20), "H2B": (1, 35), "H3": (1, 45), "H4": (1, 30)}

SITE_STATUS_VS_BASELINE = ("seen in both", "gained", "lost", "unseen in both")


def load_histone_inputs(directory: str | Path) -> tuple[ProteinSet, pd.DataFrame, pd.DataFrame]:
    """Histone set, reviewed annotation and reviewed known marks (refuses DRAFT tables)."""
    directory = Path(directory)
    proteins = load_histone_set(directory)
    annotation = load_annotation(directory / "histone_annotation.tsv")
    marks = load_known_marks(directory / "known_marks.tsv")
    if problems := validate_annotation(annotation, proteins):
        raise ValueError(f"histone annotation problems: {problems}")
    return proteins, annotation, marks


def _peptides_by_entry(an) -> dict[int, list[tuple[int, int]]]:
    """Identifiable peptides per entry index as (start0, end0) local 0-based half-open."""
    d = an.digest
    m = an.identifiable
    out: dict[int, list[tuple[int, int]]] = {}
    for p, s, e in zip(d.protein[m].tolist(), d.start0[m].tolist(), d.end0[m].tolist(),
                       strict=True):  # fmt: skip
        out.setdefault(p, []).append((s, e))
    return out


def _check_same_proteins(analyses) -> ProteinSet:
    proteins = analyses[0].proteins
    if any(an.proteins is not proteins for an in analyses[1:]):
        raise ValueError("all analyses must digest the same ProteinSet object")
    return proteins


def _status(in_base: bool, in_alt: bool) -> str:
    return SITE_STATUS_VS_BASELINE[{(1, 1): 0, (0, 1): 1, (1, 0): 2, (0, 0): 3}[
        (int(in_base), int(in_alt))]]  # fmt: skip


# ----- 1. Known-mark visibility -----------------------------------------------------------
def mark_visibility(
    analyses,
    marks: pd.DataFrame,
    annotation: pd.DataFrame,
    entries: list[str] | None = None,
    modifications: list[str] | None = None,
) -> pd.DataFrame:
    """Known-mark visibility per site and digest.

    A site is covered when its lysine lies on >=1 identifiable peptide. Rows are per site
    (all modifications of the site listed in ``modifications``).

    Args:
        analyses: Baseline first, then alternatives (same ProteinSet).
        marks: Known-marks table.
        annotation: Histone annotation.
        entries: Accessions to include (default: canonical entries, :func:`default_entries`).
        modifications: ``mod_abbrev`` values to include (default: all).

    Returns:
        Columns ``digest, accession, display_name, site, position_uniprot, position_mature,
        modifications, covered, covering_peptides, n_K, status``. ``covering_peptides`` lists
        ``SEQUENCE(start-end)`` in mature numbering; ``n_K`` gives each peptide's K count
        (same order). ``status`` compares with the baseline (``baseline`` on its own rows).
    """
    proteins = _check_same_proteins(analyses)
    entries = default_entries(annotation) if entries is None else entries
    m = marks[marks["accession"].isin(entries)]
    if modifications is not None:
        m = m[m["mod_abbrev"].isin(modifications)]
    sites = (
        m.groupby(["accession", "position_uniprot", "position_mature"], sort=False)["mod_abbrev"]
        .agg(lambda s: ",".join(dict.fromkeys(s)))
        .reset_index()
        .rename(columns={"mod_abbrev": "modifications"})
    )
    names = dict(zip(annotation["accession"], annotation["display_name"], strict=True))
    seqs = dict(zip(proteins.accessions, proteins.sequences, strict=True))
    rows = []
    base_cov = None
    for k, an in enumerate(analyses):
        peps = _peptides_by_entry(an)
        covered_now = []
        for site in sites.itertuples():
            idx = proteins.index_of(site.accession)
            shift = int(an.digest.mature_offset[idx])
            i = site.position_uniprot - 1
            cov = [(s, e) for s, e in peps.get(idx, []) if s <= i < e]
            seq = seqs[site.accession]
            covered = bool(cov)
            covered_now.append(covered)
            rows.append(
                {
                    "digest": an.label,
                    "accession": site.accession,
                    "display_name": names.get(site.accession, site.accession),
                    "site": f"K{site.position_mature}",
                    "position_uniprot": site.position_uniprot,
                    "position_mature": site.position_mature,
                    "modifications": site.modifications,
                    "covered": covered,
                    "covering_peptides": "; ".join(
                        f"{seq[s:e]}({s + 1 - shift}-{e - shift})" for s, e in cov
                    ),
                    "n_K": ";".join(str(seq[s:e].count("K")) for s, e in cov),
                    "status": "baseline"
                    if k == 0
                    else _status(base_cov[len(covered_now) - 1], covered),
                }
            )
        if k == 0:
            base_cov = covered_now
    return pd.DataFrame(rows)


def mark_grid(visibility: pd.DataFrame) -> pd.DataFrame:
    """Marks x digests grid of status (baseline column shows covered / not covered)."""
    v = visibility.copy()
    v.loc[v["status"] == "baseline", "status"] = np.where(
        v.loc[v["status"] == "baseline", "covered"], "covered", "not covered"
    )
    order = list(dict.fromkeys(v["digest"]))
    grid = v.pivot_table(index=["display_name", "site"], columns="digest", values="status",
                         aggfunc="first", sort=False)  # fmt: skip
    return grid[order]


# ----- 2. Isoform discrimination ------------------------------------------------------------
def isoform_discrimination(
    analyses,
    annotation: pd.DataFrame,
    entries: list[str] | None = None,
    method: str = "substring",
) -> pd.DataFrame:
    """Pairwise discrimination of entries within each family, per digest.

    A discriminating peptide of entry A against entry B is an identifiable peptide of A that
    does not map to B: with ``method="substring"`` (default, user decision 2026-11-14) B's
    sequence does not contain it (I = L); with ``"digest"`` it is not among B's cleavage
    products.

    Status: ``resolved`` (both directions), ``one-way``, ``unresolved``, or
    ``identical sequence`` (I/L-folded sequences equal; never resolvable).

    Returns:
        Columns ``digest, family, accession_a, entry_a, accession_b, entry_b, status,
        n_discriminating_a, n_discriminating_b, discriminating_peptides_a,
        discriminating_peptides_b`` (peptides as ``SEQUENCE(start-end)``, mature numbering).
    """
    from .metrics import peptide_matches

    proteins = _check_same_proteins(analyses)
    entries = default_entries(annotation) if entries is None else entries
    ann = annotation[annotation["accession"].isin(entries)]
    names = dict(zip(annotation["accession"], annotation["display_name"], strict=True))
    folded = [s.replace("I", "L") for s in proteins.sequences]
    rows = []
    for an in analyses:
        d = an.digest
        m = an.identifiable
        pairs = peptide_matches(d, proteins if method == "substring" else d, m, method)
        hits: dict[int, set[int]] = {}
        for q, p in zip(pairs["query"].tolist(), pairs["protein"].tolist(), strict=True):
            hits.setdefault(q, set()).add(p)
        own = d.protein[m].tolist()
        s0, e0 = d.start0[m].tolist(), d.end0[m].tolist()
        by_entry: dict[int, list[int]] = {}
        for q, p in enumerate(own):
            by_entry.setdefault(p, []).append(q)

        def disc(a: int, b: int, *, _h=hits, _be=by_entry, _s=s0, _e=e0, _d=d) -> list[str]:
            seq = proteins.sequences[a]
            shift = int(_d.mature_offset[a])
            return [
                f"{seq[_s[q] : _e[q]]}({_s[q] + 1 - shift}-{_e[q] - shift})"
                for q in _be.get(a, [])
                if b not in _h.get(q, set())
            ]

        for fam, grp in ann.groupby("family", sort=False):
            accs = grp["accession"].tolist()
            for i, acc_a in enumerate(accs):
                for acc_b in accs[i + 1 :]:
                    a, b = proteins.index_of(acc_a), proteins.index_of(acc_b)
                    if folded[a] == folded[b]:
                        da, db, status = [], [], "identical sequence"
                    else:
                        da, db = disc(a, b), disc(b, a)
                        status = (
                            "resolved" if da and db else "one-way" if da or db else "unresolved"
                        )
                    rows.append(
                        {
                            "digest": an.label,
                            "family": fam,
                            "accession_a": acc_a,
                            "entry_a": names[acc_a],
                            "accession_b": acc_b,
                            "entry_b": names[acc_b],
                            "status": status,
                            "n_discriminating_a": len(da),
                            "n_discriminating_b": len(db),
                            "discriminating_peptides_a": "; ".join(da),
                            "discriminating_peptides_b": "; ".join(db),
                        }
                    )
    return pd.DataFrame(rows)


def discrimination_matrix(table: pd.DataFrame, digest: str, family: str) -> pd.DataFrame:
    """Entry x entry matrix for one digest and family.

    Cell (row A, column B) = number of A's identifiable peptides that do not map to B;
    the diagonal is empty. Pair status is in :func:`isoform_discrimination`.
    """
    t = table[(table["digest"] == digest) & (table["family"] == family)]
    names = list(dict.fromkeys(list(t["entry_a"]) + list(t["entry_b"])))
    mat = pd.DataFrame(pd.NA, index=names, columns=names, dtype="Int64")
    for r in t.itertuples():
        mat.loc[r.entry_a, r.entry_b] = r.n_discriminating_a
        mat.loc[r.entry_b, r.entry_a] = r.n_discriminating_b
    return mat


# ----- 3. Lysine co-occurrence ---------------------------------------------------------------
def _letters(n: int) -> list[str]:
    out = []
    for i in range(n):
        s, k = "", i
        while True:
            s = chr(65 + k % 26) + s
            k = k // 26 - 1
            if k < 0:
                break
        out.append(s)
    return out


def lysine_cooccurrence(
    analyses, annotation: pd.DataFrame, entries: list[str] | None = None
) -> dict[str, pd.DataFrame]:
    """Lysine sites and K-K pairs per entry, compared with the baseline.

    Definitions (CLAUDE.md section 8): a site is covered when the K lies on >=1 identifiable
    peptide; a pair is two K on the same identifiable peptide. Versus the baseline: ``new``
    site/pair (alternative only), ``kept`` (both), ``split`` pair (baseline pair whose two
    sites the alternative covers, but never on one peptide), ``lost`` pair (the alternative
    misses one or both sites); sites can also be ``lost``.

    Args:
        entries: Accessions (default: :func:`cooccurrence_entries`, i.e. H2A type 1,
            H2B type 1-C/E/F/G/I, H3.1, H3.3 and H4).

    Returns:
        ``letters`` (alternative digest -> letter), ``sites`` (per digest and K),
        ``pairs`` (per alternative and pair, with status) and ``matrix`` (per entry and
        K_i < K_j: ``baseline_pair`` and the letters of alternatives giving ``new``, ``split``
        or ``lost`` status), all K positions in mature numbering.
    """
    proteins = _check_same_proteins(analyses)
    entries = cooccurrence_entries(annotation) if entries is None else entries
    names = dict(zip(annotation["accession"], annotation["display_name"], strict=True))
    alts = analyses[1:]
    letters = dict(zip([a.label for a in alts], _letters(len(alts)), strict=True))

    def sites_pairs(an, idx):
        seq = proteins.sequences[idx]
        sites, pairs = set(), set()
        for s, e in _peptides_by_entry(an).get(idx, []):
            ks = [i for i in range(s, e) if seq[i] == "K"]
            sites.update(ks)
            pairs.update((ks[x], ks[y]) for x in range(len(ks)) for y in range(x + 1, len(ks)))
        return sites, pairs

    site_rows, pair_rows, mat_rows = [], [], []
    for acc in entries:
        idx = proteins.index_of(acc)
        seq = proteins.sequences[idx]
        shift = int(analyses[0].digest.mature_offset[idx])
        lys = [i for i, c in enumerate(seq) if c == "K"]
        b_sites, b_pairs = sites_pairs(analyses[0], idx)
        per_alt = {an.label: sites_pairs(an, idx) for an in alts}
        for an in analyses:
            sites = b_sites if an is analyses[0] else per_alt[an.label][0]
            for i in lys:
                site_rows.append(
                    {
                        "digest": an.label,
                        "letter": letters.get(an.label, ""),
                        "accession": acc,
                        "display_name": names.get(acc, acc),
                        "K": i + 1 - shift,
                        "position_uniprot": i + 1,
                        "covered": i in sites,
                        "status": "baseline"
                        if an is analyses[0]
                        else {
                            "seen in both": "kept",
                            "gained": "new",
                            "lost": "lost",
                            "unseen in both": "unseen",
                        }[_status(i in b_sites, i in sites)],  # fmt: skip
                    }
                )
        for x in range(len(lys)):
            for y in range(x + 1, len(lys)):
                i, j = lys[x], lys[y]
                in_base = (i, j) in b_pairs
                cell = {"new": [], "split": [], "lost": []}
                for label, (a_sites, a_pairs) in per_alt.items():
                    if (i, j) in a_pairs:
                        status = "kept" if in_base else "new"
                    elif in_base:
                        status = "split" if i in a_sites and j in a_sites else "lost"
                    else:
                        continue
                    if status != "kept":
                        cell[status].append(letters[label])
                    pair_rows.append(
                        {
                            "digest": label,
                            "letter": letters[label],
                            "accession": acc,
                            "display_name": names.get(acc, acc),
                            "K_i": i + 1 - shift,
                            "K_j": j + 1 - shift,
                            "position_uniprot_i": i + 1,
                            "position_uniprot_j": j + 1,
                            "status": status,
                        }
                    )
                mat_rows.append(
                    {
                        "accession": acc,
                        "display_name": names.get(acc, acc),
                        "K_i": i + 1 - shift,
                        "K_j": j + 1 - shift,
                        "baseline_pair": in_base,
                        "new_in": ",".join(cell["new"]),
                        "split_in": ",".join(cell["split"]),
                        "lost_in": ",".join(cell["lost"]),
                    }
                )
    return {
        "letters": pd.DataFrame({"letter": list(letters.values()), "digest": list(letters)}),
        "sites": pd.DataFrame(site_rows),
        "pairs": pd.DataFrame(
            pair_rows,
            columns=[
                "digest",
                "letter",
                "accession",
                "display_name",
                "K_i",
                "K_j",
                "position_uniprot_i",
                "position_uniprot_j",
                "status",
            ],
        ),  # fmt: skip
        "matrix": pd.DataFrame(mat_rows),
    }


# ----- 4. Peptide maps for histone tails ------------------------------------------------------
def tail_peptide_map(
    analyses, annotation: pd.DataFrame, accession: str, window: tuple[int, int] | None = None
) -> dict[str, pd.DataFrame]:
    """:func:`proteasescope.metrics.peptide_map` limited to the entry's N-terminal tail
    window (:data:`TAIL_WINDOWS` by family unless ``window`` is given; mature numbering)."""
    from .metrics import peptide_map

    fam = annotation.set_index("accession").loc[accession, "family"]
    return peptide_map(analyses, accession, window or TAIL_WINDOWS[fam])
