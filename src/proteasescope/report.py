"""Settings headers and Methods text for exports.

Every exported table starts with a ``#``-prefixed header that states the results are
predictions and lists the settings that produced them.
"""

from __future__ import annotations

import datetime as dt
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from . import __version__
from .chemistry import Chemistry
from .digest import DigestSettings
from .metrics import Analysis, Background
from .rules import PSI_MS_VERSION, CleavageRule

PREDICTION_NOTICE = (
    "PREDICTED (in silico) results. These are predictions from cleavage rules, not observed data."
)


def rule_lines(rule: CleavageRule) -> list[str]:
    """Description of one cleavage rule: label, regex, source/version, citation, flags."""
    src = rule.source
    if rule.source == "PSI-MS":
        src += f" {rule.accession}, psi-ms.obo v{PSI_MS_VERSION['data_version']}"
    lines = [f"{rule.label}: regex {rule.regex}  [{src}]", f"  Citation: {rule.citation}"]
    if rule.flags:
        lines.append(f"  Flags: {', '.join(rule.flags)}")
    return lines


def chemistry_lines(chem: Chemistry) -> list[str]:
    """Chemistry settings as ``key: value`` lines."""
    return [f"{k}: {v}" for k, v in chem.describe().items()]


def digest_lines(settings: DigestSettings) -> list[str]:
    """Settings of one digest (single or sequential)."""
    kind = "sequential (union of cut sites)" if settings.is_sequential else "single"
    mz = (
        f"on, {settings.mz_min:g}-{settings.mz_max:g} m/z, >=1 charge state in "
        f"{'/'.join(f'{z}+' for z in settings.charges)}"
        if settings.mz_filter
        else "off"
    )
    lines = [f"Digest: {settings.label} [{kind}]"]
    for r in settings.rules:
        lines += ["  " + x for x in rule_lines(r)]
    lines += [
        f"  Peptide length: {settings.min_length}-{settings.max_length} aa",
        f"  Max missed cleavages: {settings.missed_cleavages}",
        f"  m/z filter: {mz}",
    ]
    lines += ["  " + x for x in chemistry_lines(settings.chemistry)]
    return lines


def header_lines(
    digests: Sequence[tuple[DigestSettings, ...]],
    proteins_name: str,
    background_description: str,
    title: str = "",
    when: dt.datetime | None = None,
) -> list[str]:
    """Header lines (without ``#``) from digest settings alone.

    Args:
        digests: Settings of each digest, baseline first; a tuple with several settings is a
            parallel digest.
        proteins_name: Description of the digested entries.
        background_description: Uniqueness background and mapping method.
        title: Table name.
        when: Timestamp (defaults to now, local time).
    """
    when = when or dt.datetime.now().astimezone()
    lines = [
        f"ProteaseScope {__version__}" + (f" - {title}" if title else ""),
        PREDICTION_NOTICE,
        f"Generated: {when.isoformat(timespec='seconds')}",
        f"Proteins digested: {proteins_name}",
        f"Uniqueness background: {background_description}",
    ]
    for i, settings in enumerate(digests):
        if len(digests) == 1:
            role = "Digest"
        else:
            role = "Baseline" if i == 0 else f"Alternative {i}"
        if len(settings) > 1:
            lines.append(f"{role}: parallel digests (union of peptide sets)")
        else:
            lines.append(f"{role}:")
        for s in settings:
            lines += ["  " + x for x in digest_lines(s)]
    lines += [
        "Numbering: UniProt (initiator Met = 1); *_mature = UniProt - 1 when Met is removed.",
        "Masses: neutral monoisotopic and average (ExPASy average residue masses); "
        "m/z = (M + z*1.007276)/z.",
        "Counts: peptides are counted by position (entry, start, end).",
    ]
    return lines


def settings_header(
    analyses: Sequence[Analysis] | Analysis,
    title: str = "",
    background: Background | None = None,
    when: dt.datetime | None = None,
) -> list[str]:
    """Header lines (without ``#``) for an export computed from analyses.

    Args:
        analyses: The analysis/analyses the table was computed from (baseline first).
        title: Table name, e.g. "Output 1: Identifiable peptides".
        background: Uniqueness background (defaults to the first analysis's).
        when: Timestamp (defaults to now, local time).
    """
    if isinstance(analyses, Analysis):
        analyses = [analyses]
    background = background or analyses[0].background
    return header_lines(
        [an.digest.settings for an in analyses],
        str(analyses[0].proteins.name or len(analyses[0].proteins)),
        background.description,
        title,
        when,
    )


def tsv_bytes(df: pd.DataFrame, header: Sequence[str]) -> bytes:
    """TSV text (UTF-8) with ``# ``-prefixed header lines."""
    head = "".join(f"# {line}\n" for line in header)
    return (head + df.to_csv(sep="\t", index=False)).encode()


def excel_bytes(sheets: dict[str, pd.DataFrame], header: Sequence[str]) -> bytes:
    """An .xlsx workbook: a "Settings" sheet with the header lines, then one sheet per table.

    Needs ``openpyxl``.
    """
    import io

    buf = io.BytesIO()
    with pd.ExcelWriter(buf, engine="openpyxl") as xw:
        pd.DataFrame({"ProteaseScope export": list(header)}).to_excel(
            xw, sheet_name="Settings", index=False
        )
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name[:31], index=False)
    return buf.getvalue()


def write_tsv(df: pd.DataFrame, path: str | Path, header: Sequence[str]) -> Path:
    """Write ``df`` as TSV preceded by ``# ``-prefixed header lines."""
    path = Path(path)
    with path.open("w", newline="") as fh:
        for line in header:
            fh.write(f"# {line}\n")
        df.to_csv(fh, sep="\t", index=False)
    return path


def methods_text(analyses: Sequence[Analysis] | Analysis) -> str:
    """A Methods paragraph: rule sources and versions, database release, settings, software."""
    if isinstance(analyses, Analysis):
        analyses = [analyses]
    return methods_text_from(
        [an.digest.settings for an in analyses],
        str(analyses[0].proteins.name or len(analyses[0].proteins)),
        analyses[0].background.description,
    )


def methods_text_from(
    digests: Sequence[tuple[DigestSettings, ...]], proteins_name: str, background_description: str
) -> str:
    """Methods paragraph from digest settings alone (baseline first; tuples = parallel)."""
    parts_all = []
    rules: dict[str, CleavageRule] = {}
    for settings in digests:
        parts = []
        for s in settings:
            for r in s.rules:
                rules[r.label] = r
            parts.append(
                f"{s.label} ({s.min_length}-{s.max_length} aa, "
                + (
                    f"up to {s.missed_cleavages} missed cleavage"
                    f"{'s' if s.missed_cleavages != 1 else ''}"
                    if s.missed_cleavages
                    else "no missed cleavages"
                )
                + (
                    f", m/z {s.mz_min:g}-{s.mz_max:g} for at least one of charges "
                    f"{', '.join(map(str, s.charges))}"
                    if s.mz_filter
                    else ""
                )
                + ")"
            )
        parts_all.append(" in parallel with ".join(parts))
    chem = digests[0][0].chemistry
    rule_txt = "; ".join(
        f"{lab}: {r.regex} ({r.source}{' ' + r.accession if r.accession else ''})"
        for lab, r in rules.items()
    )
    has_psi = any(r.source == "PSI-MS" for r in rules.values())
    psi = (
        f" PSI-MS rules were taken from psi-ms.obo version {PSI_MS_VERSION['data_version']}."
        if has_psi
        else ""
    )
    cm = (
        "carbamidomethylation of cysteine (+57.021464 Da) as a fixed modification"
        if chem.carbamidomethyl
        else "no fixed modifications"
    )
    mz = " and m/z window" if any(s.mz_filter for d in digests for s in d) else ""
    return (
        f"Proteolytic digestion was predicted in silico with ProteaseScope {__version__}. "
        f"Proteins ({proteins_name}) were digested with: {'; '.join(parts_all)}. Cleavage "
        f"rules (zero-width regular expressions): {rule_txt}.{psi} Sequential digests used the "
        f"union of the cut sites of their rules; parallel digests the union of their peptide "
        f"sets. Chemistry: {cm}; derivatization: {chem.propionylation.label.lower()}; "
        f"initiator methionine removal: {chem.met_removal}. Peptides within the length "
        f"window{mz} were considered identifiable. Uniqueness was assessed against "
        f"{background_description}. All results are predictions and were not observed "
        "experimentally."
    )
