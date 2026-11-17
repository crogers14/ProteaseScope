"""Rule equivalence with pyteomics (test-time only) and with psi-ms.obo."""

import random
import re

import pytest

from proteasescope import rules
from proteasescope.rules import available_enzymes, custom_rule, get_rule

from .conftest import PSI_MS_OBO

parser = pytest.importorskip("pyteomics.parser")

# vendored (mode, name) -> (pyteomics dict, key, exception key)
PYTEOMICS_EQUIVALENT = {
    ("psi-ms", "Trypsin"): ("psims", "Trypsin", None),
    ("psi-ms", "Trypsin/P"): ("psims", "Trypsin/P", None),
    ("psi-ms", "Arg-C"): ("psims", "Arg-C", None),
    ("psi-ms", "Lys-C"): ("psims", "Lys-C", None),
    ("psi-ms", "Lys-C/P"): ("psims", "Lys-C/P", None),
    ("psi-ms", "Asp-N"): ("psims", "Asp-N", None),
    ("psi-ms", "Asp-N ambic"): ("psims", "Asp-N ambic", None),
    ("psi-ms", "Glu-C"): ("psims", "glutamyl endopeptidase", None),
    ("psi-ms", "Chymotrypsin"): ("psims", "Chymotrypsin", None),
    ("expasy", "Trypsin"): ("expasy", "trypsin", "trypsin_exception"),
    ("expasy", "Arg-C"): ("expasy", "arg-c", None),
    ("expasy", "Lys-C"): ("expasy", "lysc", None),
    ("expasy", "Asp-N"): ("expasy", "asp-n", None),
    ("expasy", "Glu-C"): ("expasy", "glutamyl endopeptidase", None),
    ("expasy", "Chymotrypsin high specificity"): ("expasy", "chymotrypsin high specificity", None),
    ("expasy", "Chymotrypsin low specificity"): ("expasy", "chymotrypsin low specificity", None),
    ("expasy", "Proteinase K"): ("expasy", "proteinase k", None),
}
# Rules with no pyteomics counterpart (checked by explicit cases below).
NO_PYTEOMICS = {("psi-ms", "Lys-N"), ("expasy", "Lys-N"), ("literature", "ProAlanase"),
                ("literature", "Arg-C Ultra"), ("literature", "Recombinant chymotrypsin")}

ALPHABET_UNIFORM = "ACDEFGHIKLMNPQRSTVWYBZXUO"
ALPHABET_ENRICHED = "KRPDECWMHYFLKRPDEW"  # residues used by rules and their exceptions


def _random_sequences(n=10_000, seed=20261114):
    rng = random.Random(seed)
    out = []
    for i in range(n):
        alphabet = ALPHABET_ENRICHED if i % 2 else ALPHABET_UNIFORM
        out.append("".join(rng.choice(alphabet) for _ in range(rng.randint(1, 60))))
    return out


SEQUENCES = _random_sequences()


def _pyteomics_sites(seq, regex, exception):
    peps = parser.icleave(seq, regex, 0, min_length=1, regex=True, exception=exception)
    return sorted({start for start, _ in peps} - {0})


def test_every_vendored_rule_is_covered():
    assert set(rules.RULES) == set(PYTEOMICS_EQUIVALENT) | NO_PYTEOMICS


@pytest.mark.parametrize("key", sorted(PYTEOMICS_EQUIVALENT), ids=lambda k: f"{k[0]}:{k[1]}")
def test_identical_cut_positions_to_pyteomics(key):
    which, pt_key, exc_key = PYTEOMICS_EQUIVALENT[key]
    table = parser.psims_rules if which == "psims" else parser.expasy_rules
    regex = table[pt_key]
    exception = parser.expasy_rules[exc_key] if exc_key else None
    rule = rules.RULES[key]
    for seq in SEQUENCES:
        assert rule.cut_sites(seq) == _pyteomics_sites(seq, regex, exception), (key, seq)


def test_rules_are_zero_width():
    for seq in SEQUENCES[:2000]:
        for rule in rules.RULES.values():
            assert all(m.end() == m.start() for m in rule.pattern.finditer(seq)), rule.name


def test_lys_n_cuts_before_k():
    for mode in rules.MODES:
        rule = get_rule("Lys-N", mode)
        assert rule.cut_sites("AKAKKPKA") == [1, 3, 4, 6]


def _obo_terms():
    text = PSI_MS_OBO.read_text()
    terms = {}
    for block in text.split("\n\n"):
        m = re.search(r"^id: (MS:\d+)$", block, re.M)
        if m:
            terms[m.group(1)] = block
    return terms


@pytest.mark.skipif(not PSI_MS_OBO.exists(), reason="psi-ms.obo not downloaded")
def test_psi_ms_regexes_match_obo():
    terms = _obo_terms()
    assert f"data-version: {rules.PSI_MS_VERSION['data_version']}" in PSI_MS_OBO.read_text()
    for (mode, name), rule in rules.RULES.items():
        if mode != "psi-ms":
            continue
        agent = terms[rule.accession]
        rx = re.search(r"^relationship: has_regexp (MS:\d+) ! (.*)$", agent, re.M)
        if name == "Lys-N":
            # CV inconsistency: links to the Lys-C regex ID but comments "(?=K)".
            assert rx.group(1) == "MS:1001335"
            assert rx.group(2).strip() == "(?=K)" == rule.regex
            continue
        rx_name = re.search(r"^name: (.*)$", terms[rx.group(1)], re.M).group(1)
        assert rx_name.replace("\\!", "!") == rule.regex, name


def test_available_enzymes_v1_list():
    assert available_enzymes("psi-ms") == [
        "Trypsin",
        "Trypsin/P",
        "Arg-C",
        "Lys-C",
        "Lys-C/P",
        "Lys-N",
        "Asp-N",
        "Asp-N ambic",
        "Glu-C",
        "Chymotrypsin",
        "Chymotrypsin high specificity",
        "Chymotrypsin low specificity",
        "Proteinase K",
        "ProAlanase",
        "Arg-C Ultra",
        "Recombinant chymotrypsin",
    ]
    assert available_enzymes("expasy")[-3:] == ["ProAlanase", "Arg-C Ultra",
                                                "Recombinant chymotrypsin"]  # fmt: skip
    assert get_rule("Proteinase K", "psi-ms").source == "ExPASy"
    assert "low confidence" in get_rule("Proteinase K", "expasy").flags
    with pytest.raises(KeyError):
        get_rule("Trypsin/P", "expasy")


def test_known_source_differences():
    # Arg-C: PSI-MS skips R-P, ExPASy cuts it
    assert get_rule("Arg-C", "psi-ms").cut_sites("AARPAA") == []
    assert get_rule("Arg-C", "expasy").cut_sites("AARPAA") == [3]
    # Glu-C: PSI-MS does not cut after EE
    assert get_rule("Glu-C", "psi-ms").cut_sites("AEEAEA") == [2, 5]
    assert get_rule("Glu-C", "expasy").cut_sites("AEEAEA") == [2, 3, 5]
    # Asp-N: PSI-MS also before B
    assert get_rule("Asp-N", "psi-ms").cut_sites("AABAADA") == [2, 5]
    assert get_rule("Asp-N", "expasy").cut_sites("AABAADA") == [5]


def test_custom_rule():
    r = custom_rule(after="kr", except_before_p=True)
    assert r.regex == "(?<=[KR])(?!P)"
    assert r.cut_sites("AKPARA") == [5]
    r2 = custom_rule(after="E", before="D")
    assert r2.regex == "(?<=[E])|(?=[D])"
    assert r2.cut_sites("AEADA") == [2, 3]
    r3 = custom_rule(after="E", before="D", except_before_p=True)
    assert r3.regex == "(?:(?<=[E])|(?=[D]))(?!P)"
    with pytest.raises(ValueError):
        custom_rule(after="K#")
    with pytest.raises(ValueError):
        custom_rule()


# ----- ProAlanase (literature-defined) ----------------------------------------------------------
def test_proalanase_cuts_after_every_p_and_a():
    rule = get_rule("ProAlanase", "psi-ms")
    assert get_rule("ProAlanase", "expasy") is rule
    assert (rule.source, rule.mode, rule.regex) == ("Literature", "literature", "(?<=[PA])")
    for seq in SEQUENCES + ["PPPAAPA", "APKPRPPDAE"]:
        expected = [i + 1 for i, aa in enumerate(seq[:-1]) if aa in "PA"]
        assert rule.cut_sites(seq) == expected, seq
    # no restriction on the next residue: cuts before P, K, R, D, E alike
    assert rule.cut_sites("GAPGPKAKPRAD") == [2, 3, 5, 7, 9, 11]


def test_proalanase_metadata():
    rule = get_rule("ProAlanase")
    assert "Samodova" in rule.citation and "to verify" in rule.citation
    assert "Promega" in rule.citation
    assert "P-P" in rule.notes and "unverified" in rule.notes
    assert rule.ui_note.startswith("Acidic digestion (around pH 1.5–2)")
    assert "≤1–2 missed cleavages" in rule.ui_note
    assert get_rule("Trypsin").ui_note == ""


@pytest.mark.skipif(not PSI_MS_OBO.exists(), reason="psi-ms.obo not downloaded")
def test_proalanase_not_in_psi_ms_obo():
    """Documents the source check: psi-ms.obo has no ProAlanase cleavage agent, so the rule is
    literature-defined. If a later CV release adds one, this fails and the rule should move."""
    agents = [b for b in _obo_terms().values() if "is_a: MS:1001045" in b]
    assert len(agents) > 20
    names = [re.search(r"^name: (.*)$", b, re.M).group(1).lower() for b in agents]
    assert not [n for n in names if "alanase" in n or "alanine" in n]


# ----- Arg-C Ultra and recombinant chymotrypsin (Turner et al. 2026) ------------------------------
def test_arg_c_ultra_equals_psi_ms_arg_c():
    ultra = get_rule("Arg-C Ultra", "psi-ms")
    assert get_rule("Arg-C Ultra", "expasy") is ultra
    assert (ultra.source, ultra.mode) == ("Literature", "literature")
    argc = get_rule("Arg-C", "psi-ms")
    for seq in SEQUENCES:
        assert ultra.cut_sites(seq) == argc.cut_sites(seq), seq
    assert ultra.cut_sites("AARPAARAAK") == [7]  # not before P, never after K


def test_recombinant_chymotrypsin_rule():
    rule = get_rule("Recombinant chymotrypsin")
    assert rule.regex == "(?<=[FYLM])(?!P)"
    for seq in SEQUENCES + ["FPWAHAYALMKMP"]:
        expected = [i + 1 for i in range(len(seq) - 1)
                    if seq[i] in "FYLM" and seq[i + 1] != "P"]  # fmt: skip
        sites = rule.cut_sites(seq)
        assert sites == expected, seq
        assert all(seq[i - 1] not in "WH" and seq[i] != "P" for i in sites)
    assert rule.cut_sites("AWAHAFPAYA") == [9]


def test_turner_enzyme_metadata():
    ultra, chymo = get_rule("Arg-C Ultra"), get_rule("Recombinant chymotrypsin")
    assert (ultra.suggested_missed_cleavages, chymo.suggested_missed_cleavages) == (1, 3)
    for r in (ultra, chymo):
        assert "Turner et al. (2026) Anal Chem 98:18694" in r.citation
        assert "10.1021/acs.analchem.6c01147" in r.citation
        assert r.evidence.startswith("Observed in Turner et al. 2026 (HEK293T histones): ")
    assert "Promega Arg-C Ultra" in ultra.citation
    assert "rChymoselect" in chymo.citation
    assert "peptides like propionylated trypsin" in ultra.evidence
    assert "poor H3 coverage (10–20%)" in chymo.evidence
    assert "ambiguous" in ultra.notes  # 'not before P' for Arg-C Ultra to confirm
    assert get_rule("Trypsin").suggested_missed_cleavages is None
