# RIPUP (Turner et al. 2026) vs Standard histone digest: all goals (PREDICTED, in silico)

- PREDICTED (in silico) results. These are predictions from cleavage rules, not observed data.
- Generated: 2026-11-16T10:39:02-05:00
- Histone scopes: H2A type 1 (P0C0S8), residues 1-20 mature numbering (N-terminal tail); H2A type 1 (P0C0S8), residues 1-129 mature numbering (Whole protein); H2B type 1-C/E/F/G/I (P62807), residues 1-35 mature numbering (N-terminal tail); H2B type 1-C/E/F/G/I (P62807), residues 1-125 mature numbering (Whole protein); H3.1 (P68431), residues 1-40 mature numbering (N-terminal tail); H3.1 (P68431), residues 1-135 mature numbering (Whole protein); H4 (P62805), residues 1-25 mature numbering (N-terminal tail); H4 (P62805), residues 1-102 mature numbering (Whole protein)
- Digest:
-   Digest: Trypsin (PSI-MS) [single]
-     Trypsin (PSI-MS): regex (?<=[KR])(?!P)  [PSI-MS MS:1001251, psi-ms.obo v4.2.17]
-       Citation: HUPO-PSI PSI-MS controlled vocabulary (psi-ms.obo v4.2.17); Mayer G et al. (2013) Database 2013:bat009.
-     Peptide length: 4-40 aa
-     Max missed cleavages: 1
-     m/z filter: off
-     Fixed modifications: Carbamidomethyl (C)
-     Propionylation: Propionylate before and after digestion
-     Initiator Met removal: MetAP rule
- Numbering: UniProt (initiator Met = 1); *_mature = UniProt - 1 when Met is removed.
- Masses: neutral monoisotopic and average (ExPASy average residue masses); m/z = (M + z*1.007276)/z.
- Counts: peptides are counted by position (entry, start, end).
- Compared with: propionylation + trypsin · 4–40 aa · ≤1 missed cleavage
- Measurable alone: peptides with at most 1 lysine(s) (all K on the peptide count). Known-mark filter: any.
- Strategies: alone = enzyme instead of the standard digest; parallel = standard + separate enzyme digest (union of identifiable peptides); sequential = standard then enzyme in the same tube (union of cut sites, always propionylated). Unmodified = separate prep without propionylation (lysine-cleaving enzymes).
- Separable pairs: only pairs that share an identifiable peptide in the compared-against digest (Standard histone digest) are counted.
- Net = Strategy - Standard = Gained - Lost.
- Custom strategy: Standard histone digest: Trypsin, propionylate before and after digestion
- Custom strategy: RIPUP (Turner et al. 2026): Arg-C Ultra, no derivatization; Recombinant chymotrypsin, no derivatization (separate digests, results combined)
- Compared against: Standard histone digest

## H2A type 1 (P0C0S8), N-terminal tail (1–20)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 4 → 4 | 0 | 0 | +0 | — | — |
| Measure lysines one at a time | 0 → 0 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 4 → 4 | 0 | 0 | +0 | — | — |
| Read lysines together (combinations / crosstalk) | 6 → 6 | 0 | 0 | +0 | — | — |
| Cover known modification sites | 4 → 4 | 0 | 0 | +0 | — | — |
| Tell variants and isoforms apart | 13 → 13 | 0 | 0 | +0 | — | — |

## H2A type 1 (P0C0S8), Whole protein (1–129)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 7 → 14 | 7 | 0 | +7 | K95, K99, K118, K119, K125, K127, K129 | — |
| Measure lysines one at a time | 1 → 3 | 2 | 0 | +2 | K95, K99 | — |
| Separate lysines that share a peptide | 4 → 4 | 0 | 0 | +0 | — | — |
| Read lysines together (combinations / crosstalk) | 7 → 22 | 15 | 0 | +15 | K99–K118, K99–K119, K99–K125, K99–K127, K99–K129, K118–K119, K118–K125, K118–K127, K118–K129, K119–K125, K119–K127, K119–K129, K125–K127, K125–K129, K127–K129 | — |
| Cover known modification sites | 7 → 12 | 5 | 0 | +5 | K95, K99, K118, K119, K125 | — |
| Tell variants and isoforms apart | 15 → 17 | 2 | 0 | +2 | H2A type 1-H, H2A type 1-J | — |

## H2B type 1-C/E/F/G/I (P62807), N-terminal tail (1–35)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 12 → 12 | 0 | 0 | +0 | — | — |
| Measure lysines one at a time | 1 → 1 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 10 → 10 | 0 | 0 | +0 | — | — |
| Read lysines together (combinations / crosstalk) | 55 → 66 | 11 | 0 | +11 | K5–K34, K11–K34, K12–K34, K15–K34, K16–K34, K20–K34, K23–K34, K24–K34, K27–K34, K28–K34, K30–K34 | — |
| Cover known modification sites | 9 → 9 | 0 | 0 | +0 | — | — |
| Tell variants and isoforms apart | 16 → 16 | 0 | 0 | +0 | — | — |

## H2B type 1-C/E/F/G/I (P62807), Whole protein (1–125)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 20 → 20 | 0 | 0 | +0 | — | — |
| Measure lysines one at a time | 2 → 4 | 2 | 0 | +2 | K43, K125 | — |
| Separate lysines that share a peptide | 10 → 18 | 8 | 0 | +8 | K34–K43, K34–K46, K34–K57, K43–K46, K43–K57, K108–K125, K116–K125, K120–K125 | — |
| Read lysines together (combinations / crosstalk) | 67 → 78 | 11 | 0 | +11 | K5–K34, K11–K34, K12–K34, K15–K34, K16–K34, K20–K34, K23–K34, K24–K34, K27–K34, K28–K34, K30–K34 | — |
| Cover known modification sites | 16 → 16 | 0 | 0 | +0 | — | — |
| Tell variants and isoforms apart | 17 → 17 | 0 | 0 | +0 | — | — |

## H3.1 (P68431), N-terminal tail (1–40)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 8 → 8 | 0 | 0 | +0 | — | — |
| Measure lysines one at a time | 1 → 1 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 12 → 13 | 1 | 0 | +1 | K18–K23 | — |
| Read lysines together (combinations / crosstalk) | 17 → 18 | 1 | 0 | +1 | K4–K18 | — |
| Cover known modification sites | 8 → 8 | 0 | 0 | +0 | — | — |
| Tell variants and isoforms apart | 6 → 6 | 0 | 0 | +0 | — | — |

## H3.1 (P68431), Whole protein (1–135)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 13 → 13 | 0 | 0 | +0 | — | — |
| Measure lysines one at a time | 6 → 6 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 13 → 14 | 1 | 0 | +1 | K18–K23 | — |
| Read lysines together (combinations / crosstalk) | 18 → 20 | 2 | 0 | +2 | K4–K18, K115–K122 | — |
| Cover known modification sites | 13 → 13 | 0 | 0 | +0 | — | — |
| Tell variants and isoforms apart | 8 → 8 | 0 | 0 | +0 | — | — |

## H4 (P62805), N-terminal tail (1–25)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 5 → 5 | 0 | 0 | +0 | — | — |
| Measure lysines one at a time | 1 → 1 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 0 → 4 | 4 | 0 | +4 | K5–K12, K5–K16, K8–K12, K8–K16 | — |
| Read lysines together (combinations / crosstalk) | 6 → 10 | 4 | 0 | +4 | K5–K20, K8–K20, K12–K20, K16–K20 | — |
| Cover known modification sites | 5 → 5 | 0 | 0 | +0 | — | — |
| Tell variants and isoforms apart | 0 → 0 | 0 | 0 | +0 | — | — |

## H4 (P62805), Whole protein (1–102)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 11 → 11 | 0 | 0 | +0 | — | — |
| Measure lysines one at a time | 5 → 6 | 1 | 0 | +1 | K91 | — |
| Separate lysines that share a peptide | 4 → 9 | 5 | 0 | +5 | K5–K12, K5–K16, K8–K12, K8–K16, K79–K91 | — |
| Read lysines together (combinations / crosstalk) | 11 → 18 | 7 | 0 | +7 | K5–K20, K8–K20, K12–K20, K12–K31, K16–K20, K16–K31, K31–K44 | — |
| Cover known modification sites | 11 → 11 | 0 | 0 | +0 | — | — |
| Tell variants and isoforms apart | 0 → 0 | 0 | 0 | +0 | — | — |
