# RIPUP (Turner et al. 2026) vs Standard histone digest: all goals (PREDICTED, in silico)

- PREDICTED (in silico) results. These are predictions from cleavage rules, not observed data.
- Generated: 2026-11-16T10:39:02-05:00
- Histone scopes: H2A type 1 (P0C0S8), residues 1-20 mature numbering (N-terminal tail); H2A type 1 (P0C0S8), residues 1-129 mature numbering (Whole protein); H2B type 1-C/E/F/G/I (P62807), residues 1-35 mature numbering (N-terminal tail); H2B type 1-C/E/F/G/I (P62807), residues 1-125 mature numbering (Whole protein); H3.1 (P68431), residues 1-40 mature numbering (N-terminal tail); H3.1 (P68431), residues 1-135 mature numbering (Whole protein); H4 (P62805), residues 1-25 mature numbering (N-terminal tail); H4 (P62805), residues 1-102 mature numbering (Whole protein)
- Digest:
-   Digest: Trypsin (PSI-MS) [single]
-     Trypsin (PSI-MS): regex (?<=[KR])(?!P)  [PSI-MS MS:1001251, psi-ms.obo v4.2.17]
-       Citation: HUPO-PSI PSI-MS controlled vocabulary (psi-ms.obo v4.2.17); Mayer G et al. (2013) Database 2013:bat009.
-     Peptide length: 7-30 aa
-     Max missed cleavages: 0
-     m/z filter: off
-     Fixed modifications: Carbamidomethyl (C)
-     Propionylation: Propionylate before and after digestion
-     Initiator Met removal: MetAP rule
- Numbering: UniProt (initiator Met = 1); *_mature = UniProt - 1 when Met is removed.
- Masses: neutral monoisotopic and average (ExPASy average residue masses); m/z = (M + z*1.007276)/z.
- Counts: peptides are counted by position (entry, start, end).
- Compared with: propionylation + trypsin · 7–30 aa · no missed cleavages
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
| See lysines the standard digest misses | 2 → 4 | 2 | 0 | +2 | K13, K15 | — |
| Measure lysines one at a time | 0 → 0 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 0 → 0 | 0 | 0 | +0 | — | — |
| Read lysines together (combinations / crosstalk) | 1 → 6 | 5 | 0 | +5 | K5–K13, K5–K15, K9–K13, K9–K15, K13–K15 | — |
| Cover known modification sites | 2 → 4 | 2 | 0 | +2 | K13, K15 | — |
| Tell variants and isoforms apart | 7 → 13 | 6 | 0 | +6 | H2A type 1-C, H2A type 1-A (TH2A), H2A type 2-A, H2A type 2-B, H2A type 2-C, H2A type 3 | — |

## H2A type 1 (P0C0S8), Whole protein (1–129)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 3 → 13 | 10 | 0 | +10 | K13, K15, K74, K75, K99, K118, K119, K125, K127, K129 | — |
| Measure lysines one at a time | 1 → 2 | 1 | 0 | +1 | K99 | — |
| Separate lysines that share a peptide | 0 → 0 | 0 | 0 | +0 | — | — |
| Read lysines together (combinations / crosstalk) | 1 → 17 | 16 | 0 | +16 | K5–K13, K5–K15, K9–K13, K9–K15, K13–K15, K74–K75, K118–K119, K118–K125, K118–K127, K118–K129, K119–K125, K119–K127, K119–K129, K125–K127, K125–K129, K127–K129 | — |
| Cover known modification sites | 3 → 11 | 8 | 0 | +8 | K13, K15, K74, K75, K99, K118, K119, K125 | — |
| Tell variants and isoforms apart | 14 → 17 | 3 | 0 | +3 | H2A type 1-C, H2A type 1-H, H2A type 1-J | — |

## H2B type 1-C/E/F/G/I (P62807), N-terminal tail (1–35)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 10 → 10 | 0 | 0 | +0 | — | — |
| Measure lysines one at a time | 0 → 0 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 0 → 0 | 0 | 0 | +0 | — | — |
| Read lysines together (combinations / crosstalk) | 45 → 45 | 0 | 0 | +0 | — | — |
| Cover known modification sites | 8 → 8 | 0 | 0 | +0 | — | — |
| Tell variants and isoforms apart | 14 → 14 | 0 | 0 | +0 | — | — |

## H2B type 1-C/E/F/G/I (P62807), Whole protein (1–125)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 15 → 17 | 2 | 0 | +2 | K46, K57 | — |
| Measure lysines one at a time | 1 → 1 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 0 → 0 | 0 | 0 | +0 | — | — |
| Read lysines together (combinations / crosstalk) | 51 → 52 | 1 | 0 | +1 | K46–K57 | — |
| Cover known modification sites | 12 → 14 | 2 | 0 | +2 | K46, K57 | — |
| Tell variants and isoforms apart | 16 → 16 | 0 | 0 | +0 | — | — |

## H3.1 (P68431), N-terminal tail (1–40)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 7 → 8 | 1 | 0 | +1 | K4 | — |
| Measure lysines one at a time | 0 → 0 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 0 → 1 | 1 | 0 | +1 | K18–K23 | — |
| Read lysines together (combinations / crosstalk) | 5 → 13 | 8 | 0 | +8 | K4–K9, K4–K14, K4–K18, K9–K18, K14–K18, K23–K27, K23–K36, K23–K37 | — |
| Cover known modification sites | 7 → 8 | 1 | 0 | +1 | K4 | — |
| Tell variants and isoforms apart | 6 → 6 | 0 | 0 | +0 | — | — |

## H3.1 (P68431), Whole protein (1–135)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 10 → 12 | 2 | 0 | +2 | K4, K115 | — |
| Measure lysines one at a time | 3 → 3 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 0 → 1 | 1 | 0 | +1 | K18–K23 | — |
| Read lysines together (combinations / crosstalk) | 5 → 14 | 9 | 0 | +9 | K4–K9, K4–K14, K4–K18, K9–K18, K14–K18, K23–K27, K23–K36, K23–K37, K115–K122 | — |
| Cover known modification sites | 10 → 12 | 2 | 0 | +2 | K4, K115 | — |
| Tell variants and isoforms apart | 7 → 8 | 1 | 0 | +1 | H3.2 | — |

## H4 (P62805), N-terminal tail (1–25)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 4 → 5 | 1 | 0 | +1 | K20 | — |
| Measure lysines one at a time | 0 → 0 | 0 | 0 | +0 | — | — |
| Separate lysines that share a peptide | 0 → 4 | 4 | 0 | +4 | K5–K12, K5–K16, K8–K12, K8–K16 | — |
| Read lysines together (combinations / crosstalk) | 6 → 8 | 2 | 0 | +2 | K12–K20, K16–K20 | — |
| Cover known modification sites | 4 → 5 | 1 | 0 | +1 | K20 | — |
| Tell variants and isoforms apart | 0 → 0 | 0 | 0 | +0 | — | — |

## H4 (P62805), Whole protein (1–102)

| Goal | Standard → RIPUP | Gained | Lost | Net | Gained items | Lost items |
| --- | --- | --- | --- | --- | --- | --- |
| See lysines the standard digest misses | 9 → 11 | 2 | 0 | +2 | K20, K44 | — |
| Measure lysines one at a time | 3 → 5 | 2 | 0 | +2 | K44, K91 | — |
| Separate lysines that share a peptide | 0 → 5 | 5 | 0 | +5 | K5–K12, K5–K16, K8–K12, K8–K16, K79–K91 | — |
| Read lysines together (combinations / crosstalk) | 7 → 10 | 3 | 0 | +3 | K12–K20, K16–K20, K77–K79 | — |
| Cover known modification sites | 9 → 11 | 2 | 0 | +2 | K20, K44 | — |
| Tell variants and isoforms apart | 0 → 0 | 0 | 0 | +0 | — | — |
