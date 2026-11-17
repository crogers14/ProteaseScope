# ProteaseScope tutorial

All results in ProteaseScope are **predictions** from cleavage rules, not observed data.

## Enzymes

Pick the rule source in the sidebar (**Enzymes → Rule source**):

- **PSI-MS** is the default.
- **ExPASy** reproduces Protein Cleaver.

Enzymes that only one source defines are offered in both modes and keep their own source label, so every export names the rule that was actually used. Hover an enzyme's chip under **Alternative enzymes** to see its rule, source and notes.

| Enzyme | Cuts | Source | Notes |
| --- | --- | --- | --- |
| Trypsin | after K or R, not before P | PSI-MS; ExPASy version has extra exceptions | the standard bottom-up enzyme |
| Trypsin/P | after K or R, also before P | PSI-MS | |
| Arg-C | after R (PSI-MS: not before P; ExPASy: also before P) | PSI-MS / ExPASy | trypsin on propionylated protein behaves like Arg-C |
| Lys-C | after K (PSI-MS: not before P; ExPASy: also before P) | PSI-MS / ExPASy | needs unmodified lysines |
| Lys-C/P | after K, also before P | PSI-MS | needs unmodified lysines |
| Lys-N | before K | PSI-MS; Protein Cleaver in ExPASy mode | needs unmodified lysines |
| Asp-N | before D (PSI-MS: also before B) | PSI-MS / ExPASy | |
| Asp-N ambic | before D or E | PSI-MS | Asp-N in ammonium bicarbonate |
| Glu-C | after E (PSI-MS: not after EE) | PSI-MS / ExPASy | |
| Chymotrypsin | after F, Y, W or L, not before P | PSI-MS | |
| Chymotrypsin high specificity | after F or Y (not before P), after W (not before M or P) | ExPASy | |
| Chymotrypsin low specificity | adds L, M and H to the high-specificity sites, with exceptions | ExPASy | |
| Proteinase K | after A, E, F, I, L, T, V, W or Y | ExPASy | nonspecific, low confidence; left out of the Strategy picker |
| ProAlanase | after P or A, before any residue | literature (Samodova et al. 2020, *Mol Cell Proteomics*; Promega; citation to verify) | does not cut at K, so it works on propionylated protein. Acidic digestion (around pH 1.5–2); more missed cleavages expected than trypsin, so consider ≤1–2 missed cleavages. Behavior at P–P bonds is unverified |
| Arg-C Ultra | after R, not before P (same sites as PSI-MS Arg-C) | literature (Promega; Turner et al. 2026, *Anal Chem*) | does not cut at K, so no derivatization is needed for Arg-C-like peptides; suggested ≤1 missed cleavage. Whether "not before P" applies is to be confirmed |
| Recombinant chymotrypsin | after F, Y, L or M, not before P; never after W or H | literature (Promega rChymoselect; Turner et al. 2026) | does not cut at K; about 20–25% of peptides fully cleaved in the paper, so suggested ≤3 missed cleavages |

You can also define a **custom enzyme** ("cut after", "cut before", "except before P"). The app shows the regex it compiles to.

### Lysine-cleaving enzymes and propionylation

Propionylation blocks every lysine, so trypsin, Lys-C, Lys-C/P and Lys-N stop cutting at K. On propionylated protein, trypsin cuts only after R. The other enzymes, ProAlanase included, cut the same way with or without propionylation.

The Strategy picker uses this distinction:

- **K-cleaving enzymes** (trypsin, Lys-C, Lys-N) are also offered on unmodified protein. Those strategies carry an "Unmodified lysines only" badge, because acetylated or methylated lysines block the cut.
- **Every other enzyme**, ProAlanase, Arg-C Ultra and recombinant chymotrypsin included, appears only in its propionylated form: alone, as a separate digest, or after trypsin in the same tube.

## Build your own strategy

The Strategy picker (Histones tab) ranks ready-made strategies. To test a workflow of your own, such as a published one, click **Build your own strategy** at the top of the picker. A panel opens on the right; the rest of the page stays as it is.

### What a strategy is

A strategy is one to three **digests**, and how they are **combined**.

Each digest has an enzyme and one chemistry choice:

| Chemistry | Lysines during digestion | Labels after digestion | Typical use |
| --- | --- | --- | --- |
| Propionylate before and after digestion | blocked: no enzyme cuts at K | peptide N-termini (lysines are already propionylated) | the standard histone workflow with trypsin |
| Propionylate after digestion only | free: K-cleaving enzymes still cut | every free K and every peptide N-terminus | Arg-C Ultra + propionylation |
| No derivatization | free | nothing | RIPUP |

With two or three digests, choose how they are combined:

- **Separate digests, results combined:** the sample is split; each digest is run on its own with its own chemistry, and the peptides of all digests are pooled. It costs one sample prep and one LC-MS run per digest.
- **Same tube, one after another:** all enzymes cut the same sample, so it is one prep and one run. The cut sites of all enzymes add up, and every digest shares the first digest's chemistry.

The length window, missed cleavages, m/z filter, fixed modifications and Met removal come from the sidebar, as for every other strategy.

### Guard rails

Short notes under a digest flag common problems:

- **Lysine blocking:** a K-cleaving enzyme (trypsin, Lys-C, Lys-N) with "propionylate before and after" barely cuts at K.
- **No derivatization:** short, very basic peptides are often lost in practice. The prediction does not model this.
- **Same tube:** all digests share digest 1's chemistry.
- **Duplicates:** two identical digests show a red note and turn off **Compare** and **Save** until you change one.

### Worked example: RIPUP (Turner et al. 2026)

RIPUP pairs Arg-C Ultra and recombinant chymotrypsin, each without derivatization, as two separate digests.

1. Choose the **Histone PTM** preset in the sidebar (4–40 aa, ≤1 missed cleavage).
2. Open the **Histones** tab, then **Strategy picker**, and pick **H4** and its **N-terminal tail (1–25)**.
3. Click **Build your own strategy**. Under **Start from…**, choose **RIPUP (Turner et al. 2026)**. The two digests fill in. Hover **source ⓘ** for the citation and **About … ⓘ** for each enzyme's notes.
4. Leave **Compare against** on *Standard digest (current preset)* and click **Compare**.

The picker now shows one card for RIPUP, an **All goals** table, and the usual pictures:

| Goal (H4 1–25) | Standard → RIPUP | Gained | Lost |
| --- | --- | --- | --- |
| See lysines the standard digest misses | 5 → 5 | 0 | 0 |
| Separate lysines that share a peptide | 0 → 4 | 4 (K5–K12, K5–K16, K8–K12, K8–K16) | 0 |
| Read lysines together | 6 → 10 | 4 (K5–K20, K8–K20, K12–K20, K16–K20) | 0 |

**Why:**
- **The standard digest:** every peptide carrying K5 or K8 also carries K12 and K16 (4–17, 1–17 and 4–19), so these pairs can never be measured apart.
- **Separating them:** recombinant chymotrypsin adds 1–10 (K5, K8) and 11–22 (K12, K16, K20), which split K5/K8 from K12/K16.
- **Reading K20 with the others:** its missed-cleavage peptide 1–22 carries all five tail lysines.
- **The cost:** 2 sample preps, 2 enzymes and 2 LC-MS runs.

The peptide map shows the standard digest, then one row per RIPUP digest (Arg-C Ultra, recombinant chymotrypsin), then the combined row.

Switch the histone to **H2A type 1, whole protein**:
- RIPUP sees 14 lysines instead of 7 (+K95, K99, K118, K119, K125, K127, K129), nothing lost.
- It reads 22 lysine pairs together instead of 7, mostly across the C-terminal tail.

These are predictions from cleavage rules. Turner et al. report that recombinant chymotrypsin gave poor H3 coverage in practice (10–20%) despite higher predicted coverage.

### "Separate lysines that share a peptide": what counts

A pair counts only if its two lysines share at least one identifiable peptide in the digest you compare against ("tangled"), and the strategy then gives each lysine a peptide without the other. Pairs that were never on one peptide are not counted, even if they end up on different peptides. On the H4 tail with the regression settings (7–30 aa, 0 missed cleavages), K5–K20 is not counted, because the standard digest never sees K20 at all.

### Saving, comparing and sharing

- **Save:** keeps the strategy for this browser session. It is then ranked in the picker's cards, rankings and advanced table like any other strategy, and you can choose it under **Compare against**, for example to compare a same-tube version with the separate one.
- **Download / Load:** save a strategy as a small JSON file and load it again later, or share it:

  ```json
  {"name": "RIPUP (Turner et al. 2026)",
   "digests": [{"enzyme": "Arg-C Ultra", "source": "literature", "chemistry": "none"},
               {"enzyme": "Recombinant chymotrypsin", "source": "literature", "chemistry": "none"}],
   "combine": "separate"}
  ```

- **Exports:** every export header, the downloaded answer and the Methods text list each custom strategy's digests, chemistry and combination.
- **Back to the recommendations:** this link in the blue banner returns to the ranked cards.
