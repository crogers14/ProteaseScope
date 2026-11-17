# Changelog

## Repository and GitHub Pages deployment (2026-11-17, "v0.1 beta")

- **Git repository:** the project is now a git repository, pushed to `crogers14/ProteaseScope` on GitHub. The first commit is labeled "v0.1 beta"; the package version stays 0.7.0.dev2.
- **`.gitignore`:** excludes caches, build outputs (`build/`, `site/`), local assistant settings (`.claude/`), Node modules, and large raw downloads. The full UniProt FASTA (`data/reference/*.fasta`, about 29 MB) is re-created by `scripts/fetch_data.py`; its provenance stays in `data/reference/metadata.json`. `data/precomputed/` is committed (largest file 20.4 MB), because rebuilding it takes about 18 minutes and 4.7 GB of memory.
- **GitHub Actions (`.github/workflows/pages.yml`):** on every push to main, runs ruff and pytest, then builds the Shinylive site and deploys it to GitHub Pages. Build and deploy run only if the tests pass.
- **Crawler blocking:** `scripts/build_site.py` adds `<meta name="robots" content="noindex">` to every HTML page and writes a `robots.txt` that disallows everything. On a GitHub Pages project site, crawlers only read `robots.txt` at the domain root, so the meta tag is what keeps the pages out of search results.
- **Test change:** `test_downloaded_files_have_provenance` no longer requires a raw FASTA download to be present. It still requires the file's provenance record, and still checks the file when it exists. Reason: the FASTA is deliberately not committed, so it is absent in CI.

## 0.7.0.dev2 (2026-11-16): Phase 7 step 2, "Build your own strategy" panel

The builder UI, following `CUSTOM_STRATEGY.md` §3–6. It adds no new tab and no new results layout.

**Entry point:** one **Build your own strategy** button at the top of the Strategy picker (Histones tab) opens a right-hand side panel; the rest of the page is unchanged. The panel can be collapsed again with its chevron.

**Layout fix (2026-11-16, after review):** the panel's container had made the picker narrower and indented even with the panel closed: 24 px of padding at the top, 48 px on each side, plus a collapse arrow over the content. The picker now has the same width and alignment as in 0.6.0. The arrow is hidden while the panel is closed, since the button opens it, and appears as the panel's close control once it is open. Nothing else in the picker changed.

**Panel, top to bottom:**
- **Name:** default "My strategy 1".
- **Start from…:** Standard histone digest, RIPUP (Turner et al. 2026), or Arg-C Ultra + propionylation (Turner et al. 2026). Choosing one fills the rows, which stay editable; a "source ⓘ" chip shows the citation.
- **Digest rows (1–3):** each has an enzyme dropdown (Arg-C Ultra and recombinant chymotrypsin first), a chemistry dropdown with the three plain-language choices, and a × button. "+ Add digest" disappears at 3 rows. New rows start with the current preset's chemistry and the first enzyme not already in the strategy.
- **Enzyme tooltip:** an "About <enzyme> ⓘ" chip under each row shows the enzyme's tooltip: bench note, the Turner et al. evidence line, suggested missed cleavages, the regex and the notes. The sidebar's enzyme chips show the same tooltip.
- **How they're combined:** shown only with 2–3 digests: "Separate digests, results combined" (default) or "Same tube, one after another", each with its tooltip.
- **Compare against:** "Standard digest (current preset)" or any saved strategy.
- **Buttons:** Compare, Save, Download (JSON) and Load… (JSON).

**Guard rails** appear as inline notes under each row:
- **Same tube:** the other rows' chemistry dropdowns are replaced by "Same chemistry as digest 1", with the note "In one tube, all digests share the same chemistry." The model uses the first row's chemistry for every row.
- **K-cleaving enzyme with "before and after":** "Lysines are blocked, so this enzyme will barely cut at K."
- **"No derivatization":** "Short, very basic peptides are often lost without derivatization." No realism flag exists yet to link to.
- **Two identical rows:** "This digest is already in the strategy", shown in red. This note disables Compare and Save.

**Results** (Compare) reuse the picker's results view for the current histone, region and goal:
- **Banner:** the strategy and what it is compared against, with "Back to the recommendations".
- **Answer card:** one card for the strategy: the headline gain for the goal ("+0 nothing new for this goal" when it gains nothing), "What you give up", and **Work** (e.g. "2 sample preps, 2 enzymes, 2 LC-MS runs"; same tube: "1 sample prep, N enzymes, 1 LC-MS run").
- **All goals:** a table with one row per goal: Goal, "<compared-against> → <strategy>", Gained, Lost, Net. Hover a number to list its lysines, pairs or variants.
- **Lysine picture, arc view and variant list:** as in the picker. Labels name the compared-against strategy when it is not the standard digest.
- **Peptide map:** the compared-against row, one row per digest, then the combined row.
- **Text download:** the card plus the all-goals table.

**Saving:**
- **Save** adds the strategy under its name for this browser session. It is then scored in the picker's cards, rankings and advanced table like any other strategy, and appears under "Compare against". The name then advances to the next free "My strategy N".
- **Download / Load:** a strategy as a small JSON file in the documented format. Loading fills the rows; a strategy with an enzyme not available under the current rules gets an error message.
- **Exports:** the picker's answer text and table headers, and the Overview's Methods text, list each saved strategy's digests, chemistry and combination.

**Not done:**
- **Builder presets and missed cleavages:** the builder presets do not change missed cleavages. The spec's preset table defines none, so the sidebar setting applies.
- **Tutorial:** step 3 is done (below).

**Step 3, docs:** TUTORIAL.md has a new "Build your own strategy" section:
- **Explained:** what a strategy is (digests, the three chemistry choices, separate vs same tube), the guard rails, and what "Separate lysines that share a peptide" counts.
- **Worked example:** RIPUP under the Histone PTM preset, with H4 1–25 (0 → 4 pairs newly separable, 6 → 10 pairs read together, and the peptides that cause them) and H2A type 1 (7 → 14 lysines).
- **Sharing:** Save, Compare against, Download / Load, with a JSON example that loads as the RIPUP preset.
- **Caveat:** Turner et al. observed poor H3 coverage with recombinant chymotrypsin.

**Checked:**
- **Native app:** RIPUP on H4 1–25 and on H2A type 1 (whole protein) under Histone PTM settings. The checks covered the guard rails (duplicate rows, same-tube sync), Save, and the saved strategy appearing in the rankings and under "Compare against".
- **Browser (Shinylive) build:** the same RIPUP comparison on H2A type 1 gives 7 → 14 lysines and 7 → 22 pairs read together, matching `results/phase7/ripup_vs_standard_histone_ptm.md`.

**Tests:** 197 pass. New tests check:
- the builder controls are in the UI
- the "My strategy N" naming
- the split peptide-map rows: the combined row is the union of the digest rows
- strict-JSON serialization of the comparison figures

The precomputed data are unchanged from 0.7.0.dev1, and the manifest still says dev1.

## 0.7.0.dev1 (2026-11-16): Phase 7 step 1, custom strategies core (`CUSTOM_STRATEGY.md`)

Core only; the builder panel (step 2) and the tutorial section (step 3) are not built yet.

**Two new enzymes**, literature-defined (`source="Literature"`) and offered in both rule-source modes:
- **Arg-C Ultra:** `(?<=R)(?!P)`, after R, not before P. Same cut sites as the PSI-MS Arg-C rule; kept under the product name. Suggested ≤1 missed cleavage.
- **Recombinant chymotrypsin** (Promega rChymoselect): `(?<=[FYLM])(?!P)`, after F, Y, L or M, not before P, never after W or H. Suggested ≤3 missed cleavages.
- **New `CleavageRule` fields:** `suggested_missed_cleavages` and `evidence`. The evidence is a one-line note labeled "Observed in Turner et al. 2026 (HEK293T histones)".
- **Citation:** Turner et al. 2026, *Anal Chem* 98:18694, doi:10.1021/acs.analchem.6c01147, as provided by the user, plus the Promega product information. The Ryzhaya et al. citation is incomplete and marked "to verify".
- **To confirm (user):** the paper's text is ambiguous about whether "not before P" also applies to Arg-C Ultra. The rule assumes it does, as `CUSTOM_STRATEGY.md` §1 says. If Arg-C Ultra does cut before P, drop `(?!P)`; the Arg-C equivalence test will then need the ExPASy Arg-C rule instead.
- **Every preset gains six alternatives** (each new enzyme alone, in parallel and sequential): 16 enzymes for proteome, histone PTM and middle-down, 12 for Protein Cleaver compatible. The precomputed data were rebuilt.
- **Strategy picker:** both enzymes are propionylated-only candidates (neither cleaves at K). On H4 K16 + K20 (ExPASy, 7–30 aa, "Read lysines together"), recombinant chymotrypsin 11–22 (GKGGAKRHRKVL) reads both lysines and ties Asp-N alone; the full-set test now expects Asp-N, ProAlanase, recombinant chymotrypsin.

**Per-digest chemistry** (`CUSTOM_STRATEGY.md` §2). The old propionylation setting is replaced by one plain-language choice per digest:
- **The choices** (`Propionylation` values):
  - `before_and_after`, "Propionylate before and after digestion": K blocked, K and peptide N-termini labeled
  - `after_only`, "Propionylate after digestion only": K not blocked, every free K and peptide N-terminus labeled
  - `none`, "No derivatization"
- **Migration:** the old names still load. `protein + peptide N-term` → `before_and_after` (same cleavage and masses). `none` → `none`.
- **Behavior change:** the old `protein-level` option (K blocked, no N-terminal label) also maps to `before_and_after`, as the spec says. Masses for that option now include +56.026215 Da on each peptide N-terminus. No preset used it.
- **UI:** the sidebar select is now labeled "Derivatization" and shows the three plain-language choices. Export headers and the Methods text use the same wording.

**Strategy model** (`CUSTOM_STRATEGY.md` §7; new module `workflow.py`):
- **Structure:** `Workflow(name, digests, combine)`, where each `DigestStep` holds an enzyme and its chemistry, and `combine` is `separate` or `same_tube`.
  - Separate: each digest runs with its own chemistry, then the union of identifiable peptides is taken.
  - Same tube: one chemistry (the first row's), the union of cut sites, one digest.
  - JSON save/load uses the documented format.
- **One code path:** the standard digest, the picker's `alone`, `parallel` and `sequential`, and the proteome presets' `DigestRecipe` are all expressed with it.
  - Tests check they produce the same settings as before.
  - "Separate [standard, X]" equals parallel + X, and "same tube [trypsin, X], before and after" equals sequential + X, on goals and peptides.
- **Guard rails:** same-tube chemistries are synchronized, with "In one tube, all digests share the same chemistry." Each row also gets plain notes:
  - K-cleaving enzymes with "before and after": "Lysines are blocked, so this enzyme will barely cut at K."
  - "No derivatization" rows: "Short, very basic peptides are often lost without derivatization."
  - identical rows: "This digest is already in the strategy". Only this note blocks Compare.
- **Builder presets ("Start from…"):** Standard histone digest; RIPUP (Turner et al. 2026); Arg-C Ultra + propionylation (Turner et al. 2026). Each carries its citation.
- **Picker integration:**
  - custom strategies are scored and ranked like any other strategy (`pick_strategies(custom=…)`)
  - any strategy can be compared against a saved one (`against=…`, `StrategyResults.reference`)
  - `compare_workflows` is the shortcut for one custom strategy
  - `card_for` gives the one answer card, with absolute work, e.g. "2 sample preps, 2 enzymes, 2 LC-MS runs"
  - `all_goals_summary` gives the "All goals" table (Compared-against → Custom, Gained, Lost, Net)
  - export headers list each custom strategy's digests, chemistry and combination

**Separable-pairs correction** (`CUSTOM_STRATEGY.md` §8, also written into `STRATEGY_PICKER_2.md` §1, §4, §5 and §6):
- **Why:** the old definition counted any two lysines that end up on different peptides, including far-apart pairs that were never a problem.
- **New definition:** a pair counts only if its two lysines share an identifiable peptide in the compared-against digest ("tangled") and, with the strategy, each lysine also sits on a peptide without the other. The goal's total is now the number of tangled pairs, and its sentences say "pairs that share a peptide".
- **Arc view:** the newly-separable arcs use the same definition.
- **Tests:** the H4 parallel + low-specificity chymotrypsin case still gives exactly K5–K12, K5–K16, K8–K12 and K8–K16. In RIPUP vs standard on H4 1–25, K5–K20 and K8–K20 are separable but not counted.

**RIPUP regression** (exact, PSI-MS, 7–30 aa, 0 missed cleavages, Met removed):
- H4 1–25: 4 → 5 lysines (+K20); new pairs K12–K20 and K16–K20; newly separable K5–K12, K5–K16, K8–K12 and K8–K16.
- H2A type 1 whole protein: 3 → 13 lysines (+K13, K15, K74, K75, K99, K118, K119, K125, K127, K129), nothing lost.

Summaries for H2A type 1, H2B type 1-C, H3.1 and H4 (N-terminal tail and whole protein) are in `results/phase7/` (`scripts/export_custom_summary.py`). Two settings are exported: the regression settings, and the Histone PTM preset's 4–40 aa with ≤1 missed cleavage.

## 0.6.0.dev4 (2026-11-15): Overview leads with gained vs lost

Requested by the user: the four baseline value boxes (identifiable proteins, peptides, unique peptides, median coverage) never changed with the alternatives, so they are gone. The gained/lost chart now opens the Overview tab with exact numbers and a median-coverage panel.

**"Gained and lost versus the baseline: <baseline>"** (top card):
- **Layout:** a 2 × 2 grid of panels:
  - identifiable peptides
  - K sites
  - S/T/Y sites
  - median coverage

  Each panel title gives the baseline total, e.g. "(baseline 2,844,568)".
- **Bars:** red to the left = lost, blue to the right = gained, with the exact count printed on every bar (+1,872,078 / −929,366). A digest that changes nothing says "no change".
- **Axis sizing:** each panel's axis leaves room for its labels, so numbers never overlap bars or row names. Below 820 px the figure scrolls sideways instead of squeezing.
- **Hover:** full digest name, gained, lost, net, net as a % of the baseline, and baseline total → new total.
- **Median coverage:** the change in median per-protein coverage (alternative − baseline, percentage points; the same quantity as `median_coverage_pct_diff` in the summary table). Bars are blue up, red down, labeled e.g. "+10.2 pts"; hover shows "79.9% → 90.1%".
- **Row labels:** short and in words ("Asp-N · parallel", "Asp-N · sequential", "Asp-N · alone") instead of "Trypsin (PSI-MS) | Asp-N (PSI-MS)". (No note under the chart, at the user's request.)
- **Reconciliation:** checked on all 39 Histone PTM digests: baseline + gained − lost equals the alternative's total for peptides, K sites and each S/T/Y residue. The net and new-total hover values rely on this.
- **Summary table:** still follows, unchanged.

**Tests:** 162 pass. New: `test_gained_lost_figure_labels_every_bar` checks exact labels, "no change", the coverage-point labels, hover net/% and totals, the baseline totals in the titles, and strict JSON.

## 0.6.0.dev3 (2026-11-15): ProAlanase

**New enzyme: ProAlanase**, requested by the user.
- **Rule:** cuts after P or A, with no restriction on the next residue: `(?<=[PA])`. The rule's notes say its behavior at P–P bonds is unverified; it cuts between every P–P pair.
- **Source:** `psi-ms.obo` v4.2.17 (2026-11-11) has no ProAlanase cleavage-agent term, and neither ExPASy nor `pyteomics` has one. It is therefore a literature-defined rule: `source="Literature"`, `mode="literature"`, offered in both rule-source modes like the ExPASy-only enzymes.
- **Citation**, marked "[to verify]": Samodova et al. 2020, *Mol Cell Proteomics*, "ProAlanase is an effective alternative to trypsin for proteomics applications and disulfide bond mapping", plus the Promega product information.
- **Chemistry:** it does not cut at K, so it works on propionylated protein. The Strategy picker includes it as a candidate in `alone`, `parallel` and `sequential`, propionylated only; `cleaves_at_lysine` is False, so no unmodified variant is generated.
- **UI:** `CleavageRule` has a new `ui_note` field. Under "Alternative enzymes" the sidebar shows one chip per selected enzyme; hovering shows its rule, source and notes. ProAlanase's chip is tinted blue and its tooltip leads with "Acidic digestion (around pH 1.5–2); more missed cleavages expected than trypsin, so consider ≤1–2 missed cleavages."
- **Every preset gains three alternatives** (ProAlanase alone, in parallel and sequential): the proteome, histone PTM and middle-down presets go from 13 to 14 enzymes, and Protein Cleaver compatible from 9 to 10. The precomputed data were rebuilt.

**Effect on the Strategy picker:** ProAlanase can now take answer cards. For the STRATEGY_PICKER_2.md specific-lysines case (H4 K16 + K20, "Read lysines together", ExPASy 7–30 aa), its peptide 16–32 (KRHRKVLRDNIQGITKP) reads both lysines. It ties "Use Asp-N instead" on gain, losses and runs, ranks second by name, and moves "Add a separate Asp-N digest" to third. The spec test still runs on the v2 candidate set, and a new assertion covers the full set.

**Tests:**
- cuts after every P and A on 10,000 random sequences and on hand-made cases, including before P, K, R and D
- metadata: the citation is marked "to verify", the notes flag P–P, and the UI note text
- `psi-ms.obo` has no ProAlanase or alanine cleavage-agent term (this test fails if a later CV release adds one)
- picker candidates: only propionylated `alone`, `parallel` and `sequential`
- Histone PTM preset, ProAlanase alone:
  - H3.1 gives 1–7 ARTKQTA (K4), 8–15 RKSTGGKA (K9, K14), 22–25 TKAA (K23) and 31–38 ATGGVKKP (K36, K37)
  - H2B type 1-C gives 11–17 KKGSKKA (K11, K12, K15, K16) and 18–21 VTKA (K20)

**Precomputed data rebuilt** (`scripts/precompute.py`: 17.9 min natively, peak memory 4.5 GB):

| Preset | Digests | Build time | `.npz` size |
| --- | --- | --- | --- |
| Proteome | 40 | 355 s | 18.0 MB |
| Histone PTM | 40 | 357 s | 16.0 MB |
| Middle-down | 40 | 193 s | 18.5 MB |
| Protein Cleaver compatible | 28 | 124 s | 13.9 MB |

- Every file is still under 25 MB; `data/precomputed/` totals 85 MB and the site 156 MB.
- The run started before the version bump, so `manifest.json` says "ProteaseScope 0.6.0.dev2". The code that ran is this release's.
- `tests/test_results.py`: precomputed-vs-live digest count updated from 37 to 40.

**Verified:**
- tests: 161 pass
- natively and in the Shinylive build, ProAlanase loads from precomputed data on the proteome preset as parallel and sequential rows, and its tooltip shows the UI note
- `results/phase6/` regenerated. ProAlanase takes one top-3 card in the three examples: H4 1–25 "Read lysines together" (#3, +4 combinations, tied with Asp-N)

**Docs:**
- CLAUDE.md §3: rules-table row, plus an "Added after v1" line
- new `TUTORIAL.md` with the enzyme list (the "Choosing a strategy" section is still Phase 6 step 3)

## 0.6.0.dev2 (2026-11-14): Phase 6 step 2, Strategy picker UI

**Histones tab → Strategy picker** (the first sub-tab; the other histone sub-tabs keep their own header):
- **Three choices:**
  1. Histone, grouped by family with canonical entries first.
  2. Where to look: N-/C-terminal tail presets for the entry, the whole protein, a custom range, or specific lysines.
  3. What are you looking for: the six goals. An N = 1 / 2 toggle appears only for "Measure lysines one at a time"; a modification-type filter appears only for "Cover known modification sites".
- **"Compared with:" line.** It shows the standard digest (propionylation + trypsin with the sidebar's rule source, lengths, missed cleavages, m/z filter, fixed mods and Met removal). When the lengths or missed cleavages differ from the Histone PTM preset, a "Use the Histone PTM preset settings" link appears.
- **Answer cards (top 3)** with title, big headline number and sentence, the goal-loss line when there is one, "What you give up", "Extra work" and warning badges with tooltips. Hover shows full lists, mark lists and discriminating peptides. Clicking a card selects it; card 1 is selected by default. When nothing improves the goal, the template message is shown instead.
- **Visuals:**
  - lysine picture: standard digest + top 3, with known-mark dots
  - arc view: standard digest and the selected card; arcs above = read together, dashed arcs below = newly separable
  - variant list, with hover peptides
  - region peptide map: standard digest + the selected card
  Each visual has a one-line "How to read this" note.
- **"Compare across all histones"** (button): one row per canonical histone region for the current goal. Clicking a row loads that histone and region.
- **"Show all strategies (advanced)"** (collapsed): the full table, with TSV and Excel downloads that carry the settings header.
- **"Download this answer (text)".**
- A settings strip above the figures says "Predicted (in silico)".

**Core additions (`strategy.py`):**
- `lysine_states`, `pair_states`, `variant_states`, `map_window`, `map_data` (the existing peptide-map format) and `answer_text`.
- Messages now read correctly for a single item ("sees the one lysine").

**Figures:**
- `lysine_picture` and `arc_view` are new.
- The peptide-map legend can say "vs standard digest".
- Colours follow the app palette: dark grey = as with the standard digest, blue = new, red = lost, orange dashed = newly separable.

**Verified:**
- natively and in the Shinylive build: goal switch about 0.6 s, entry switch (new digests) 0.62 s, all within the 1 s target
- card selection, grid click (loads the region across an entry change), specific lysines, custom range, the downloads and the preset link
- tests: 155 pass; the new figures are serialized as strict JSON for several scopes, including no-lysine and one-lysine edge cases

## 0.6.0.dev1 (2026-11-14): Phase 6 step 1 redone for STRATEGY_PICKER_2.md

Spec: `STRATEGY_PICKER_2.md`, which replaces the earlier addendum. The user said STRATEGY_PICKER.md had been replaced with v2; the file attached was `STRATEGY_PICKER_2.md`, and STRATEGY_PICKER.md on disk is still the old draft. CLAUDE.md §11/§12 now point to `STRATEGY_PICKER_2.md`.

**`strategy.py` rework** (strategy generation, chemistry and the Net = Gained − Lost structure are kept):
- **Goals** replace M1–M5 in the API. There are six: see lysines, measure one at a time (N = 1 by default, or 2), separate pairs, read together, known modification sites (filter: any / acetylation / methylation / ubiquitination) and variants. N and the mark filter are query-time options, so changing them needs no re-digest.
- **Separable pairs** (new): some peptide carries the first lysine but not the second, and some peptide carries the second but not the first.
- **Specific-lysines scope** (new).
- **Ranking:**
  - strategies with the same in-scope peptides as the standard digest are removed
  - strategies of one kind with identical in-scope peptides are merged into one recommendation with a combined title
  - only goal net > 0 is shown
  - sort order: net ↓, lysines lost ↑, LC-MS runs ↑, sample preps ↑, name
- **Unmodified variants** are now generated only for enzymes that cleave at lysine (trypsin, trypsin/P, Lys-C, Lys-C/P, Lys-N).
- **Text:** card titles and short names, headline and loss sentences, "what you give up", extra work, warning badges and tooltips, the no-card message, and the "Compared with:" line all come from template functions.
- **`compare_all_histones()`** builds the grid.
- **`advanced_table()`** lists every strategy with each goal's Standard / Strategy / Gained / Lost / Net.
- Region presets now come per family: N-terminal tail, C-terminal tail and whole protein.
- **Removed (v1 API):** `ranked(sort=...)`, `detail()`, `compare_top()`, `full_table()`, `SORTS`, `METRICS`.

**Decisions where v2 left room** (please review):
- **Specific lysines:** a peptide is in scope only if it carries ≥1 selected lysine. A peptide that merely lies between them is not. This also applies to "Tell variants apart".
- **Merging** happens only within one strategy kind. Strategies of different kinds stay separate even when their peptides are identical, because the bench work differs. Warnings on a merged card are the ones every member shares.
- **Headline number** = items gained, so it matches the list. When the goal also loses items (net is still > 0), a second line says so: "but 1 lysine no longer measurable on its own (K30)". For "See lysines" this information is already in "What you give up", so the line is omitted.
- **Mark filters:** acetylation = ac, acme; methylation = me, me1, me2, me3, acme; ubiquitination = ub (SUMO and UFM1 are not included).
- **C-terminal tail** preset (H2A 110–129) is offered only for H2A entries within 10 residues of canonical length. It is not offered for macroH2A, where residues 110–129 are not a tail.
- **Grid cells:** a fourth status, "nothing to measure", is used when the scope has nothing for the goal, e.g. H4 has no other H4 entry to tell apart, or there are no ubiquitination sites.
- **Grid rows:** H2A type 1 N- and C-tail, H2B type 1-C/E/F/G/I tail, H3.1 tail, H3.3 tail and H4 tail.

**User decision (2026-11-14):** strategies without propionylation stay eligible for answer cards by default and always carry the "Unmodified lysines only" badge. There is no toggle to hide them. This matches the current code, so no behaviour changes.

**Tests:** `tests/test_strategy.py` (20) covers every exact case in v2 §6. The parts of the v1 regression that still apply are kept. The tests were mutation-checked; two mutations survived at first, and stronger tests were added for both:
- treating a pair as separable when only one lysine is seen
- treating selected lysines as a span

## 0.6.0.dev0 (2026-11-14): Phase 6 step 1, histone strategy picker core

Spec: `STRATEGY_PICKER_1.md` (the newer of the two addendum files; it supersedes `STRATEGY_PICKER.md`, which had a proteome-wide scope).

**CLAUDE.md changes** (requested in the addendum):
- §8 Output 2: a discriminating peptide is now defined explicitly as one whose sequence (I = L) is not contained anywhere in the other entry's full sequence, not one that is absent from the other entry's digest. The code already did this since 0.1.1 (substring mapping), so no results change.
- §12: the "enzyme optimizer or recommender" exclusion is narrowed. The histone strategy picker is in scope; it enumerates and sorts and does not optimize. A proteome-wide picker is not built.
- §11: added Phase 6.

**`strategy.py`** (core, no UI):
- Strategies: `alone`, `parallel` and `sequential` per candidate enzyme, with costs of 1/1/1, 2/2/2 and 1/2/1 (preps / enzyme steps / LC-MS runs).
- Metrics M1–M5 over a region in mature numbering, each reported as Baseline, Strategy, Gained, Lost, Net (M4 also "of which split"), with named items.
- Outputs: ranked table with sort shortcuts and stable tie-breaks, detail tables, compare-top-3, region peptide list, settings header.
- Speed: H2B 1–35 with 57 strategies takes 0.03 s.

**Decisions where the addendum left room** (please review):
- **Skipped no-op strategies:** a propionylated `alone`/`parallel`/`sequential` strategy whose enzyme is already the baseline enzyme reproduces the baseline, so it is skipped. Trypsin under the Histone PTM preset therefore appears only as `Trypsin alone (unmodified)` and `parallel + Trypsin (unmodified)`.
- **Variants:** both variants (propionylated, unmodified) are generated for every enzyme. The UI toggle will choose which are shown.
- **K-dependent warning:** shown only for unmodified variants of enzymes that cleave at K. "Cleaves at K" means some cut site disappears when K is blocked, tested on all xKy contexts; this is true for Trypsin, Trypsin/P, Lys-C, Lys-C/P and Lys-N.
- **Propionylated chemistry:** if the baseline itself is not propionylated, the propionylated variant and `sequential` use `protein + peptide N-term`, the Histone PTM option.
- **Other warnings:** "middle-down length (>40 aa)" is shown when an in-scope identifiable peptide is longer than 40 aa. Rule flags (e.g. Proteinase K "nonspecific, low confidence") are added as warnings when such a rule is included.
- **`total_lost`** (tie-break) = sum of Lost over M1–M5 minus split pairs, since split pairs are not a loss of information.
- **Mature numbering** follows the baseline's Met-removal setting, as in the other histone views.
- **M5** compares against every other entry of the same family in `histone_annotation.tsv`, not only canonical entries.
- **Regions:** the region end is clipped to the entry length.

**Tests** (`tests/test_strategy.py`, 18):
- the exact H2B type 1-C regression from the addendum
- Net = Strategy − Baseline = Gained − Lost for every metric and strategy
- `parallel` never loses
- split pairs
- region rules on a synthetic family
- enumeration and chemistry
- K-dependence
- sequential Lys-C/Lys-N = baseline
- warnings
- stable ranking for every sort
- under 1 s per whole-protein table

The region-rule tests were mutation-checked.

## 0.5.1 (2026-11-14): scatter click opens the protein

**Fix:** clicking a point in the baseline-vs-alternative coverage scatter switched to the Protein viewer but left the protein box empty, so nothing was drawn. The protein box is a server-side selectize: the browser holds only the first 1,000 of 42,545 entries, and `update_selectize(selected=...)` without the choices silently drops any other entry. The click handler now resends the choices with `server=True`, so Shiny sends the selected entry to the browser. Clicking the same point twice now also reopens it. Verified natively and in the Shinylive build (Q9BUQ8-2: peptide map, peptides and the K/S/T/Y gained/lost table all load).

## 0.5.0 (2026-11-14): Phase 5, the app

**App** (`app/`, Shiny for Python 1.9 + shinywidgets/Plotly 5.23; runs natively or fully in the browser with Shinylive 0.10.16 / Pyodide 0.27.7):
- **Sidebar:** proteins (reviewed human proteome, or an uploaded FASTA with background "reference + upload" or "upload only"); preset; rule source; baseline and alternative enzymes; compare as alone / parallel / sequential; custom enzyme ("cut after", "cut before", "except before P", with the compiled regex shown); length; missed cleavages; m/z filter; chemistry; an Apply button; PNG/SVG figure download format.
- **Tabs:**
  - **Overview:** baseline value boxes; summary table; gained/lost diverging bars; baseline-vs-alternative coverage scatter (clicking a protein opens it in the Protein viewer); seen-vs-unseen lists; coverage distributions; all-pairs ranking table; Methods.
  - **Peptides:** live for the chosen proteins; Output 1, or Output 2 with the unique-only switch.
  - **Proteins:** Output 3.
  - **Digest stats:** length histogram with the identifiable window shaded; coverage distribution; amino-acid detectability with K/S/T/Y/C highlighted.
  - **Histones:** known-mark grid; isoform discrimination matrix and table; K × K co-occurrence matrix with letter legend; tail peptide maps.
  - **Protein viewer:** peptide map per digest, peptides, K/S/T/Y sites gained and lost.
- "Predicted (in silico)" badge in the header, and a settings strip above every figure.
- **Exports:** TSV and Excel for every table, each with the settings header listing the digests the table came from; Methods text download; figures via the Plotly toolbar in PNG or SVG.
- **Modes:** whole proteome with unchanged preset settings opens precomputed (instant). Uploads or changed settings run live with a progress bar; the live all-pairs ranking table is computed on demand.

**Core additions**
- `presets.DigestRecipe` / `recipes_for()`: one definition of each alternative digest and its label, shared by precompute and live mode.
- `results.py`: `PrecomputedResults` / `LiveResults`, one interface for the app. They give identical tables for the same input: tested on all 37 Histone PTM digests and on 3 proteome digests.
- `report.header_lines()`, `methods_text_from()`, `tsv_bytes()`, `excel_bytes()`.
- `precomputed.load_reference_bytes()` / `load_preset_bytes()`, so the browser can fetch files over HTTP.
- `scripts/build_site.py` builds `site/`: about 150 MB, no file over 25 MB, with `.nojekyll` for GitHub Pages.

**Fixed (found by running the app in the browser and the test suite inside Pyodide 0.27.7)**
- numpy 2.0.2 in wasm32 also rejects int64 input to `bincount` → `digest.bincount()`.
- Figures with NaN broke the widget session, because shinywidgets serializes with `allow_nan=False`. Heatmap gaps are now `None`, and a test checks every figure.
- A live whole-proteome run ran out of browser memory. Fixed with lean `LiveResults`, int32 digest positions and coverage counters, and numpy cut-site collection. The heap peak fell from 2.43 to 1.84 GB, and 7 live proteome digests now take 39 s in the browser.
- Progress bars now update in the browser: the live loops yield to the event loop.

**Deviation from CLAUDE.md, with reason**
- **Gained/lost colours are blue/red, not green/red** (§7 peptide map strips, and the gained/lost charts and mark grid for consistency). The palette validator fails green #0ca30c vs red #d03b3b for colour-blind readers (deutan ΔE 4.1, floor 8); blue #2a78d6 vs red passes (ΔE 23.8). Co-occurrence uses grey/blue/orange, as specified.

**Not done here:** deployment to GitHub Pages; it needs the user's repository and approval to publish.

## 0.4.0 (2026-11-14): Phase 4, precompute

**Uniqueness strategy (decided from Pyodide measurements; see `bench/pyodide/RESULTS.md`).** Live substring mapping of the whole proteome takes 16.5 s in Pyodide (k-mer index), and a plain scan takes 26 s for titin alone. Instead, `metrics.SharedPrefixIndex` stores, for every reference residue, the longest prefix that also occurs (I = L) in another entry, and in an entry of another gene (uint8, capped at 255). A peptide of length L at residue s is unique iff the stored value is < L. This is exactly the substring definition for **any** enzyme or setting:
- It was verified natively on all 5.6 M identifiable peptides of three proteome digests, against brute force on random proteins, and on the histones.
- It is built by prefix doubling in 31–35 s, once per reference release, and ships as a 6.5 MB file. Lookup takes 0.14 s in Pyodide for all 2.84 M proteome peptides.
- `Background(..., index=...)` uses it automatically. Entries uploaded by the user that are not in the reference fall back to live mapping, and they also correctly make reference peptides non-unique (tested). With the index, `n_proteins`/`n_genes` are 1 for unique peptides and <NA> otherwise.

**Precomputed data** (`scripts/precompute.py` → `data/precomputed/`, 16.4 min natively, peak memory 4.5 GB):

| File | MB |
|---|---|
| `preset_proteome.npz` / `.json` | 16.4 / 0.1 |
| `preset_histone_ptm.npz` / `.json` | 14.5 / 0.1 |
| `preset_middle_down.npz` / `.json` | 17.1 / 0.1 |
| `preset_protein_cleaver.npz` / `.json` | 12.3 / 0.07 |
| `reference.fasta.gz` | 12.8 |
| `uniqueness_index.npz` | 6.5 |
| **Total** (every file < 25 MB) | **79.9** |

- Each preset holds its baseline plus every enzyme of its mode alone, in parallel and sequentially: 37 digests, or 25 for the ExPASy Protein Cleaver preset.
- Per-entry arrays (uint16): identifiable, unique and total peptides; covered and unique-covered residues; residues and K/S/T/Y sites gained and lost versus the baseline.
- Per-digest JSON: summary metrics, comparison counts (the ranking table), length histogram (1–100 aa plus overflow) and amino-acid detectability.
- `manifest.json` records versions, the reference release, preset settings, sizes, sha256 and timings.
- Peptide-level detail stays live per protein. In Pyodide, opening a preset takes 0.09–0.14 s and the protein viewer takes 0.03–0.05 s.

**Other changes**
- `presets.py`: the four §5 presets as core data. **Protein Cleaver compatible uses Met removal "never"**, matching Protein Cleaver and the golden tests; §5 does not specify it.
- `precomputed.py`: build, write and load app data.
- `metrics.comparison_counts()` / `comparison_per_protein()`: vectorized comparison counts, now used by `ranking_table()`. `compare()` uses numpy keys for peptides gained and lost.
- **Fixed for wasm32:** `np.repeat`/`reduceat` with int64 counts or indices failed in Pyodide, where `np.intp` is int32. Both are now cast to `np.intp`.
- `bench/pyodide/`: a reproducible Pyodide (314.0.7) benchmark harness.

## 0.3.0 (2026-11-14): Phase 3, histone module

- **Review:** both tables are marked `REVIEWED`. The user approved them in the session ("Done, build outputs"); the files had not been edited, so the status lines were changed on the user's behalf. Re-running `scripts/build_known_marks.py` resets `known_marks.tsv` to DRAFT, so a regenerated table needs a new review.
- `histones.load_histone_inputs()` loads the 55-entry histone set together with the reviewed annotation and marks, and refuses tables marked DRAFT.
- **Output 1, known-mark visibility** (`mark_visibility`, `mark_grid`): one row per site and digest. Each row gives the modifications at the site, whether it is covered, the covering peptides (mature numbering), the K count on each peptide, and the status versus the baseline (seen in both / gained / lost / unseen in both). Optional filter by `mod_abbrev`.
- **Output 2, isoform discrimination** (`isoform_discrimination`, `discrimination_matrix`): for each pair of entries in the same family and each digest, the status is resolved / one-way / unresolved / identical sequence (identical after I = L folding), with the discriminating peptides listed. A peptide counts as discriminating if it does not map to the other entry under substring mapping (user decision 2026-11-14); `method="digest"` remains available. Matrix cell (A, B) = number of A's peptides that do not map to B.
- **Output 3, lysine co-occurrence** (`lysine_cooccurrence`): results per site and per pair (kept / new / split / lost; sites can also be lost). Matrix data per entry has `baseline_pair` and the comma-separated letters of the alternatives that make a pair new, split or lost; `letters` gives the legend.
- **Output 4, peptide maps:** `metrics.peptide_map()` is generic and also serves the §7 protein viewer. It returns peptides, lysine ticks, gained/lost strips versus the baseline, and coverage. `histones.tail_peptide_map()` limits it to an N-terminal window per family (`TAIL_WINDOWS`, mature numbering: H1 1–40, H2A 1–20, H2B 1–35, H3 1–45, H4 1–30). These windows are display defaults only.
- **Defaults:** visibility and discrimination use canonical entries. Co-occurrence uses H2A type 1, H2B type 1-C/E/F/G/I, H3.1, H3.3 and H4.
- `scripts/export_histone_results.py` writes `results/phase3/histones/`: the Histone PTM baseline plus 36 alternatives.
- **Tests:** each output is checked against a brute-force recomputation from the peptide tables, plus synthetic cases for identical-sequence, one-way and split statuses. A live histone digest plus all outputs runs in under 1 s.

## 0.3.0-draft, revision 1 (2026-11-14): user decisions on the histone annotation

**`category` replaces `is_canonical`.** The values are `canonical`, `canonical_divergent`, `tissue_specific`, `replacement_variant` and `putative`; `validate_annotation()` rejects any other value.

**Decisions (user, 2026-11-14)**
- **H2AC1** (Q96QV6) = `tissue_specific` (TH2A).
- **H2AC25** (Q7L7L0) = `canonical`, or `canonical_divergent` if it differs from canonical H2A type 1. It differs from P0C0S8 at UniProt 17 (T>S), 41 (A>S) and 100 (K>R), so it is **`canonical_divergent`**. Its note says "H2AW" is a legacy alias unrelated to plant H2A.W.
- **H2BC26** (Q8N257) = `canonical_divergent`.
- **H2BC12L** (P57053) = `replacement_variant`. Its note says K28R (mature 27) creates an Arg-C cleavage site after UniProt residue 28 (…DGR|KRKR…). A test checks this.
- **H3-7** (Q5TEC6) = `putative`.
- **macroH2A.1: both splice forms are separate rows.** The UniProt O75367 isoform section (release 2026_04, saved as `data/reference/uniprot_O75367.json`) gives:
  - isoform 1 = mH2A1.1 = isoform ID `O75367-2`, the displayed sequence → row `O75367`, "macroH2A1.1"
  - isoform 2 = mH2A1.2 = `O75367-1` → new row `O75367-1`, "macroH2A1.2"

  Isoform 3 (`O75367-3`, no common name) is not included. The `O75367-1` record is copied verbatim from the reference FASTA into `data/histones/histone_isoforms.fasta`; provenance is in `data/histones/metadata.json`. The user FASTA is unchanged.
- **Default histone views use `canonical` entries only** (`default_entries()`). Lysine co-occurrence defaults to one canonical representative per core family (H2A type 1, H2B type 1-C/E/F/G/I, H3.1, H4) **plus H3.3 as an explicit exception** (user decision 2026-11-14, consistent with CLAUDE.md §8). H3.3 stays a `replacement_variant` and is not in the other default views. The exception is listed in `histones.COOCCURRENCE_EXCEPTIONS`; `validate_annotation()` rejects any other non-canonical co-occurrence default.

**Other changes**
- `load_histone_set()` returns the user FASTA plus `histone_isoforms.fasta`: 55 entries.
- `transfer_marks()` copies UniProt marks from O75367 onto O75367-1 where the site lies in identical sequence. All 9 sites transfer; K320 maps to K323. Transferred rows carry `evidence_type = "transferred from O75367 (identical sequence)"`. `known_marks.tsv` now has 2,764 rows, 632 sites and 51 entries.
- I assigned the remaining categories with my own judgment and marked four of them `REVIEW:`: H1.1 (`canonical`; expression is tissue-restricted), H2A type 2-B / H2AC21 (`canonical`; its C-terminal tail differs from H2A type 1 at 6 positions), H3.Y and H3.X (`replacement_variant`; expression is limited).
- Both tables are still marked DRAFT, pending final review.

## 0.3.0-draft (2026-11-14): Phase 3 part 1, histone tables for review

- `histones.py` adds a UniProt feature parser, `known_marks_from_features()`, loaders and `validate_annotation()`. The loaders **refuse files whose `# status:` line says DRAFT** until the user marks them reviewed (CLAUDE.md section 8).
- `data/histones/known_marks.tsv` (DRAFT) is generated by `scripts/build_known_marks.py` from UniProt 2026_04 `ft_mod_res` and `ft_crosslnk`. It has 2,752 rows: 623 lysine sites on 50 entries; H2A.B, CENP-A, H1.7 and H1.9 have no lysine features. It keeps every lysine modification type in UniProt, including recent acylations (lactyl, isonicotinyl, benzoyl, …) and ub/SUMO2/UFM1 isopeptides, with UniProt evidence. 675 rows are "by similarity" only.
- Added columns beyond the spec, for convenience: `mod_abbrev`, `feature`, `evidence_type` and `uniprot_note` in `known_marks.tsv`; `gene`, `pe` and `cooccurrence_default` in the annotation.
- `data/histones/histone_annotation.tsv` (DRAFT) has 54 entries drafted from UniProt names using the Talbert et al. 2012 / HistoneDB 2.0 nomenclature. `is_canonical` = replication-dependent histone. Six judgment calls are marked `REVIEW:`.
- The histone FASTA keeps its file name `human_histones.fasta` (CLAUDE.md says `histones.fasta`). All 54 entries are kept; the three PE=5 entries are marked non-canonical.

## 0.2.0 (2026-11-14): Phase 2, data and outputs

**Added**
- `scripts/fetch_data.py` downloads two files from UniProt release **2026_04** (04-Nov-2026):
  - the reviewed human proteome with isoforms: 42,545 entries (20,414 canonical + 22,131 isoforms), 24.5 M residues
  - UniProt features (`ft_mod_res`, `ft_crosslnk`) for the 54 histone accessions

  Source URL, release, download date, SHA-256 and entry counts are recorded in `data/reference/metadata.json`. All 54 sequences in `human_histones.fasta` are identical to this release.
- `metrics.py`:
  - `Background`: uniqueness background plus mapping method
  - `Analysis`: Outputs 1–4
  - `compare()`: Output 5 lists and counts
  - `summary_table()`
  - `ranking_table()`: all-pairs table of alone, parallel with baseline, and sequential baseline + enzyme
- `report.py`: a settings header for every export, with the "PREDICTED (in silico)" notice, rule regexes, sources, versions and citations, and all settings; `write_tsv()`; `methods_text()`.
- `read_reference_fasta()` names a ProteinSet with its provenance, and `ProteinSet.name` is printed in exports.

**Definitions chosen where CLAUDE.md is not explicit** (flag any you want changed)
- Peptide counts (`n_identifiable_peptides`, `n_unique_peptides`, `n_total_peptides`, "confident" = ≥2 unique peptides) count peptides **by position** (entry, start, end), as Protein Cleaver does. With missed cleavages, overlapping peptides each count.
- `median_coverage_pct` is the median over **all** entries; non-identifiable entries count as 0 %. `overall_coverage_pct` = covered mature residues / all mature residues.
- "Proteins gaining unique peptides" = entries with more unique peptides in the alternative than in the baseline.
- Residues and K/S/T/Y sites gained or lost compare residues covered by identifiable peptides. Peptides gained or lost compare identifiable peptides by position.
- Output 3 lists identifiable proteins by default; `protein_table(identifiable_only=False)` returns every entry.
- In a parallel digest, the `digest` column of Output 1 names the component digest(s) that produce each peptide.
- Mapped entries are counted by distinct accession. When building the background, an uploaded entry is dropped only if both its accession and its sequence match an existing entry.

**Fixed**
- `parallel()` kept the missed-cleavage count of the first digest. It now keeps the smallest, so "total peptides" (0 missed cleavages) is right for parallel digests.

**Performance** (natively, whole reference, Proteome preset: Trypsin, 7–30 aa, ≤1 missed cleavage)
- Digest 0.9 s; uniqueness 5.7 s, including the ~3 s index build; Outputs 1–4 about 6 s.
- One comparison takes 6–11 s.
- Histone ranking table (Histone PTM preset, 13 enzymes, proteome background): 37 rows in 3.8 s.

## 0.1.1 (2026-11-14): peptide-to-protein mapping definition

**Decision (user, 2026-11-14).** A peptide maps to **every entry in the reference database whose sequence contains it**, with I and L treated as identical. This holds whether or not that entry's digest would produce the peptide. The same definition applies to **uniqueness** (`unique_protein`, `unique_gene`) and to **histone isoform discrimination** (Phase 3): a discriminating peptide of entry A must not occur as a substring of entry B.

**Effect.** This is the stricter of the two definitions: a peptide that occurs anywhere in another entry counts as shared with it. It replaces CLAUDE.md §8's wording "absent from the other's digest" for isoform discrimination.

**Changes**
- Added `metrics.py` with `peptide_matches()` (all peptide–entry pairs) and `peptide_mapping()` (counts and the `unique_protein` / `unique_gene` flags). `method="substring"` is the default. The previous definition is still available as `method="digest"`: a peptide maps only to entries whose cleavage products, under the same settings, include it. Mapping code moved here from `digest.py`.
- Substring search searches the full UniProt sequences; a removed initiator Met is kept. The search uses an exact index of the 8 residues starting at every reference position (packed into a uint64, so there are no hash collisions). Peptides of ≤8 residues are range lookups. Longer peptides are checked residue by residue. Matches never span two entries.
- Tests: a case where the two definitions disagree (SAMPLER is a tryptic product only of A, but B contains SAMPIER; substring → 2 entries, digest → 1). Brute-force `in` checks on the histone FASTA for both methods. A test that matches never cross entry junctions.
- Timing on a proteome-sized synthetic set (42,500 entries, 24 M residues): the index is built once per reference in ~3 s. Mapping takes 0.4 s (trypsin, 0 missed cleavages, 1.2 M peptides), 1.4 s (≤2 missed cleavages, 4.6 M peptides) and 3.9 s (Arg-C, 4–40 aa, ≤1 missed cleavage, 1.75 M peptides). Peak memory is a few GB natively; Pyodide feasibility is a Phase 4 question.
- Consequence: under this definition, very short peptides are rarely unique against a whole proteome. In the synthetic test only 13% of 4–6-aa peptides were protein-unique, against 99.5% of ≥7-aa peptides. A 4-mer occurs by chance in about 150 proteins.

## 0.1.0 (2026-11-14): Phase 1, core engine

Added `rules.py`, `fasta.py`, `chemistry.py` and `digest.py`, with rule-equivalence, unit, golden (Protein Cleaver), regression and performance tests.

### Decisions and their reasons

- **Rule source versions.** The PSI-MS regexes match `psi-ms.obo` data-version 4.2.17 (commit `4d2a234`, downloaded 2026-11-14; provenance is in `data/reference/metadata.json`). Every vendored rule gives the same cut positions as `pyteomics` 5.0.2 on 10,000 random sequences. ExPASy rules are stored in zero-width form; ExPASy trypsin includes its `trypsin_exception` sites.
- **Lys-N (PSI-MS mode) = `(?=K)`.** `psi-ms.obo` now contains Lys-N (MS:1003093). Its definition says it cleaves on the amino side of K, and its `has_regexp` comment reads `(?=K)`. The linked ID, however, is MS:1001335, the Lys-C regex `(?<=K)(?!P)`. This is a CV inconsistency. We follow the definition, the comment and Taouatas et al. 2008.
- **Proteinase K comes from ExPASy only.** `psi-ms.obo` now has "peptidase K" (MS:1003981), but its regex `(?=[AEFILTVWY])` cuts *before* these residues, which contradicts ExPASy and the enzyme's P1 specificity. Per the spec, Proteinase K uses the ExPASy rule (cut after) and is flagged "nonspecific; low confidence".
- **ExPASy-only enzymes in PSI-MS mode.** Chymotrypsin high/low specificity and Proteinase K can be chosen in PSI-MS mode, but they keep `source = "ExPASy"`, so exports always name the rule actually used.
- **Average masses use the ExPASy average residue table** (A 71.0788 … Y 163.1760, water 18.01524). A least-squares fit to all 824 distinct Protein Cleaver peptides recovers this table exactly. Our masses match Protein Cleaver within 4e-5 Da. Average masses of the mods use Unimod values (carbamidomethyl 57.0513, propionyl 56.0633).
- **Blocked lysines.** A cut site counts as K-derived if it disappears when every K is replaced by an unreactive placeholder. Blocking therefore only removes sites and never creates them; for example, the ExPASy trypsin exception CR|K stays uncut.
- **Met removal.** Rules are applied to the mature sequence (residue 2 onward) when Met is removed. `prev_aa` then reports `M`.
- **Ambiguous residues.** B, Z, J and X have no mass, so peptides containing them get NaN masses and fail the m/z filter. U and O have defined masses.
- **Uniqueness:** see the 0.1.1 entry above. In 0.1.0 a provisional digest-product definition was the default.
- **Directory names** were lowercased to match the spec (`Data/Histones` → `data/histones`, `Tests/Data/Protein_Cleaver` → `tests/data/protein_cleaver`). File names are unchanged.
