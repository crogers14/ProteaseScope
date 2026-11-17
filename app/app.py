"""ProteaseScope app (Shiny for Python; runs natively or in the browser via Shinylive).

The app only calls the core library (``proteasescope``). All results are predictions.
"""

from __future__ import annotations

import asyncio
import gc
import html
import itertools
import re
import sys
from dataclasses import dataclass, replace
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(APP_DIR))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import plotly.graph_objects as go  # noqa: E402
from shiny import App, reactive, render, req, ui  # noqa: E402
from shinywidgets import output_widget, render_plotly  # noqa: E402

import data_source as DS  # noqa: E402
import figures as F  # noqa: E402
from proteasescope import __version__  # noqa: E402
from proteasescope import histones as H  # noqa: E402
from proteasescope import metrics as M  # noqa: E402
from proteasescope import precomputed as P  # noqa: E402
from proteasescope import report as R  # noqa: E402
from proteasescope import strategy as S  # noqa: E402
from proteasescope import workflow as W  # noqa: E402
from proteasescope.chemistry import Chemistry, MetRemoval, Propionylation, met_removed  # noqa: E402
from proteasescope.digest import DigestSettings  # noqa: E402
from proteasescope.fasta import ProteinSet, parse_fasta  # noqa: E402
from proteasescope.presets import PRESETS, preset_recipes, recipes_for  # noqa: E402
from proteasescope.results import LiveResults, PrecomputedResults, Results  # noqa: E402
from proteasescope.rules import CleavageRule, available_enzymes, custom_rule, get_rule  # noqa: E402

PRE = "data/precomputed"
KINDS = {"alone": "Alone", "parallel": "Parallel with baseline",
         "sequential": "Sequential after baseline"}  # fmt: skip
DEFAULT_ALTS = {
    "proteome": ["Lys-C", "Glu-C", "Asp-N", "Chymotrypsin"],
    "histone_ptm": ["Asp-N", "Glu-C", "Chymotrypsin"],
    "middle_down": ["Asp-N", "Lys-C"],
    "protein_cleaver": ["Asp-N", "Lys-N", "Chymotrypsin high specificity"],
}
DEFAULT_KINDS = ["parallel", "sequential"]
LIVE_PROTEOME_WARN = 6  # live whole-proteome digests beyond this get a warning
OUTPUT3_COLS = M.OUTPUT3_COLUMNS


@dataclass(frozen=True)
class Snapshot:
    """The settings the results are computed from (taken on Apply or preset change)."""

    source: str
    preset: str
    mode: str
    baseline: DigestSettings
    rules: tuple[CleavageRule, ...]
    kinds: tuple[str, ...]
    background: str


def preset_snapshot(key: str, source: str, background: str) -> Snapshot:
    p = PRESETS[key]
    rules = tuple(get_rule(n, p.mode) for n in DEFAULT_ALTS[key])
    return Snapshot(source, key, p.mode, p.baseline, rules, tuple(DEFAULT_KINDS), background)


def strip_text(res: Results, indices: list[int]) -> str:
    """One-line settings strip shown above figures."""
    b = res.recipes[0].baseline
    c = b.chemistry
    alts = len([i for i in indices if i != 0])
    mz = f"m/z {b.mz_min:g}-{b.mz_max:g}" if b.mz_filter else "m/z filter off"
    return (
        f"Predicted (in silico) · {res.mode} · Baseline {b.label} · {alts} alternative "
        f"digest(s) · {b.min_length}-{b.max_length} aa · <={b.missed_cleavages} missed "
        f"cleavage(s) · {mz} · {'Carbamidomethyl C' if c.carbamidomethyl else 'no fixed mods'}"
        f" · {c.propionylation.label.lower()} · Met: {c.met_removal} · "
        f"{res.proteins.name or len(res.proteins)}"
    )


def round_df(df: pd.DataFrame, digits: int = 2) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        if pd.api.types.is_float_dtype(out[c]):
            out[c] = out[c].round(digits)
    return out


def dl_buttons(name: str) -> ui.Tag:
    return ui.div(
        ui.download_button(f"dl_{name}_tsv", "TSV", class_="btn-sm btn-outline-secondary"),
        ui.download_button(f"dl_{name}_xlsx", "Excel", class_="btn-sm btn-outline-secondary"),
        class_="d-flex gap-2 mb-2",
    )


def strip(id_: str) -> ui.Tag:
    return ui.div(ui.output_text(id_), class_="settings-strip")


# =========================================================================================
# UI
# =========================================================================================
CSS = """
.settings-strip { font-size: 0.78rem; color: #52514e; background: #f0efec;
  border-radius: 6px; padding: 4px 8px; margin: 4px 0 8px; }
.predicted-badge { font-size: 0.72rem; font-weight: 600; background: #fff3cd; color: #5c4500;
  border: 1px solid #e0c46c; border-radius: 999px; padding: 2px 10px; margin-left: 10px;
  vertical-align: middle; }
.note { color: #52514e; font-size: 0.85rem; }
.enz-info { display: flex; flex-wrap: wrap; gap: 4px; margin: -6px 0 10px; }
.enz-chip { font-size: 0.75rem; color: #52514e; background: #f0efec; border-radius: 999px;
  padding: 1px 8px; cursor: help; }
.enz-chip.has-note { background: #e6f0fb; color: #1d5ea8; }
.tooltip-inner .enz-tip { text-align: left; }
.wide-fig { overflow-x: auto; }
.wide-fig > * { min-width: 820px; }
.bslib-value-box .value-box-value { font-size: 1.5rem; }
.sp-standard { font-size: 0.9rem; color: #0b0b0b; margin: 6px 0 10px; }
.sp-card { border: 1px solid #e1e0d9; border-radius: 10px; padding: 12px 14px; height: 100%;
  background: #fff; cursor: pointer; transition: border-color .15s, box-shadow .15s; }
.sp-card:hover { border-color: #9ec5f4; }
.sp-card.selected { border: 2px solid #2a78d6; box-shadow: 0 1px 6px rgba(42,120,214,.18); }
.sp-rank { font-size: 0.75rem; color: #898781; font-weight: 600; }
.sp-title { font-size: 1rem; margin: 2px 0 8px; }
.sp-big { font-size: 2rem; font-weight: 700; color: #2a78d6; line-height: 1; margin-right: 8px; }
.sp-headline { margin-bottom: 8px; }
.sp-line { font-size: 0.86rem; margin: 3px 0; color: #0b0b0b; }
.sp-muted { color: #52514e; }
.sp-badge { display: inline-block; font-size: 0.72rem; background: #fff3cd; color: #5c4500;
  border: 1px solid #e0c46c; border-radius: 999px; padding: 1px 8px; margin: 6px 4px 0 0;
  cursor: help; }
.sp-message { background: #f0efec; border-radius: 8px; padding: 12px 14px; margin: 6px 0; }
.sp-howto { font-size: 0.8rem; color: #52514e; margin: 2px 0 14px; }
.sp-hover { text-decoration: underline dotted #898781; cursor: help; }
.sp-chip { display: inline-block; border-radius: 999px; padding: 1px 10px; font-size: 0.78rem;
  border: 1px solid #c3c2b7; }
.sp-grid td, .sp-grid th { padding: 6px 10px; font-size: 0.86rem; vertical-align: top; }
.sp-grid tr.sp-click { cursor: pointer; }
.sp-grid tr.sp-click:hover { background: #eef5fd; }
.bld-open { margin: 8px 0 10px; }
/* The builder's sidebar wrapper must not change the picker's width or indent: no padding
   around the main area, and no collapse arrow while closed (the button opens it). The arrow
   stays rendered (visibility, not display): bslib finishes opening the panel when the
   arrow's rotation transition ends. */
.bslib-sidebar-layout > .main.sp-layout,
.bslib-sidebar-layout.sidebar-collapsed > .main.sp-layout { padding: 0; }
.bslib-sidebar-layout:not(.sidebar-collapsed) > .main.sp-layout { padding-right: 1rem; }
.bslib-sidebar-layout.sidebar-collapsed:has(> .main.sp-layout) > .collapse-toggle {
  visibility: hidden; }
.bld-row { border: 1px solid #e1e0d9; border-radius: 8px; padding: 8px 10px 4px;
  margin-bottom: 8px; background: #fff; }
.bld-row-head { display: flex; justify-content: space-between; align-items: center;
  min-height: 1.6rem; }
.bld-label { font-weight: 600; font-size: 0.85rem; margin: 0; }
.bld-rm { border: none; background: none; color: #898781; font-size: 1.15rem; line-height: 1;
  padding: 0 4px; }
.bld-rm:hover { color: #d03b3b; background: none; }
.bld-row .shiny-input-container { margin-bottom: 6px; }
.bld-note { font-size: 0.8rem; color: #5c4500; background: #fff3cd; border-radius: 6px;
  padding: 3px 8px; margin: 2px 0 6px; }
.bld-note.blocking { color: #fff; background: #d03b3b; }
.bld-same { font-size: 0.82rem; color: #52514e; padding: 2px 0 8px; }
.bld-start { display: flex; gap: 6px; align-items: baseline; }
.bld-start .shiny-input-container { flex: 1; margin-bottom: 8px; }
.bld-buttons { display: flex; gap: 8px; margin: 4px 0 8px; }
.bld-files { display: flex; gap: 8px; align-items: flex-start; }
.bld-files .shiny-input-container { margin-bottom: 0; flex: 1; }
.bld-saved { font-size: 0.8rem; color: #52514e; margin-top: 8px; }
.sp-banner { background: #eef5fd; border: 1px solid #9ec5f4; border-radius: 8px;
  padding: 8px 12px; margin: 6px 0 10px; font-size: 0.9rem; }
.sp-allgoals td.num { text-align: right; font-variant-numeric: tabular-nums; }
.sp-allgoals td.gain { color: #2a78d6; font-weight: 600; }
.sp-allgoals td.loss { color: #d03b3b; font-weight: 600; }
"""

sidebar = ui.sidebar(
    ui.input_radio_buttons(
        "source", "Proteins",
        {"proteome": "Reviewed human proteome", "upload": "Upload FASTA"},
    ),  # fmt: skip
    ui.panel_conditional(
        "input.source === 'upload'",
        ui.input_file("fasta", "FASTA file", accept=[".fasta", ".fa", ".faa", ".txt"]),
        ui.input_radio_buttons(
            "background", "Uniqueness background",
            {"ref_upload": "Reference proteome + upload", "upload": "Upload only"},
        ),  # fmt: skip
    ),
    ui.input_select("preset", "Preset", {k: p.name for k, p in PRESETS.items()}),
    ui.accordion(
        ui.accordion_panel(
            "Enzymes",
            ui.input_radio_buttons("mode", "Rule source",
                                   {"psi-ms": "PSI-MS", "expasy": "ExPASy"}, inline=True),
            ui.input_select("baseline", "Baseline enzyme", available_enzymes("psi-ms")),
            ui.input_selectize("alts", "Alternative enzymes", available_enzymes("psi-ms"),
                               multiple=True, selected=DEFAULT_ALTS["proteome"]),
            ui.output_ui("enzyme_info"),
            ui.input_checkbox_group("kinds", "Compare as", KINDS, selected=DEFAULT_KINDS),
            ui.input_checkbox("use_custom", "Add a custom enzyme"),
            ui.panel_conditional(
                "input.use_custom",
                ui.input_text("cut_after", "Cut after (residues)", placeholder="e.g. KR"),
                ui.input_text("cut_before", "Cut before (residues)", placeholder="e.g. D"),
                ui.input_checkbox("except_p", "Except before P", True),
                ui.output_code("custom_regex"),
            ),
        ),  # fmt: skip
        ui.accordion_panel(
            "Peptides",
            ui.layout_columns(ui.input_numeric("min_len", "Min length", 7, min=1, max=200),
                              ui.input_numeric("max_len", "Max length", 30, min=1, max=255)),
            ui.input_numeric("mc", "Max missed cleavages", 1, min=0, max=5),
            ui.input_switch("mz_filter", "m/z filter (charges 2-4)", False),
            ui.panel_conditional(
                "input.mz_filter",
                ui.layout_columns(ui.input_numeric("mz_min", "m/z min", 350),
                                  ui.input_numeric("mz_max", "m/z max", 1400)),
            ),
        ),  # fmt: skip
        ui.accordion_panel(
            "Chemistry",
            ui.input_checkbox("cam", "Carbamidomethyl C (fixed)", True),
            ui.input_select("prop", "Derivatization",
                            {p.value: p.label for p in Propionylation}),
            ui.input_select("met", "Initiator Met removal", [m.value for m in MetRemoval],
                            selected=MetRemoval.METAP.value),
        ),  # fmt: skip
        open=["Enzymes"],
    ),
    ui.input_action_button("apply", "Apply settings", class_="btn-primary w-100"),
    ui.output_ui("mode_note"),
    ui.input_radio_buttons("fig_fmt", "Figure download format", {"png": "PNG", "svg": "SVG"},
                           inline=True),  # fmt: skip
    width=340,
)


def tab_overview():
    return ui.nav_panel(
        "Overview",
        ui.card(
            ui.card_header(ui.output_text("div_header", inline=True)),
            strip("strip_div"),
            ui.div(output_widget("div_fig"), class_="wide-fig"),
        ),
        ui.card(ui.card_header("Summary versus the baseline"), strip("strip_summary"),
                dl_buttons("summary"), ui.output_data_frame("summary_tbl")),
        ui.card(
            ui.card_header("Seen vs unseen"),
            ui.input_select("ov_alt", "Alternative digest", []),
            strip("strip_scatter"),
            output_widget("scatter_fig"),
            ui.p("Click a protein to open it in the Protein viewer.", class_="note"),
            ui.navset_underline(
                ui.nav_panel("Rescued proteins", dl_buttons("rescued"),
                             ui.output_data_frame("rescued_tbl")),
                ui.nav_panel("Lost proteins", dl_buttons("lost"),
                             ui.output_data_frame("lost_tbl")),
                ui.nav_panel("Gaining unique peptides", dl_buttons("gaining"),
                             ui.output_data_frame("gaining_tbl")),
                ui.nav_panel("Residues and K/S/T/Y sites", dl_buttons("changes"),
                             ui.output_data_frame("changes_tbl")),
            ),  # fmt: skip
        ),
        ui.card(strip("strip_dist"), output_widget("dist_fig")),
        ui.card(
            ui.card_header("All-pairs ranking table (every enzyme alone, parallel and "
                           "sequential with the baseline)"),
            ui.output_ui("rank_controls"),
            dl_buttons("ranking"),
            ui.output_data_frame("ranking_tbl"),
        ),  # fmt: skip
        ui.card(ui.card_header("Methods"), ui.output_code("methods"),
                ui.download_button("dl_methods", "Download methods (.txt)",
                                   class_="btn-sm btn-outline-secondary")),  # fmt: skip
        value="overview",
    )


def tab_peptides():
    return ui.nav_panel(
        "Peptides",
        ui.layout_columns(
            ui.input_select("pep_digest", "Digest", []),
            ui.input_text("pep_query", "Proteins (accessions or genes)", placeholder="H3C1 P62805"),
            ui.input_switch("pep_unique", "Unique peptides only (Output 2)", False),
        ),  # fmt: skip
        ui.output_ui("pep_note"),
        dl_buttons("peptides"),
        ui.output_data_frame("pep_tbl"),
        value="peptides",
    )


def tab_proteins():
    return ui.nav_panel(
        "Proteins",
        ui.layout_columns(ui.input_select("prot_digest", "Digest", []),
                          ui.input_switch("prot_all", "Include non-identifiable entries")),
        dl_buttons("proteins"),
        ui.output_data_frame("prot_tbl"),
        value="proteins",
    )  # fmt: skip


def tab_stats():
    return ui.nav_panel(
        "Digest stats",
        ui.input_select("stats_digest", "Digest", []),
        strip("strip_stats"),
        ui.output_ui("cov_summary"),
        ui.layout_columns(ui.card(output_widget("len_fig"), dl_buttons("lenhist")),
                          ui.card(output_widget("covhist_fig"))),
        ui.card(output_widget("aa_fig"), dl_buttons("aa")),
        value="stats",
    )  # fmt: skip


GOAL_CHOICES = {k: g.title for k, g in S.GOALS.items()}
MARK_FILTER_CHOICES = {
    "any": "Any modification",
    "acetylation": "Acetylation",
    "methylation": "Methylation",
    "ubiquitination": "Ubiquitination",
}
HOWTO = {
    "cards": "Each card is one way to change the digest. The big number is what you gain for "
             "your goal; \"What you give up\" lists lysines here that the standard digest "
             "sees and this strategy does not. Click a card to show it below.",
    "lysine": "Each column is a lysine, each row a digest. Dark grey: as with the standard "
              "digest; blue: new with that strategy; red: lost; white: not achieved by either. "
              "A dot above a column marks a lysine with known modifications (hover to list them).",
    "pair": "Lysines sit on the line. An arc above joins two lysines found together on one "
            "peptide (grey: also with the standard digest; blue: new; red dotted: no longer). "
            "An orange dashed arc below joins two lysines that this strategy lets you measure "
            "separately for the first time.",
    "variant": "Each row is another entry of the same histone family. \"Told apart\" means at "
               "least one predicted peptide from here exists in your histone but nowhere in that "
               "entry (I and L count as the same). Hover a row to see the peptides.",
    "map": "Bars are the predicted peptides that touch this region, one row per digest; "
           "triangles mark lysines. Blue strips mark residues only the strategy covers, red "
           "strips residues only the standard digest covers. The % is coverage of the region.",
    "grid": "Each row answers the selected goal for one histone region: the best strategy and "
            "what it gains, or why there is nothing to gain. Click a row to open it above.",
}


CHEM_CHOICES = {p.value: p.label for p in Propionylation}
COMBINE_CHOICES = {k: ui.span(W.COMBINE_LABELS[k], title=W.COMBINE_TOOLTIPS[k])
                   for k in W.COMBINES}  # fmt: skip
START_CHOICES = {"": "Start from…", **{k: k for k in W.BUILDER_PRESETS}}
STANDARD_CHOICE = {"standard": "Standard digest (current preset)"}


def builder_row(i: int) -> ui.Tag:
    """One digest row of the builder (shown while ``bld_n > i``)."""
    chem = ui.input_select(f"bld_chem{i}", None, CHEM_CHOICES, width="100%")
    if i > 0:
        chem = ui.div(
            ui.panel_conditional("input.bld_combine !== 'same_tube'", chem),
            ui.panel_conditional("input.bld_combine === 'same_tube'",
                                 ui.div("Same chemistry as digest 1", class_="bld-same")),
        )  # fmt: skip
    head = ui.div(
        ui.tags.label(f"Digest {i + 1}", class_="bld-label"),
        ui.panel_conditional("input.bld_n > 1", ui.input_action_button(
            f"bld_rm{i}", "×", class_="bld-rm", title="Remove this digest")),
        class_="bld-row-head",
    )  # fmt: skip
    row = ui.div(head, ui.input_select(f"bld_enz{i}", None, [], width="100%"), chem,
                 ui.output_ui(f"bld_note{i}"), class_="bld-row")  # fmt: skip
    return ui.panel_conditional(f"input.bld_n > {i}", row)


def builder_panel() -> list:
    """Contents of the "Build your own strategy" side panel (CUSTOM_STRATEGY.md section 3)."""
    return [
        ui.input_text("bld_name", "Name", "My strategy 1", width="100%"),
        ui.div(ui.input_select("bld_start", None, START_CHOICES), ui.output_ui("bld_cite"),
               class_="bld-start"),
        ui.div(ui.input_numeric("bld_n", None, 1, min=1, max=W.MAX_DIGESTS),
               style="display:none"),
        *[builder_row(i) for i in range(W.MAX_DIGESTS)],
        ui.panel_conditional(
            f"input.bld_n < {W.MAX_DIGESTS}",
            ui.input_action_button("bld_add", "+ Add digest",
                                   class_="btn-sm btn-outline-secondary mb-2"),
        ),
        ui.panel_conditional(
            "input.bld_n > 1",
            ui.input_radio_buttons("bld_combine", "How they're combined", COMBINE_CHOICES),
            ui.panel_conditional("input.bld_combine === 'same_tube'",
                                 ui.div(W.NOTE_SAME_TUBE, class_="bld-note")),
        ),
        ui.input_select("bld_against", "Compare against", STANDARD_CHOICE, width="100%"),
        ui.div(ui.input_action_button("bld_compare", "Compare", class_="btn-primary btn-sm"),
               ui.input_action_button("bld_save", "Save", class_="btn-outline-primary btn-sm"),
               class_="bld-buttons"),
        ui.div(ui.download_button("bld_download", "Download", class_="btn-sm btn-link"),
               ui.input_file("bld_load", None, accept=[".json"], button_label="Load…",
                             placeholder="strategy .json"),
               class_="bld-files"),
        ui.output_ui("bld_saved"),
    ]  # fmt: skip


def tab_picker():
    return ui.nav_panel(
        "Strategy picker",
        ui.layout_sidebar(
            ui.sidebar(*builder_panel(), position="right", open="closed", id="sp_builder",
                       width=380, title="Build your own strategy"),
            *picker_main(),
            border=False, fillable=False, class_="sp-layout",
        ),
        value="picker",
    )  # fmt: skip


def picker_main() -> list:
    return [
        ui.input_action_button("bld_open", "Build your own strategy",
                               class_="btn-outline-primary btn-sm bld-open"),
        ui.layout_columns(
            ui.input_selectize("sp_entry", "1. Histone", choices={}),
            ui.div(
                ui.input_select("sp_where", "2. Where to look", {}),
                ui.panel_conditional(
                    "input.sp_where === 'custom'",
                    ui.layout_columns(ui.input_numeric("sp_from", "From (mature)", 1, min=1),
                                      ui.input_numeric("sp_to", "To (mature)", 40, min=1)),
                ),
                ui.panel_conditional(
                    "input.sp_where === 'lysines'",
                    ui.input_selectize("sp_lys", "Lysines (mature numbering)", [],
                                       multiple=True),
                ),
            ),
            ui.div(
                ui.input_select("sp_goal", "3. What are you looking for?", GOAL_CHOICES),
                ui.panel_conditional(
                    "input.sp_goal === 'alone'",
                    ui.input_radio_buttons("sp_n", "Lysines allowed on the peptide",
                                           {"1": "N = 1", "2": "N = 2"}, inline=True),
                ),
                ui.panel_conditional(
                    "input.sp_goal === 'marks'",
                    ui.input_select("sp_marks", "Modification type", MARK_FILTER_CHOICES),
                ),
            ),
            col_widths=[4, 4, 4],
        ),
        ui.output_ui("sp_standard"),
        ui.output_ui("sp_cards"),
        ui.output_ui("sp_all_goals"),
        ui.p(HOWTO["cards"], class_="sp-howto"),
        strip("strip_sp"),
        ui.panel_conditional("input.sp_goal !== 'variants'", output_widget("sp_picture")),
        ui.output_ui("sp_variants"),
        ui.output_ui("sp_picture_note"),
        output_widget("sp_map"),
        ui.p(HOWTO["map"], class_="sp-howto"),
        ui.download_button("dl_sp_answer", "Download this answer (text)",
                           class_="btn-sm btn-outline-secondary mb-3"),
        ui.hr(),
        ui.input_action_button("sp_grid_go", "Compare across all histones",
                               class_="btn-outline-primary btn-sm"),
        ui.output_ui("sp_grid"),
        ui.accordion(
            ui.accordion_panel(
                "Show all strategies (advanced)",
                ui.p("Every strategy, including those identical to the standard digest, with "
                     "each goal as Standard, Strategy, Gained, Lost and Net "
                     "(Net = Strategy − Standard = Gained − Lost).", class_="note"),
                dl_buttons("sp_adv"),
                ui.output_data_frame("sp_adv"),
            ),
            open=False,
            class_="mt-3",
        ),
    ]  # fmt: skip


def tab_histones():
    return ui.nav_panel(
        "Histones",
        ui.panel_conditional(
            "input.hist_tabs !== 'picker'",
            ui.p("Histone outputs use the bundled histone set (54 user-provided entries + "
                 "macroH2A1.2) with the current baseline and alternative digests.",
                 class_="note"),
            strip("strip_hist"),
            ui.output_data_frame("hist_letters"),
        ),  # fmt: skip
        ui.navset_underline(
            tab_picker(),
            ui.nav_panel(
                "Known marks",
                ui.layout_columns(
                    ui.input_selectize("hm_entries", "Entries", [], multiple=True),
                    ui.input_selectize("hm_mods", "Modifications (empty = all)", [],
                                       multiple=True),
                ),
                output_widget("hm_fig"), dl_buttons("marks"), ui.output_data_frame("hm_tbl"),
            ),
            ui.nav_panel(
                "Isoform discrimination",
                ui.layout_columns(
                    ui.input_select("iso_family", "Family", ["H3", "H2A", "H2B", "H1"]),
                    ui.input_select("iso_digest", "Digest", []),
                    ui.input_radio_buttons("iso_scope", "Entries",
                                           {"canonical": "Canonical", "all": "All in family"},
                                           inline=True),
                ),
                output_widget("iso_fig"), dl_buttons("iso"), ui.output_data_frame("iso_tbl"),
            ),
            ui.nav_panel(
                "Lysine co-occurrence",
                ui.input_select("co_entry", "Entry", []),
                output_widget("co_fig"), dl_buttons("cooc"), ui.output_data_frame("co_tbl"),
            ),
            ui.nav_panel(
                "Tail peptide maps",
                ui.layout_columns(ui.input_select("tail_entry", "Entry", []),
                                  ui.input_numeric("tail_from", "From (mature)", 1, min=1),
                                  ui.input_numeric("tail_to", "To (mature)", 45, min=1)),
                output_widget("tail_fig"),
            ),
            id="hist_tabs",
        ),  # fmt: skip
        value="histones",
    )


def tab_viewer():
    return ui.nav_panel(
        "Protein viewer",
        ui.input_selectize("protein", "Protein", choices=[], width="100%"),
        strip("strip_viewer"),
        output_widget("map_fig"),
        ui.navset_underline(
            ui.nav_panel("Peptides", dl_buttons("viewer_peps"), ui.output_data_frame("viewer_peps")),
            ui.nav_panel("K/S/T/Y sites gained and lost", dl_buttons("viewer_sites"),
                         ui.output_data_frame("viewer_sites")),
        ),  # fmt: skip
        value="viewer",
    )


app_ui = ui.page_sidebar(
    sidebar,
    ui.head_content(ui.tags.style(CSS)),
    ui.navset_card_underline(
        tab_overview(), tab_peptides(), tab_proteins(), tab_stats(), tab_histones(),
        tab_viewer(), id="main_tabs",
    ),  # fmt: skip
    title=ui.tags.span("ProteaseScope", ui.tags.span("Predicted (in silico)",
                                                     class_="predicted-badge")),
    fillable=False,
)


def enzyme_tip(r: CleavageRule) -> list:
    """Tooltip body for an enzyme: bench note, literature evidence, suggested missed
    cleavages, rule and notes."""
    tip: list = []
    if r.ui_note:
        tip += [ui.tags.b(r.ui_note), ui.tags.br()]
    if r.evidence:
        tip += [ui.tags.i(r.evidence), ui.tags.br()]
    if r.suggested_missed_cleavages is not None and f"≤{r.suggested_missed_cleavages}" not in (
            r.ui_note or ""):
        tip += [f"Suggested max missed cleavages: {r.suggested_missed_cleavages}", ui.tags.br()]
    return tip + [f"{r.label}: regex {r.regex}", ui.tags.br(), r.notes]


def next_strategy_name(taken) -> str:
    """``My strategy N`` with the smallest N not in ``taken``."""
    k = 1
    while f"My strategy {k}" in taken:
        k += 1
    return f"My strategy {k}"


# =========================================================================================
# Server
# =========================================================================================
def server(input, output, session):  # noqa: C901 - one server function by Shiny convention
    snapshot = reactive.value(preset_snapshot("proteome", "proteome", "ref_upload"))
    clicked = reactive.value(None)  # (accession, click number)
    clicks = itertools.count(1)
    protein_choices: reactive.Value[dict[str, str]] = reactive.value({})
    live_ranking = reactive.value(None)
    _presets: dict[str, P.PrecomputedPreset] = {}

    # ----- settings ------------------------------------------------------------------------
    def form_baseline() -> DigestSettings:
        rule = get_rule(input.baseline(), input.mode())
        chem = Chemistry(carbamidomethyl=input.cam(), propionylation=input.prop(),
                         met_removal=input.met())  # fmt: skip
        return DigestSettings(
            (rule,), missed_cleavages=int(input.mc() or 0), min_length=int(input.min_len() or 1),
            max_length=int(input.max_len() or 1), mz_filter=input.mz_filter(),
            mz_min=float(input.mz_min() or 350), mz_max=float(input.mz_max() or 1400),
            chemistry=chem,
        )  # fmt: skip

    def custom() -> CleavageRule | None:
        if not input.use_custom():
            return None
        return custom_rule(input.cut_after() or "", input.cut_before() or "", input.except_p())

    @render.ui
    def enzyme_info():
        """One chip per selected enzyme; hover shows its rule, source, notes and any bench note
        (e.g. ProAlanase's acidic digestion)."""
        mode = input.mode()
        names = list(dict.fromkeys([input.baseline(), *(input.alts() or [])]))
        chips = []
        for name in names:
            try:
                r = get_rule(name, mode)
            except KeyError:
                continue
            noted = bool(r.ui_note or r.evidence)
            chip = ui.span(f"{name} ⓘ", class_="enz-chip has-note" if noted else "enz-chip")
            chips.append(ui.tooltip(chip, ui.div(*enzyme_tip(r), class_="enz-tip"),
                                    placement="right"))  # fmt: skip
        return ui.div(*chips, class_="enz-info") if chips else None

    @render.code
    def custom_regex():
        try:
            r = custom()
            return f"regex: {r.regex}" if r else ""
        except ValueError as e:
            return str(e)

    @reactive.effect
    @reactive.event(input.preset)
    def _apply_preset():
        apply_preset(input.preset())

    def apply_preset(key: str) -> None:
        """Set every sidebar control to the preset and recompute (also used by the picker's
        "Use the Histone PTM settings" link)."""
        p = PRESETS[key]
        b = p.baseline
        enz = available_enzymes(p.mode)
        ui.update_radio_buttons("mode", selected=p.mode)
        ui.update_select("baseline", choices=enz, selected=b.rules[0].name)
        ui.update_selectize("alts", choices=enz, selected=DEFAULT_ALTS[key])
        ui.update_checkbox_group("kinds", selected=DEFAULT_KINDS)
        ui.update_numeric("min_len", value=b.min_length)
        ui.update_numeric("max_len", value=b.max_length)
        ui.update_numeric("mc", value=b.missed_cleavages)
        ui.update_switch("mz_filter", value=b.mz_filter)
        ui.update_checkbox("cam", value=b.chemistry.carbamidomethyl)
        ui.update_select("prop", selected=b.chemistry.propionylation.value)
        ui.update_select("met", selected=b.chemistry.met_removal.value)
        ui.update_checkbox("use_custom", value=False)
        snapshot.set(preset_snapshot(key, input.source(), input.background()))

    @reactive.effect
    @reactive.event(input.mode, ignore_init=True)
    def _mode_changed():
        enz = available_enzymes(input.mode())
        with reactive.isolate():
            base = input.baseline() if input.baseline() in enz else enz[0]
            alts = [a for a in (input.alts() or []) if a in enz]
        ui.update_select("baseline", choices=enz, selected=base)
        ui.update_selectize("alts", choices=enz, selected=alts)

    @reactive.effect
    @reactive.event(input.apply)
    def _apply():
        try:
            rules = [get_rule(n, input.mode()) for n in (input.alts() or [])]
            c = custom()
            if c is not None:
                rules.append(c)
            snap = Snapshot(input.source(), input.preset(), input.mode(), form_baseline(),
                            tuple(rules), tuple(input.kinds() or ("parallel",)),
                            input.background())  # fmt: skip
        except ValueError as e:
            ui.notification_show(f"Invalid settings: {e}", type="error")
            return
        live_ranking.set(None)
        snapshot.set(snap)

    # ----- data --------------------------------------------------------------------------
    @reactive.calc
    async def reference():
        with ui.Progress() as p:
            p.set(message="Loading the reference proteome and uniqueness index…")
            man = await DS.read_bytes(f"{PRE}/manifest.json")
            fasta = await DS.read_bytes(f"{PRE}/reference.fasta.gz")
            idx = await DS.read_bytes(f"{PRE}/uniqueness_index.npz")
            return P.load_reference_bytes(man, fasta, idx)

    @reactive.calc
    def upload() -> ProteinSet | None:
        files = input.fasta()
        if not files:
            return None
        f = files[0]
        try:
            ps = parse_fasta(Path(f["datapath"]).read_text())
        except ValueError as e:
            ui.notification_show(f"FASTA error: {e}", type="error", duration=10)
            return None
        ps.name = f"{f['name']} ({len(ps)} entries)"
        return ps

    async def load_preset(key: str, ref: ProteinSet) -> P.PrecomputedPreset:
        if key not in _presets:
            with ui.Progress() as p:
                p.set(message=f"Loading precomputed {PRESETS[key].name} results…")
                js = await DS.read_bytes(f"{PRE}/preset_{key}.json")
                npz = await DS.read_bytes(f"{PRE}/preset_{key}.npz")
                _presets[key] = P.load_preset_bytes(key, js, npz, ref.accessions)
        return _presets[key]

    @reactive.calc
    async def results() -> tuple[Results, list[int]]:
        snap = snapshot()
        if snap.source == "upload":
            up = upload()
            req(up is not None)
            if snap.background == "upload":
                bg = M.Background(up)
            else:
                ref, idx, _ = await reference()
                bg = M.Background(ProteinSet.concat([ref, up]), index=idx)
            recipes = recipes_for(snap.baseline, list(snap.rules), snap.kinds)
            res = LiveResults(up, bg, recipes)
            return await _compute_all(res), list(range(len(recipes)))
        ref, idx, _ = await reference()
        bg = M.Background(ref, index=idx)
        key = next((k for k, p in PRESETS.items() if p.baseline == snap.baseline
                    and all(r.source != "Custom" for r in snap.rules)), None)  # fmt: skip
        if key is not None:
            res = PrecomputedResults(ref, bg, preset_recipes(PRESETS[key]),
                                     await load_preset(key, ref))  # fmt: skip
            sel = [i for i, r in enumerate(res.recipes)
                   if i == 0 or (r.rule in snap.rules and r.kind in snap.kinds)]  # fmt: skip
            return res, sel
        recipes = recipes_for(snap.baseline, list(snap.rules), snap.kinds)
        if len(recipes) > LIVE_PROTEOME_WARN:
            ui.notification_show(
                f"Live whole-proteome digest of {len(recipes)} digests; this can take a "
                "minute in the browser.", duration=8)  # fmt: skip
        res = LiveResults(ref, bg, recipes)
        return await _compute_all(res), list(range(len(recipes)))

    async def _compute_all(res: LiveResults) -> LiveResults:
        with ui.Progress(min=0, max=len(res.recipes)) as p:
            for i, r in enumerate(res.recipes):
                p.set(i, message=f"Digesting ({i + 1}/{len(res.recipes)})", detail=r.label)
                await asyncio.sleep(0.05)  # let the browser show progress (single thread)
                res.record(i)
                gc.collect()
        return res

    @render.ui
    async def mode_note():
        res, sel = await results()
        if res.mode == "precomputed":
            msg = "Precomputed whole-proteome results (instant)."
        else:
            msg = f"Live results: {len(res.proteins):,} entries, {len(res.recipes)} digests."
        return ui.p(msg, class_="note mt-2")

    # ----- digest choice lists -------------------------------------------------------------
    @reactive.effect
    async def _digest_choices():
        res, sel = await results()
        labels = {str(i): res.recipes[i].label for i in sel}
        alts = {k: v for k, v in labels.items() if k != "0"}
        ui.update_select("ov_alt", choices=alts, selected=next(iter(alts), None))
        for id_ in ("pep_digest", "prot_digest", "stats_digest"):
            ui.update_select(id_, choices=labels, selected="0")
        hist_labels = {str(k): lab for k, lab in enumerate(labels.values())}
        ui.update_select("iso_digest", choices=hist_labels, selected="0")
        choices = {a: f"{g or '—'} · {a} · {n or ''}"[:90] for a, g, n in zip(
            res.proteins.meta["accession"], res.proteins.meta["gene"],
            res.proteins.meta["protein_name"], strict=True)}  # fmt: skip
        default = "P68431" if "P68431" in choices else next(iter(choices))
        protein_choices.set(choices)
        ui.update_selectize("protein", choices=choices, selected=default, server=True)

    # ----- settings strips -----------------------------------------------------------------
    async def _strip() -> str:
        res, sel = await results()
        return strip_text(res, sel)

    for _sid in ("strip_summary", "strip_div", "strip_scatter", "strip_dist", "strip_stats",
                 "strip_hist", "strip_viewer"):  # fmt: skip
        output(id=_sid)(render.text(_strip))

    # ----- downloads -----------------------------------------------------------------------
    def register_download(name: str, frame_fn, title: str, indices_fn=None):
        async def header_and_frame():
            res, sel = await results()
            df = await frame_fn()
            idx = await indices_fn() if indices_fn else sel
            return res.header(title, idx), df

        @output(id=f"dl_{name}_tsv")
        @render.download(filename=f"proteasescope_{name}.tsv")
        async def _tsv():
            head, df = await header_and_frame()
            yield R.tsv_bytes(df, head)

        @output(id=f"dl_{name}_xlsx")
        @render.download(filename=f"proteasescope_{name}.xlsx")
        async def _xlsx():
            head, df = await header_and_frame()
            yield R.excel_bytes({name: df}, head)

    def fig_widget(fig: go.Figure, name: str) -> go.FigureWidget:
        w = go.FigureWidget(fig)
        w._config = {"displaylogo": False, "toImageButtonOptions": {
            "format": input.fig_fmt(), "filename": f"proteasescope_{name}", "scale": 2}}
        return w

    # ----- Overview ------------------------------------------------------------------------
    async def summary_frame():
        res, sel = await results()
        return round_df(res.summary_table(sel))

    @render.data_frame
    async def summary_tbl():
        return render.DataGrid(await summary_frame(), width="100%")

    register_download("summary", summary_frame, "Output 5: Summary versus the baseline")

    @render.text
    async def div_header():
        res, _ = await results()
        return f"Gained and lost versus the baseline: {res.recipes[0].label}"

    @render_plotly
    async def div_fig():
        res, sel = await results()
        s0 = res.summary(0)
        aa0 = res.aa_detectability(0).set_index("residue")["on_identifiable_peptides"]
        base = {"peptides": s0["n_identifiable_peptides"], "K_sites": aa0.get("K", 0),
                "STY_sites": sum(aa0.get(r, 0) for r in "STY"),
                "median_coverage_pct": s0["median_coverage_pct"]}  # fmt: skip
        rows = []
        for i in sel[1:]:
            c = res.comparison_counts(i)
            rec = res.recipes[i]
            rows.append({
                "digest": f"{rec.rule.name} · {rec.kind}", "full_label": rec.label,
                "peptides_gained": c["peptides_gained"], "peptides_lost": c["peptides_lost"],
                "K_sites_gained": c["K_sites_gained"], "K_sites_lost": c["K_sites_lost"],
                "STY_sites_gained": sum(c[f"{r}_sites_gained"] for r in "STY"),
                "STY_sites_lost": sum(c[f"{r}_sites_lost"] for r in "STY"),
                "median_coverage_pct": res.summary(i)["median_coverage_pct"],
            })  # fmt: skip
        return fig_widget(F.diverging_changes(pd.DataFrame(rows), base), "gained_lost")

    async def _pair_indices() -> list[int]:
        return [0, await alt_index()]

    async def alt_index() -> int:
        res, sel = await results()
        v = input.ov_alt()
        req(v)
        i = int(v)
        req(i in sel)
        return i

    @render_plotly
    async def scatter_fig():
        res, _ = await results()
        i = await alt_index()
        fig = F.coverage_scatter(res.coverage_pair(i), res.recipes[0].label, res.recipes[i].label)
        w = fig_widget(fig, "coverage_scatter")

        def on_click(trace, points, state):
            if points.point_inds:
                # the counter makes a repeated click on the same point fire again
                clicked.set((str(trace.customdata[points.point_inds[0]][0]), next(clicks)))

        w.data[1].on_click(on_click)
        return w

    @reactive.effect
    @reactive.event(clicked)
    def _open_clicked():
        acc = clicked()[0] if clicked() else None
        choices = protein_choices()
        if acc and acc in choices:
            # Server-side selectize: the browser holds only the first 1000 options, so
            # `selected` alone is dropped for any other entry. Resending the choices with
            # server=True makes Shiny ship the selected option to the client.
            ui.update_selectize("protein", choices=choices, selected=acc, server=True)
            ui.update_navs("main_tabs", selected="viewer")

    def seen_frame(key: str):
        async def f():
            res, _ = await results()
            return round_df(res.seen_vs_unseen(await alt_index())[key])

        return f

    for _key, _name in (("rescued_proteins", "rescued"), ("lost_proteins", "lost"),
                        ("proteins_gaining_unique", "gaining"),
                        ("residue_and_site_changes", "changes")):  # fmt: skip
        _fn = seen_frame(_key)

        def _make(fn=_fn):
            async def _render():
                return render.DataGrid(await fn(), width="100%", filters=True)

            return _render

        output(id=f"{_name}_tbl")(render.data_frame(_make()))
        register_download(_name, _fn, f"Output 5: {_key.replace('_', ' ')}",
                          lambda: _pair_indices())  # fmt: skip

    @render_plotly
    async def dist_fig():
        res, sel = await results()
        vals = {res.recipes[i].label: res.entry_stats(i)["max_coverage_pct"].to_numpy()
                for i in sel}  # fmt: skip
        return fig_widget(F.coverage_distributions(vals), "coverage_distributions")

    @render.ui
    async def rank_controls():
        res, _ = await results()
        if res.mode == "precomputed":
            return ui.p("Precomputed: all enzymes of the preset's rule source.", class_="note")
        return ui.div(
            ui.p("Live mode: the table needs every enzyme alone, in parallel and sequentially "
                 "(may take minutes for the whole proteome).", class_="note"),
            ui.input_action_button("rank_go", "Compute all-pairs table",
                                   class_="btn-sm btn-outline-primary"),
        )  # fmt: skip

    @reactive.effect
    @reactive.event(input.rank_go)
    async def _rank_live():
        res, _ = await results()
        snap = snapshot()
        rules = [get_rule(n, snap.mode) for n in available_enzymes(snap.mode)]
        rules += [r for r in snap.rules if r.source == "Custom"]
        full = LiveResults(res.proteins, res.background, recipes_for(snap.baseline, rules))
        n = len(full.recipes)
        with ui.Progress(min=0, max=n) as p:
            for i, r in enumerate(full.recipes):
                p.set(i, message=f"Ranking {i + 1}/{n}", detail=r.label)
                await asyncio.sleep(0.05)
                full.record(i)
                gc.collect()
        live_ranking.set(full.ranking_table())

    async def ranking_frame():
        res, _ = await results()
        if res.mode == "precomputed":
            return round_df(res.ranking_table())
        t = live_ranking()
        req(t is not None)
        return round_df(t)

    @render.data_frame
    async def ranking_tbl():
        return render.DataGrid(await ranking_frame(), width="100%", filters=True)

    register_download("ranking", ranking_frame, "Output 5: All-pairs ranking table")

    @render.code
    async def methods():
        res, sel = await results()
        text = R.methods_text_from([res.recipes[i].settings for i in sel],
                                   str(res.proteins.name), res.background.description)
        if saved():
            text += (" Custom strategies (Strategy picker): "
                     + "; ".join(w.description for w in saved().values()) + ".")  # fmt: skip
        return text

    @render.download(filename="proteasescope_methods.txt")
    async def dl_methods():
        res, sel = await results()
        text = R.methods_text_from([res.recipes[i].settings for i in sel],
                                   str(res.proteins.name), res.background.description)
        yield "\n".join(res.header("Methods", sel)) + "\n\n" + text + "\n"

    # ----- Peptides --------------------------------------------------------------------------
    async def peptide_frame():
        res, sel = await results()
        v = input.pep_digest()
        req(v)
        i = int(v)
        q = (input.pep_query() or "").replace(",", " ").split()
        meta = res.proteins.meta
        if q:
            ql = {t.lower() for t in q}
            accs = meta.loc[meta["accession"].str.lower().isin(ql)
                            | meta["gene"].fillna("").str.lower().isin(ql), "accession"]
            accs = accs.tolist()[:200]  # fmt: skip
        elif len(res.proteins) <= 500:
            accs = res.proteins.accessions
        else:
            accs = []
        req(accs)
        df = res.peptide_table(i, accs)
        if input.pep_unique():
            df = df[df["unique_protein"]].reset_index(drop=True)
        return round_df(df, 4)

    @render.ui
    async def pep_note():
        res, _ = await results()
        if len(res.proteins) > 500 and not (input.pep_query() or "").strip():
            return ui.p("Enter accessions or gene names (up to 200 proteins); peptides are "
                        "computed live for the proteins you choose.", class_="note")  # fmt: skip
        return None

    @render.data_frame
    async def pep_tbl():
        return render.DataGrid(await peptide_frame(), width="100%", filters=True)

    async def _pep_idx() -> list[int]:
        return [int(input.pep_digest())]

    register_download("peptides", peptide_frame, "Output 1/2: Identifiable peptides", _pep_idx)

    # ----- Proteins -------------------------------------------------------------------------
    async def protein_frame():
        res, _ = await results()
        v = input.prot_digest()
        req(v)
        st = res.entry_stats(int(v))
        if not input.prot_all():
            st = st[st["identifiable"]]
        return round_df(st[OUTPUT3_COLS].reset_index(drop=True))

    @render.data_frame
    async def prot_tbl():
        return render.DataGrid(await protein_frame(), width="100%", filters=True)

    async def _prot_idx() -> list[int]:
        return [int(input.prot_digest())]

    register_download("proteins", protein_frame, "Output 3: Identifiable proteins", _prot_idx)

    # ----- Digest stats -----------------------------------------------------------------------
    async def stats_i() -> int:
        v = input.stats_digest()
        req(v)
        return int(v)

    @render.ui
    async def cov_summary():
        res, _ = await results()
        i = await stats_i()
        s = res.summary(i)
        st = res.entry_stats(i)
        return ui.p(
            f"{res.recipes[i].label}: overall coverage {s['overall_coverage_pct']:.1f}% of "
            f"mature residues; median per protein {s['median_coverage_pct']:.1f}%, mean "
            f"{st['max_coverage_pct'].mean():.1f}%; {s['n_total_peptides']:,} cleavage products "
            "(0 missed cleavages).", class_="note")  # fmt: skip

    async def len_frame():
        res, _ = await results()
        return res.length_histogram(await stats_i())

    async def aa_frame():
        res, _ = await results()
        return round_df(res.aa_detectability(await stats_i()), 4)

    @render_plotly
    async def len_fig():
        return fig_widget(F.length_histogram(await len_frame()), "length_histogram")

    @render_plotly
    async def covhist_fig():
        res, _ = await results()
        cov = res.entry_stats(await stats_i())["max_coverage_pct"].to_numpy()
        return fig_widget(F.coverage_histogram(cov), "coverage_histogram")

    @render_plotly
    async def aa_fig():
        return fig_widget(F.aa_detectability(await aa_frame()), "aa_detectability")

    async def _stats_idx() -> list[int]:
        return [await stats_i()]

    register_download("lenhist", len_frame, "Output 4: Cleaved peptide frequency", _stats_idx)
    register_download("aa", aa_frame, "Output 4: Amino acid detectability", _stats_idx)

    # ----- Histones -------------------------------------------------------------------------
    @reactive.calc
    def histone_inputs():
        return H.load_histone_inputs(DS.histone_dir())

    @reactive.calc
    async def histone_analyses() -> list[M.Analysis]:
        res, sel = await results()
        hset, _, _ = histone_inputs()
        cache: dict = {}
        return [M.Analysis(res.recipes[i].run(hset, cache)) for i in sel]

    @reactive.effect
    def _histone_choices():
        hset, ann, marks = histone_inputs()
        names = dict(zip(ann["accession"], ann["display_name"], strict=True))
        cats = dict(zip(ann["accession"], ann["category"], strict=True))
        entries = {a: f"{names[a]} ({a}, {cats[a]})" for a in hset.accessions}
        ui.update_selectize("hm_entries", choices=entries, selected=["P68431", "P62805"])
        ui.update_selectize("hm_mods", choices=sorted(set(marks["mod_abbrev"])), selected=[])
        ui.update_select("co_entry", choices=entries, selected="P68431")
        ui.update_select("tail_entry", choices=entries, selected="P68431")

    @render.data_frame
    async def hist_letters():
        an = await histone_analyses()
        letters = ["Baseline"] + H._letters(len(an) - 1)
        return render.DataGrid(pd.DataFrame({"letter": letters, "digest": [a.label for a in an]}),
                               width="100%")  # fmt: skip

    async def marks_frame():
        an = await histone_analyses()
        _, ann, marks = histone_inputs()
        entries = list(input.hm_entries() or [])
        req(entries)
        mods = list(input.hm_mods() or []) or None
        return H.mark_visibility(an, marks, ann, entries=entries, modifications=mods)

    @render_plotly
    async def hm_fig():
        an = await histone_analyses()
        letters = dict(zip([a.label for a in an], ["Base"] + H._letters(len(an) - 1),
                           strict=True))  # fmt: skip
        return fig_widget(F.mark_grid(await marks_frame(), letters), "known_marks")

    @render.data_frame
    async def hm_tbl():
        return render.DataGrid(await marks_frame(), width="100%", filters=True)

    register_download("marks", marks_frame, "Histones: known-mark visibility")

    async def iso_frame():
        an = await histone_analyses()
        _, ann, _ = histone_inputs()
        fam = ann[ann["family"] == input.iso_family()]
        entries = (fam["accession"] if input.iso_scope() == "all"
                   else fam.loc[fam["category"] == "canonical", "accession"]).tolist()  # fmt: skip
        req(len(entries) >= 2)
        return H.isoform_discrimination(an, ann, entries=entries)

    @render_plotly
    async def iso_fig():
        an = await histone_analyses()
        tab = await iso_frame()
        v = input.iso_digest()
        req(v)
        label = an[int(v)].label
        mat = H.discrimination_matrix(tab, label, input.iso_family())
        return fig_widget(F.isoform_matrix(mat, f"{input.iso_family()}: peptides of the row "
                                           f"entry absent from the column entry — {label}"),
                          "isoform_matrix")  # fmt: skip

    @render.data_frame
    async def iso_tbl():
        return render.DataGrid(await iso_frame(), width="100%", filters=True)

    register_download("iso", iso_frame, "Histones: isoform discrimination")

    async def cooc():
        an = await histone_analyses()
        _, ann, _ = histone_inputs()
        acc = input.co_entry()
        req(acc)
        return H.lysine_cooccurrence(an, ann, entries=[acc])

    @render_plotly
    async def co_fig():
        co = await cooc()
        sites = co["sites"]
        base = sites[sites["status"] == "baseline"]
        lys = base["K"].tolist()
        name = base["display_name"].iloc[0] if len(base) else input.co_entry()
        return fig_widget(F.cooccurrence_matrix(co["matrix"], lys,
                                                f"Lysine co-occurrence — {name}"),
                          "cooccurrence")  # fmt: skip

    async def cooc_frame():
        return (await cooc())["pairs"]

    @render.data_frame
    async def co_tbl():
        return render.DataGrid(await cooc_frame(), width="100%", filters=True)

    register_download("cooc", cooc_frame, "Histones: lysine co-occurrence pairs")

    @reactive.effect
    @reactive.event(input.tail_entry)
    def _tail_window():
        _, ann, _ = histone_inputs()
        fam = ann.set_index("accession").loc[input.tail_entry(), "family"]
        lo, hi = H.TAIL_WINDOWS[fam]
        ui.update_numeric("tail_from", value=lo)
        ui.update_numeric("tail_to", value=hi)

    @render_plotly
    async def tail_fig():
        an = await histone_analyses()
        _, ann, _ = histone_inputs()
        acc = input.tail_entry()
        req(acc, input.tail_from(), input.tail_to())
        pm = M.peptide_map(an, acc, (int(input.tail_from()), int(input.tail_to())))
        name = ann.set_index("accession").loc[acc, "display_name"]
        return fig_widget(F.peptide_map(pm, [a.label for a in an], f"{name} tail peptide map"),
                          "tail_map")  # fmt: skip

    # ----- Strategy picker (Histones tab) --------------------------------------------------
    sp_pending_where: reactive.Value[str | None] = reactive.value(None)
    #: Saved builder strategies (this browser session), by name.
    saved: reactive.Value[dict[str, W.Workflow]] = reactive.value({})
    #: (strategy, compared-against strategy or None) shown after Compare; None = picker.
    bld_shown: reactive.Value[tuple[W.Workflow, W.Workflow | None] | None] = \
        reactive.value(None)  # fmt: skip
    sp_selected = reactive.value(0)
    sp_grid_on = reactive.value(False)

    @reactive.calc
    def sp_baseline() -> DigestSettings:
        """The standard digest: propionylation + trypsin with the sidebar's settings."""
        snap = snapshot()
        b = snap.baseline
        return replace(b, rules=(get_rule("Trypsin", snap.mode),), name=None,
                       chemistry=S.propionylated(b.chemistry))  # fmt: skip

    @reactive.calc
    def sp_candidates() -> list[CleavageRule]:
        snap = snapshot()
        rules = S.default_candidates([get_rule(n, snap.mode) for n in available_enzymes(snap.mode)])
        return rules + [r for r in snap.rules if r.mode == "custom"]

    def sp_shift(seq: str) -> int:
        with reactive.isolate():
            return int(met_removed(seq, sp_baseline().chemistry.met_removal))

    @reactive.effect
    def _sp_entries():
        hset, ann, _ = histone_inputs()
        grouped: dict[str, dict[str, str]] = {}
        a = ann[ann["accession"].isin(hset.accessions)].copy()
        a["first"] = a["category"] != "canonical"
        for fam in H.FAMILIES:
            sub = a[a["family"] == fam].sort_values("first", kind="mergesort")
            grouped[fam] = {
                r.accession: f"{r.display_name} ({r.accession})"
                + ("" if r.category == "canonical" else f" · {r.category.replace('_', ' ')}")
                for r in sub.itertuples()
            }
        ui.update_selectize("sp_entry", choices=grouped, selected="P68431")

    @reactive.effect
    @reactive.event(input.sp_entry)
    def _sp_where_choices():
        hset, ann, _ = histone_inputs()
        acc = input.sp_entry()
        req(acc in hset.accessions)
        seq = hset.sequences[hset.index_of(acc)]
        shift = sp_shift(seq)
        length = len(seq) - shift
        fam = ann.set_index("accession").loc[acc, "family"]
        choices = {}
        for name, rng in S.region_presets(fam, length).items():
            a, b = rng if rng else (1, length)
            choices[name] = f"{name} ({a}–{min(b, length)})"
        choices |= {"custom": "Custom range", "lysines": "Specific lysines"}
        with reactive.isolate():
            pending = sp_pending_where()
        sp_pending_where.set(None)
        ui.update_select("sp_where", choices=choices,
                         selected=pending if pending in choices else next(iter(choices)))
        lys = [i + 1 - shift for i, c in enumerate(seq) if c == "K" and i >= shift]
        ui.update_selectize("sp_lys", choices={str(k): f"K{k}" for k in lys}, selected=[])
        ui.update_numeric("sp_to", value=min(40, length))

    @reactive.calc
    def sp_state() -> tuple[S.StrategyResults | None, str | None]:
        """Picker results for the current entry and scope, or an explanation."""
        hset, ann, marks = histone_inputs()
        acc, where = input.sp_entry(), input.sp_where()
        req(acc in hset.accessions, where)
        region, lysines = where, None
        if where == "custom":
            a, b = input.sp_from(), input.sp_to()
            if not a or not b or a < 1 or b < a:
                return None, "Enter a residue range with From ≤ To (mature numbering)."
            region = (int(a), int(b))
        elif where == "lysines":
            lysines = [int(x) for x in (input.sp_lys() or [])]
            if not lysines:
                return None, "Pick one or more lysines."
            region = None
        shown = bld_shown()
        try:
            if shown is not None:
                wf, against = shown
                res = S.compare_workflows(hset, ann, marks, acc, region, sp_baseline(), wf,
                                          against, lysines=lysines)  # fmt: skip
            else:
                res = S.pick_strategies(hset, ann, marks, acc, region, sp_baseline(),
                                        sp_candidates(), lysines=lysines,
                                        custom=list(saved().values()))  # fmt: skip
        except KeyError:
            req(False)  # region list is being refreshed for a new entry
        except ValueError as e:
            return None, str(e)
        return res, None

    def sp_res() -> S.StrategyResults:
        res, _ = sp_state()
        req(res is not None)
        return res

    def sp_goal() -> tuple[str, int, str]:
        return input.sp_goal() or "see", int(input.sp_n() or 1), input.sp_marks() or "any"

    @reactive.calc
    def sp_cards_list() -> list[S.Card]:
        goal, n, mk = sp_goal()
        shown = bld_shown()
        if shown is not None:
            return [sp_res().card_for(S.Strategy.custom(shown[0]), goal, n, mk)]
        return sp_res().cards(goal, n, mk)

    @reactive.effect
    def _sp_reset_selection():
        sp_cards_list()
        sp_selected.set(0)

    @reactive.effect
    @reactive.event(input.sp_card_click)
    def _sp_card_click():
        sp_selected.set(int(input.sp_card_click()))

    def sp_option() -> S.Option | None:
        cards = sp_cards_list()
        i = sp_selected()
        return cards[i].option if cards and i < len(cards) else None

    def _md(text: str) -> ui.HTML:
        return ui.HTML(re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", html.escape(text)))

    def _hover(text: str, full: list[str], details: dict[str, list[str]] | None = None):
        tip = ", ".join(full)
        if details:
            tip = "\n".join(f"{k}: {', '.join(v)}" for k, v in details.items())
        return ui.span(text, title=tip, class_="sp-hover") if tip else ui.span(text)

    def _card(c: S.Card, i: int, selected: bool) -> ui.Tag:
        sentence = c.headline.split("** ", 1)[1]
        custom_view = bld_shown() is not None
        parts = [
            ui.div("Your strategy" if custom_view else f"Option {c.rank}", class_="sp-rank"),
            ui.div(_md(c.title), class_="sp-title"),
            ui.div(ui.span(f"+{c.gain_n}", class_="sp-big"),
                   _hover(sentence, c.headline_items, c.hover), class_="sp-headline"),
        ]  # fmt: skip
        if c.goal_loss:
            parts.append(ui.p(_hover(c.goal_loss.capitalize(), []), class_="sp-line sp-muted"))
        parts += [
            ui.p(ui.strong("What you give up: "), _hover(c.gives_up, c.gives_up_items),
                 class_="sp-line"),
            (ui.p(ui.strong("Work: "), c.work, class_="sp-line") if c.option.kind == "custom"
             else ui.p(ui.strong("Extra work: "), c.extra_work, class_="sp-line")),
            ui.div(*[ui.span(b, title=tip, class_="sp-badge") for b, tip in c.warnings]),
        ]  # fmt: skip
        return ui.div(
            *parts, class_="sp-card" + (" selected" if selected else ""), role="button",
            tabindex="0",
            onclick=f"Shiny.setInputValue('sp_card_click', {i}, {{priority: 'event'}})",
        )  # fmt: skip

    @render.ui
    def sp_standard():
        shown = bld_shown()
        if shown is not None:
            wf, against = shown
            ref = (_md(f"**{against.name}** ({against.description.split(': ', 1)[1]})")
                   if against else _md(S.standard_line(sp_baseline())[len("Compared with: "):]))
            return ui.div(
                ui.strong(wf.name), ": ", wf.description.split(": ", 1)[1], ui.br(),
                ui.span("Compared against: ", class_="sp-muted"), ref, " · ",
                ui.input_action_link("sp_back", "Back to the recommendations"),
                class_="sp-banner")  # fmt: skip
        b = sp_baseline()
        ptm = PRESETS["histone_ptm"].baseline
        same = (b.min_length, b.max_length, b.missed_cleavages) == (
            ptm.min_length, ptm.max_length, ptm.missed_cleavages)
        hint = [] if same else [
            " ", ui.span("(These lengths come from the sidebar.", class_="sp-muted"), " ",
            ui.input_action_link("sp_use_ptm", "Use the Histone PTM preset settings"),
            ui.span(")", class_="sp-muted")]  # fmt: skip
        return ui.p(_md(S.standard_line(b)), *hint, class_="sp-standard")

    @reactive.effect
    @reactive.event(input.sp_use_ptm)
    def _sp_use_ptm():
        if input.preset() == "histone_ptm":
            apply_preset("histone_ptm")
        else:
            ui.update_select("preset", selected="histone_ptm")

    @render.ui
    def sp_cards():
        res, err = sp_state()
        if err:
            return ui.div(err, class_="sp-message")
        req(res is not None)
        goal, n, mk = sp_goal()
        cards = sp_cards_list()
        if not cards:
            return ui.div(_md(res.message(goal, n, mk) or ""), class_="sp-message")
        sel = sp_selected()
        return ui.layout_columns(*[_card(c, i, i == sel) for i, c in enumerate(cards)],
                                 col_widths=[12 // len(cards)] * len(cards))  # fmt: skip

    @render.text
    def strip_sp():
        res = sp_res()
        c = res.baseline.chemistry
        shown = bld_shown()
        what = (f"{shown[0].description} vs {res.ref_name}" if shown is not None
                else f"{len(sp_candidates())} candidate enzymes"
                + (f" + {len(saved())} saved strateg{'y' if len(saved()) == 1 else 'ies'}"
                   if saved() else ""))  # fmt: skip
        return (f"Predicted (in silico) · {S.standard_line(res.baseline)[len('Compared with: '):]}"
                f" · {res.scope.description} · Met: {c.met_removal} · {what}")  # fmt: skip

    ACHIEVED = {"see": "seen", "alone": "measurable on its own", "marks": "seen"}

    @render_plotly
    def sp_picture():
        res = sp_res()
        goal, n, mk = sp_goal()
        req(goal != "variants")
        where = res.scope.where
        if S.GOALS[goal].unit == "lysine":
            opts = [c.option for c in sp_cards_list()]
            states = res.lysine_states(goal, opts, n, mk)
            numbered = bld_shown() is None
            labels = {o.label: f"{i}. {o.label}" if numbered else o.label
                      for i, o in enumerate(opts, start=1)}  # fmt: skip
            states["row"] = states["row"].map(lambda r: labels.get(r, r))
            ref = "the standard digest" if res.reference.kind == "baseline" else res.ref_name
            fig = F.lysine_picture(states, ACHIEVED[goal], f"Lysines in {where}", ref)
        else:
            opt = sp_option()
            label = (opt.label if bld_shown() is not None
                     else f"{sp_selected() + 1}. {opt.label}") if opt else None  # fmt: skip
            fig = F.arc_view(res.pair_states(opt),
                             [res.scope.k_name(k) for k in res.scope.lysines], label,
                             f"Lysine pairs in {where}", res.ref_label)  # fmt: skip
        return fig_widget(fig, "strategy_lysines")

    @render.ui
    def sp_picture_note():
        goal, _, _ = sp_goal()
        return ui.p(HOWTO[S.GOALS[goal].unit], class_="sp-howto")

    VARIANT_STYLE = {
        "told apart": ("#6b6a65", "#fff"), "newly told apart": ("#2a78d6", "#fff"),
        "no longer told apart": ("#d03b3b", "#fff"), "not told apart": ("#fff", "#0b0b0b"),
        "identical sequence": ("#f0efec", "#52514e"),
    }  # fmt: skip

    @render.ui
    def sp_variants():
        goal, _, _ = sp_goal()
        req(goal == "variants")
        res = sp_res()
        opt = sp_option()
        v = res.variant_states(opt)
        if v.empty:
            return ui.div(f"No other {res.scope.family} entries in the histone set.",
                          class_="sp-message")  # fmt: skip
        who = f"Option {sp_selected() + 1}: {opt.label}" if opt else "The standard digest"
        rows = []
        for r in v.itertuples():
            bg, fg = VARIANT_STYLE[r.status]
            rows.append(ui.tags.tr(
                ui.tags.td(r.entry), ui.tags.td(ui.span(r.status, class_="sp-chip",
                                                        style=f"background:{bg};color:{fg}")),
                title=r.peptides.replace("; ", "\n") or "no discriminating peptide here",
            ))  # fmt: skip
        return ui.div(
            ui.p(ui.strong(f"{res.scope.display_name} vs the other {res.scope.family} entries"),
                 ui.span(f" — {who}", class_="sp-muted")),
            ui.tags.table(ui.tags.thead(ui.tags.tr(ui.tags.th("Entry"), ui.tags.th("Status"))),
                          ui.tags.tbody(*rows), class_="table table-sm sp-grid"),
        )  # fmt: skip

    @render_plotly
    def sp_map():
        res = sp_res()
        opt = sp_option()
        custom_view = bld_shown() is not None
        label = (f"{opt.label} (combined)" if custom_view else f"{sp_selected() + 1}. {opt.label}"
                 ) if opt else ""  # fmt: skip
        pm = res.map_data(opt, label, split=custom_view)
        w0, w1 = res.map_window()
        fig = F.peptide_map(pm, pm["labels"], f"Peptides in {res.scope.where}", res.ref_name)
        fig.update_xaxes(range=[w0 - 1, w1 + 1])
        return fig_widget(fig, "strategy_peptide_map")

    @render.download(filename="proteasescope_strategy_answer.txt")
    def dl_sp_answer():
        goal, n, mk = sp_goal()
        res = sp_res()
        text = res.answer_text(goal, n, mk)
        shown = bld_shown()
        if shown is not None:
            card = res.card_for(S.Strategy.custom(shown[0]), goal, n, mk)
            summary = res.all_goals_summary(S.Strategy.custom(shown[0]), n, mk)
            text = (text.split("\n\n", 1)[0] + "\n\n" + card.text().replace("**", "")
                    + "\n\nAll goals:\n" + summary.to_csv(sep="\t", index=False))  # fmt: skip
        yield text

    # -- compare across all histones
    @reactive.effect
    @reactive.event(input.sp_grid_go)
    def _sp_grid_on():
        sp_grid_on.set(True)

    @reactive.calc
    def sp_grid_frame() -> pd.DataFrame:
        req(sp_grid_on())
        hset, ann, marks = histone_inputs()
        goal, n, mk = sp_goal()
        return S.compare_all_histones(hset, ann, marks, sp_baseline(), sp_candidates(),
                                      goals=[goal], n=n, mark_filter=mk)  # fmt: skip

    @render.ui
    def sp_grid():
        g = sp_grid_frame()
        goal, _, _ = sp_goal()
        rows = []
        for r in g.itertuples():
            cell = ([ui.strong(r.best), ui.br(), ui.span(r.gain)] if r.status == "best"
                    else [ui.span(r.status, class_="sp-muted")])  # fmt: skip
            key = f"{r.accession}|{r.region_name}"
            rows.append(ui.tags.tr(
                ui.tags.td(r.histone), ui.tags.td(f"{r.region_name} ({r.region})"),
                ui.tags.td(*cell), class_="sp-click",
                onclick=f"Shiny.setInputValue('sp_grid_click', '{key}', "
                        "{priority: 'event'})"))  # fmt: skip
        return ui.div(
            ui.p(ui.strong(f"All canonical histones — {S.GOALS[goal].title}"), class_="mt-3"),
            ui.tags.table(ui.tags.thead(ui.tags.tr(ui.tags.th("Histone"), ui.tags.th("Region"),
                                                   ui.tags.th("Best strategy"))),
                          ui.tags.tbody(*rows), class_="table table-sm sp-grid"),
            ui.p(HOWTO["grid"], class_="sp-howto"),
        )  # fmt: skip

    @reactive.effect
    @reactive.event(input.sp_grid_click)
    def _sp_grid_click():
        acc, region = input.sp_grid_click().split("|", 1)
        if input.sp_entry() == acc:
            ui.update_select("sp_where", selected=region)
        else:
            sp_pending_where.set(region)
            ui.update_selectize("sp_entry", selected=acc)

    # -- custom strategy view
    @reactive.effect
    @reactive.event(input.sp_back)
    def _sp_back():
        bld_shown.set(None)

    @render.ui
    def sp_all_goals():
        shown = bld_shown()
        req(shown is not None)
        res = sp_res()
        _, n, mk = sp_goal()
        s = res.all_goals_summary(S.Strategy.custom(shown[0]), n, mk)
        head = f"{res.ref_label} → {shown[0].name}"

        def cell(v, cls, items):
            return ui.tags.td(str(v), class_=f"num {cls if v else ''}", title=items or None)

        rows = [ui.tags.tr(
            ui.tags.td(r["goal"]), ui.tags.td(r["Compared-against → Custom"], class_="num"),
            cell(r["Gained"], "gain", r["gained_items"]), cell(r["Lost"], "loss", r["lost_items"]),
            ui.tags.td(f"{r['Net']:+d}", class_="num")) for _, r in s.iterrows()]  # fmt: skip
        return ui.div(
            ui.p(ui.strong("All goals"), ui.span(f" — {res.scope.where}; hover a number to "
                                                  "list the lysines, pairs or variants",
                                                  class_="sp-muted"), class_="mt-3 mb-1"),
            ui.tags.table(
                ui.tags.thead(ui.tags.tr(*[ui.tags.th(h) for h in
                                           ("Goal", head, "Gained", "Lost", "Net")])),
                ui.tags.tbody(*rows), class_="table table-sm sp-grid sp-allgoals"),
        )  # fmt: skip

    # -- builder ("Build your own strategy" side panel)
    @reactive.effect
    @reactive.event(input.bld_open)
    def _bld_open():
        ui.update_sidebar("sp_builder", show=True)

    def bld_n() -> int:
        return max(1, min(W.MAX_DIGESTS, int(input.bld_n() or 1)))

    def bld_default_chem() -> str:
        """New digests start with the chemistry of the current settings (preset)."""
        return snapshot().baseline.chemistry.propionylation.value

    @reactive.effect
    def _bld_enzymes():
        enz = W.builder_enzymes(available_enzymes(snapshot().mode))
        with reactive.isolate():
            cur = [input[f"bld_enz{i}"]() for i in range(W.MAX_DIGESTS)]
        for i in range(W.MAX_DIGESTS):
            ui.update_select(f"bld_enz{i}", choices=enz,
                             selected=cur[i] if cur[i] in enz else enz[min(i, len(enz) - 1)])

    @reactive.effect
    @reactive.event(input.preset)
    def _bld_chem_default():
        chem = PRESETS[input.preset()].baseline.chemistry.propionylation.value
        with reactive.isolate():
            n, started = bld_n(), bool(input.bld_start())
        for i in range(W.MAX_DIGESTS):
            if i >= n or (n == 1 and not started):
                ui.update_select(f"bld_chem{i}", selected=chem)

    @reactive.calc
    def bld_workflow() -> W.Workflow:
        mode = snapshot().mode
        steps = []
        for i in range(bld_n()):
            name = input[f"bld_enz{i}"]()
            req(name)
            try:
                rule = get_rule(name, mode)
            except KeyError:
                req(False)
            steps.append(W.DigestStep(rule, input[f"bld_chem{i}"]() or "none"))
        combine = (input.bld_combine() or "separate") if len(steps) > 1 else "separate"
        name = (input.bld_name() or "").strip() or next_strategy_name(saved())
        return W.Workflow(name, tuple(steps), combine)

    def bld_fill(wf: W.Workflow) -> None:
        """Fill the panel's rows from a strategy."""
        enz = W.builder_enzymes(available_enzymes(snapshot().mode))
        missing = [d.enzyme.name for d in wf.digests if d.enzyme.name not in enz]
        if len(wf.digests) > W.MAX_DIGESTS or missing or any(len(d.enzymes) > 1
                                                             for d in wf.digests):
            ui.notification_show(
                "This strategy cannot be shown in the builder"
                + (f": {', '.join(missing)} is not available with these rules" if missing
                   else f" (at most {W.MAX_DIGESTS} digests of one enzyme each)"),
                type="error")  # fmt: skip
            return
        ui.update_text("bld_name", value=wf.name)
        ui.update_numeric("bld_n", value=len(wf.digests))
        for i, d in enumerate(wf.digests):
            ui.update_select(f"bld_enz{i}", selected=d.enzyme.name)
            ui.update_select(f"bld_chem{i}", selected=d.chemistry.value)
        ui.update_radio_buttons("bld_combine", selected=wf.combine)

    @reactive.effect
    @reactive.event(input.bld_start)
    def _bld_start():
        key = input.bld_start()
        req(key)
        bld_fill(W.builder_preset(key, snapshot().mode))

    @render.ui
    def bld_cite():
        key = input.bld_start()
        if not key:
            return None
        return ui.tooltip(ui.span("source ⓘ", class_="enz-chip"),
                          W.BUILDER_PRESETS[key].citation, placement="left")  # fmt: skip

    @reactive.effect
    @reactive.event(input.bld_add)
    def _bld_add():
        n = bld_n()
        req(n < W.MAX_DIGESTS)
        chem = bld_default_chem()
        used = {(input[f"bld_enz{i}"](), input[f"bld_chem{i}"]()) for i in range(n)}
        enz = W.builder_enzymes(available_enzymes(snapshot().mode))
        pick = next((e for e in enz if (e, chem) not in used), enz[0])
        ui.update_select(f"bld_enz{n}", selected=pick)
        ui.update_select(f"bld_chem{n}", selected=chem)
        ui.update_numeric("bld_n", value=n + 1)

    def _bld_remover(i: int):
        @reactive.effect
        @reactive.event(input[f"bld_rm{i}"])
        def _remove():
            n = bld_n()
            req(n > 1 and i < n)
            enz = [input[f"bld_enz{j}"]() for j in range(n)]
            chem = [input[f"bld_chem{j}"]() for j in range(n)]
            del enz[i], chem[i]
            for j in range(i, n - 1):
                ui.update_select(f"bld_enz{j}", selected=enz[j])
                ui.update_select(f"bld_chem{j}", selected=chem[j])
            ui.update_numeric("bld_n", value=n - 1)

        return _remove

    def _bld_noter(i: int):
        def _note():
            req(i < bld_n())
            wf = bld_workflow()
            rule = wf.digests[i].enzyme
            notes = [ui.div(text, class_="bld-note blocking" if block else "bld-note")
                     for j, text, block in wf.notes() if j == i]  # fmt: skip
            info = ui.tooltip(ui.span(f"About {rule.name} ⓘ", class_="enz-chip"),
                              ui.div(*enzyme_tip(rule), class_="enz-tip"), placement="left")
            return ui.div(*notes, ui.div(info, class_="mb-1"))

        return _note

    for _i in range(W.MAX_DIGESTS):
        _bld_remover(_i)
        output(id=f"bld_note{_i}")(render.ui(_bld_noter(_i)))

    @reactive.effect
    def _bld_block():
        ok = bld_workflow().can_compare
        ui.update_action_button("bld_compare", disabled=not ok)
        ui.update_action_button("bld_save", disabled=not ok)

    @reactive.effect
    def _bld_against_choices():
        names = list(saved())
        with reactive.isolate():
            cur = input.bld_against()
        ui.update_select("bld_against", choices={**STANDARD_CHOICE, **{x: x for x in names}},
                         selected=cur if cur in names else "standard")  # fmt: skip

    @reactive.effect
    @reactive.event(input.bld_compare)
    def _bld_compare():
        wf = bld_workflow()
        if not wf.can_compare:
            ui.notification_show(W.NOTE_DUPLICATE, type="error")
            return
        key = input.bld_against()
        against = saved().get(key) if key and key != "standard" else None
        if against is not None and against.name == wf.name and against != wf:
            wf = replace(wf, name=f"{wf.name} (edited)")
        bld_shown.set((wf, against))

    @reactive.effect
    @reactive.event(input.bld_save)
    def _bld_save():
        wf = bld_workflow()
        if not wf.can_compare:
            ui.notification_show(W.NOTE_DUPLICATE, type="error")
            return
        new = dict(saved())
        new[wf.name] = wf
        saved.set(new)
        ui.notification_show(f"Saved “{wf.name}”: it is now ranked with the other strategies "
                             "and can be chosen under “Compare against”.", duration=6)
        ui.update_text("bld_name", value=next_strategy_name(new))

    @render.ui
    def bld_saved():
        names = list(saved())
        if not names:
            return None
        return ui.div("Saved this session: " + ", ".join(names), class_="bld-saved")

    def _bld_filename() -> str:
        try:
            name = bld_workflow().name
        except Exception:  # noqa: BLE001 - any failure falls back to a generic name
            name = "strategy"
        return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") + ".json"

    @render.download(filename=_bld_filename)
    def bld_download():
        yield bld_workflow().to_json()

    @reactive.effect
    @reactive.event(input.bld_load)
    def _bld_load():
        files = input.bld_load()
        req(files)
        try:
            text = Path(files[0]["datapath"]).read_text()
            wf = W.Workflow.from_json(text, snapshot().mode)
        except (OSError, ValueError, KeyError, TypeError) as e:
            ui.notification_show(f"Could not load this strategy file: {e}", type="error")
            return
        bld_fill(wf)

    # -- advanced table
    def sp_adv_frame() -> pd.DataFrame:
        goal, n, mk = sp_goal()
        return sp_res().advanced_table(n, mk)

    @render.data_frame
    def sp_adv():
        return render.DataGrid(sp_adv_frame(), width="100%", filters=True)

    @render.download(filename="proteasescope_strategy_all.tsv")
    def dl_sp_adv_tsv():
        goal, n, mk = sp_goal()
        res = sp_res()
        yield R.tsv_bytes(sp_adv_frame(),
                          res.header_lines("Strategy picker: all strategies", n, mk))  # fmt: skip

    @render.download(filename="proteasescope_strategy_all.xlsx")
    def dl_sp_adv_xlsx():
        goal, n, mk = sp_goal()
        res = sp_res()
        yield R.excel_bytes({"strategies": sp_adv_frame()},
                            res.header_lines("Strategy picker: all strategies", n, mk))  # fmt: skip

    # ----- Protein viewer -------------------------------------------------------------------
    @reactive.calc
    async def viewer():
        res, sel = await results()
        acc = input.protein()
        req(acc and acc in set(res.proteins.accessions))
        return acc, res.protein_analyses(acc, sel)

    @render_plotly
    async def map_fig():
        acc, an = await viewer()
        pm = M.peptide_map(an, acc)
        res, _ = await results()
        meta = res.proteins.meta.set_index("accession").loc[acc]
        title = f"{meta['gene'] or ''} {acc} — {meta['protein_name'] or ''}".strip()
        return fig_widget(F.peptide_map(pm, [a.label for a in an], title), f"peptide_map_{acc}")

    async def viewer_peps_frame():
        _, an = await viewer()
        return round_df(pd.concat([a.peptide_table() for a in an], ignore_index=True), 4)

    async def viewer_sites_frame():
        _, an = await viewer()
        frames = []
        for a in an[1:]:
            s = M.compare(an[0], a).sites
            s.insert(0, "digest", a.label)
            frames.append(s)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    @render.data_frame
    async def viewer_peps():
        return render.DataGrid(await viewer_peps_frame(), width="100%", filters=True)

    @render.data_frame
    async def viewer_sites():
        return render.DataGrid(await viewer_sites_frame(), width="100%", filters=True)

    register_download("viewer_peps", viewer_peps_frame, "Protein viewer: peptides")
    register_download("viewer_sites", viewer_sites_frame, "Protein viewer: sites gained/lost")


app = App(app_ui, server)
_ = (np, __version__)  # used in templates / kept for introspection
