"""Plotly figure builders for the app (pure functions: DataFrames in, figures out).

Colours follow the validated data-viz palette: one series colour (blue) for single-series
charts; gained/lost use the diverging pair blue/red (green/red fails colour-blind checks);
co-occurrence uses grey (baseline) / blue (new) / orange (split) as CLAUDE.md specifies.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
NEUTRAL = "#f0efec"
BLUE = "#2a78d6"  # series 1, "gained", "new"
ORANGE = "#eb6834"  # "split"
RED = "#d03b3b"  # "lost"
GREY = "#b9b8b0"  # baseline / seen in both
FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif"


def _base(fig: go.Figure, height: int, title: str = "") -> go.Figure:
    fig.update_layout(
        template="plotly_white",
        height=height,
        title={"text": title, "font": {"size": 15, "color": INK}, "x": 0.01, "xanchor": "left",
               "y": 0.995, "yref": "container", "yanchor": "top"},
        font={"family": FONT, "size": 12, "color": INK2},
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        margin={"l": 10, "r": 10, "t": 84 if title else 40, "b": 10},
        hoverlabel={"bgcolor": "white", "font": {"family": FONT, "color": INK}},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.0, "x": 0, "xanchor": "left"},
        barcornerradius=4,
    )
    fig.update_xaxes(gridcolor=GRID, linecolor=AXIS, zerolinecolor=AXIS, automargin=True)
    fig.update_yaxes(gridcolor=GRID, linecolor=AXIS, zerolinecolor=AXIS, automargin=True)
    return fig


def empty(message: str, height: int = 200) -> go.Figure:
    """A blank figure carrying a message."""
    fig = go.Figure()
    fig.add_annotation(text=message, showarrow=False, font={"color": MUTED, "size": 13})
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return _base(fig, height)


def _gaps(z: np.ndarray) -> list:
    """2-D array -> nested lists with None for NaN (widgets reject NaN in JSON)."""
    return [[None if (isinstance(v, float) and np.isnan(v)) else v for v in row]
            for row in np.asarray(z, dtype=object).tolist()]


def _wrap(label: str, width: int = 34) -> str:
    out, line = [], ""
    for word in label.replace(" | ", " | ").split(" "):
        if len(line) + len(word) + 1 > width and line:
            out.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    out.append(line)
    return "<br>".join(out)


# ----- Overview --------------------------------------------------------------------------
#: Gained/lost panels: column prefix, panel title, unit for hover text.
CHANGE_PANELS = [("peptides", "Identifiable peptides", "peptides"),
                 ("K_sites", "K sites", "K sites"),
                 ("STY_sites", "S/T/Y sites", "S/T/Y sites")]  # fmt: skip


def _signed(v: float, digits: int = 0) -> str:
    """'+1,234' / '−56' / '0' (true minus sign)."""
    if round(v, digits) == 0:
        return "0"
    txt = f"{abs(v):,.{digits}f}"
    return ("+" if v > 0 else "−") + txt


def _fit_range(lost: np.ndarray, gained: np.ndarray, lost_txt: list[str],
               gained_txt: list[str], width_px: float = 300.0) -> list[float]:
    """x range [-left, right] so every value label fits beside its bar (outside, on its side)
    in a panel about ``width_px`` wide (11 px digits are ~6.6 px; 10 px padding)."""

    def w(t: str) -> float:
        return 6.6 * len(t) + 10

    lost = np.abs(np.asarray(lost, dtype=float))
    gained = np.asarray(gained, dtype=float)
    span = max(float(lost.max(initial=0) + gained.max(initial=0)), 1e-9)
    left = right = 0.0
    for _ in range(60):
        left = max((v + w(t) * span / width_px for v, t in zip(lost, lost_txt, strict=True)),
                   default=0.0)  # fmt: skip
        right = max((v + w(t) * span / width_px for v, t in zip(gained, gained_txt, strict=True)),
                    default=0.0)  # fmt: skip
        new = left + right
        if abs(new - span) <= 1e-6 * span:
            break
        span = new
    return [-left, right]


def diverging_changes(counts: pd.DataFrame, base: dict[str, float]) -> go.Figure:
    """Gained (right, blue) and lost (left, red) per alternative digest, with exact counts on
    every bar, plus the change in median coverage; a 2 x 2 grid of panels.

    ``counts`` has one row per alternative digest: ``digest`` (short row label), optional
    ``full_label`` (shown on hover; defaults to ``digest``), ``{panel}_gained`` and
    ``{panel}_lost`` for ``peptides``, ``K_sites`` and ``STY_sites``, and
    ``median_coverage_pct`` (the alternative's median per-protein coverage, %). ``base`` has
    the baseline totals ``peptides``, ``K_sites``, ``STY_sites`` (identifiable peptides and
    sites on identifiable peptides) and ``median_coverage_pct``. Coverage change is in
    percentage points.
    """
    if counts.empty:
        return empty("Select at least one alternative digest.")
    sub = f"<span style='font-size:11px;color:{MUTED}'>(baseline {{}})</span>"
    titles = [f"<b>{t}</b> {sub.format(f'{base[k]:,.0f}')}" for k, t, _ in CHANGE_PANELS]
    cov0 = f"{base['median_coverage_pct']:.1f}%"
    titles.append(f"<b>Median coverage</b> {sub.format(cov0)}")
    fig = make_subplots(rows=2, cols=2, shared_yaxes=True, horizontal_spacing=0.06,
                        vertical_spacing=0.12, subplot_titles=titles)  # fmt: skip
    cells = [(1, 1), (1, 2), (2, 1), (2, 2)]
    y = [_wrap(d) for d in counts["digest"]]
    full = counts["full_label"] if "full_label" in counts else counts["digest"]
    for (key, _, unit), (r, c) in zip(CHANGE_PANELS, cells, strict=False):
        g = counts[f"{key}_gained"].to_numpy(dtype=float)
        lo = counts[f"{key}_lost"].to_numpy(dtype=float)
        net = g - lo
        b = float(base[key])
        pct = [f" ({_signed(100 * n / b, 1)}% of baseline)" if b else "" for n in net]
        hover = [f"<b>{d}</b><br>Gained: {_signed(gi)} {unit}<br>Lost: {_signed(-li)}<br>"
                 f"Net: {_signed(n)}{p}<br>Total: {b:,.0f} → {b + n:,.0f}"
                 for d, gi, li, n, p in zip(full, g, lo, net, pct, strict=True)]
        # zero bars carry no label, except "no change" when nothing is gained or lost
        texts = {"lost": [_signed(-v) if v else "" for v in lo],
                 "gained": [_signed(v) if v else ("" if li else "no change")
                            for v, li in zip(g, lo, strict=True)]}  # fmt: skip
        for status, vals, color in (("lost", -lo, RED), ("gained", g, BLUE)):
            fig.add_trace(
                go.Bar(
                    y=y, x=vals, orientation="h", marker_color=color,
                    name=f"{status.capitalize()} vs baseline", legendgroup=status,
                    showlegend=(r, c) == (1, 1), text=texts[status],
                    textposition="outside", cliponaxis=False,
                    textfont={"size": 11, "color": [INK2 if v else MUTED for v in vals]},
                    hovertext=hover, hovertemplate="%{hovertext}<extra></extra>",
                ),
                row=r, col=c,
            )  # fmt: skip
        fig.update_xaxes(range=_fit_range(lo, g, texts["lost"], texts["gained"]), row=r, col=c)
    cov = counts["median_coverage_pct"].to_numpy(dtype=float)
    d_cov = np.round(cov - float(base["median_coverage_pct"]), 1)
    cov_txt = [f"{_signed(v, 1)} pts" if v else "no change" for v in d_cov]
    fig.add_trace(
        go.Bar(
            y=y, x=d_cov, orientation="h", showlegend=False,
            marker_color=[BLUE if v > 0 else RED if v < 0 else GREY for v in d_cov],
            text=cov_txt, textposition="outside", cliponaxis=False,
            textfont={"size": 11, "color": [INK2 if v else MUTED for v in d_cov]},
            hovertext=[f"<b>{d}</b><br>Median coverage {float(base['median_coverage_pct']):.1f}%"
                       f" → {cv:.1f}%<br>Change: {_signed(v, 1)} percentage points"
                       for d, cv, v in zip(full, cov, d_cov, strict=True)],
            hovertemplate="%{hovertext}<extra></extra>",
        ),
        row=2, col=2,
    )  # fmt: skip
    neg = np.where(d_cov < 0, -d_cov, 0.0)
    pos = np.where(d_cov > 0, d_cov, 0.0)
    # labels sit on the side the bar points to; zero bars put "0 pts" on the right
    fig.update_xaxes(range=_fit_range(neg, np.maximum(pos, 0.5 * (d_cov == 0)),
                                      [t if v < 0 else "" for t, v in zip(cov_txt, d_cov,
                                                                           strict=True)],
                                      [t if v >= 0 else "" for t, v in zip(cov_txt, d_cov,
                                                                            strict=True)]),
                     row=2, col=2)  # fmt: skip
    fig.update_yaxes(autorange="reversed")
    fig.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor=AXIS, showticklabels=False,
                     showgrid=False)  # fmt: skip
    for ann in fig.layout.annotations:
        ann.font = {"size": 13, "color": INK}
    fig = _base(fig, 150 + 2 * 30 * len(counts))
    fig.update_layout(barmode="relative", bargap=0.3, margin={"t": 40, "b": 40},
                      legend={"orientation": "h", "y": -0.02, "yanchor": "top", "x": 0})
    return fig


def coverage_scatter(df: pd.DataFrame, base_label: str, alt_label: str) -> go.Figure:
    """Per-protein coverage, baseline (x) vs alternative (y). ``customdata[0]`` = accession."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, 100], y=[0, 100], mode="lines", hoverinfo="skip",
                             line={"color": AXIS, "width": 1}, showlegend=False))  # fmt: skip
    fig.add_trace(
        go.Scattergl(
            x=df["coverage_pct_baseline"], y=df["coverage_pct_alternative"], mode="markers",
            marker={"size": 6, "color": BLUE, "opacity": 0.45,
                    "line": {"width": 0}},
            customdata=np.stack([df["accession"], df["gene"].fillna(""),
                                 df["protein_name"].fillna("")], axis=1),
            hovertemplate="%{customdata[1]} %{customdata[0]}<br>%{customdata[2]}<br>"
            "baseline %{x:.1f}% → alternative %{y:.1f}%<extra>click to open</extra>",
            showlegend=False,
        )
    )  # fmt: skip
    fig.update_xaxes(title=_wrap(f"Coverage, {base_label} (%)", 60), range=[-2, 102])
    fig.update_yaxes(title=_wrap(f"Coverage, {alt_label} (%)", 50), range=[-2, 102])
    return _base(fig, 460, "Per-protein coverage: baseline vs alternative")


def coverage_distributions(values: dict[str, np.ndarray]) -> go.Figure:
    """Per-protein coverage distribution (box) per digest."""
    fig = go.Figure()
    for label, v in values.items():
        fig.add_trace(go.Box(x=v, name=_wrap(label), marker_color=BLUE, boxpoints=False,
                             line={"width": 1.5}, fillcolor="rgba(42,120,214,0.10)",
                             showlegend=False, hoverinfo="x"))  # fmt: skip
    fig.update_xaxes(title="Coverage per protein (%)", range=[-2, 102])
    fig.update_yaxes(autorange="reversed")
    return _base(fig, 110 + 46 * len(values), "Coverage distribution per digest")


# ----- Digest stats --------------------------------------------------------------------------
def length_histogram(df: pd.DataFrame) -> go.Figure:
    """Cleavage-product lengths with the identifiable window shaded."""
    fig = go.Figure(go.Bar(
        x=df["length_value"], y=df["n_peptides"], marker_color=BLUE,
        customdata=df["length"], hovertemplate="%{customdata} aa: %{y:,}<extra></extra>",
    ))  # fmt: skip
    win = df.loc[df["in_window"], "length_value"]
    if len(win):
        fig.add_vrect(x0=win.min() - 0.5, x1=win.max() + 0.5, fillcolor=BLUE, opacity=0.08,
                      line_width=0, annotation_text="identifiable window",
                      annotation_position="top left",
                      annotation_font={"color": INK2, "size": 11})  # fmt: skip
    last = int(df["length_value"].max())
    ticks = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, last]
    fig.update_xaxes(title="Peptide length (aa)", tickvals=ticks,
                     ticktext=[str(t) for t in ticks[:-1]] + [f">{last - 1}"])  # fmt: skip
    fig.update_yaxes(title="Cleavage products (0 missed cleavages)")
    fig.update_layout(bargap=0.15)
    return _base(fig, 340, "Cleaved peptide frequency")


def coverage_histogram(coverage_pct: np.ndarray) -> go.Figure:
    """Distribution of per-protein coverage by identifiable peptides."""
    fig = go.Figure(go.Histogram(x=coverage_pct, xbins={"start": 0, "end": 100.001, "size": 2},
                                 marker_color=BLUE,
                                 hovertemplate="%{x}%: %{y:,} proteins<extra></extra>"))  # fmt: skip
    fig.update_xaxes(title="Coverage by identifiable peptides (%)")
    fig.update_yaxes(title="Proteins")
    fig.update_layout(bargap=0.1)
    return _base(fig, 320, "Cleaved peptide coverage per protein")


def aa_detectability(df: pd.DataFrame) -> go.Figure:
    """Fraction of each residue type on identifiable peptides; K, S, T, Y, C highlighted."""
    d = df.sort_values("residue")
    colors = [BLUE if h else AXIS for h in d["highlight"]]
    fig = go.Figure(go.Bar(
        x=d["residue"], y=100 * d["fraction_detectable"], marker_color=colors,
        customdata=np.stack([d["on_identifiable_peptides"], d["occurrences"]], axis=1),
        hovertemplate="%{x}: %{y:.1f}% (%{customdata[0]:,} of %{customdata[1]:,})<extra></extra>",
    ))  # fmt: skip
    fig.update_yaxes(title="Occurrences on identifiable peptides (%)", range=[0, 100])
    fig.update_xaxes(title="Residue (K, S, T, Y and C highlighted)")
    fig.update_layout(bargap=0.35)
    return _base(fig, 320, "Amino acid detectability")


# ----- Peptide map ---------------------------------------------------------------------------
def peptide_map(pm: dict[str, pd.DataFrame], labels: list[str], title: str,
                reference: str = "baseline") -> go.Figure:
    """One row per digest: peptide bars (lanes for overlaps), K ticks, gained/lost strips,
    coverage % at the right. ``pm`` is the output of ``metrics.peptide_map`` (mature numbering).
    """
    peps, lys, strips, cov = pm["peptides"], pm["lysines"], pm["strips"], pm["coverage"]
    lane_h, gap = 1.0, 1.4
    rows, y0 = [], 0.0
    fig = go.Figure()
    xs, ys, hov = [], [], []
    for label in labels:
        p = peps[peps["digest"] == label].sort_values(["start_mature", "end_mature"]) \
            if len(peps) else peps  # fmt: skip
        lane_end: list[int] = []
        lanes = []
        for s, e in zip(p["start_mature"], p["end_mature"], strict=True) if len(p) else []:
            for k, end in enumerate(lane_end):
                if end < s:
                    lane_end[k] = e
                    lanes.append(k)
                    break
            else:
                lane_end.append(e)
                lanes.append(len(lane_end) - 1)
        n_lanes = max(1, len(lane_end))
        top = y0
        rows.append((label, top, top + n_lanes * lane_h))
        if len(p):
            y = top + (np.asarray(lanes) + 0.5) * lane_h
            fig.add_trace(go.Bar(
                y=y, x=(p["end_mature"] - p["start_mature"] + 1) - 0.3,
                base=p["start_mature"] - 0.35,
                orientation="h", width=0.7 * lane_h, marker_color=BLUE, showlegend=False,
                customdata=np.stack([p["peptide"], p["start_mature"], p["end_mature"],
                                     p["length"], p["n_K"]], axis=1),
                hovertemplate="%{customdata[0]}<br>%{customdata[1]}–%{customdata[2]} "
                "(%{customdata[3]} aa), %{customdata[4]} K<extra>" + label + "</extra>",
            ))  # fmt: skip
        s_rows = strips[strips["digest"] == label] if len(strips) else strips
        for st in s_rows.itertuples() if len(s_rows) else []:
            xs += [st.start_mature - 0.5, st.end_mature + 0.5, None]
            ys += [top - 0.35] * 2 + [None]
            hov += [st.status] * 3
        y0 = top + n_lanes * lane_h + gap
    for status, color in (("gained", BLUE), ("lost", RED)):
        sx = [x for x, h in zip(xs, hov, strict=True) if h == status]
        sy = [y for y, h in zip(ys, hov, strict=True) if h == status]
        if sx:
            fig.add_trace(go.Scatter(x=sx, y=sy, mode="lines", line={"color": color, "width": 5},
                                     name=f"Residues {status} vs {reference}",
                                     hoverinfo="name"))  # fmt: skip
    if len(lys):
        kx, ky = [], []
        for k in lys["position_mature"]:
            kx += [k, k, None]
            ky += [-0.8, y0, None]
        fig.add_trace(go.Scatter(x=kx, y=ky, mode="lines", line={"color": AXIS, "width": 1},
                                 name="Lysine", hoverinfo="skip"))  # fmt: skip
        fig.add_trace(go.Scatter(x=lys["position_mature"], y=[-0.8] * len(lys), mode="markers",
                                 marker={"symbol": "triangle-down", "size": 8, "color": INK2},
                                 showlegend=False,
                                 hovertemplate="K%{x}<extra></extra>"))  # fmt: skip
    covmap = dict(zip(cov["digest"], cov["coverage_window_pct"], strict=True)) if len(cov) else {}
    for label, top, bottom in rows:
        fig.add_annotation(x=1.0, xref="paper", xanchor="left", y=(top + bottom) / 2,
                           text=f"{covmap.get(label, 0):.0f}%", showarrow=False,
                           font={"color": INK2, "size": 11})  # fmt: skip
    fig.update_yaxes(tickvals=[(t + b) / 2 for _, t, b in rows],
                     ticktext=[_wrap(lbl, 30) for lbl, _, _ in rows], autorange="reversed",
                     showgrid=False, zeroline=False)  # fmt: skip
    fig.update_xaxes(title="Residue (mature numbering)", showgrid=False)
    fig = _base(fig, int(110 + 18 * y0), title)
    fig.update_layout(bargap=0, margin={"r": 52})
    return fig


# ----- Histone figures ---------------------------------------------------------------------
STATUS_COLORS = {
    "covered": GREY, "seen in both": GREY, "gained": BLUE, "lost": RED,
    "not covered": NEUTRAL, "unseen in both": NEUTRAL,
}  # fmt: skip


def mark_grid(vis: pd.DataFrame, letters: dict[str, str]) -> go.Figure:
    """Known marks x digests: status per cell (columns lettered; legend in ``letters``)."""
    if vis.empty:
        return empty("No known marks for the selected entries.")
    v = vis.copy()
    base = v["status"] == "baseline"
    v.loc[base, "status"] = np.where(v.loc[base, "covered"], "covered", "not covered")
    v["row"] = v["display_name"] + " " + v["site"]
    rows = list(dict.fromkeys(v["row"]))
    cols = list(dict.fromkeys(v["digest"]))
    order = list(STATUS_COLORS)
    code = {s: i for i, s in enumerate(order)}
    z = np.full((len(rows), len(cols)), np.nan)
    text = np.full((len(rows), len(cols)), "", dtype=object)
    ri = {r: i for i, r in enumerate(rows)}
    ci = {c: i for i, c in enumerate(cols)}
    for r in v.itertuples():
        z[ri[r.row], ci[r.digest]] = code[r.status]
        text[ri[r.row], ci[r.digest]] = (f"{r.row} ({r.modifications})<br>{r.digest}<br>"
                                          f"{r.status}<br>{r.covering_peptides or '—'}")  # fmt: skip
    n = len(order)
    scale = []
    for i, s in enumerate(order):
        scale += [[i / n, STATUS_COLORS[s]], [(i + 1) / n, STATUS_COLORS[s]]]
    fig = go.Figure(go.Heatmap(z=_gaps(z), x=[letters.get(c, c) for c in cols], y=rows, zmin=0,
                               zmax=n, colorscale=scale, showscale=False, xgap=2, ygap=2,
                               text=text.tolist(),
                               hovertemplate="%{text}<extra></extra>"))  # fmt: skip
    for s in ("covered", "gained", "lost", "not covered"):
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=s.capitalize() if s
                                 != "covered" else "Covered / seen in both",
                                 marker={"symbol": "square", "size": 12,
                                         "color": STATUS_COLORS[s]}))  # fmt: skip
    fig.update_yaxes(autorange="reversed", showgrid=False)
    fig.update_xaxes(side="top", showgrid=False, title="Digest (see legend table)")
    return _base(fig, 140 + 16 * len(rows), "Known-mark visibility")


def isoform_matrix(mat: pd.DataFrame, title: str) -> go.Figure:
    """Entry x entry: number of the row entry's peptides that do not map to the column entry."""
    z = mat.astype("float").to_numpy()
    fig = go.Figure(go.Heatmap(
        z=_gaps(z), x=list(mat.columns), y=list(mat.index),
        colorscale=[[0, "#cde2fb"], [1, "#104281"]],
        text=np.where(np.isnan(z), "", np.nan_to_num(z).astype(int).astype(str)).tolist(),
        texttemplate="%{text}",
        xgap=2, ygap=2, colorbar={"title": "peptides", "thickness": 10},
        hovertemplate="%{y} peptides absent from %{x}: %{z}<extra></extra>",
    ))  # fmt: skip
    fig.update_yaxes(autorange="reversed", showgrid=False)
    fig.update_xaxes(showgrid=False, tickangle=-35)
    return _base(fig, 160 + 30 * len(mat), title)


def cooccurrence_matrix(m: pd.DataFrame, lysines: list[int], title: str) -> go.Figure:
    """Upper-triangle K x K: grey = baseline pair, blue = new, orange = split; letters name the
    alternative digests responsible; the diagonal shows the lysine number."""
    n = len(lysines)
    if n == 0:
        return empty("This entry has no lysines.")
    idx = {k: i for i, k in enumerate(lysines)}
    z = np.full((n, n), np.nan)
    text = np.full((n, n), "", dtype=object)
    hover = np.full((n, n), "", dtype=object)
    for r in m.itertuples():
        i, j = idx[r.K_i], idx[r.K_j]
        state = 3 if r.split_in else 2 if r.new_in else 1 if r.baseline_pair else 0
        z[i, j] = state
        text[i, j] = (r.split_in or r.new_in)
        hover[i, j] = (f"K{r.K_i}–K{r.K_j}<br>baseline pair: {'yes' if r.baseline_pair else 'no'}"
                       f"<br>new in: {r.new_in or '—'}<br>split in: {r.split_in or '—'}"
                       f"<br>lost in: {r.lost_in or '—'}")  # fmt: skip
    for k, i in idx.items():
        text[i, i] = f"K{k}"
        hover[i, i] = f"K{k}"
    scale = [[0, SURFACE], [0.25, SURFACE], [0.25, GREY], [0.5, GREY], [0.5, BLUE], [0.75, BLUE],
             [0.75, ORANGE], [1, ORANGE]]  # fmt: skip
    labels = [f"K{k}" for k in lysines]
    fig = go.Figure(go.Heatmap(z=_gaps(z), x=labels, y=labels, zmin=0, zmax=4, colorscale=scale,
                               showscale=False, xgap=2, ygap=2, text=text.tolist(),
                               texttemplate="%{text}", textfont={"size": 10, "color": INK},
                               customdata=hover.tolist(),
                               hovertemplate="%{customdata}<extra></extra>"))  # fmt: skip
    for name, color in (("Baseline pair", GREY), ("New pair", BLUE), ("Split pair", ORANGE)):
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=name,
                                 marker={"symbol": "square", "size": 12, "color": color}))
    fig.update_yaxes(autorange="reversed", showgrid=False, showticklabels=False)
    fig.update_xaxes(showgrid=False, showticklabels=False)
    side = max(360, 26 * n + 120)
    return _base(fig, side, title)


# ----- Strategy picker ------------------------------------------------------------------------
DARK_GREY = "#6b6a65"  # "seen, as with the standard digest"
WHITE = "#ffffff"
LYSINE_STATES = {
    "standard": (DARK_GREY, "as with the standard digest"),
    "gained": (BLUE, "new with this strategy"),
    "lost": (RED, "lost compared with the standard digest"),
    "neither": (WHITE, "not achieved by either"),
    "not counted": (NEUTRAL, "no known modification (not counted)"),
}


def lysine_picture(states: pd.DataFrame, achieved: str, title: str = "",
                   reference: str = "the standard digest") -> go.Figure:  # fmt: skip
    """Lysine picture: one row per digest, one column per lysine in scope; cell colour =
    state (see :data:`LYSINE_STATES`). Known-mark lysines get a dot above the column.

    ``achieved`` names what the goal counts, e.g. "seen" or "measurable on its own";
    ``reference`` names the compared-against digest in hover and legend text.
    """
    if states.empty:
        return empty("There are no lysines here.")
    rows = list(dict.fromkeys(states["row"]))
    cols = list(dict.fromkeys(states["lysine"]))
    order = list(LYSINE_STATES)
    code = {s: i for i, s in enumerate(order)}
    z = np.full((len(rows), len(cols)), np.nan)
    hover = np.full((len(rows), len(cols)), "", dtype=object)
    ri, ci = {r: i for i, r in enumerate(rows)}, {c: i for i, c in enumerate(cols)}
    for r in states.itertuples():
        z[ri[r.row], ci[r.lysine]] = code[r.state]
        what = {"standard": f"{achieved}, as with {reference}",
                "gained": f"newly {achieved}", "lost": f"no longer {achieved}",
                "neither": f"not {achieved}",
                "not counted": "no known modification of this type (not counted)"}[r.state]
        hover[ri[r.row], ci[r.lysine]] = f"{r.lysine}: {what}<br>{r.row}"
    n = len(order)
    scale = []
    for i, s in enumerate(order):
        scale += [[i / n, LYSINE_STATES[s][0]], [(i + 1) / n, LYSINE_STATES[s][0]]]
    fig = go.Figure(go.Heatmap(z=_gaps(z), x=cols, y=list(range(len(rows))), zmin=0,
                               zmax=n, colorscale=scale, showscale=False, xgap=3, ygap=3,
                               customdata=hover.tolist(),
                               hovertemplate="%{customdata}<extra></extra>"))  # fmt: skip
    marked = states.drop_duplicates("lysine")
    marked = marked[marked["marks"] != ""]
    if len(marked):
        fig.add_trace(go.Scatter(
            x=marked["lysine"], y=[-0.75] * len(marked), mode="markers", name="Known modifications",
            marker={"size": 7, "color": INK2}, customdata=marked["marks"],
            hovertemplate="%{x}: %{customdata}<extra>known modifications</extra>"))  # fmt: skip
    present = set(states["state"])
    labels = {"standard": f"{achieved.capitalize()}, as with {reference}",
              "gained": f"Newly {achieved}", "lost": "Lost", "neither": f"Not {achieved}",
              "not counted": "No known modification"}  # fmt: skip
    for s in order:
        if s in present:
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", name=labels[s],
                                     marker={"symbol": "square", "size": 12,
                                             "color": LYSINE_STATES[s][0],
                                             "line": {"color": AXIS, "width": 1}}))  # fmt: skip
    fig.update_yaxes(showgrid=False, zeroline=False, range=[len(rows) - 0.5, -1.1],
                     tickvals=list(range(len(rows))), ticktext=[_wrap(r, 40) for r in rows])
    fig.update_xaxes(side="bottom", showgrid=False, type="category")
    fig = _base(fig, 190 + 40 * len(rows), title)
    fig.update_layout(plot_bgcolor=GRID, margin={"b": 70},
                      legend={"yref": "container", "y": 0.01, "yanchor": "bottom"})  # fmt: skip
    return fig


ARC_STYLES = {
    "together": (GREY, "solid", "Read together (also by the standard digest)"),
    "new together": (BLUE, "solid", "Newly read together"),
    "lost together": (RED, "dot", "No longer read together"),
    "new separable": (ORANGE, "dash", "Newly separable"),
    "lost separable": (RED, "dash", "No longer separable"),
}


def _arc(i: int, j: int, below: bool, n_pts: int = 24) -> tuple[list, list]:
    t = np.linspace(0, np.pi, n_pts)
    h = 0.35 + 0.18 * (j - i)
    xs = (i + j) / 2 - (j - i) / 2 * np.cos(t)
    ys = h * np.sin(t) * (-1 if below else 1)
    return xs.tolist() + [None], ys.tolist() + [None]


def arc_view(pairs: pd.DataFrame, lysines: list[str], strategy_label: str | None,
             title: str = "", reference: str = "Standard digest") -> go.Figure:  # fmt: skip
    """Arc view for pair goals: lysines on a line; arcs above join pairs read together
    (grey/blue/red dotted); dashed arcs below join pairs newly separable (orange) or no
    longer separable (red). One panel for the standard digest, one for the strategy."""
    n = len(lysines)
    if n < 2:
        return empty("Fewer than two lysines here: no pairs to show.")
    panels = [("standard", reference)]
    if strategy_label:
        panels.append(("strategy", strategy_label))
    fig = make_subplots(rows=len(panels), cols=1, shared_xaxes=True, vertical_spacing=0.12,
                        subplot_titles=[_wrap(lab, 70) for _, lab in panels])  # fmt: skip
    shown = set()
    top = 0.35 + 0.18 * (n - 1)
    for row, (panel, _) in enumerate(panels, start=1):
        p = pairs[pairs["panel"] == panel]
        for cat, (color, dash, name) in ARC_STYLES.items():
            if reference != "Standard digest":
                name = name.replace("the standard digest", reference)
            sub = p[p["category"] == cat]
            if sub.empty:
                continue
            xs, ys, hov = [], [], []
            for r in sub.itertuples():
                x, y = _arc(r.i, r.j, below="separable" in cat)
                xs += x
                ys += y
                hov += [r.pair] * (len(x) - 1) + [None]
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name=name, legendgroup=cat,
                                     showlegend=cat not in shown,
                                     line={"color": color, "dash": dash,
                                           "width": 1.4 if cat == "together" else 2},
                                     text=hov, hovertemplate="%{text}<extra>" + name
                                     + "</extra>"), row=row, col=1)  # fmt: skip
            shown.add(cat)
        fig.add_trace(go.Scatter(x=[-0.5, n - 0.5], y=[0, 0], mode="lines",
                                 line={"color": AXIS, "width": 1}, showlegend=False,
                                 hoverinfo="skip"), row=row, col=1)  # fmt: skip
        for x, name in enumerate(lysines):
            fig.add_annotation(x=x, y=0, text=name, showarrow=False, row=row, col=1,
                               font={"size": 10, "color": INK}, bgcolor=SURFACE,
                               bordercolor=INK2, borderwidth=1, borderpad=2)  # fmt: skip
        has_below = (p["category"].str.contains("separable")).any()
        fig.update_yaxes(visible=False, range=[-(top + 0.4) if has_below else -0.6, top + 0.2],
                         row=row, col=1)  # fmt: skip
    fig.update_xaxes(visible=False, range=[-0.7, n - 0.3])
    fig = _base(fig, 160 + 190 * len(panels), title)
    fig.update_layout(legend={"y": -0.02, "yanchor": "top"}, margin={"b": 40})
    for a in fig.layout.annotations[: len(panels)]:  # subplot titles
        a.update(font={"size": 12, "color": INK2}, x=0, xanchor="left")
    return fig
