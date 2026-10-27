"""PEC-vs-MEC box-and-whisker figure: chemicals on the y axis, concentration (log) on the x axis, one PEC box (turquoise) and one MEC box (white) per chemical,
every field plotted as a dot.

box = interquartile range, line = median, whiskers = minimum and maximum over the fields, dots = the individual field values; all statistics are taken on log10
concentrations (the axis is logarithmic).  PEC dots are the central estimate per field; the thin line through each PEC dot is the property range (PEC_low to PEC_high),
which is a point when the chemical's properties were measured (clipped at the axis edge).  PEC colour: solid = measured properties, light + hatched = estimated properties.
MEC dots: filled = MS2 confirmed, open = MS1 only / weak or conflicting MS2.  Chemicals run top to bottom in order of PEC/MEC (largest first)."""
from __future__ import annotations
import numpy as np
import pandas as pd

PEC_SOLID, PEC_LIGHT, PEC_EDGE = "#1FB5A6", "#8FE3D9", "#0B6E66"
INK, GRID = "#222222", "#DDDDDD"


def _stats(v: np.ndarray):
    lv = np.log10(v); q1, med, q3 = np.percentile(lv, [25, 50, 75])
    return dict(med=10 ** med, q1=10 ** q1, q3=10 ** q3, whislo=10 ** lv.min(), whishi=10 ** lv.max(), fliers=[])


def _nice_floor(v):
    d = 10 ** np.floor(np.log10(v)); m = v / d
    return d * (5 if m >= 5 else 2 if m >= 2 else 1)


def _nice_ceil(v):
    d = 10 ** np.floor(np.log10(v)); m = v / d
    return d * (1 if m <= 1 else 2 if m <= 2 else 5 if m <= 5 else 10)


def _ticks(lo, hi):
    """1-2-5 ticks when the axis spans <= 3 decades, whole decades otherwise."""
    span = np.log10(hi / lo); k0, k1 = int(np.floor(np.log10(lo))) - 1, int(np.ceil(np.log10(hi))) + 1
    mult = (1, 2, 5) if span <= 3.2 else (1,)
    t = [m * 10.0 ** k for k in range(k0, k1 + 1) for m in mult]
    return [x for x in t if lo * 0.999 <= x <= hi * 1.001]


def _plain(v):
    return np.format_float_positional(v, trim="-", precision=8) if v < 1 else f"{v:g}"


def plot_box_chart(S: pd.DataFrame, P: pd.DataFrame, path: str, title: str, matrix_label="maize grain", order_col="PEC_over_MEC_geomean", ascending=False, row_height=0.62):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    S = S[S.cand_id.isin(P.cand_id)].copy()
    S = S.sort_values(order_col, ascending=ascending, na_position="last")
    rows = []
    for r in S.itertuples():
        g = P[P.cand_id == r.cand_id]
        pec = g.dropna(subset=["PEC_central"]); mec = g.dropna(subset=["MEC_ng_g"])
        pec = pec[pec.PEC_central > 0]; mec = mec[mec.MEC_ng_g > 0]
        if pec.empty and mec.empty: continue
        rows.append((r, pec, mec))
    if not rows: return None
    n = len(rows)
    allv = np.concatenate([np.r_[p.PEC_central.values, m.MEC_ng_g.values] for _, p, m in rows])      # axis fitted to the data; property-range lines are clipped at the axis edge
    allv = allv[np.isfinite(allv) & (allv > 0)]
    lo, hi = _nice_floor(allv.min() / 1.4), _nice_ceil(allv.max() * 1.4)
    if n < 10: row_height = row_height + (10 - n) * 0.07            # few chemicals: keep the rows from looking squashed
    H = row_height * n + 3.3
    fig, ax = plt.subplots(figsize=(10.5 if n >= 10 else 9.2, H), dpi=150)
    for i, (r, pec, mec) in enumerate(rows):
        y0 = n - 1 - i                                            # first chemical at the top
        solid = str(getattr(r, "prop_status", "")) == "table"
        for vals, y, kind in ((pec.PEC_central.values, y0 + 0.19, "pec"), (mec.MEC_ng_g.values, y0 - 0.19, "mec")):
            if len(vals) == 0: continue
            st = _stats(vals)
            face = (PEC_SOLID if solid else PEC_LIGHT) if kind == "pec" else "white"
            kw = dict(positions=[y], widths=0.3, showfliers=False, patch_artist=True, manage_ticks=False,
                      boxprops=dict(facecolor=face, edgecolor=INK, linewidth=0.9, hatch=("///" if (kind == "pec" and not solid) else None)),
                      medianprops=dict(color=INK, linewidth=1.6), whiskerprops=dict(color=INK, linewidth=0.9), capprops=dict(color=INK, linewidth=0.9))
            try: ax.bxp([st], orientation="horizontal", **kw)
            except TypeError: ax.bxp([st], vert=False, **kw)
            rng = np.random.default_rng(i * 7 + (0 if kind == "pec" else 1)); jit = (rng.random(len(vals)) - 0.5) * 0.12
            if kind == "pec":
                cd = (pec.src.values == "CD") if "src" in pec.columns else np.zeros(len(vals), bool)
                ax.scatter(vals[~cd], (y + jit)[~cd], s=22, facecolor=PEC_EDGE, edgecolor="white", linewidth=0.5, zorder=4)
                ax.scatter(vals[cd], (y + jit)[cd], s=34, marker="D", facecolor=PEC_EDGE, edgecolor="white", linewidth=0.5, zorder=4)
            else:
                cd = (mec.src.values == "CD") if "src" in mec.columns else np.zeros(len(vals), bool)
                filled = (mec.grain_id.fillna("").values == "MS2 confirmed") & ~cd
                ax.scatter(vals[filled], (y + jit)[filled], s=24, facecolor=INK, edgecolor="white", linewidth=0.5, zorder=4)
                ax.scatter(vals[~filled & ~cd], (y + jit)[~filled & ~cd], s=24, facecolor="white", edgecolor=INK, linewidth=1.1, zorder=4)
                ax.scatter(vals[cd], (y + jit)[cd], s=34, marker="D", facecolor="white", edgecolor=INK, linewidth=1.1, zorder=4)
        if i > 0: ax.axhline(y0 + 0.5, color="#888888", linestyle=(0, (3, 3)), linewidth=0.8, zorder=0)
    ax.set_xscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(-0.55, n - 0.45)
    ax.set_yticks(range(n))
    lab = lambda r: (str(r.name)[:44] + ("..." if len(str(r.name)) > 44 else "")) + ((f"  [{r.method}]" if r.level == "CD" else f"  [{r.level} | {r.method}]") if hasattr(r, "method") else f"  [{r.level}]")
    ax.set_yticklabels([lab(r) for r, _, _ in rows][::-1], fontsize=9, color=INK)
    ax.set_xlabel(f"{matrix_label} concentration (ng/g dw, semi-quantified)", fontsize=10, color=INK)
    ticks = _ticks(lo, hi); ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: _plain(v))); ax.xaxis.set_minor_locator(NullLocator())
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", labelsize=9, colors=INK); ax.tick_params(axis="y", length=0)
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    handles = [Patch(facecolor=PEC_SOLID, edgecolor=INK, label="PEC (measured properties)"), Patch(facecolor=PEC_LIGHT, edgecolor=INK, hatch="///", label="PEC (estimated properties)"),
               Patch(facecolor="white", edgecolor=INK, label="MEC (semi-quantified)"),
               Line2D([], [], marker="o", linestyle="", markerfacecolor=PEC_EDGE, markeredgecolor="white", label="PEC per field"),
               Line2D([], [], marker="o", linestyle="", markerfacecolor=INK, markeredgecolor="white", label="MEC per field, MS2 confirmed"),
               Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor=INK, label="MEC per field, MS1 only / weak MS2")]
    if any("src" in p.columns and (p.src == "CD").any() for _, p, _ in rows) or any("src" in m.columns and (m.src == "CD").any() for _, _, m in rows):
        handles += [Line2D([], [], marker="D", linestyle="", markerfacecolor=PEC_EDGE, markeredgecolor="white", label="PEC from Compound Discoverer water value"),
                    Line2D([], [], marker="D", linestyle="", markerfacecolor="white", markeredgecolor=INK, label="MEC from Compound Discoverer (no MS2 level)")]
    nleg = int(np.ceil(len(handles) / 3)); cap_y = 0.06 + 0.27 * nleg + 0.12
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.06 / H))
    cap = "box = interquartile range, line = median, whiskers = min-max over fields, dots = individual fields"
    has_method = any(hasattr(r, "method") for r, _, _ in rows)
    if has_method: cap += "\nidentification method in brackets: level | toolkit = mzML library / in-silico match, CD = Compound Discoverer annotation, toolkit + CD = both (values pooled)"
    fig.text(0.5, cap_y / H, cap, ha="center", va="bottom", fontsize=8, color="#555555")
    fig.tight_layout(rect=(0, (cap_y + (0.42 if has_method else 0.2)) / H, 1, 1))
    fig.savefig(path); plt.close(fig)
    return path
