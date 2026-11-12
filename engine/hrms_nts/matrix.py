"""Matrix effect calculated after the run (post-extraction spike method, Matuszewski et al. 2003):

    ME (%) = 100 x (area in spiked control extract / area in neat solvent standard - 1)
             < 0 = ion suppression, > 0 = enhancement; |ME| <= 20 % is commonly regarded as negligible.

Inputs: targeted results (Project.targeted) containing injections of type 'solvent_std' (mix in e.g. 10 % MeOH / water)
and 'matrix_spike' (the same mix spiked into a control / blank matrix extract). Replicates are averaged per type.
Outputs: one row per compound (ME %, n, RT) and a map of ME across the chromatogram (rolling median over RT), used to
flag suppression zones and, optionally, to correct semi-quantitative concentrations at a target's RT.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def matrix_effect(assigned: pd.DataFrame, sheet: pd.DataFrame, solvent_type="solvent_std", matrix_type="matrix_spike",
                  ratings=("Good", "Fair")) -> pd.DataFrame:
    files_s = set(sheet.loc[sheet.type == solvent_type, "file"])
    files_m = set(sheet.loc[sheet.type == matrix_type, "file"])
    if not files_s or not files_m:
        raise ValueError(f"matrix effect needs sample-sheet injections of type '{solvent_type}' and '{matrix_type}'")
    a = assigned[assigned.int_rating.astype(str).str.startswith(tuple(ratings)) & (assigned.int_area > 0)]
    rows = []
    for name, g in a.groupby("name"):
        s = g[g.file.isin(files_s)].int_area
        m = g[g.file.isin(files_m)].int_area
        if s.empty or m.empty:
            continue
        me = 100 * (m.mean() / s.mean() - 1)
        rows.append(dict(name=name, rt=float(g.int_apex_rt.median()), area_solvent=float(s.mean()), area_matrix=float(m.mean()),
                         n_solvent=len(s), n_matrix=len(m), ME_pct=me,
                         effect=("suppression" if me < -20 else "enhancement" if me > 20 else "negligible (|ME| <= 20 %)")))
    cols = ["name", "rt", "area_solvent", "area_matrix", "n_solvent", "n_matrix", "ME_pct", "effect"]
    return pd.DataFrame(rows, columns=cols).sort_values("rt").reset_index(drop=True)


def me_map(me: pd.DataFrame, window_min=1.0, step_min=0.1) -> pd.DataFrame:
    """Rolling median of ME % across the chromatogram (window +-window_min/2)."""
    if me.empty:
        return pd.DataFrame(columns=["rt", "ME_pct_median", "n"])
    grid = np.arange(np.floor(me.rt.min()), np.ceil(me.rt.max()) + step_min, step_min)
    out = []
    for t in grid:
        w = me[(me.rt - t).abs() <= window_min / 2]
        if len(w):
            out.append(dict(rt=round(float(t), 2), ME_pct_median=float(w.ME_pct.median()), n=len(w)))
    return pd.DataFrame(out)


def me_at_rt(mp: pd.DataFrame, rt: float) -> float:
    """Matrix effect (%) at a retention time from the map (nearest grid point; NaN outside the mapped range)."""
    if mp is None or mp.empty or not np.isfinite(rt) or rt < mp.rt.min() - 0.5 or rt > mp.rt.max() + 0.5:
        return np.nan
    return float(mp.iloc[(mp.rt - rt).abs().argmin()].ME_pct_median)


def plot_me(me: pd.DataFrame, mp: pd.DataFrame, path: str, title="Matrix effect across the chromatogram"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cap = 200.0                                     # enhancement above this is drawn at the cap (triangle + value)
    fig, ax = plt.subplots(figsize=(8, 3.8), dpi=130)
    ax.axhspan(-20, 20, color="#e8efe9", zorder=0, label="|ME| ≤ 20 % (negligible)")
    ax.axhline(0, color="#7a8580", lw=0.8)
    y = me.ME_pct.clip(upper=cap)
    over = me.ME_pct > cap
    col = np.where(me.ME_pct < -20, "#c0392b", np.where(me.ME_pct > 20, "#2a78d6", "#5a6663"))
    ax.scatter(me.rt[~over], y[~over], s=24, c=col[~over], edgecolor="white", linewidth=0.6, zorder=3)
    ax.scatter(me.rt[over], y[over], s=46, c="#2a78d6", marker="^", edgecolor="white", linewidth=0.6, zorder=3,
               label=f"off scale (> +{cap:.0f} %)")
    if len(mp):
        ax.plot(mp.rt, mp.ME_pct_median.clip(upper=cap), color="#1b2422", lw=1.6, zorder=4, label="running median")
    # in-range extremes labelled on the plot; off-scale compounds listed in a box (they cluster and would overlap)
    inr = me[~over]
    top = inr.reindex(inr.ME_pct.abs().sort_values(ascending=False).index).head(6)
    for _, r in top[top.ME_pct.abs() > 50].iterrows():
        ax.annotate(r["name"].split(" ")[0], (r.rt, r.ME_pct), fontsize=6.5, xytext=(4, 3), textcoords="offset points", color="#1b2422")
    if over.any():
        lst = "\n".join(f"{r['name'].split(' ')[0]}  {r.ME_pct:+.0f} %  ({r.rt:.1f} min)"
                        for _, r in me[over].sort_values("ME_pct", ascending=False).iterrows())
        ax.text(0.01, 0.97, "off scale:\n" + lst, transform=ax.transAxes, fontsize=6.5, va="top", color="#1b2422",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#d9dfdc"))
    ax.set_ylim(-105, cap + 25)
    ax.set_xlabel("retention time (min)")
    ax.set_ylabel("matrix effect (%)  ← suppression · enhancement →")
    ax.set_title(title, loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path
