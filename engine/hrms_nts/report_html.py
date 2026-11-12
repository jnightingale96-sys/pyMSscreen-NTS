"""Self-contained HTML report (NTS_report.html) for a workflow.Project - one file, figures embedded, opens in any browser.

Sections: summary, injections, QC (IS response per injection, IS reproducibility), calibration (table + curves), matrix
effect (map + table), identified compounds (level, evidence, confidence), semi-quantification (compound x sample with
uncertainty), methods and definitions.
"""
from __future__ import annotations
import base64
import datetime
import html
import io
import os
import numpy as np
import pandas as pd

LEVEL_ORDER = ["1", "2a", "2b", "3", "4", "5"]
CSS = """
:root{--ink:#1b2422;--muted:#5a6663;--rule:#d9dfdc;--bg:#f7f8f6;--card:#ffffff;--accent:#1f6f78;--accent-soft:#e3f0f1;
--good:#2e7d4f;--warn:#b7791f;--bad:#b03a2e}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.55 "Segoe UI",Arial,sans-serif}
header{background:var(--accent);color:#fff;padding:28px max(16px,4vw)}
header h1{margin:0 0 4px;font-size:26px;font-weight:600}
header p{margin:0;opacity:.85;font-size:13px}
main{max-width:1180px;margin:0 auto;padding:24px max(16px,3vw) 64px}
section{margin-top:32px}
h2{font-size:18px;margin:0 0 10px;padding-bottom:6px;border-bottom:2px solid var(--accent);color:var(--accent)}
h3{font-size:14px;margin:18px 0 6px}
p.note{color:var(--muted);font-size:12.5px;margin:4px 0 10px;max-width:90ch}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}
.card{background:var(--card);border:1px solid var(--rule);border-radius:6px;padding:12px 14px}
.card b{display:block;font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted)}
.card span{font-size:22px;font-weight:600}
.card small{display:block;color:var(--muted);font-size:12px}
.tw{overflow-x:auto;background:var(--card);border:1px solid var(--rule);border-radius:6px}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th{background:var(--accent-soft);text-align:left;padding:6px 8px;position:sticky;top:0;font-weight:600;white-space:nowrap}
td{padding:5px 8px;border-top:1px solid var(--rule);vertical-align:top}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.badge{display:inline-block;padding:1px 7px;border-radius:10px;font-size:11.5px;font-weight:600;color:#fff}
.L1,.L2a{background:var(--good)}.L2b{background:var(--accent)}.L3{background:var(--warn)}.L4,.L5{background:var(--bad)}
.ok{color:var(--good);font-weight:600}.warn{color:var(--warn);font-weight:600}.bad{color:var(--bad);font-weight:600}
img{max-width:100%;background:#fff;border:1px solid var(--rule);border-radius:6px}
dl{display:grid;grid-template-columns:200px 1fr;gap:6px 16px;margin:0}dt{font-weight:600}dd{margin:0;color:var(--muted)}
footer{color:var(--muted);font-size:12px;text-align:center;padding:20px}
"""


def _img(fig_or_path):
    if isinstance(fig_or_path, str):
        if not os.path.exists(fig_or_path):
            return ""
        data = open(fig_or_path, "rb").read()
    else:
        buf = io.BytesIO()
        fig_or_path.savefig(buf, format="png", dpi=130, bbox_inches="tight")
        data = buf.getvalue()
    return f'<img alt="" src="data:image/png;base64,{base64.b64encode(data).decode()}">'


def _fmt(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (float, np.floating)):
        a = abs(v)
        return f"{v:.3g}" if a and (a < 0.01 or a >= 1e5) else f"{v:,.2f}".rstrip("0").rstrip(".")
    return html.escape(str(v))


def _table(df: pd.DataFrame, max_rows=400, level_col=None):
    if df is None or df.empty:
        return '<p class="note">No data.</p>'
    d = df.head(max_rows)
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in d.columns)
    rows = []
    for _, r in d.iterrows():
        cells = []
        for c in d.columns:
            v = r[c]
            if c == level_col:
                cells.append(f'<td><span class="badge L{html.escape(str(v))}">{html.escape(str(v))}</span></td>')
            else:
                num = isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))
                cells.append(f'<td class="num">{_fmt(v)}</td>' if num else f"<td>{_fmt(v)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    more = f'<p class="note">Showing {max_rows} of {len(df)} rows - the Excel report has all of them.</p>' if len(df) > max_rows else ""
    return f'<div class="tw"><table><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>{more}'


def _cal_figure(fits, units, max_n=12):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .calibration import predict
    names = list(fits)[:max_n]
    if not names:
        return None
    n = len(names); cols = min(4, n); rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(3.2 * cols, 2.6 * rows), squeeze=False)
    for ax, nm in zip(axes.flat, names):
        f = fits[nm]["chosen"]
        x, y = np.asarray(f["x"], float), np.asarray(f["y"], float)
        xx = np.linspace(0, x.max() * 1.05, 100)
        ax.plot(xx, predict(f, xx), color="#1f6f78", lw=1.5)
        ax.scatter(x, y, s=16, color="#1b2422", zorder=3)
        ax.set_title(f"{nm}\n{f['model']}, R² {f['r2']:.4f}", fontsize=8, loc="left")
        ax.tick_params(labelsize=7)
        ax.set_xlabel(units, fontsize=7)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.tight_layout()
    out = _img(fig)
    plt.close(fig)
    return out


def write_html_report(P, path, notes=None):
    from . import __version__
    S = P._load("screen"); Q = P._load("qa"); C = P._load("calibration"); TG = P._load("targeted")
    SQ = P._load("semiquant"); ME = P._load("matrix_effect")
    parts = []
    # ---- summary cards
    cards = []
    if P.sheet is not None:
        vc = P.sheet.type.value_counts()
        cards.append(("Injections", str(len(P.sheet)), ", ".join(f"{v} {k}" for k, v in vc.items())))
    T = Q["identifications"] if Q is not None else (S["top"] if S is not None else None)
    if T is not None and len(T):
        best = T.assign(lv=T.level.map({k: i for i, k in enumerate(LEVEL_ORDER)})).sort_values("lv").drop_duplicates("name")
        vc = best.level.value_counts()
        cards.append(("Compounds identified", str(int(best.level.isin(["1", "2a", "2b", "3"]).sum())),
                      " · ".join(f"L{k}: {vc[k]}" for k in LEVEL_ORDER if k in vc.index)))
    if TG is not None:
        cards.append(("Internal standards found", str(TG["assigned"].name.nunique()), "Q1 / Q3 checked in every injection"))
    if ME is not None and len(ME["table"]):
        m = ME["table"].ME_pct
        cards.append(("Matrix effect (median)", f"{m.median():+.0f} %", f"{(m < -20).sum()} suppressed · {(m > 20).sum()} enhanced of {len(m)}"))
    if SQ is not None:
        cards.append(("Semi-quantification", "IS surrogate" if SQ.get("method") == "is" else "Elucidation distance",
                      f"{SQ['table']['sq_conc'].notna().sum()} detections quantified"))
    parts.append('<section><h2>Summary</h2><div class="cards">' + "".join(
        f"<div class='card'><b>{html.escape(a)}</b><span>{html.escape(b)}</span><small>{html.escape(c)}</small></div>"
        for a, b, c in cards) + "</div></section>")
    # ---- injections
    if P.sheet is not None:
        parts.append("<section><h2>Injections</h2>" + _table(P.sheet[[c for c in ("file", "type", "group", "replicate", "cal_level", "conc", "notes")
                                                                         if c in P.sheet.columns]]) + "</section>")
    # ---- QC
    if TG is not None:
        q = ["<section><h2>Quality control</h2>",
             '<p class="note">Internal-standard response per injection relative to the batch median (values below 0.5 flag '
             "suppression or an injection problem), and the reproducibility of each IS across injections.</p>",
             "<h3>Injections</h3>", _table(TG["injection_qc"]), "<h3>Internal standards</h3>",
             _table(TG["is_qc"].sort_values("area_rsd_pct", ascending=False) if "area_rsd_pct" in TG["is_qc"] else TG["is_qc"])]
        parts.append("".join(q) + "</section>")
    # ---- calibration
    if C is not None and len(C["summary"]):
        cs = C["summary"][[c for c in ("name", "model", "n_levels", "range_low", "range_high", "r2", "lod", "loq", "accuracy_pct", "warning")
                           if c in C["summary"].columns]]
        parts.append("<section><h2>Calibration</h2><p class='note'>Weighted linear and quadratic fits; the quadratic is used only when "
                     "the curvature is significant (F-test), improves AICc and is monotonic over the range.</p>" + _table(cs)
                     + (_cal_figure(C["fits"], C.get("units", "")) or "") + "</section>")
    # ---- matrix effect
    if ME is not None and len(ME["table"]):
        parts.append("<section><h2>Matrix effect</h2><p class='note'>ME % = 100 × (area in spiked control extract / area in solvent "
                     "standard − 1). Below −20 % = ion suppression, above +20 % = enhancement. The line is the running median "
                     "across the chromatogram.</p>" + _img(ME["png"])
                     + _table(ME["table"][["name", "rt", "ME_pct", "effect", "n_solvent", "n_matrix"]].sort_values("ME_pct")) + "</section>")
    # ---- identified compounds
    if T is not None and len(T):
        cols = {"name": "Compound", "level": "Level", "formula": "Formula", "adduct": "Adduct", "rt_apex": "RT (min)",
                "ppm_error": "Q1 error (ppm)", "lib_score": "Library score", "bi_expl_adj_n": "Explained product ions",
                "presence_confidence_pct": "Presence confidence (%)", "peak_rating": "Peak rating (0-10)"}
        ok = T[T.level.isin(["1", "2a", "2b", "3"])].copy()
        ok["lv"] = ok.level.map({k: i for i, k in enumerate(LEVEL_ORDER)})
        b = ok.sort_values(["lv", "presence_confidence_pct" if "presence_confidence_pct" in ok else "lv"], ascending=[True, False])
        n_s = b.groupby("name")["sample"].nunique() if "sample" in b else None
        b = b.drop_duplicates("name")
        tab = b[[c for c in cols if c in b.columns]].rename(columns=cols)
        if n_s is not None:
            tab.insert(2, "Samples detected", b.name.map(n_s).values)
        parts.append("<section><h2>Identified compounds</h2><p class='note'>One row per compound (best evidence across samples), "
                     "Levels 1–3 after Schymanski et al. (2014). Every detection per sample is in the Excel report.</p>"
                     + _table(tab, level_col="Level") + "</section>")
    # ---- semi-quant
    if SQ is not None and len(SQ["table"]):
        st = SQ["table"].dropna(subset=["sq_conc"])
        unit = st["sq_unit"].iloc[0] if len(st) else ""
        if len(st) and "sample" in st:
            pv = st.assign(v=[f"{c:.3g} ({lo:.2g}–{hi:.2g})" for c, lo, hi in zip(st.sq_conc, st.sq_low, st.sq_high)]) \
                   .pivot_table(index="name", columns="sample", values="v", aggfunc="first").reset_index().fillna("")
            meth = ("Each compound takes the molar response of the most similar internal standards (inverse-distance weighted kNN) "
                    "measured in the IS-mix injection, corrected for matrix effect; the range is the reported uncertainty."
                    if SQ.get("method") == "is" else
                    "Each compound is matched to the nearest calibrant by elucidation distance; the range is the reported uncertainty.")
            parts.append(f"<section><h2>Semi-quantification ({html.escape(unit)})</h2><p class='note'>{meth} "
                         "Values are estimates, not calibrated concentrations.</p>" + _table(pv) + "</section>")
    # ---- methods
    if notes:
        parts.append("<section><h2>Methods and definitions</h2><dl>" + "".join(
            f"<dt>{html.escape(str(a))}</dt><dd>{html.escape(str(b))}</dd>" for a, b in notes) + "</dl></section>")
    doc = (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
           f"<title>HRMS NTS report</title><style>{CSS}</style></head><body><header><h1>HRMS non-target screening report</h1>"
           f"<p>pyMSscreen-NTS engine {__version__} · generated {datetime.datetime.now():%d %b %Y %H:%M} · data: "
           f"{html.escape(P.data_dir)}</p></header><main>{''.join(parts)}</main>"
           f"<footer>Results folder: {html.escape(P.out_dir)}</footer></body></html>")
    open(path, "w", encoding="utf-8").write(doc)
    return path
