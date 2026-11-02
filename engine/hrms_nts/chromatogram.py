"""Q1 / Q3 chromatograms for one compound in one injection: precursor XIC (Q1) and product-ion XICs (Q3), integrated, rated.

Product ions come from (in this order) the caller (e.g. library-matched / annotated ions of an identification), the
correlated DIA pseudo-MS2 or the DDA spectrum at the peak, or in-silico prediction from the SMILES (flagged 'predicted').
DIA: each product ion is traced in the SWATH window that contains the precursor. DDA: product-ion intensities are read
from every MS2 scan whose precursor matches (sparse trace)."""
from __future__ import annotations
import numpy as np
import pandas as pd

from .targeted import integrate, integration_rating
from .quality import peak_rating, presence_confidence


def _dda_q3(run, mz0, frag_mz, rt_lo, rt_hi, tol_da=0.02, ppm=20.0):
    rts, vals = [], []
    for s in run.ms2:
        if rt_lo <= s["rt"] <= rt_hi and abs(s["precursor_mz"] - mz0) <= max(tol_da, s.get("iso_width", 1.0) / 2):
            m = s["mz"]; i = s["inten"]
            k = np.abs(m - frag_mz) <= frag_mz * ppm * 1e-6
            rts.append(s["rt"]); vals.append(float(i[k].max()) if k.any() else 0.0)
    o = np.argsort(rts)
    return np.asarray(rts)[o], np.asarray(vals)[o]


def product_ions_at_peak(run, mz0, rt_start, rt_end, cfg=None, n=6, smiles=None, adduct="[M+H]+"):
    """Observed product ions (DIA correlated / DDA spectrum); falls back to in-silico predicted ions if none observed."""
    from .targeted import ms2_for_peak
    cfg = cfg or {"dia": {}, "ms2": dict(precursor_tol_da=0.02, rt_pad_min=0.1, merge_scans=3, rel_min=0.01, top_n=100)}
    desc, fr = ms2_for_peak(run, mz0, dict(rt=(rt_start + rt_end) / 2, rt_start=rt_start, rt_end=rt_end), cfg)
    ions = [dict(mz=f["mz"], label=f"{f['mz']:.4f}", source="observed (" + desc.split(" ")[0] + ")") for f in fr if f["mz"] < mz0 - 1.5][:n]
    if not ions and isinstance(smiles, str):
        from .insilico import predict_product_ions
        p = predict_product_ions(smiles, adduct, top=n)
        ions = [dict(mz=float(r.mz), label=f"{r.mz:.4f} {r.label}", source="predicted (in-silico)") for r in p.itertuples()]
    return ions


def extract(run, mz0, rt0, product_ions=None, cfg=None, q1_ppm=10.0, q3_ppm=20.0, window=1.0, smiles=None, adduct="[M+H]+",
            n_q3=6, noise_floor=150.0):
    """-> dict(q1=(rt, I), q3=[dict(mz, label, source, rt, I, area, apex_rt)], integration, rating, ion_ratios)."""
    rt, I, MZ = run.xic_max(mz0, q1_ppm)
    m = np.abs(rt - rt0) <= window
    g = integrate(rt, I, rt0, noise_floor)
    pr = peak_rating(rt, I, g["start_rt"], g["end_rt"], g["sn"], noise_floor) if np.isfinite(g["start_rt"]) else peak_rating([], [], 0, 0)
    sel = (rt >= g["start_rt"]) & (rt <= g["end_rt"]) & np.isfinite(MZ)
    ppm = float(np.average((MZ[sel] - mz0) / mz0 * 1e6, weights=I[sel])) if sel.any() and I[sel].sum() > 0 else np.nan
    if product_ions is None:
        product_ions = product_ions_at_peak(run, mz0, g["start_rt"] if np.isfinite(g["start_rt"]) else rt0 - 0.1,
                                            g["end_rt"] if np.isfinite(g["end_rt"]) else rt0 + 0.1, cfg, n_q3, smiles, adduct)
    q3 = []
    win = run.dia_window(mz0) if run.acquisition == "DIA" else None
    for ion in product_ions[:n_q3]:
        if win is not None:
            r3, I3, _ = run.dia_xic(win[0], ion["mz"], q3_ppm)
        else:
            r3, I3 = _dda_q3(run, mz0, ion["mz"], rt0 - window, rt0 + window, ppm=q3_ppm)
        k = np.abs(r3 - rt0) <= window
        r3, I3 = r3[k], I3[k]
        s3 = (r3 >= g["start_rt"]) & (r3 <= g["end_rt"])
        area3 = float(np.trapezoid(I3[s3], r3[s3]) * 60) if s3.sum() > 1 else 0.0
        q3.append(dict(**ion, rt=r3, I=I3, area=area3, apex_rt=float(r3[s3][np.argmax(I3[s3])]) if s3.any() and I3[s3].max() > 0 else np.nan))
    ratios = {}
    if q3 and q3[0]["area"] > 0:
        ratios = {f"{x['mz']:.4f}": round(x["area"] / q3[0]["area"], 3) for x in q3[1:]}
    return dict(mz=mz0, rt0=rt0, ppm=ppm, q1=(rt[m], I[m]), integration=g, rating=pr, int_rating=integration_rating(g), q3=q3,
                ion_ratios=ratios, dia_window=f"{win[1]:.1f}-{win[2]:.1f} ({win[3]:g} eV)" if win else "DDA")


def summary(tr, lib_score=None, lib_matched=0, n_explained=None, replicate_frac=None, blank_ratio=None, iso_score=None, ppm_tol=5.0, weights=None):
    g = tr["integration"]
    conf, comp = presence_confidence(tr["rating"]["peak_rating"], g["sn"], tr["ppm"], ppm_tol, iso_score, lib_score, lib_matched,
                                     n_explained if n_explained is not None else sum(1 for x in tr["q3"] if x["area"] > 0),
                                     replicate_frac, blank_ratio, weights)
    return dict(apex_rt=round(g["apex_rt"], 3), start_rt=round(g["start_rt"], 3), end_rt=round(g["end_rt"], 3), area=g["area"],
                height=g["height"], sn=round(g["sn"], 1), q1_ppm=round(tr["ppm"], 2) if np.isfinite(tr["ppm"]) else None,
                peak_rating=tr["rating"]["peak_rating"], integration=tr["int_rating"], presence_confidence_pct=conf,
                n_q3_with_signal=sum(1 for x in tr["q3"] if x["area"] > 0), ion_ratios=tr["ion_ratios"], components=comp)


def plot(tr, title="", show_integration=True):
    """Plotly figure: Q1 (top) with the integrated region shaded, Q3 product-ion traces (bottom)."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08, row_heights=[0.5, 0.5],
                        subplot_titles=(f"Q1 precursor m/z {tr['mz']:.4f}", f"Q3 product ions ({tr['dia_window']})"))
    rt, I = tr["q1"]; g = tr["integration"]
    fig.add_trace(go.Scatter(x=rt, y=I, mode="lines", name="Q1", line=dict(color="#2a78d6", width=2)), row=1, col=1)
    if show_integration and np.isfinite(g["start_rt"]):
        s = (rt >= g["start_rt"]) & (rt <= g["end_rt"])
        fig.add_trace(go.Scatter(x=rt[s], y=I[s], fill="tozeroy", mode="none", name="integrated", fillcolor="rgba(42,120,214,0.25)"), row=1, col=1)
    pal = ["#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#7d4bc2", "#c43c3c", "#555555"]
    for k, x in enumerate(tr["q3"]):
        fig.add_trace(go.Scatter(x=x["rt"], y=x["I"], mode="lines", name=f"{x['label']}  [{x['source']}]",
                                 line=dict(color=pal[k % len(pal)], width=1.5, dash="dot" if "predicted" in x["source"] else "solid")), row=2, col=1)
    rs = tr["rating"]
    fig.update_layout(title=f"{title}<br><sup>RT {g['apex_rt']:.2f} min | area {g['area']:,.0f} | S/N {g['sn']:.0f} | "
                            f"peak rating {rs['peak_rating']}/10 | Q1 {tr['ppm']:+.1f} ppm</sup>", height=620,
                      legend=dict(font=dict(size=10)), margin=dict(l=60, r=20, t=90, b=40), template="plotly_white")
    fig.update_xaxes(title_text="RT (min)", row=2, col=1); fig.update_yaxes(title_text="counts", row=1, col=1); fig.update_yaxes(title_text="counts", row=2, col=1)
    return fig
