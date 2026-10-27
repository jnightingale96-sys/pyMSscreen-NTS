"""DIA (SWATH / AIF / MSE) MS2: correlation-based pseudo-MS2 spectra.

In DIA every precursor in a 25-200 Da window is fragmented together, so the MS2 scan cannot be assigned to a feature by
precursor m/z the way DDA scans are.  A fragment is taken to belong to a feature when, inside the feature's elution window,
its XIC in the feature's DIA window correlates with the precursor's MS1 XIC (Gaussian-smoothed Pearson r >= corr_min) and
apexes within apex_tol_scans cycles of it.  The resulting pseudo-MS2 fills the same ms2_* columns as the DDA path, so the
library search, in-silico scoring and the Level rules are unchanged.  Validated on SCIEX X500R SWATH data (2026-10)."""
from __future__ import annotations
import logging
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

from .ms2 import peaks_to_string

log = logging.getLogger("hrms_nts")


def pseudo_ms2(run, mz0, rt_start, rt_end, ms1_rt=None, ms1_int=None, corr_min=0.7, apex_tol=3, frag_ppm=20.0,
               max_candidates=150, min_frag_int=30.0, pad=0.05):
    """-> (window tuple or None, list of dict(mz, inten, r)) sorted by intensity.  ms1_rt/ms1_int: precursor MS1 XIC."""
    w = run.dia_window(mz0)
    if w is None:
        return None, []
    d = run.dia
    key = w[0]
    sc_all = np.flatnonzero(d["win"] == key)
    rts = d["rt"][sc_all]
    inpk = sc_all[(rts >= rt_start - pad) & (rts <= rt_end + pad)]
    if len(inpk) < 3:
        return w, []
    if ms1_rt is None:
        ms1_rt, ms1_int, _ = run.xic_max(mz0, 10.0)
    rt_in = d["rt"][inpk]
    prof_s = gaussian_filter1d(np.interp(rt_in, ms1_rt, ms1_int), 1.0)
    apex_j = int(np.argmax(prof_s))
    cand = []
    for j in range(max(0, apex_j - 1), min(len(inpk), apex_j + 2)):
        s = inpk[j]; a, b = d["off"][s], d["off"][s + 1]
        cand.append(np.c_[d["mz"][a:b], d["inten"][a:b]])
    cand = np.vstack(cand) if cand else np.empty((0, 2))
    cand = cand[(cand[:, 0] < mz0 + 1.5) & (cand[:, 1] >= min_frag_int)]
    cand = cand[np.argsort(-cand[:, 1])]
    taken = []
    for m, _ in cand:
        if all(abs(m - t) / t * 1e6 > 30 for t in taken):
            taken.append(m)
        if len(taken) >= max_candidates:
            break
    frags = []
    for m in taken:
        _, I, MZ = run.dia_xic(key, m, frag_ppm, scans=inpk)
        if I.max() <= 0 or np.count_nonzero(I) < 3:
            continue
        Is = gaussian_filter1d(I, 1.0)
        r = np.corrcoef(Is, prof_s)[0, 1] if np.std(Is) > 0 and np.std(prof_s) > 0 else 0.0
        if r >= corr_min and abs(int(np.argmax(Is)) - apex_j) <= apex_tol:
            top = I >= 0.3 * I.max()
            mzw = float(np.nansum(MZ[top] * I[top]) / I[top].sum())
            frags.append(dict(mz=mzw, inten=float(I[max(0, apex_j - 1):apex_j + 2].sum()), r=float(r)))
    frags.sort(key=lambda f: -f["inten"])
    return w, frags


def assign_dia_ms2_to_features(feat: pd.DataFrame, run, cfg: dict) -> pd.DataFrame:
    """DIA counterpart of ms2.assign_ms2_to_features: one correlation pseudo-MS2 per feature (most intense features first)."""
    dc = cfg.get("dia", {}); m2 = cfg["ms2"]
    rec = []
    order = feat.sort_values("max_int", ascending=False)
    elig = order[order["max_int"] >= dc.get("min_feature_int", 2000)].head(dc.get("max_features", 6000))
    elig_ids = set(elig["feature_id"])
    for r in feat.itertuples():
        empty = dict(feature_id=r.feature_id, ms2_n_scans=0, ms2_peaks="", ms2_n_peaks=0, ms2_rt=np.nan, ms2_max_int=0.0,
                     ms2_good_int=False, ms2_near_apex=False, ms2_ce=np.nan, ms2_dia_window="")
        if r.feature_id not in elig_ids:
            rec.append(empty); continue
        rt, I, _ = run.xic_max(r.mz_apex, cfg["ms1"]["feature_ppm"])
        w, fr = pseudo_ms2(run, r.mz_apex, r.rt_start, r.rt_end, rt, I, dc.get("corr_min", 0.7), dc.get("apex_tol_scans", 3),
                           dc.get("frag_ppm", 20.0), dc.get("max_candidates", 150), dc.get("min_frag_int", 30.0))
        if not fr:
            rec.append(empty); continue
        mz = np.array([f["mz"] for f in fr]); it = np.array([f["inten"] for f in fr])
        keep = it >= m2["rel_min"] * it.max()
        mz, it = mz[keep], it[keep]
        o = np.argsort(mz)
        rec.append(dict(feature_id=r.feature_id, ms2_n_scans=1, ms2_peaks=peaks_to_string(mz[o], it[o]), ms2_n_peaks=int(len(mz)),
                        ms2_rt=float(r.rt_apex), ms2_max_int=float(it.max()), ms2_good_int=bool(it.max() > m2["good_int"]),
                        ms2_near_apex=True, ms2_ce=float(w[3]), ms2_dia_window=f"{w[1]:.1f}-{w[2]:.1f}"))
    log.info("%s: DIA pseudo-MS2 built for %d of %d features (%d with >= 3 correlated fragments)", run.name, len(elig_ids), len(feat),
             sum(1 for x in rec if x["ms2_n_peaks"] >= 3))
    return feat.merge(pd.DataFrame(rec), on="feature_id", how="left")
