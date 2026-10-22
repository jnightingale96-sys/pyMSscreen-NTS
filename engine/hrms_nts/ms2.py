"""DDA MS2 assignment to MS1 features / candidate precursors + pyMSscreen-style QA flags."""
from __future__ import annotations
import numpy as np
import pandas as pd


def clean_spectrum(mz, it, precursor_mz=None, rel_min=0.01, top_n=100, max_above_precursor=1.5):
    mz = np.asarray(mz, float); it = np.asarray(it, float)
    if len(mz) == 0:
        return mz, it
    if precursor_mz is not None:
        k = mz <= precursor_mz + max_above_precursor
        mz, it = mz[k], it[k]
    if len(it) == 0:
        return mz, it
    k = it >= rel_min * it.max()
    mz, it = mz[k], it[k]
    if len(it) > top_n:
        idx = np.argsort(it)[::-1][:top_n]
        idx.sort(); mz, it = mz[idx], it[idx]
    return mz, it


def merge_spectra(spectra, ppm=10.0):
    """Sum several MS2 scans of the same precursor (bin by ppm) -> consensus spectrum."""
    allmz = np.concatenate([s[0] for s in spectra]); allit = np.concatenate([s[1] for s in spectra])
    if len(allmz) == 0:
        return allmz, allit
    o = np.argsort(allmz); allmz, allit = allmz[o], allit[o]
    grp = np.concatenate([[0], np.cumsum(np.diff(allmz) / allmz[:-1] * 1e6 > ppm)])
    s = pd.DataFrame({"g": grp, "mz": allmz, "i": allit, "w": allmz * allit}).groupby("g").agg(i=("i", "sum"), w=("w", "sum"), n=("i", "size"))
    return (s["w"] / s["i"]).values, (s["i"] / len(spectra)).values


class MS2Index:
    def __init__(self, run):
        self.run = run
        m = run.ms2
        self.pmz = np.array([x["precursor_mz"] for x in m]) if m else np.array([])
        self.rt = np.array([x["rt"] for x in m]) if m else np.array([])
        self.o = np.argsort(self.pmz)
        self.pmz_s = self.pmz[self.o]

    def find(self, mz, rt_start, rt_end, tol_da=0.02, rt_pad=0.1):
        if len(self.pmz_s) == 0:
            return []
        lo, hi = np.searchsorted(self.pmz_s, [mz - tol_da, mz + tol_da])
        idx = self.o[lo:hi]
        return [int(i) for i in idx if rt_start - rt_pad <= self.rt[i] <= rt_end + rt_pad]

    def best_spectrum(self, mz, rt_apex, rt_start, rt_end, tol_da=0.02, rt_pad=0.1, merge_n=3, rel_min=0.01, top_n=100):
        ids = self.find(mz, rt_start, rt_end, tol_da, rt_pad)
        if not ids:
            return None
        ms2 = self.run.ms2
        # rank by closeness to apex, then total ion current
        ids = sorted(ids, key=lambda i: (abs(ms2[i]["rt"] - rt_apex), -ms2[i]["tic"]))
        take = ids[:merge_n]
        specs = [clean_spectrum(ms2[i]["mz"], ms2[i]["inten"], mz, rel_min, top_n) for i in take]
        mzs, its = merge_spectra(specs) if len(specs) > 1 else specs[0]
        mzs, its = clean_spectrum(mzs, its, mz, rel_min, top_n)
        return dict(mz=mzs, inten=its, n_scans=len(ids), rt=ms2[take[0]]["rt"], precursor_mz=ms2[take[0]]["precursor_mz"],
                    ce=ms2[take[0]]["ce"], scan_id=ms2[take[0]]["scan_id"], rt_delta=abs(ms2[take[0]]["rt"] - rt_apex))


def peaks_to_string(mz, it):
    return ";".join(f"{a:.4f}:{b:.0f}" for a, b in zip(mz, it))


def string_to_peaks(s):
    if not isinstance(s, str) or not s:
        return np.array([]), np.array([])
    a = [p.split(":") for p in s.split(";") if ":" in p]
    return np.array([float(x[0]) for x in a]), np.array([float(x[1]) for x in a])


def assign_ms2_to_features(feat: pd.DataFrame, run, tol_da=0.02, rt_pad=0.1, merge_n=3, rel_min=0.01, top_n=100,
                           ms2_int_thresh=2500.0) -> pd.DataFrame:
    """One consensus MS2 per feature.  Adds ms2_* columns (peak list as 'mz:int;...')."""
    idx = MS2Index(run)
    rec = []
    for r in feat.itertuples():
        sp = idx.best_spectrum(r.mz_apex, r.rt_apex, r.rt_start, r.rt_end, tol_da, rt_pad, merge_n, rel_min, top_n)
        if sp is None:
            rec.append(dict(feature_id=r.feature_id, ms2_n_scans=0, ms2_peaks="", ms2_n_peaks=0, ms2_rt=np.nan,
                            ms2_max_int=0.0, ms2_good_int=False, ms2_near_apex=False, ms2_ce=np.nan))
        else:
            rec.append(dict(feature_id=r.feature_id, ms2_n_scans=sp["n_scans"], ms2_peaks=peaks_to_string(sp["mz"], sp["inten"]),
                            ms2_n_peaks=len(sp["mz"]), ms2_rt=sp["rt"], ms2_max_int=float(sp["inten"].max()) if len(sp["inten"]) else 0.0,
                            ms2_good_int=bool(len(sp["inten"]) and sp["inten"].max() > ms2_int_thresh),
                            ms2_near_apex=sp["rt_delta"] <= max(0.1, rt_pad), ms2_ce=sp["ce"]))
    return feat.merge(pd.DataFrame(rec), on="feature_id", how="left")
