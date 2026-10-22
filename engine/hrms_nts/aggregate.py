"""Cross-sample alignment, blank subtraction and detection-frequency summaries."""
from __future__ import annotations
import numpy as np
import pandas as pd


def blank_filter(features: pd.DataFrame, blank_samples: list[str], ratio=3.0, ppm=10.0, rt_tol=0.2) -> pd.DataFrame:
    """Flag features present in procedural/solvent blanks (sample area < ratio x blank area)."""
    f = features.copy(); f["blank_flag"] = False; f["blank_ratio"] = np.nan
    if not blank_samples:
        return f
    b = f[f["sample"].isin(blank_samples)]
    if b.empty:
        return f
    bm, bt, ba = b["mz"].values, b["rt_apex"].values, b["area"].values
    o = np.argsort(bm); bm, bt, ba = bm[o], bt[o], ba[o]
    for i in f.index[~f["sample"].isin(blank_samples)]:
        mz, rt, ar = f.at[i, "mz"], f.at[i, "rt_apex"], f.at[i, "area"]
        lo, hi = np.searchsorted(bm, [mz * (1 - ppm * 1e-6), mz * (1 + ppm * 1e-6)])
        if hi > lo:
            m = np.abs(bt[lo:hi] - rt) <= rt_tol
            if m.any():
                bmax = ba[lo:hi][m].max()
                f.at[i, "blank_ratio"] = ar / bmax if bmax > 0 else np.inf
                f.at[i, "blank_flag"] = ar < ratio * bmax
    return f


def align_features(features: pd.DataFrame, ppm=10.0, rt_tol=0.15) -> pd.DataFrame:
    """Greedy m/z-RT alignment -> consensus feature id per row (column 'consensus_id') + presence table."""
    f = features.sort_values(["polarity", "mz"]).reset_index(drop=True)
    cid = np.full(len(f), -1, dtype=np.int64); n = 0
    for pol, g in f.groupby("polarity"):
        idx = g.index.values; mz = g["mz"].values; rt = g["rt_apex"].values
        used = np.zeros(len(idx), bool)
        order = np.argsort(-g["max_int"].values)
        for k in order:
            if used[k]:
                continue
            lo, hi = np.searchsorted(mz, [mz[k] * (1 - ppm * 1e-6), mz[k] * (1 + ppm * 1e-6)])
            cand = [j for j in range(lo, hi) if not used[j] and abs(rt[j] - rt[k]) <= rt_tol]
            # one feature per sample per consensus: keep most intense per sample
            sm = g["sample"].values
            seen = {}
            for j in cand:
                if sm[j] not in seen or g["max_int"].values[j] > g["max_int"].values[seen[sm[j]]]:
                    seen[sm[j]] = j
            for j in seen.values():
                cid[idx[j]] = n; used[j] = True
            n += 1
    f["consensus_id"] = [f"C{c:06d}" for c in cid]
    return f


def presence_table(aligned: pd.DataFrame, value="area") -> pd.DataFrame:
    p = aligned.pivot_table(index="consensus_id", columns="sample", values=value, aggfunc="max")
    p["n_samples"] = p.notna().sum(axis=1)
    p["detection_frequency"] = p["n_samples"] / (p.shape[1] - 1)
    return p
