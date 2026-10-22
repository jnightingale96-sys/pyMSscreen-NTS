"""MS1 feature detection.

Unlike a base-peak-per-scan approach, ALL centroided peaks above the intensity threshold are used, so
co-eluting ions in a busy scan are not lost.  Peaks are clustered in m/z by ppm adjacency, split in RT
on gaps, optionally split again at chromatographic valleys (continuous-background ions otherwise form
one run-length 'feature'), then summarised (apex, width, area, continuity)."""
from __future__ import annotations
import logging
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from .common import ISO_SPACING

log = logging.getLogger("hrms_nts")


def _mz_clusters(mz_sorted: np.ndarray, ppm: float) -> np.ndarray:
    """Cluster id per sorted m/z: gap-based, with greedy re-windowing for over-wide (chained) clusters."""
    if len(mz_sorted) == 0:
        return np.array([], dtype=np.int64)
    gap_ppm = np.diff(mz_sorted) / mz_sorted[:-1] * 1e6
    cl = np.concatenate([[0], np.cumsum(gap_ppm > ppm)]).astype(np.int64)
    # chained clusters: span > 3*ppm -> greedy windows of 2*ppm anchored on the lowest peak
    starts = np.flatnonzero(np.concatenate([[True], np.diff(cl) != 0]))
    ends = np.concatenate([starts[1:], [len(cl)]])
    span = (mz_sorted[ends - 1] - mz_sorted[starts]) / mz_sorted[starts] * 1e6
    bad = np.flatnonzero(span > 3 * ppm)
    if len(bad):
        cl = cl.copy()
        next_id = cl.max() + 1
        for b in bad:
            i, e = starts[b], ends[b]
            while i < e:
                j = np.searchsorted(mz_sorted[i:e], mz_sorted[i] * (1 + 2 * ppm * 1e-6), side="right") + i
                cl[i:j] = next_id; next_id += 1; i = j
    return cl


def detect_features(run, ppm_tol=10.0, rt_gap_tol=0.3, min_scans=3, min_intensity=1e5,
                    rt_range=None, split_broad_min=0.6, split_prominence=0.2, max_width_flag=2.0,
                    bg_occupancy=0.30, bg_peak_to_median=15.0) -> pd.DataFrame:
    sel = run.inten >= min_intensity
    mz = run.mz[sel]; it = run.inten[sel].astype(np.float64)
    sc = run.scan_index_of_peaks()[sel]
    rt = run.ms1_rt[sc]
    if rt_range:
        m = (rt >= rt_range[0]) & (rt <= rt_range[1])
        mz, it, sc, rt = mz[m], it[m], sc[m], rt[m]
    if len(mz) == 0:
        return pd.DataFrame()
    o = np.argsort(mz, kind="mergesort")
    mz, it, sc, rt = mz[o], it[o], sc[o], rt[o]
    cl = _mz_clusters(mz, ppm_tol)
    df = pd.DataFrame({"cl": cl, "mz": mz, "inten": it, "scan": sc, "rt": rt})
    df = df.sort_values(["cl", "rt"], kind="mergesort").reset_index(drop=True)
    newg = (df["cl"].values[1:] != df["cl"].values[:-1]) | (np.diff(df["rt"].values) > rt_gap_tol)
    df["g"] = np.concatenate([[0], np.cumsum(newg)])
    df["w"] = df["inten"] * df["mz"]

    # per (group, scan): max intensity (a scan may hold several peaks in the same ppm cluster)
    ps = (df.sort_values(["g", "scan", "inten"], ascending=[True, True, False])
            .drop_duplicates(["g", "scan"])[["g", "scan", "rt", "inten", "mz"]]
            .sort_values(["g", "rt"]).reset_index(drop=True))

    # valley-splitting of broad groups
    width = ps.groupby("g")["rt"].agg(lambda s: s.max() - s.min())
    broad = width[width > split_broad_min].index.values
    if len(broad):
        gmax = ps["g"].max() + 1
        ps = ps.set_index("g")
        add = []
        for g in broad:
            sub = ps.loc[g]
            y = sub["inten"].values
            if len(y) < 9:
                continue
            ys = np.convolve(y, np.ones(3) / 3, mode="same")
            pk, _ = find_peaks(ys, prominence=split_prominence * ys.max())
            if len(pk) < 2:
                continue
            cuts = [pk[i] + int(np.argmin(ys[pk[i]:pk[i + 1] + 1])) for i in range(len(pk) - 1)]
            bounds = [0] + [c + 1 for c in cuts] + [len(y)]
            lab = np.zeros(len(y), dtype=np.int64)
            for k in range(len(bounds) - 1):
                lab[bounds[k]:bounds[k + 1]] = k
            newids = np.where(lab == 0, g, gmax + lab - 1)
            gmax += lab.max()
            add.append((g, sub.index.values, newids))
        ps = ps.reset_index()
        if add:
            # assign by position within each original group
            pos = ps.groupby("g").cumcount().values
            ps_g = ps["g"].values.copy()
            for g, _, newids in add:
                idx = np.flatnonzero(ps["g"].values == g)
                ps_g[idx] = newids[:len(idx)]
            ps["g"] = ps_g
            ps = ps.sort_values(["g", "rt"]).reset_index(drop=True)

    ps["trap"] = 0.0
    same = ps["g"].values[1:] == ps["g"].values[:-1]
    seg = 0.5 * (ps["inten"].values[1:] + ps["inten"].values[:-1]) * np.diff(ps["rt"].values)
    ps.loc[1:, "trap"] = np.where(same, seg, 0.0)

    agg = ps.groupby("g").agg(n_scans=("scan", "size"), rt_start=("rt", "min"), rt_end=("rt", "max"),
                              max_int=("inten", "max"), sum_int=("inten", "sum"), area=("trap", "sum"))
    apex = ps.loc[ps.groupby("g")["inten"].idxmax(), ["g", "rt", "mz", "scan"]].set_index("g")
    agg["rt_apex"] = apex["rt"]; agg["mz_apex"] = apex["mz"]; agg["apex_scan"] = apex["scan"]
    # intensity-weighted m/z over the whole feature (use per-scan maxima)
    wmz = ps.assign(w=ps["inten"] * ps["mz"]).groupby("g")[["w", "inten"]].sum()
    agg["mz"] = wmz["w"] / wmz["inten"]
    agg["rt_width"] = agg["rt_end"] - agg["rt_start"]
    agg = agg[(agg["n_scans"] >= min_scans) & (agg["max_int"] >= min_intensity)].copy()
    # continuity: fraction of MS1 scans inside the feature window that contain the ion
    lo = np.searchsorted(run.ms1_rt, agg["rt_start"].values, side="left")
    hi = np.searchsorted(run.ms1_rt, agg["rt_end"].values, side="right")
    agg["continuity"] = agg["n_scans"].values / np.maximum(hi - lo, 1)
    agg["broad_flag"] = agg["rt_width"] > max_width_flag
    # background / chemical-noise test: how many MS1 scans of the whole run contain this m/z (any intensity >= store floor)?
    allmz = run.mz; allsc = run.scan_index_of_peaks(); allit = run.inten
    go = np.argsort(allmz, kind="mergesort"); gmz, gsc, git = allmz[go], allsc[go], allit[go]
    glo = np.searchsorted(gmz, agg["mz"].values * (1 - ppm_tol * 0.5e-6), side="left")
    ghi = np.searchsorted(gmz, agg["mz"].values * (1 + ppm_tol * 0.5e-6), side="right")
    occ = np.zeros(len(agg)); p2m = np.zeros(len(agg))
    for k in range(len(agg)):
        if ghi[k] > glo[k]:
            occ[k] = len(np.unique(gsc[glo[k]:ghi[k]])) / max(run.n_ms1, 1)
            p2m[k] = agg["max_int"].values[k] / max(float(np.median(git[glo[k]:ghi[k]])), 1.0)
    agg["mz_occupancy"] = occ; agg["peak_to_median"] = p2m
    agg["background_flag"] = (occ > bg_occupancy) & (p2m < bg_peak_to_median)
    # single-scan-integrated area (trapz) can be 0 for 2-pt features; guard
    agg["area"] = agg["area"].where(agg["area"] > 0, agg["sum_int"] * np.median(np.diff(run.ms1_rt)) if run.n_ms1 > 1 else 0)
    agg = agg.reset_index(drop=True)
    agg.insert(0, "feature_id", [f"{run.name}_F{i:06d}" for i in range(len(agg))])
    agg["sample"] = run.name
    agg["polarity"] = run.polarity
    log.info("%s: %d features (>= %d scans, >= %g counts)", run.name, len(agg), min_scans, min_intensity)
    return agg[["feature_id", "sample", "polarity", "mz", "mz_apex", "rt_apex", "rt_start", "rt_end", "rt_width",
                "n_scans", "max_int", "sum_int", "area", "continuity", "broad_flag", "mz_occupancy", "peak_to_median",
                "background_flag", "apex_scan"]]


def annotate_relations(feat: pd.DataFrame, ppm=7.0, rt_tol=0.06) -> pd.DataFrame:
    """Flag isotopologue peaks (M+1/M+2 of a more intense co-eluting feature) and adduct partners.
    Adds: iso_of, iso_k, adduct_partner_of, adduct_delta."""
    f = feat.sort_values("rt_apex").reset_index(drop=True)
    n = len(f)
    iso_of = np.full(n, "", dtype=object); iso_k = np.zeros(n, dtype=int)
    add_of = np.full(n, "", dtype=object); add_d = np.full(n, "", dtype=object)
    rt = f["rt_apex"].values; mz = f["mz"].values; mi = f["max_int"].values; ids = f["feature_id"].values
    deltas = {"Na-H": 21.981942, "K-H": 37.955882, "NH4-H": 17.026549, "HCOO-H(neg)": 46.005479,
              "Cl-H(neg)": 35.976678, "H2O-loss": -18.010565}
    for i in range(n):
        lo = np.searchsorted(rt, rt[i] - rt_tol, side="left"); hi = np.searchsorted(rt, rt[i] + rt_tol, side="right")
        if hi - lo < 2:
            continue
        for j in range(lo, hi):
            if j == i:
                continue
            d = mz[i] - mz[j]
            # i is isotopologue of j?
            if mi[j] > mi[i]:
                for k in (1, 2):
                    if abs(d - k * ISO_SPACING) / mz[j] * 1e6 <= ppm and 0.003 <= mi[i] / mi[j] <= 1.3:
                        iso_of[i] = ids[j]; iso_k[i] = k
            for nm, dd in deltas.items():
                if abs(d - dd) / mz[j] * 1e6 <= ppm and mi[j] >= mi[i] * 0.2 and not add_of[i]:
                    add_of[i] = ids[j]; add_d[i] = nm
    f["iso_of"] = iso_of; f["iso_k"] = iso_k; f["adduct_partner_of"] = add_of; f["adduct_delta"] = add_d
    return f


def eic_area(run, mz0, rt_apex, ppm=5.0, rt_half_window=0.25):
    """Re-integrate an EIC (trapezoid, minutes) around an RT apex and return (area, apex_rt, apex_int, snr)."""
    rt, y = run.eic(mz0, ppm, rt_apex - 3 * rt_half_window, rt_apex + 3 * rt_half_window)
    if len(y) < 3 or y.max() <= 0:
        return 0.0, rt_apex, 0.0, 0.0
    c = np.argmax(np.where(np.abs(rt - rt_apex) <= rt_half_window, y, 0))
    apex_int = y[c]
    # walk down to <5% of apex or minimum
    l = c
    while l > 0 and y[l - 1] <= y[l] and y[l - 1] > 0.05 * apex_int:
        l -= 1
    r = c
    while r < len(y) - 1 and y[r + 1] <= y[r] and y[r + 1] > 0.05 * apex_int:
        r += 1
    area = float(np.trapezoid(y[l:r + 1], rt[l:r + 1])) if r > l else 0.0
    base = np.concatenate([y[:l], y[r + 1:]])
    noise = np.median(base[base > 0]) if (base > 0).any() else 0.0
    snr = apex_int / noise if noise > 0 else np.inf
    return area, float(rt[c]), float(apex_int), float(snr)
