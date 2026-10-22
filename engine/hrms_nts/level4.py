"""Level 4 (and 5): accurate-mass / formula matching of MS1 features against candidate databases.

Every compound is searched against every adduct of the polarity uniformly (no adduct bias).  Isotope-pattern
agreement (M+1, M+2 from the molecular formula) is scored where the isotopologue peaks are measurable."""
from __future__ import annotations
import logging
import numpy as np
import pandas as pd
from .common import isotope_pattern, ISO_SPACING, ppm_err

log = logging.getLogger("hrms_nts")


def observed_isotope_ratios(run, feat: pd.DataFrame, ppm=6.0) -> pd.DataFrame:
    """Observed M+1/M0 and M+2/M0 at the apex scan (max over apex +-1 scans)."""
    r1 = np.full(len(feat), np.nan); r2 = np.full(len(feat), np.nan); m0i = np.zeros(len(feat))
    for n, (mz0, sc) in enumerate(zip(feat["mz_apex"].values, feat["apex_scan"].values.astype(int))):
        best0 = best1 = best2 = 0.0
        for s in range(max(sc - 1, 0), min(sc + 2, run.n_ms1)):
            mz, it = run.ms1_scan(s)
            if len(mz) == 0:
                continue
            def inten(target, tol_ppm):
                lo, hi = np.searchsorted(mz, [target * (1 - tol_ppm * 1e-6), target * (1 + tol_ppm * 1e-6)])
                return float(it[lo:hi].max()) if hi > lo else 0.0
            i0 = inten(mz0, ppm)
            if i0 > best0:
                best0 = i0
                best1 = inten(mz0 + ISO_SPACING, ppm * 1.5)
                lo, hi = np.searchsorted(mz, [mz0 + 1.9958 - 0.006, mz0 + 2.0067 + 0.006])
                best2 = float(it[lo:hi].max()) if hi > lo else 0.0
        m0i[n] = best0
        if best0 > 0:
            r1[n] = best1 / best0; r2[n] = best2 / best0
    out = feat[["feature_id"]].copy()
    out["obs_m1"] = r1; out["obs_m2"] = r2; out["m0_scan_int"] = m0i
    return out


def _iso_score(formula, obs_m1, obs_m2, m0_int, floor):
    try:
        pat = isotope_pattern(formula)
    except Exception:
        return np.nan, np.nan, np.nan
    e1, e2 = pat[1], pat[2]
    if not np.isfinite(obs_m1):
        return np.nan, e1, e2
    scores = []
    # M+1 informative only if the expected isotopologue is above the stored-peak floor
    if m0_int * e1 >= 3 * floor:
        tol = max(0.35 * e1, 0.03)
        scores.append(max(0.0, 1 - abs(obs_m1 - e1) / tol))
    if e2 >= 0.08 and m0_int * e2 >= 3 * floor and np.isfinite(obs_m2):   # Cl / Br / S signatures
        tol2 = max(0.4 * e2, 0.04)
        scores.append(max(0.0, 1 - abs(obs_m2 - e2) / tol2))
    if not scores:
        return np.nan, e1, e2
    return float(np.mean(scores)), e1, e2


def match_features(feat: pd.DataFrame, db, adducts: dict, ppm: float, run=None, store_floor=1e4,
                   drop_isotopes=True, min_iso_score=0.0, use_isotopes=True) -> pd.DataFrame:
    if feat.empty:
        return pd.DataFrame()
    f = feat
    if drop_isotopes and "iso_of" in f.columns:
        f = f[f["iso_of"].astype(str) == ""]
    f = f.reset_index(drop=True)
    iso = observed_isotope_ratios(run, f) if (run is not None and use_isotopes) else None
    rows = []
    for adduct, shift in adducts.items():
        neutral = f["mz"].values - shift
        qi, di = db.search(neutral, ppm)
        if len(qi) == 0:
            continue
        h = f.iloc[qi].reset_index(drop=True)
        c = db.df.iloc[di].reset_index(drop=True)
        calc = c["exact_mass"].values + shift
        h = pd.concat([h, c], axis=1)
        h["adduct"] = adduct; h["calc_mz"] = calc; h["ppm_error"] = ppm_err(h["mz"].values, calc)
        rows.append(h)
    # user compounds given as direct precursor m/z (no adduct)
    dp = db.df[db.df["direct_precursor_mz"].notna()]
    if len(dp):
        for _, c in dp.iterrows():
            e = ppm_err(f["mz"].values, c["direct_precursor_mz"])
            m = np.abs(e) <= ppm
            if m.any():
                h = f[m].reset_index(drop=True)
                for col in c.index: h[col] = c[col]
                h["adduct"] = "direct"; h["calc_mz"] = c["direct_precursor_mz"]; h["ppm_error"] = e[m]
                rows.append(h)
    if not rows:
        log.info("no Level 4 matches")
        return pd.DataFrame()
    hits = pd.concat(rows, ignore_index=True)
    if iso is not None:
        hits = hits.merge(iso, on="feature_id", how="left")
        sc = [_iso_score(fo, a, b, c, store_floor) for fo, a, b, c in
              zip(hits["formula"], hits["obs_m1"], hits["obs_m2"], hits["m0_scan_int"])]
        hits["iso_score"] = [s[0] for s in sc]; hits["exp_m1"] = [s[1] for s in sc]; hits["exp_m2"] = [s[2] for s in sc]
        hits["iso_contradicts"] = hits["iso_score"].notna() & (hits["iso_score"] <= min_iso_score) & (min_iso_score > 0)
        if min_iso_score > 0:
            hits = hits[~hits["iso_contradicts"]]
    else:
        hits["iso_score"] = np.nan
    # corroborating adducts: same compound seen at several adducts co-eluting within the sample
    hits = hits.reset_index(drop=True)
    # co-eluting adducts of the same compound: count distinct adducts within +-0.1 min (two overlapping RT bins)
    hits["n_adducts_seen"] = 1
    has_ik = hits["ik1"].astype(str).str.len() > 3
    sub = hits[has_ik]
    best = pd.Series(1, index=sub.index)
    for shift in (0.0, 0.1):
        key = (sub["rt_apex"] + shift) // 0.2
        cnt = sub.groupby([sub["sample"], sub["ik1"], key])["adduct"].transform("nunique")
        best = np.maximum(best, cnt)
    hits.loc[sub.index, "n_adducts_seen"] = best.values
    hits["n_candidates_for_feature"] = hits.groupby("feature_id")["feature_id"].transform("size")
    hits["level"] = "4"
    log.info("Level 4: %d hits on %d features (%d unique compounds)", len(hits), hits["feature_id"].nunique(),
             hits["ik1"].nunique())
    return hits
