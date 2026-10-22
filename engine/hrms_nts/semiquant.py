"""Semi-quantification of non-target / suspect features.

Modes
  single : conc = ref_conc * area / ref_area   (one universal reference, e.g. carbamazepine 0.5 ug/mL)
  logie  : ionisation-efficiency-corrected ratio to the nearest calibrant in 'elucidation space'
           (z-scored logIE, RT, MW, logP distance; conc = cal_conc * (area/cal_area) * 10^(logIE_cal - logIE_cmpd))
  curve  : external multi-level calibration curves (calibrants.csv with several conc levels per compound) - used for
           compounds that have their own curve (targeted quantification), others fall back to --fallback mode
  none

IMPORTANT: area units must be the same for samples and calibrants.  Use `calibrate` (measures calibrant areas with this
toolkit from a standard mzML) rather than pasting areas from another package (e.g. Compound Discoverer) into config.yaml."""
from __future__ import annotations
import logging
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")


def logie_estimate(smiles):
    """Crude logIE proxy (positive ESI) from RDKit descriptors, as used previously in this project:
    0.379 + 0.653*logP - 0.206*HBA - 0.264*HBD + 0.00098*MW.  Replace with a trained model (e.g. Liigand 2020) if available."""
    from rdkit import Chem
    from rdkit.Chem import Crippen, Lipinski, Descriptors
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return None
    lp = Crippen.MolLogP(m); mw = Descriptors.MolWt(m)
    return dict(logie=0.379 + 0.653 * lp - 0.206 * Lipinski.NumHAcceptors(m) - 0.264 * Lipinski.NumHDonors(m) + 0.00098 * mw,
                logp=lp, mw=mw)


def matrix_factor(cfg: dict, matrix: str) -> float:
    M = cfg["semiquant"]["matrices"]
    if matrix not in M:
        raise KeyError(f"matrix '{matrix}' not defined in config.yaml semiquant.matrices ({list(M)})")
    return float(M[matrix]["factor"])


def calibrants_table(path):
    from .common import read_table
    c = read_table(path)
    c.columns = [x.lower().strip() for x in c.columns]
    return c


def semiquantify(top: pd.DataFrame, cfg: dict, matrix: str, calibrants: pd.DataFrame | None = None, polarity="positive",
                 area_col="area") -> pd.DataFrame:
    """top: one row per feature/candidate (needs smiles, rt_apex, area, max_int).  Adds conc_ug_mL_extract, conc_final, bounds."""
    sq = cfg["semiquant"]
    mode = sq.get("mode", "single")
    h = top.copy()
    if mode == "none" or h.empty:
        return h
    h["_area"] = h[area_col] if area_col in h.columns else h["max_int"]
    conc = np.full(len(h), np.nan); how = np.array([""] * len(h), dtype=object); calname = np.array([""] * len(h), dtype=object)

    if mode == "single":
        ref = sq["single_reference"]
        conc = ref["conc_ug_mL"] * h["_area"].values / ref["area"]
        how[:] = f"single ref ({ref['name']})"
    elif mode in ("logie", "curve"):
        if calibrants is None or calibrants.empty:
            raise ValueError("calibrants.csv required for logie/curve modes")
        cal = calibrants.copy()
        if "polarity" in cal.columns:
            cal = cal[(cal["polarity"].astype(str).str.lower() == polarity) | cal["polarity"].isna()]
        # collapse to single-level calibrants (conc, area, rt, smiles) for logIE
        cal = cal.dropna(subset=["area", "conc_ug_ml"]) if "conc_ug_ml" in cal.columns else cal
        feats = []
        for _, r in cal.iterrows():
            e = logie_estimate(r.get("smiles"));
            if e: feats.append(dict(name=r["name"], rt=float(r["rt"]), area=float(r["area"]), conc=float(r["conc_ug_ml"]), **e))
        C = pd.DataFrame(feats)
        if C.empty:
            raise ValueError("no usable calibrants (need name, smiles, rt, area, conc_ug_ml)")
        # --- curve mode first: compounds with >= 3 levels get their own regression ---
        curve_fit = {}
        if mode == "curve" and "name" in calibrants.columns:
            for nm, g in calibrants.dropna(subset=["area", "conc_ug_ml"]).groupby("name"):
                if g["conc_ug_ml"].nunique() >= 3:
                    k, b = np.polyfit(g["conc_ug_ml"].astype(float), g["area"].astype(float), 1)
                    curve_fit[str(nm).lower()] = (k, b)
        sd = {k: (C[k].std(ddof=0) or 1.0) for k in ("logie", "rt", "mw", "logp")}
        mu = {k: C[k].mean() for k in ("logie", "rt", "mw", "logp")}
        for i, r in enumerate(h.itertuples()):
            nm = str(getattr(r, "name", "")).lower()
            if nm in curve_fit:
                k, b = curve_fit[nm]; conc[i] = max((r._area - b) / k, 0); how[i] = "own curve"; calname[i] = nm; continue
            e = logie_estimate(getattr(r, "smiles", None))
            if e is None:
                continue
            keys = ("logie", "rt", "mw", "logp")
            zc = np.array([((e[k] if k != "rt" else r.rt_apex) - mu[k]) / sd[k] for k in keys])
            Z = np.column_stack([(C[k].values - mu[k]) / sd[k] for k in keys])
            d = np.sqrt(((Z - zc) ** 2).sum(axis=1))          # 'elucidation distance'
            j = int(np.argmin(d)); c = C.iloc[j]
            conc[i] = c["conc"] * (r._area / c["area"]) * 10 ** (c["logie"] - e["logie"])
            how[i] = f"logIE-corrected (elucidation dist {d[j]:.2f})"; calname[i] = c["name"]
    else:
        raise ValueError(f"unknown semiquant mode {mode}")

    f = matrix_factor(cfg, matrix)
    unc = float(sq.get("uncertainty_factor", 10.0))
    h["sq_method"] = how; h["sq_calibrant"] = calname
    h["conc_ug_mL_extract"] = conc
    h[f"conc_{matrix}"] = conc * f
    h[f"conc_{matrix}_low"] = conc * f / unc
    h[f"conc_{matrix}_high"] = conc * f * unc
    h["sq_unit"] = cfg["semiquant"]["matrices"][matrix].get("unit", "")
    loq = sq.get("loq_ug_mL_extract")
    if loq:
        h["below_loq"] = h["conc_ug_mL_extract"] < loq
    return h.drop(columns=["_area"])


def measure_calibrants(run, calibrants: pd.DataFrame, db_masses=None, ppm=5.0, rt_tol=0.3, adduct_shift=None):
    """Fill 'area' for calibrants (name, smiles, rt, conc_ug_ml, [mz]) from a standard mzML run (EIC integration)."""
    from .features import eic_area
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    shift = adduct_shift if adduct_shift is not None else (1.007276 if run.polarity == "positive" else -1.007276)
    out = []
    for r in calibrants.itertuples():
        mz = getattr(r, "mz", np.nan)
        if not np.isfinite(mz):
            m = Chem.MolFromSmiles(r.smiles)
            mz = Descriptors.ExactMolWt(m) + shift
        a, rt, apex, snr = eic_area(run, float(mz), float(r.rt), ppm, rt_tol)
        out.append(dict(name=r.name, mz=mz, rt_found=rt, area=a, apex=apex, snr=snr))
    return pd.DataFrame(out)
