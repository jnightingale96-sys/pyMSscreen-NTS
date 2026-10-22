"""Identification confidence (Schymanski et al. 2014 scheme, with the project's 2a/2b split).

  Level 1   authentic standard on this method: standard RT within tolerance AND experimental MS2 matches a library spectrum
  Level 2a  experimental library spectrum match (cosine + >= N matched fragments), RT not contradicting
  Level 2b  diagnostic in-silico fragmentation evidence: >= N explained product ions, decoy-corrected score, CLEAR top
            candidate for the feature, RT not contradicting
  Level 3   tentative: MS2 present with >= 2 explained product ions (or weak library match) AND RT plausible
  Level 4   unequivocal molecular formula: accurate mass (+ isotope pattern not contradicting)
  Level 5   exact mass of interest only (feature without any candidate)
Evidence columns are kept so every assignment is auditable."""
from __future__ import annotations
import numpy as np
import pandas as pd

ORDER = {"1": 0, "2a": 1, "2b": 2, "3": 3, "4": 4, "5": 5}


def mark_clear_top(h: pd.DataFrame, margin=0.15, min_score=0.2) -> pd.DataFrame:
    """Per feature: is the best in-silico-scored planar structure clearly ahead of the runner-up?"""
    h = h.copy()
    # MetFrag 'Score' is rank-normalised (the top of N candidates is always 1.0, also when N == 1), so it is only usable to
    # compare >= 2 scored candidates of one feature; otherwise use the fraction of used peaks it explained (damped).
    bi = h["bi_score"].fillna(0.0) if "bi_score" in h.columns else pd.Series(0.0, index=h.index)
    if "mf_score" in h.columns:
        n_mf = h.groupby("feature_id")["mf_score"].transform("count")
        frac = (h["mf_expl"] / h["mf_used"].replace(0, np.nan)) if "mf_used" in h.columns else pd.Series(np.nan, index=h.index)
        mf = np.where(n_mf >= 2, h["mf_score"], 0.6 * frac.fillna(0.0))
        mf = pd.Series(mf, index=h.index).fillna(0.0)
        h["ins_score"] = np.maximum(bi, mf)
    else:
        h["ins_score"] = bi
    h["clear_top"] = False; h["stereo_ambiguous"] = False
    single_min = 0.5
    for fid, g in h.groupby("feature_id"):
        by_ik1 = g.groupby("ik1")["ins_score"].max().sort_values(ascending=False)
        top = by_ik1.iloc[0]
        second = by_ik1.iloc[1] if len(by_ik1) > 1 else 0.0
        clear = (top >= min_score) and ((top - second) >= margin) if len(by_ik1) > 1 else (top >= single_min)
        topk = by_ik1.index[0]
        idx = g.index[g["ik1"] == topk]
        h.loc[idx, "clear_top"] = bool(clear)
        h.loc[idx, "stereo_ambiguous"] = g.loc[idx, "inchikey"].nunique() > 1
    return h


def assign_levels(h: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    L = cfg["levels"]
    h = h.copy()
    for c, d in (("lib_score", 0.0), ("lib_matched", 0), ("bi_expl_n", 0), ("bi_score", 0.0), ("std_rt_diff", np.nan),
                 ("rt_agree", False), ("rt_conflict", False), ("rt_available", False), ("ms2_n_peaks", 0), ("ms2_good_int", False),
                 ("clear_top", False), ("ins_score", 0.0), ("iso_score", np.nan), ("n_adducts_seen", 1), ("broad_flag", False),
                 ("lib_inhouse", False), ("mf_expl", np.nan), ("bi_expl_adj_n", 0.0), ("background_flag", False)):
        if c not in h.columns:
            h[c] = d
    lv, why = [], []
    for r in h.itertuples():
        has_ms2 = (r.ms2_n_peaks or 0) >= L["min_ms2_peaks"] and bool(r.ms2_good_int or L.get("ignore_ms2_intensity", False))
        lib_ok = (r.lib_score or 0) >= L["lib_cosine"] and (r.lib_matched or 0) >= L["lib_min_matched"]
        std_ok = pd.notna(r.std_rt_diff) and abs(r.std_rt_diff) <= L["level1_rt_tol"]
        rt_bad = bool(r.rt_conflict)
        rt_good = bool(r.rt_agree) or (not r.rt_available)
        expl = max(getattr(r, "bi_expl_adj_n", 0) or 0, r.mf_expl if pd.notna(r.mf_expl) else 0)
        ins_ok = has_ms2 and expl >= L["min_explained_peaks_2b"] and (r.ins_score or 0) >= L["insilico_min_score"] and bool(r.clear_top)
        n = []
        if std_ok and lib_ok and (r.lib_inhouse or not L["level1_require_inhouse_ms2"]):
            level = "1"; n.append("standard RT + MS2 match")
        elif lib_ok and not rt_bad:
            level = "2a"; n.append(f"library cosine {r.lib_score:.2f}/{int(r.lib_matched)} peaks")
            if r.rt_agree: n.append("RT agrees")
        elif ins_ok and not rt_bad:
            level = "2b"; n.append(f"{int(expl)} explained ions, clear top")
            if r.rt_agree: n.append("RT agrees")
        elif has_ms2 and (expl >= L["min_explained_peaks_3"] or (r.lib_score or 0) >= L["lib_cosine_weak"]) and rt_good:
            level = "3"; n.append("MS2 supports tentatively")
        else:
            level = "4"; n.append("mass/formula only")
        lowmass = pd.notna(getattr(r, "exact_mass", np.nan)) and r.exact_mass < L.get("lowmass_da", 120)
        bg = bool(getattr(r, "background_flag", False))
        if level in ("2b", "3") and (lowmass or bg):
            n.append("capped at 4: " + ("background-like ion" if bg else f"neutral mass < {L.get('lowmass_da', 120)} Da (too few diagnostic fragments)"))
            level = "4"
        if rt_bad and level in ("4", "3"):
            n.append("RT CONFLICT with prediction")
        if r.broad_flag: n.append("broad/continuous signal")
        if std_ok and level != "1": n.append("standard RT matches (no MS2 proof)")
        lv.append(level); why.append("; ".join(n))
    h["level"] = lv; h["level_basis"] = why
    # numeric confidence 0-100 (ranking aid, not a statistical probability)
    ppm_tol = cfg["level4"]["ppm"]
    s = 10 * (1 - np.clip(np.abs(h["ppm_error"].values) / ppm_tol, 0, 1))
    s += np.where(h["iso_score"].notna(), 10 * h["iso_score"].fillna(0), 5)
    s += 5 * np.clip(h["n_adducts_seen"].values - 1, 0, 2)
    s += np.where(h["rt_agree"], 15, np.where(h["rt_conflict"], -15, 5))
    s += 40 * np.where(h["lib_matched"] >= 2, h["lib_score"].fillna(0), 0)
    s += 20 * np.clip(h["ins_score"].fillna(0) / 0.5, 0, 1)
    s += np.where(h["std_rt_diff"].abs() <= L["level1_rt_tol"], 20, 0)
    s -= np.where(h["broad_flag"], 10, 0)
    s -= 2 * np.log2(np.maximum(h.get("n_candidates_for_feature", 1), 1))
    h["confidence_score"] = np.clip(s, 0, 100).round(1)
    h["level_rank"] = h["level"].map(ORDER)
    h = h.sort_values(["feature_id", "level_rank", "confidence_score"], ascending=[True, True, False]).reset_index(drop=True)
    h["is_top_candidate"] = ~h.duplicated("feature_id")
    return h
