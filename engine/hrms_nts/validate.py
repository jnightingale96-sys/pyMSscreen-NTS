"""Validation of a results folder against ground truth.

expected.csv  (compounds KNOWN to be present):  name, [smiles | inchikey], [sample], [adduct], [rt], [min_level], [true_conc_ug_ml]
                 true_conc_ug_ml = concentration in the measured EXTRACT (ug/mL) -> semi-quant accuracy (fold error) is reported
absent.csv    (compounds KNOWN to be absent, e.g. not spiked / blank matrix): name, [smiles | inchikey], [sample]
  sample = file stem ('*' or empty = any sample); rt = expected RT in min (matching window --rt-tol); min_level default 4.

Reports per expected compound: detected, best level, ppm error, RT difference, whether it was the top-ranked candidate of its feature,
and summary metrics: recovery (>= Level N), median |ppm|, median |dRT|, top-1 rate, and - with absent.csv - the false identification
rate at each level.  Level order: 1 < 2a < 2b < 3 < 4."""
from __future__ import annotations
import os, logging
import numpy as np
import pandas as pd
from .common import read_table

log = logging.getLogger("hrms_nts")
ORDER = {"1": 0, "2a": 1, "2b": 2, "3": 3, "4": 4, "5": 5}


def _ik1(row):
    ik = str(row.get("inchikey", "") or "")
    if len(ik) > 8 and ik != "nan":
        return ik.split("-")[0]
    smi = row.get("smiles")
    if isinstance(smi, str) and smi:
        try:
            from rdkit import Chem
            from rdkit.Chem.inchi import MolToInchiKey
            return MolToInchiKey(Chem.MolFromSmiles(smi)).split("-")[0]
        except Exception:
            return None
    return None


def _match(cands: pd.DataFrame, row, rt_tol):
    ik1 = _ik1(row)
    c = cands[cands["ik1"].astype(str) == ik1] if ik1 else cands[cands["name"].astype(str).str.lower() == str(row["name"]).lower()]
    s = row.get("sample")
    if isinstance(s, str) and s.strip() not in ("", "*"):
        c = c[c["sample"].astype(str) == s.strip()]
    a = row.get("adduct")
    if isinstance(a, str) and a.strip():
        c = c[c["adduct"].astype(str) == a.strip()]
    rt = row.get("rt")
    if pd.notna(rt) and str(rt) != "":
        c = c[(c["rt_apex"] - float(rt)).abs() <= rt_tol]
    return c


def validate(results_dir: str, expected: str, absent: str | None = None, rt_tol: float = 0.3, out_dir: str | None = None):
    path = os.path.join(results_dir, "all_candidates.csv")
    cands = pd.read_csv(path, low_memory=False)
    cands["level"] = cands["level"].astype(str)
    cands["level_rank"] = cands["level"].map(ORDER)
    exp = read_table(expected)
    rows = []
    for _, r in exp.iterrows():
        m = _match(cands, r, rt_tol)
        if m.empty:
            rows.append(dict(name=r["name"], sample=r.get("sample", "*"), detected=False, best_level=None, ppm_error=np.nan, rt_apex=np.nan,
                             rt_diff=np.nan, top_candidate=None, max_int=np.nan, lib_score=np.nan, min_level_required=str(r.get("min_level", "4")),
                             meets_min_level=False, n_candidate_rows=0))
            continue
        b = m.sort_values(["level_rank", "confidence_score"], ascending=[True, False]).iloc[0]
        need = str(r.get("min_level", "4")) if pd.notna(r.get("min_level", np.nan)) else "4"
        rows.append(dict(name=r["name"], sample=b["sample"], detected=True, feature_id=b["feature_id"], best_level=b["level"], ppm_error=b["ppm_error"], rt_apex=b["rt_apex"],
                         rt_diff=(b["rt_apex"] - float(r["rt"])) if pd.notna(r.get("rt", np.nan)) else np.nan, top_candidate=bool(b.get("is_top_candidate", False)),
                         max_int=b["max_int"], lib_score=b.get("lib_score", np.nan), min_level_required=need,
                         meets_min_level=ORDER.get(b["level"], 9) <= ORDER.get(need, 4), n_candidate_rows=len(m)))
    rep = pd.DataFrame(rows)
    # ---- semi-quant accuracy against spiked concentrations ------------------------------------------------------------------
    tp = os.path.join(results_dir, "top_candidates_per_feature.csv")
    rep["true_conc_ug_ml"] = pd.to_numeric(exp.get("true_conc_ug_ml", pd.Series([np.nan] * len(exp))), errors="coerce").values
    rep["pred_conc_ug_ml"] = np.nan; rep["fold_error"] = np.nan
    if os.path.exists(tp) and rep["true_conc_ug_ml"].notna().any():
        top = pd.read_csv(tp, low_memory=False)
        if "conc_ug_mL_extract" in top.columns and "feature_id" in rep.columns:
            m = top.drop_duplicates("feature_id").set_index("feature_id")["conc_ug_mL_extract"]
            rep["pred_conc_ug_ml"] = rep["feature_id"].map(m)
            ok = rep["pred_conc_ug_ml"].gt(0) & rep["true_conc_ug_ml"].gt(0)
            rep.loc[ok, "fold_error"] = rep.loc[ok, "pred_conc_ug_ml"] / rep.loc[ok, "true_conc_ug_ml"]
    n = len(rep); det = rep[rep.detected]
    summ = {"n_expected": n, "detected_any_level(>=4)": int(rep.detected.sum()), "recovery_pct": round(100 * rep.detected.mean(), 1) if n else np.nan}
    for lv in ("3", "2b", "2a", "1"):
        summ[f"reached_level<={lv}"] = int((det["best_level"].map(ORDER) <= ORDER[lv]).sum())
    summ["meets_required_min_level"] = int(rep.meets_min_level.sum())
    summ["median_abs_ppm"] = round(det["ppm_error"].abs().median(), 2) if len(det) else np.nan
    summ["median_abs_rt_diff_min"] = round(det["rt_diff"].abs().median(), 3) if det["rt_diff"].notna().any() else np.nan
    summ["top1_rate_pct"] = round(100 * det["top_candidate"].astype(float).mean(), 1) if len(det) else np.nan
    fe = rep["fold_error"].dropna()
    if len(fe):
        lf = np.abs(np.log10(fe))
        summ["semiquant_n"] = int(len(fe)); summ["semiquant_median_fold_error"] = round(float(10 ** lf.median()), 2)
        summ["semiquant_geomean_bias(pred/true)"] = round(float(10 ** np.log10(fe).mean()), 2)
        summ["semiquant_within_3x_pct"] = round(100 * float((lf <= np.log10(3)).mean()), 1)
        summ["semiquant_within_10x_pct"] = round(100 * float((lf <= 1).mean()), 1)
    if absent:
        ab = read_table(absent); fp = []
        for _, r in ab.iterrows():
            m = _match(cands, r, 99)
            fp.append(dict(name=r["name"], present_as=None if m.empty else m.sort_values("level_rank").iloc[0]["level"]))
        fp = pd.DataFrame(fp); na = len(fp)
        for lv in ("4", "3", "2b", "2a", "1"):
            k = int((fp["present_as"].map(ORDER) <= ORDER[lv]).sum()) if na else 0
            summ[f"absent_flagged_at_<=L{lv}"] = f"{k}/{na}"
        rep_absent = fp
    else:
        rep_absent = None
    out_dir = out_dir or results_dir
    rep.to_csv(os.path.join(out_dir, "validation_report.csv"), index=False)
    pd.Series(summ).to_csv(os.path.join(out_dir, "validation_summary.csv"), header=["value"])
    if rep_absent is not None:
        rep_absent.to_csv(os.path.join(out_dir, "validation_absent.csv"), index=False)
    return rep, summ
