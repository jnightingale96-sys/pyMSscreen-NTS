"""Identification-confidence filter for the screening candidates:  keep Level 1/2a, keep a Level 2b only when it is retention-time plausible OR supported by measured
MassBank-EU fragments.

RT plausibility without authentic standards: the run's own robust Level 1/2a identifications (library cosine >= `anchor_min_cos`, homologous series such as PEG/PPG/surfactants
excluded because one library entry matches every chain length) anchor a robust (Huber) map  NORMAN Pred_RTI -> retention time; the tolerance is `rt_tol_k` x the leave-one-out
RMSE of the anchors (floor `rt_tol_floor_min`).  A model is only built for a polarity with >= `min_anchors` anchors; without one, RT cannot veto a candidate.
MassBank-EU support: >= `mb_min_matched` of the fragments MassBank recorded for that structure are present in the feature's MS2 (and >= `mb_min_frac` of the reference list)."""
from __future__ import annotations
import json, logging, os, re
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")
SERIES = re.compile(r"glycol|PEG|PPG|ethoxy|polyoxy|\bEO\b|oxa|ether|alkyl|dodecyl|hexadecyl|tetradecyl|octadecyl|oligomer", re.I)


def rti_lookup(cfg: dict) -> pd.DataFrame:
    """ik1 -> NORMAN Pred_RTI (pos/neg) and the RTI model's applicability-domain flag ('Covered by Model' -> ad_* True; Aalizadeh 2021).  Sources: the optional updated NORMAN
    workbook (paths.norman_workbook / id_filter.norman_workbook, takes precedence) then SusDat.  Cached next to the spectral library."""
    from .common import resolve
    cache = os.path.join(os.path.dirname(resolve(cfg, cfg["paths"]["library_db"])), "rti_lookup_v2.pkl")
    if os.path.exists(cache): return pd.read_pickle(cache)
    def frame(df):
        ik = df["StdInChIKey"].astype(str).str.split("-").str[0]
        return pd.DataFrame({"ik1": ik, "pred_rti_pos": pd.to_numeric(df["Pred_RTI_Positive_ESI"], errors="coerce"), "pred_rti_neg": pd.to_numeric(df["Pred_RTI_Negative_ESI"], errors="coerce"),
                             "ad_pos": df["Uncertainty_RTI_pos"].fillna("").astype(str).str.contains("Covered"), "ad_neg": df["Uncertainty_RTI_neg"].fillna("").astype(str).str.contains("Covered")})
    cols = ["StdInChIKey", "Pred_RTI_Positive_ESI", "Pred_RTI_Negative_ESI", "Uncertainty_RTI_pos", "Uncertainty_RTI_neg"]
    parts = []
    wb = resolve(cfg, cfg["pec_screen"]["id_filter"].get("norman_workbook"))
    if wb and os.path.exists(wb): parts.append(frame(pd.read_excel(wb, dtype=str, usecols=cols)))
    parts.append(frame(pd.read_csv(resolve(cfg, cfg["paths"]["susdat"]), sep="\t", engine="python", dtype=str, usecols=cols)))
    d = pd.concat(parts).drop_duplicates("ik1"); d.to_pickle(cache); return d


def fit_rt_map(cands: pd.DataFrame, rti: pd.DataFrame, polarity: str, ic: dict) -> dict | None:
    from sklearn.linear_model import HuberRegressor
    col = "pred_rti_pos" if polarity == "positive" else "pred_rti_neg"; adc = "ad_pos" if polarity == "positive" else "ad_neg"
    c = cands[cands.polarity == polarity].merge(rti[["ik1", col] + ([adc] if adc in rti.columns else [])], on="ik1", how="left")
    a = c[c.level.isin(["1", "2a"]) & (c.lib_score.astype(float) >= float(ic["anchor_min_cos"])) & c[col].notna() & ~c.name.map(lambda n: bool(SERIES.search(str(n))))]
    if ic.get("use_applicability_domain", True) and adc in a.columns:
        ina = a[a[adc].fillna(False).astype(bool)]
        if len(ina) >= int(ic["min_anchors"]): a = ina                       # anchors inside the RTI model's applicability domain give the tighter, honest band
    if len(a) < int(ic["min_anchors"]): return None
    X, y = a[[col]].values, a["rt"].astype(float).values
    loo = []
    for i in range(len(a)):
        m = np.ones(len(a), bool); m[i] = False
        loo.append(float(HuberRegressor().fit(X[m], y[m]).predict(X[[i]])[0]) - y[i])
    h = HuberRegressor().fit(X, y); rmse = float(np.sqrt(np.mean(np.square(loo))))
    return dict(polarity=polarity, col=col, intercept=float(h.intercept_), slope=float(h.coef_[0]), rmse=rmse, n_anchors=len(a),
                tol=max(float(ic["rt_tol_floor_min"]), float(ic["rt_tol_k"]) * rmse), frac_within_1min=float(np.mean(np.abs(loo) <= 1.0)))


def apply_id_filter(cands: pd.DataFrame, cfg: dict, out: str | None = None, rti: pd.DataFrame | None = None, mbeu: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame, list]:
    """Returns (kept, removed, models).  Adds rt_err_min, rt_status, mb_eu_listed, mb_eu_matched, mb_eu_ref, mb_eu_support, id_filter to every candidate."""
    ic = cfg["pec_screen"]["id_filter"]
    C = cands.copy()
    rti = rti if rti is not None else rti_lookup(cfg)
    models = {}
    C["rt_err_min"] = np.nan; C["rt_status"] = "no model"; C["rt_in_ad"] = True
    for pol in sorted(C.polarity.unique()):
        mdl = fit_rt_map(C, rti, pol, ic)
        if mdl is None: log.warning("id filter: < %s RT anchors in %s mode - RT cannot veto candidates there", ic["min_anchors"], pol); continue
        models[pol] = mdl
        log.info("RT map (%s): RT = %.2f + %.4f x RTI  | %d anchors, leave-one-out RMSE %.2f min (%.0f%% within 1 min) -> tolerance +-%.1f min", pol, mdl["intercept"], mdl["slope"],
                 mdl["n_anchors"], mdl["rmse"], 100 * mdl["frac_within_1min"], mdl["tol"])
        m = C.polarity == pol
        adc = "ad_pos" if pol == "positive" else "ad_neg"
        r = C[m].merge(rti[["ik1", mdl["col"]] + ([adc] if adc in rti.columns else [])], on="ik1", how="left")
        pred = mdl["intercept"] + mdl["slope"] * r[mdl["col"]].values
        err = pred - r["rt"].astype(float).values
        inad = r[adc].fillna(False).astype(bool).values if (adc in r.columns and ic.get("use_applicability_domain", True)) else np.ones(len(r), bool)
        C.loc[m, "rt_err_min"] = err; C.loc[m, "rt_in_ad"] = inad
        C.loc[m, "rt_status"] = np.where(np.isnan(err), "no RTI", np.where(~inad, "outside AD", np.where(np.abs(err) <= mdl["tol"], "plausible", "IMPLAUSIBLE")))
    # ---- MassBank-EU fragment support -----------------------------------------------------------------------------------------------
    C["mb_eu_listed"] = False; C["mb_eu_matched"] = 0; C["mb_eu_ref"] = 0; C["mb_eu_support"] = False
    mp = ic.get("massbank_eu")
    if mbeu is None and mp:
        from .common import resolve
        from .massbank_eu import load_massbank_eu
        p = resolve(cfg, mp)
        mbeu = load_massbank_eu(p) if os.path.exists(p) else None
    if mbeu is not None:
        from .massbank_eu import fragment_support
        from .ms2 import string_to_peaks
        for i, r in C.iterrows():
            if r.ik1 not in mbeu.index: continue
            ref = mbeu.at[r.ik1, "pos_frags" if r.polarity == "positive" else "neg_frags"]
            C.at[i, "mb_eu_listed"] = True
            mz, _ = string_to_peaks(r.ms2_peaks)
            n, nref = fragment_support(mz, ref, r.mz, float(cfg["library"]["frag_tol_da"]), float(cfg["library"]["frag_tol_ppm"]))
            C.at[i, "mb_eu_matched"] = n; C.at[i, "mb_eu_ref"] = nref
            C.at[i, "mb_eu_support"] = bool(nref and n >= min(int(ic["mb_min_matched"]), nref) and n >= 2 and n / nref >= float(ic["mb_min_frac"]))
    # ---- decision --------------------------------------------------------------------------------------------------------------------
    keep, why = [], []
    for r in C.itertuples():
        if r.level in ("1", "2a"): keep.append(True); why.append("Level 1/2a")
        elif r.mb_eu_support: keep.append(True); why.append("%s + MassBank-EU fragment support" % r.level)
        elif r.rt_status == "plausible": keep.append(True); why.append("%s + RT plausible" % r.level)
        elif r.rt_status == "IMPLAUSIBLE": keep.append(False); why.append("%s removed: RT implausible (%.1f min off)" % (r.level, r.rt_err_min))
        elif r.rt_status in ("no RTI", "outside AD") and not ic.get("keep_2b_without_rti", False): keep.append(False); why.append("%s removed: RT not predictable (%s) and no MassBank-EU support" % (r.level, r.rt_status))
        else: keep.append(True); why.append("%s kept: RT cannot be judged (no model)" % r.level)
    C["id_filter"] = why
    kept, removed = C[keep].reset_index(drop=True), C[[not k for k in keep]].reset_index(drop=True)
    if out:
        removed.drop(columns=["ms2_peaks"], errors="ignore").to_csv(os.path.join(out, "candidates_removed_by_id_filter.csv"), index=False)
        json.dump(models, open(os.path.join(out, "rt_map.json"), "w"), indent=1)
    log.info("id filter: %d candidates -> %d kept (%s)", len(C), len(kept), kept.id_filter.value_counts().to_dict())
    return kept, removed, models
