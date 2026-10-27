"""PEC-vs-MEC screening of every non-target identification (Level >= `min_level`).

  water mzML (field-resolved)  --identify-->  candidates (Level 1/2a/2b)
  candidates --EIC re-integration--> semi-quant in every water injection (blank-checked)  -> C_water per field (ng/L)
  candidates --EIC + MS2 match---->  the SAME compound in every grain injection (RT-shift corrected)  -> MEC per field (ng/g dw, semi-quant)
  properties (table / estimates) --fate model--> k_fruit (range if properties missing)  -> PEC = k_fruit x C_water x 0.005 (ng/g dw)
  PEC vs MEC per field and per chemical -> class: water-explained | MEC >> PEC | PEC >> MEC | not found | not applicable

Semi-quantification is the project's single-reference method everywhere (water AND grain, no literature values):
    conc(ug/mL extract) = ref_conc x area / area_ref ;  water ng/L = conc x 1000 x 120 ;  grain ng/g dw = conc x 0.5/3 x 1000
with area_ref measured by THIS integrator on the carbamazepine standard (units therefore cancel)."""
from __future__ import annotations
import logging, os, pickle, time
import numpy as np
import pandas as pd

from . import mzml_io, common
from .common import resolve, guess_polarity
from .library import cosine
from .ms2 import MS2Index, string_to_peaks
from .levels import ORDER

log = logging.getLogger("hrms_nts")
LEVEL_SET = {"1": ["1"], "2a": ["1", "2a"], "2b": ["1", "2a", "2b"], "3": ["1", "2a", "2b", "3"], "4": ["1", "2a", "2b", "3", "4"]}
STATUS_RANK = {"MS2 confirmed": 3, "MS2 weak": 2, "MS1 only": 1, "MS2 conflict": 0, "": -1}


# ======================================================================================================================================
# integration (the project's smoothed, baseline-corrected EIC integration - identical for standard and samples)
# ======================================================================================================================================
def integrate_eic(rt, y, rt0, half=0.2, search=0.15, floor=3e4):
    """Returns (area counts.min, apex rt, raw apex in window).  area = 0 when the smoothed apex is below `floor`."""
    rt = np.asarray(rt); y = np.asarray(y, float); half, search, floor = float(half), float(search), float(floor)
    m = (rt >= rt0 - half - search) & (rt <= rt0 + half + search)
    if m.sum() < 5 or y[m].max() <= 0:
        return 0.0, np.nan, 0.0
    r, yy = rt[m], y[m]; sm = np.convolve(yy, np.ones(3) / 3, mode="same")
    c = int(np.argmax(np.where(np.abs(r - rt0) <= search, sm, -1)))
    if sm[c] < floor:
        return 0.0, float(r[c]), float(yy[c])
    l = c
    while l > 0 and sm[l - 1] <= sm[l] and sm[l - 1] > 0.02 * sm[c]: l -= 1
    rr = c
    while rr < len(sm) - 1 and sm[rr + 1] <= sm[rr] and sm[rr + 1] > 0.02 * sm[c]: rr += 1
    base = np.interp(r[l:rr + 1], [r[l], r[rr]], [yy[l], yy[rr]])
    return float(np.trapezoid(np.clip(yy[l:rr + 1] - base, 0, None), r[l:rr + 1])), float(r[c]), float(yy[c])


def reference_response(cfg: dict) -> dict:
    """Area per ug/mL of the single reference (carbamazepine) measured on the standard mzML with the same integrator."""
    from scipy.signal import find_peaks
    rc = cfg["pec_screen"]["reference"]
    path = resolve(cfg, rc.get("mzml"))
    if path and os.path.exists(path):
        run = mzml_io.read_for_config(path, common.config_for_file(cfg, path), "positive")
        rt, y = run.eic(float(rc["mz"]), cfg["level4"]["ppm"])
        rt0 = rc.get("rt")
        if not rt0:
            pk, _ = find_peaks(np.convolve(y, np.ones(3) / 3, mode="same"), prominence=0.2 * y.max())
            rt0 = float(rt[pk[np.argmax(y[pk])]]) if len(pk) else float(rt[np.argmax(y)])
        area, rtf, apex = integrate_eic(rt, y, float(rt0), half=0.25, search=0.10, floor=0.0)
        if area <= 0:
            raise RuntimeError(f"reference {rc['name']} not found in {path}")
        return dict(name=rc["name"], conc_ug_mL=float(rc["conc_ug_mL"]), area=area, area_per_ug_mL=area / float(rc["conc_ug_mL"]), rt=rtf, apex=apex, source=path)
    if rc.get("area_per_ug_mL"):
        log.warning("reference standard mzML not found - using pec_screen.reference.area_per_ug_mL = %s", rc["area_per_ug_mL"])
        return dict(name=rc["name"], conc_ug_mL=float(rc["conc_ug_mL"]), area=np.nan, area_per_ug_mL=float(rc["area_per_ug_mL"]), rt=np.nan, apex=np.nan, source="config")
    raise RuntimeError("no reference: set pec_screen.reference.mzml (standard injection) or .area_per_ug_mL")


# ======================================================================================================================================
# sample map
# ======================================================================================================================================
def read_sample_map(path: str, cfg: dict | None = None) -> pd.DataFrame:
    """CSV columns: file (path, or stem if --dir given) | matrix (water|grain) | field | role (sample|blank|exclude) | polarity (optional) | note"""
    from .common import read_table
    m = read_table(path)
    m.columns = [c.strip().lower() for c in m.columns]
    for c in ("file", "matrix", "role"):
        if c not in m.columns: raise ValueError(f"sample map needs a '{c}' column")
    for c, d in (("field", ""), ("polarity", ""), ("note", ""), ("scale", 1.0)):
        if c not in m.columns: m[c] = d
    m["scale"] = pd.to_numeric(m["scale"], errors="coerce").fillna(1.0)      # optional per-injection multiplier (extract-preparation correction)
    base = os.path.dirname(os.path.abspath(path))
    m["path"] = [p if os.path.isabs(p) else os.path.join(base, p) for p in m["file"].astype(str)]
    m["stem"] = [os.path.splitext(os.path.basename(p))[0] for p in m["path"]]
    m["matrix"] = m["matrix"].str.lower().str.strip(); m["role"] = m["role"].str.lower().str.strip()
    m["field"] = m["field"].fillna("").astype(str).str.strip(); m["polarity"] = m["polarity"].fillna("").astype(str).str.lower().str.strip()
    m["note"] = m["note"].fillna("").astype(str)
    bad = m[~m["matrix"].isin(["water", "grain", "soil"]) | ~m["role"].isin(["sample", "blank", "exclude"])]
    if len(bad): raise ValueError(f"bad matrix/role in sample map rows: {bad.index.tolist()}")
    miss = [p for p in m.loc[m.role != "exclude", "path"] if not os.path.exists(p)]
    if miss: raise FileNotFoundError("mzML not found: " + "; ".join(miss[:5]))
    return m


# ======================================================================================================================================
# step 1 - identification of the water samples (cached) and the candidate list
# ======================================================================================================================================
def identify_water(smap: pd.DataFrame, cfg: dict, out: str, min_level: str, reuse=True, ctx=None) -> dict:
    from .pipeline import Context, process_sample
    levels = LEVEL_SET[min_level]
    d = os.path.join(out, "_ident"); os.makedirs(d, exist_ok=True)
    res = {}
    todo = smap[(smap.matrix == "water") & (smap.role == "sample")]
    for r in todo.itertuples():
        pk = os.path.join(d, r.stem + ".pkl")
        if reuse and os.path.exists(pk) and os.path.getmtime(pk) > os.path.getmtime(r.path):
            cached = pickle.load(open(pk, "rb"))
            if cached.get("has_level4"):                      # older caches dropped Level 4 rows
                res[r.stem] = cached; log.info("identification reused: %s", r.stem); continue
        if ctx is None:
            ctx = Context(cfg)
        t = process_sample(r.path, ctx, "water", r.polarity or "auto", out, None, False)
        top = t["top"]
        keep = top[top["level"].isin(["1", "2a", "2b", "3", "4"])].copy() if not top.empty else top
        keep = keep[keep["smiles"].fillna("").astype(str) != ""] if not keep.empty else keep
        res[r.stem] = dict(polarity=t["polarity"], top=keep, n_features=len(t["features"]), has_level4=True)
        pickle.dump(res[r.stem], open(pk, "wb"))
    return res


def build_candidates(ident: dict, cfg: dict, min_level: str) -> pd.DataFrame:
    """One row per chemical (InChIKey block 1 x polarity): best level, RT/m/z, consensus MS2 of the most intense identification."""
    levels = LEVEL_SET[min_level]; rows = []
    for stem, r in ident.items():
        t = r["top"]
        if t is None or t.empty: continue
        t = t[t["level"].isin(levels)].copy()
        t = t[t["smiles"].fillna("").astype(str) != ""]
        if "background_flag" in t.columns: t = t[~t["background_flag"].fillna(False).astype(bool)]
        t["sample_id"] = stem; t["polarity"] = r["polarity"]
        rows.append(t)
    if not rows: return pd.DataFrame()
    A = pd.concat(rows, ignore_index=True)
    A["level_rank"] = A["level"].map(ORDER)
    out = []
    for (ik1, pol), g in A.groupby(["ik1", "polarity"]):
        g = g.sort_values(["level_rank", "confidence_score"], ascending=[True, False]); b = g.iloc[0]
        near = g[(g["rt_apex"] - b["rt_apex"]).abs() <= 0.3]
        ms2 = near.sort_values("ms2_max_int", ascending=False).iloc[0]
        out.append(dict(cand_id=f"{ik1}_{pol[:3]}", name=b["name"], formula=b["formula"], smiles=b["smiles"], ik1=ik1, inchikey=b["inchikey"], cas=b.get("cas", ""),
                        polarity=pol, adduct=b["adduct"], mz=float(b["calc_mz"]) if pd.notna(b["calc_mz"]) else float(b["mz_apex"]),
                        rt=float(near["rt_apex"].median()), level=b["level"], confidence=float(b["confidence_score"]), level_basis=b["level_basis"],
                        source=b["source"], lib_name=b.get("lib_name", ""), lib_score=float(b.get("lib_score", 0) or 0),
                        n_ident_samples=int(near["sample_id"].nunique()), ms2_peaks=ms2["ms2_peaks"]))
    C = pd.DataFrame(out).sort_values(["polarity", "mz", "rt"]).reset_index(drop=True)
    # collapse different structures assigned to the same feature in different samples (same m/z, same RT): keep the best-level one
    C["_r"] = C["level"].map(ORDER); keep = []
    order = C.sort_values(["_r", "confidence"], ascending=[True, False]).index.values
    mzv, rtv, polv = C["mz"].values.astype(float), C["rt"].values.astype(float), C["polarity"].values
    kept_idx = np.zeros(0, dtype=int)
    for i in order:
        if len(kept_idx):
            k = kept_idx
            if np.any((polv[k] == polv[i]) & (np.abs(mzv[k] - mzv[i]) / mzv[i] * 1e6 <= 5) & (np.abs(rtv[k] - rtv[i]) <= 0.15)): continue
        kept_idx = np.append(kept_idx, i); keep.append(i)
    return C.loc[sorted(keep)].drop(columns="_r").reset_index(drop=True)


# ======================================================================================================================================
# step 2 - quantify every candidate in every injection (water, blanks, grain)
# ======================================================================================================================================
def _ms2_cos(run_idx: MS2Index, c, rt_found, cfg, window=0.10):
    cq, cl = string_to_peaks(c.ms2_peaks)
    sp = run_idx.best_spectrum(c.mz, rt_found, rt_found - window, rt_found + window, cfg["ms2"]["precursor_tol_da"], 0.0, cfg["ms2"]["merge_scans"],
                               cfg["ms2"]["rel_min"], cfg["ms2"]["top_n"]) if rt_found == rt_found else None
    if sp is None or len(sp["mz"]) == 0: return 0.0, 0, 0
    lc = cfg["library"]
    ordr = np.argsort(cq); score, n = cosine(sp["mz"], sp["inten"], cq[ordr], cl[ordr], lc["frag_tol_da"], lc["frag_tol_ppm"])
    return float(score), int(n), int(sp["n_scans"])


def normalise_config(cfg: dict) -> dict:
    """YAML 1.1 reads '3.0e4' (no exponent sign) as a string - coerce every numeric pec_screen setting."""
    pc = cfg["pec_screen"]
    for k in ("water_factor", "grain_factor", "soil_factor", "soil_pec_to_ng_g", "fruit_to_grain_dw", "grain_floor_width_min", "blank_ratio", "rt_search_wide_min", "rt_shift_min_anchors",
              "rt_shift_min_apex", "ms2_cosine_min", "ms2_cosine_conflict", "within_factor", "max_interval_width"):
        pc[k] = float(pc[k])
    for k in ("half", "search", "apex_floor"): pc["integrate"][k] = float(pc["integrate"][k])
    if str(pc.get("rt_shift", "auto")) != "auto": pc["rt_shift"] = float(pc["rt_shift"] or 0.0)
    return cfg


def quantify_injections(cands: pd.DataFrame, smap: pd.DataFrame, cfg: dict, ref: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    normalise_config(cfg)
    pc = cfg["pec_screen"]; ig = pc["integrate"]; ppm = cfg["level4"]["ppm"]
    fw, fg = float(pc["water_factor"]), float(pc["grain_factor"])
    cos_min, cos_bad = float(pc["ms2_cosine_min"]), float(pc["ms2_cosine_conflict"])
    rows, shifts = [], []
    todo = smap[smap.role != "exclude"]
    for s in todo.itertuples():
        pol = (s.polarity or guess_polarity(s.path, cfg.get("polarity_patterns")) or "positive")
        cs = cands[cands.polarity == pol]
        if cs.empty: continue
        t0 = time.time()
        run = mzml_io.read_for_config(s.path, common.config_for_file(cfg, s.path), pol)
        idx = MS2Index(run); store = {}
        wide = float(pc["rt_search_wide_min"])
        for c in cs.itertuples():
            rt, y = run.eic(c.mz, ppm, c.rt - wide - 0.6, c.rt + wide + 0.6)
            store[c.cand_id] = (rt, y)
        shift = 0.0
        if s.matrix == "grain" and s.role == "sample" and str(pc.get("rt_shift", "auto")) == "auto":
            d = []
            for c in cs.itertuples():
                rt, y = store[c.cand_id]; a, rtf, apex = integrate_eic(rt, y, c.rt, 0.0, wide, ig["apex_floor"])
                if a > 0 and apex >= pc["rt_shift_min_apex"]:
                    cs_, n_, _ = _ms2_cos(idx, c, rtf, cfg)
                    if cs_ >= cos_min and n_ >= 3: d.append(rtf - c.rt)
            if len(d) >= int(pc["rt_shift_min_anchors"]):
                shift = float(np.median(d))
            shifts.append(dict(sample=s.stem, rt_shift_min=shift, n_anchors=len(d)))
        elif s.matrix == "grain" and s.role == "sample":
            shift = float(pc.get("rt_shift") or 0.0)
        for c in cs.itertuples():
            rt, y = store[c.cand_id]
            area, rtf, apex = integrate_eic(rt, y, c.rt + shift, ig["half"], ig["search"], ig["apex_floor"])
            cos_s, n_m, n_ms2 = _ms2_cos(idx, c, rtf, cfg) if area > 0 else (0.0, 0, 0)
            if area <= 0: status = ""
            elif n_ms2 == 0: status = "MS1 only"
            elif cos_s >= cos_min and n_m >= 2: status = "MS2 confirmed"
            elif cos_s < cos_bad: status = "MS2 conflict"
            else: status = "MS2 weak"
            ug = (ref["conc_ug_mL"] * area / ref["area"] if np.isfinite(ref["area"]) else area / ref["area_per_ug_mL"]) * float(s.scale)
            rows.append(dict(cand_id=c.cand_id, sample=s.stem, matrix=s.matrix, field=s.field, role=s.role, polarity=pol, rt_found=rtf, apex=apex, area=area,
                             ms2_cos=cos_s, ms2_matched=n_m, n_ms2=n_ms2, id_status=status, rt_shift=shift,
                             conc_ug_mL_extract=ug, conc=ug * (fw if s.matrix == "water" else fg)))
        log.info("quantified %d candidates in %s (%s, shift %+.2f min) %.0fs", len(cs), s.stem, s.matrix, shift, time.time() - t0)
    return pd.DataFrame(rows), pd.DataFrame(shifts)


def apply_blank_check(D: pd.DataFrame, ratio: float) -> pd.DataFrame:
    """Per candidate x matrix: apex of field injections vs the highest blank apex; below ratio x blank = not distinguishable (detected False)."""
    D = D.copy(); D["blank_apex"] = np.nan; D["blank_ratio"] = np.nan
    bl = D[D.role == "blank"].groupby(["cand_id", "matrix"])["apex"].max()
    for (cid, mat), bmax in bl.items():
        m = (D.cand_id == cid) & (D.matrix == mat) & (D.role == "sample")
        D.loc[m, "blank_apex"] = bmax
        D.loc[m, "blank_ratio"] = np.where(bmax > 0, D.loc[m, "apex"] / bmax, np.inf)
    smp = D.role == "sample"
    blank_level = smp & D.blank_apex.notna() & (D.blank_apex > 0) & (D.apex < ratio * D.blank_apex)
    D["detected"] = smp & (D.area > 0) & ~blank_level & (D.id_status != "MS2 conflict")
    D["blank_level"] = blank_level & (D.area > 0)
    return D


def aggregate_fields(D: pd.DataFrame) -> pd.DataFrame:
    """One row per candidate x matrix x field: mean over detected injections (arithmetic), n_detected / n_injections, best identity status."""
    S = D[(D.role == "sample") & (D.field != "")]
    rows = []
    for (cid, mat, fld), g in S.groupby(["cand_id", "matrix", "field"]):
        det = g[g.detected]
        rows.append(dict(cand_id=cid, matrix=mat, field=fld, n_inj=len(g), n_det=len(det), conc=float(det.conc.mean()) if len(det) else np.nan,
                         conc_min=float(det.conc.min()) if len(det) else np.nan, conc_max=float(det.conc.max()) if len(det) else np.nan,
                         id_status=max(det.id_status, key=lambda x: STATUS_RANK[x]) if len(det) else "", blank_level_n=int(g.blank_level.sum())))
    return pd.DataFrame(rows)


# ======================================================================================================================================
# step 3 - PEC and the comparison
# ======================================================================================================================================
def predict_and_compare(cands, F, props, kr, cfg, cf: dict | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    pc = cfg["pec_screen"]; fac = float(pc["fruit_to_grain_dw"]); within = float(pc["within_factor"])
    allow_ms1 = str(pc.get("grain_identity", "ms2_or_ms1")) == "ms2_or_ms1"
    floor_ng_g = pc["_grain_floor_ng_g"]; lab = pc.get("_sink_label", "grain")
    W = F[F.matrix == "water"].set_index(["cand_id", "field"]); G = F[F.matrix == "grain"].set_index(["cand_id", "field"])
    gfields = set(F[F.matrix == "grain"].field); wfields = set(F[F.matrix == "water"].field)
    fields = sorted((wfields - {"*"}) | gfields) or sorted(wfields)       # a water row with field '*' (mean of all irrigation-water injections) applies to every solid-matrix field
    out = []
    for c in cands.itertuples():
        p = props.loc[props.cand_id == c.cand_id].iloc[0]; k = kr[c.cand_id]; cfv = (cf or {}).get(c.ik1, 1.0)
        for fld in fields:
            w = W.loc[(c.cand_id, fld)] if (c.cand_id, fld) in W.index else (W.loc[(c.cand_id, "*")] if (c.cand_id, "*") in W.index else None)
            g = G.loc[(c.cand_id, fld)] if (c.cand_id, fld) in G.index else None
            cw = float(w.conc) * cfv if (w is not None and w.n_det > 0) else np.nan
            mec = float(g.conc) if (g is not None and g.n_det > 0) else np.nan
            ids = g.id_status if g is not None else ""
            if np.isfinite(mec) and ids in ("MS1 only", "MS2 weak") and not allow_ms1: mec = np.nan
            pl, pcn, ph = (k["k_low"] * cw * fac, k["k_central"] * cw * fac, k["k_high"] * cw * fac) if np.isfinite(cw) else (np.nan,) * 3
            out.append(dict(cand_id=c.cand_id, field=fld, water_ng_L=cw, water_n_det=(w.n_det if w is not None else 0), water_n_inj=(w.n_inj if w is not None else 0),
                            MEC_ng_g=mec, grain_n_det=(g.n_det if g is not None else 0), grain_n_inj=(g.n_inj if g is not None else 0), grain_id=ids,
                            PEC_low=pl, PEC_central=pcn, PEC_high=ph,
                            ratio_central=pcn / mec if np.isfinite(mec) and np.isfinite(pcn) else np.nan,
                            ratio_low=pl / mec if np.isfinite(mec) and np.isfinite(pl) else np.nan,
                            ratio_high=ph / mec if np.isfinite(mec) and np.isfinite(ph) else np.nan))
    P = pd.DataFrame(out)
    S = []
    for c in cands.itertuples():
        p = props.loc[props.cand_id == c.cand_id].iloc[0]; k = kr[c.cand_id]; g = P[P.cand_id == c.cand_id]
        pair = g[g.MEC_ng_g.notna() & g.PEC_central.notna()]
        wdet = g[g.PEC_central.notna()]; gdet = g[g.MEC_ng_g.notna()]
        gsamples_for_w = wdet[wdet.grain_n_inj > 0]
        gm = lambda v: float(np.exp(np.mean(np.log(v)))) if len(v) else np.nan
        r_c = gm(pair.ratio_central.values) if len(pair) else np.nan
        r_lo = gm(pair.ratio_low.values) if len(pair) else np.nan; r_hi = gm(pair.ratio_high.values) if len(pair) else np.nan
        if p.not_applicable: cls = "fate model not applicable"
        elif len(pair):
            # ratio interval [PEC_low, PEC_high]/MEC (a point when properties are complete) against the x`within` window
            width = r_hi / r_lo if (r_lo and np.isfinite(r_lo) and r_lo > 0) else np.inf
            if r_hi < 1 / within: cls = "MEC >> PEC (other source / contamination?)"
            elif r_lo > within: cls = "PEC >> MEC (over-predicted / low uptake)"
            elif width > float(pc["max_interval_width"]): cls = "unconstrained (properties missing)"      # window overlap is meaningless for a very wide PEC range
            else: cls = "water-explained (within x%g)" % within
            if k["k_mode"] == "grid" and not cls.startswith("unconstrained"): cls += " [bounded: properties missing]"
        elif len(gsamples_for_w):
            cls = ("predicted but not found in %s" % lab if gsamples_for_w.PEC_low.max() >= floor_ng_g else "PEC below %s detection floor (consistent)" % lab)
        elif len(gdet): cls = "%s only (no paired water)" % lab
        elif len(wdet): cls = "water only (no %s data for these fields)" % lab
        else: cls = "not quantified"
        S.append(dict(cand_id=c.cand_id, n_pairs=len(pair), n_water_fields=len(wdet), n_grain_fields=len(gdet), PEC_over_MEC_geomean=r_c, PEC_over_MEC_low=r_lo,
                      PEC_over_MEC_high=r_hi, pec_interval_width=(r_hi / r_lo if (np.isfinite(r_lo) and r_lo > 0) else np.nan), class_=cls, prop_status=p.prop_status, k_mode=k["k_mode"], k_low=k["k_low"], k_high=k["k_high"],
                      grain_identity=(max(gdet.grain_id, key=lambda x: STATUS_RANK.get(x, -1)) if len(gdet) else ""),
                      water_ng_L_max=float(g.water_ng_L.max()) if g.water_ng_L.notna().any() else np.nan,
                      MEC_ng_g_max=float(g.MEC_ng_g.max()) if g.MEC_ng_g.notna().any() else np.nan))
    return P, pd.DataFrame(S).rename(columns={"class_": "class"})


def grain_floor_ng_g(cfg: dict, ref: dict) -> float:
    """ng/g dw that corresponds to the smallest peak the integrator accepts (smoothed apex floor x typical peak width)."""
    normalise_config(cfg); pc = cfg["pec_screen"]; area = float(pc["integrate"]["apex_floor"]) * float(pc["grain_floor_width_min"])
    ug = ref["conc_ug_mL"] * area / ref["area"] if np.isfinite(ref["area"]) else area / ref["area_per_ug_mL"]
    return ug * float(pc["grain_factor"])
