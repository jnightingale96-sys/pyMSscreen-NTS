"""Compound Discoverer semi-quant tables (the project workbook) -> the same candidate / field tables that `pec-screen` builds from mzML.

Workbook layout (Israel_Modelling_Data): sheet 'Semi quant_Modelling_WW'  one row per annotated chemical x polarity, ONE value per chemical
(max peak area over the injections -> ug/mL -> ng/L, project formula x1000 x120): no per-field resolution, so that single value is used as C_water for
every grain field.  Sheet 'Semi quant modelling_Plants': Name, RT, average peak area, ug/mL, ng/g (x0.5/3x1000), Sample triple, Field (1-9 or 'FIELD n').
The CompTox annotation block of the Plants sheet is misaligned with 'Name' on many rows, so structures are taken only from rows where INPUT == Name
(water + plants + soil sheets) and otherwise from the toolkit suspect database by exact name."""
from __future__ import annotations
import logging, os, pickle, re
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")
WW, PL = "Semi quant_Modelling_WW", "Semi quant modelling_Plants"


def _f(x):
    try: return float(x)
    except (TypeError, ValueError): return np.nan


def read_cd_sheets(path: str, cache: str | None = None, max_rows=600) -> dict:
    if cache and os.path.exists(cache) and os.path.getmtime(cache) > os.path.getmtime(path):
        raw = pickle.load(open(cache, "rb"))
        if WW in raw and PL in raw and len(raw[PL]) > 150: return raw
    import openpyxl, warnings
    warnings.filterwarnings("ignore")
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    raw = {}
    for n in (WW, PL, "Semi quant modelling_Soil"):
        rows, blank = [], 0
        for r in wb[n].iter_rows(values_only=True, max_row=max_rows):
            rows.append(list(r[:30])); blank = 0 if any(c is not None for c in r[:12]) else blank + 1
            if blank > 25: break
        raw[n] = rows
    if cache: pickle.dump(raw, open(cache, "wb"))
    return raw


def _struct(r, i_in, i_cas, i_ik, i_smi, i_qsar):
    """CompTox block of one row -> dict (only used when INPUT == Name)."""
    smi = r[i_qsar] if r[i_qsar] not in (None, "", "#N/A") else r[i_smi]
    ik = str(r[i_ik] or ""); return dict(smiles=str(smi) if smi not in (None, "#N/A") else "", ik1=ik.split("-")[0] if ik not in ("", "#N/A") else "",
                                         cas=str(r[i_cas] or ""), inchikey=ik)


def build_cd_tables(raw: dict, cfg: dict | None = None, sink: str = "grain") -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Returns (water_rows, solid_rows, structures) with one row per sheet row.  sink='soil' reads 'Semi quant modelling_Soil' (columns: Name, Formula, ..., RT=6,
    sample 'Soil 8' = 10, peak area = 11, 'Normalised' = 12, ug/kg dw = 13; ug/kg dw = Normalised x 4500 = ng/g dw) instead of the grain sheet."""
    ww = [r for r in raw[WW][2:] if r[1]]
    pl = [r for r in raw[PL][1:] if r[1]]
    water = pd.DataFrame([dict(pol="positive" if str(r[0]).strip().startswith("+") else "negative", name=str(r[1]).strip(),
                               formula=str(r[2] or "").replace(" ", ""), rt=_f(r[4]), best_match=_f(r[5]), cls=str(r[7] or ""), area=_f(r[9]), ug_mL=_f(r[10]), ng_L=_f(r[11])) for r in ww])
    grain = pd.DataFrame([dict(pol="positive" if str(r[0]).strip().startswith("+") else "negative", name=str(r[1]).strip(), ppm=_f(r[2]), cls=str(r[3] or ""), rt=_f(r[7]),
                               area=_f(r[8]), ug_mL=_f(r[9]), ng_g=_f(r[10]), sample=str(r[13] or ""),
                               field=re.sub(r"(?i)^field\s*", "", str(r[14] or "")).strip()) for r in pl])
    if sink == "soil":
        so = [r for r in raw["Semi quant modelling_Soil"][2:] if r[1]]
        grain = pd.DataFrame([dict(pol="positive" if str(r[0]).strip().startswith("+") else "negative", name=str(r[1]).strip(), ppm=_f(r[3]), cls="", rt=_f(r[6]),
                                   area=_f(r[11]), ug_mL=_f(r[12]), ng_g=_f(r[13]), sample=str(r[10] or ""),
                                   field=(re.search(r"(\d+)", str(r[10] or "")) or [None, ""])[1]) for r in so])
    st = {}
    for r in ww:                                    # water block is aligned on every row
        if r[13] and str(r[13]).strip().lower() == str(r[1]).strip().lower():
            st.setdefault(str(r[1]).strip().lower(), dict(_struct(r, 13, 17, 18, 20, 23), formula=str(r[2] or "").replace(" ", "")))
    for r in (pl if sink != "soil" else []):         # plants block only where it matches the name
        if r[15] and str(r[15]).strip().lower() == str(r[1]).strip().lower():
            st.setdefault(str(r[1]).strip().lower(), _struct(r, 15, 19, 20, 22, 25))
    for r in raw.get("Semi quant modelling_Soil", [])[2:]:
        if len(r) > 25 and r[1] and r[15] and str(r[15]).strip().lower() == str(r[1]).strip().lower():
            st.setdefault(str(r[1]).strip().lower(), _struct(r, 15, 19, 20, 22, 25))
    S = pd.DataFrame([dict(name_l=k, **v) for k, v in st.items()])
    return water, grain, S


def cd_candidates(water: pd.DataFrame, grain: pd.DataFrame, S: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Candidate table, field table (same schema as pecscreen.aggregate_fields) and the long 'detections' table."""
    names = sorted(set(water.name) | set(grain.name), key=str.lower)
    sd = S.set_index("name_l").to_dict("index") if len(S) else {}
    missing = [n for n in names if n.lower() not in sd or not sd[n.lower()].get("smiles")]
    if missing:                                      # toolkit suspect database (SusDat / PubChemLite), exact name
        try:
            from .database import build_database
            db = build_database(cfg).df; dn = db.assign(_n=db["name"].astype(str).str.lower().str.strip()).drop_duplicates("_n").set_index("_n")
            for n in missing:
                if n.lower() in dn.index:
                    d = dn.loc[n.lower()]; sd[n.lower()] = dict(smiles=str(d["smiles"] or ""), ik1=str(d["ik1"] or ""), cas=str(d["cas"] or ""), inchikey=str(d["inchikey"] or ""))
        except Exception as e:
            log.warning("suspect database lookup failed (%s) - %d names stay without structure", e, len(missing))
    def cid_of(n):
        ik1 = (sd.get(n.lower(), {}) or {}).get("ik1") or ""
        return ik1 if ik1 else "N_" + re.sub(r"\W+", "_", n.lower())[:40]
    water = water.assign(cid=water.name.map(cid_of)); grain = grain.assign(cid=grain.name.map(cid_of))   # several names can share one structure
    rows = []
    for cid in sorted(set(water.cid) | set(grain.cid), key=lambda c: min([n.lower() for n in names if cid_of(n) == c])):
        w = water[water.cid == cid]; g = grain[grain.cid == cid]
        ns = list(dict.fromkeys(list(w.name) + list(g.name)))
        n = ns[0]; s_ = sd.get(n.lower(), {}) or {}
        pol = (w.pol.iloc[0] if len(w) else g.pol.iloc[0])
        rts = np.r_[w.rt.values, g.rt.values]
        rows.append(dict(cand_id=cid, name=n if len(ns) == 1 else f"{n} (= {'; '.join(ns[1:])})", formula=(w.formula.iloc[0] if len(w) else s_.get("formula", "")), smiles=s_.get("smiles", ""),
                         ik1=s_.get("ik1", ""), inchikey=s_.get("inchikey", ""), cas=s_.get("cas", ""), polarity=pol, adduct="", mz=np.nan,
                         rt=float(np.nanmedian(rts)) if np.isfinite(rts).any() else np.nan,
                         rt_water=float(np.nanmedian(w.rt.values)) if (len(w) and w.rt.notna().any()) else np.nan,
                         rt_grain=float(np.nanmedian(g.rt.values)) if (len(g) and g.rt.notna().any()) else np.nan, level="CD", cls=next((str(x) for x in list(w.cls) + list(g.cls) if str(x) not in ("", "nan", "None")), ""),
                         confidence=float(w.best_match.max()) if len(w) and w.best_match.notna().any() else np.nan, level_basis="Compound Discoverer annotation (no MS2 level assigned)",
                         source="CompoundDiscoverer", lib_name="", lib_score=0.0, n_ident_samples=0, ms2_peaks="",
                         in_water=bool(len(w)), in_grain=bool(len(g)), structure_found=bool(s_.get("smiles")), names=" | ".join(ns)))
    C = pd.DataFrame(rows)
    fields = sorted({f for f in grain.field if f}, key=lambda x: (len(x), x))
    F, D = [], []
    for cid in C.cand_id:
        w = water[(water.cid == cid) & water.ng_L.notna()]
        if len(w):
            row = w.loc[w.ng_L.idxmax()]
            for f in fields:                          # one max value, valid for every field
                F.append(dict(cand_id=cid, matrix="water", field=f, n_inj=1, n_det=1, conc=float(row.ng_L), conc_min=float(row.ng_L), conc_max=float(row.ng_L), id_status="CD annotation", blank_level_n=0))
            D.append(dict(cand_id=cid, sample="CD water max", matrix="water", field="all", role="sample", area=row.area, conc=row.ng_L, detected=True, id_status="CD annotation"))
        for f, gg in grain[(grain.cid == cid) & grain.ng_g.notna()].groupby("field"):
            v = float(gg.ng_g.max())
            F.append(dict(cand_id=cid, matrix="grain", field=f, n_inj=len(gg), n_det=len(gg), conc=v, conc_min=float(gg.ng_g.min()), conc_max=v, id_status="CD annotation", blank_level_n=0))
            D.append(dict(cand_id=cid, sample=";".join(sorted(set(gg["sample"]))), matrix="grain", field=f, role="sample", area=float(gg.area.max()), conc=v, detected=True, id_status="CD annotation"))
    F = pd.DataFrame(F).drop_duplicates(["cand_id", "matrix", "field"])
    return C, F, pd.DataFrame(D)
