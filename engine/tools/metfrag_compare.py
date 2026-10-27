"""MetFrag product-ion annotation / identification check with different candidate databases:
    NORMAN workbook (the curated pharmaceutical-type subset with RTI/TPSA)  vs  NORMAN SusDat (full)  vs  PubChemLite.
For every screening candidate (water-identified chemical with a consensus MS2) the structures of each database within 5 ppm of the neutral mass are scored by MetFrag
(FragmenterScore + offline MassBank score).  Reported per database: is OUR assigned structure in it, its MetFrag rank/score, the top structure, and the share of MS2 peaks
MetFrag explains (assigned vs top).  Product-ion table: every observed MS2 peak with the fragment formula MetFrag assigns to our structure (and the built-in predictor's label).

    python tools/metfrag_compare.py work/pec_screen_israel_pos --norman Downloads/Updated_NORMAN_with_TPSA.xlsx [--workers 4] [--limit 20]"""
import argparse, glob, os, pickle, re, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hrms_nts.common import load_config, resolve, adduct_table
from hrms_nts.database import build_database
from hrms_nts.ms2 import string_to_peaks
from hrms_nts import insilico as IS, pecscreen as PS

ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--norman", required=True); ap.add_argument("--ppm", type=float, default=5.0); ap.add_argument("--max-cands", type=int, default=3000)
ap.add_argument("--workers", type=int, default=4); ap.add_argument("--limit", type=int); ap.add_argument("--out"); ap.add_argument("--all-candidates", action="store_true", help="also candidates removed by the id filter")
a = ap.parse_args()
cfg = load_config(); ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); W = os.path.join(ROOT, "work")
ident = {os.path.basename(f)[:-4]: pickle.load(open(f, "rb")) for f in glob.glob(os.path.join(a.run, "_ident", "*.pkl"))}
C = PS.build_candidates(ident, cfg, "2b").fillna(""); C = C[C.ms2_peaks.astype(str) != ""]
if not a.all_candidates and os.path.exists(os.path.join(a.run, "pec_screen_summary.csv")):
    C = C[C.cand_id.isin(pd.read_csv(os.path.join(a.run, "pec_screen_summary.csv")).cand_id)]
if a.limit: C = C.head(a.limit)
print(f"{len(C)} candidates")

full = build_database(cfg).df.drop_duplicates(["ik1", "source"]); full = full[full.smiles.astype(str).str.len() > 1]
nw = pd.read_pickle(os.path.join(W, "norman_tpsa.pkl")) if os.path.exists(os.path.join(W, "norman_tpsa.pkl")) else pd.read_excel(a.norman, dtype=str)
nw = pd.DataFrame({"db_id": "NW" + nw["Norman_SusDat_ID"].astype(str), "name": nw["Name"], "formula": nw["Molecular_Formula"], "exact_mass": pd.to_numeric(nw["Monoiso_Mass"], errors="coerce"),
                   "smiles": nw["SMILES"], "inchikey": nw["StdInChIKey"], "source": "NORMAN-workbook"}).dropna(subset=["exact_mass", "smiles", "inchikey"])
nw["ik1"] = nw.inchikey.str.split("-").str[0]; nw = nw.drop_duplicates("ik1")
DBS = {"norman_workbook": nw, "susdat_full": full[full.source.str.contains("NORMAN")], "pubchemlite": full[full.source.str.contains("PubChem")]}
print({k: len(v) for k, v in DBS.items()})
ok, jar, java = IS.metfrag_available(cfg); assert ok, "MetFrag/Java not available"
mbf = resolve(cfg, cfg["paths"].get("metfrag_mb_folder")); work = os.path.join(W, "metfrag_compare"); os.makedirs(work, exist_ok=True)


def task(args):
    c, dbname = args
    db = DBS[dbname]
    ad = adduct_table(c.polarity).get(c.adduct, 1.007276 if c.polarity == "positive" else -1.007276); M = c.mz - ad
    cs = db[(db.exact_mass - M).abs() / M * 1e6 <= a.ppm]
    rec = dict(cand_id=c.cand_id, db=dbname, n=len(cs), in_db=bool((cs.ik1 == c.ik1).any()))
    if cs.empty: return rec, None
    if len(cs) > a.max_cands:
        keep = cs[cs.ik1 == c.ik1]; cs = pd.concat([keep, cs[cs.ik1 != c.ik1].sample(a.max_cands - len(keep), random_state=0)]); rec["capped"] = True
    mzp, itp = string_to_peaks(c.ms2_peaks); o = np.argsort(mzp); mzp, itp = mzp[o], itp[o]
    r = IS.run_metfrag(cs, mzp, itp, c.mz, c.polarity, c.adduct or "[M+H]+", work, jar, java, cfg["insilico"]["frag_tol_ppm"], cfg["insilico"]["frag_tol_da"], cfg["insilico"]["max_depth"], mbf,
                       tag=re.sub(r"\W", "_", f"{c.cand_id}_{dbname}"))
    if r is None or r.empty: rec["note"] = "no MetFrag result"; return rec, None
    r = r.merge(cs[["db_id", "ik1", "name"]].drop_duplicates("db_id"), on="db_id", how="left").sort_values("mf_score", ascending=False).reset_index(drop=True)
    mine = r.index[r.ik1 == c.ik1]
    rec.update(top1=str(r.at[0, "name"])[:50], top1_ik1=r.at[0, "ik1"], top1_score=float(r.at[0, "mf_score"]), top1_expl_frac=float(r.at[0, "mf_expl"] / max(r.at[0, "mf_used"], 1)),
               second=float(r.at[1, "mf_score"]) if len(r) > 1 else np.nan)
    if len(mine):
        i = mine[0]; rec.update(rank=int(i) + 1, score=float(r.at[i, "mf_score"]), expl_frac=float(r.at[i, "mf_expl"] / max(r.at[i, "mf_used"], 1)))
        return rec, (r.at[i, "mf_expl_formulas"], r.at[i, "mf_expl_peaks"])
    return rec, None


jobs = [(c, d) for c in C.itertuples() for d in DBS]
t0 = time.time(); recs, ions = [], {}
with ThreadPoolExecutor(max_workers=a.workers) as ex:
    for n, (rec, ann) in enumerate(ex.map(task, jobs), 1):
        recs.append(rec)
        if ann and rec["cand_id"] not in ions: ions[rec["cand_id"]] = (rec["db"], ann)
        if n % 25 == 0: print(f"{n}/{len(jobs)} MetFrag runs ({time.time() - t0:.0f}s)", flush=True)
R = pd.DataFrame(recs).merge(C[["cand_id", "name", "level"]], on="cand_id")
out = a.out or os.path.join(W, "metfrag_compare.csv"); R.to_csv(out, index=False)

# ---- product-ion table: observed peaks vs MetFrag formula for our assigned structure ---------------------------------------------------------
from hrms_nts.insilico import predict_product_ions
rows = []
for c in C.itertuples():
    mzp, itp = string_to_peaks(c.ms2_peaks); keep = itp >= 0.01 * itp.max(); mzp, itp = mzp[keep], itp[keep]
    ann = {}
    if c.cand_id in ions:
        for t in str(ions[c.cand_id][1][0]).split(";"):
            if ":" in t: mz_, f_ = t.split(":", 1); ann[round(float(mz_), 3)] = f_
    try:
        bp = predict_product_ions(c.smiles, c.adduct or "[M+H]+", top=80); bmz = bp["mz"].values; blab = bp["label"].values
    except Exception: bmz, blab = np.array([]), np.array([])
    for m, i in zip(mzp, itp):
        mf = next((v for k, v in ann.items() if abs(k - m) <= 0.003), "")
        j = np.where(np.abs(bmz - m) <= 0.005)[0]
        rows.append(dict(cand_id=c.cand_id, name=c.name, level=c.level, mz=round(float(m), 4), rel_int=round(float(i / itp.max()), 4), metfrag_db=ions.get(c.cand_id, ("", None))[0], metfrag_formula=mf,
                         builtin_label=str(blab[j[0]]) if len(j) else "", explained_by_metfrag=bool(mf), explained_by_builtin=bool(len(j))))
P = pd.DataFrame(rows); P.to_csv(os.path.join(W, "metfrag_product_ions.csv"), index=False)

# ---- summary ---------------------------------------------------------------------------------------------------------------------------------------
pd.set_option("display.width", 200); pd.set_option("display.float_format", lambda v: f"{v:,.2f}")
print("\n=== per database (candidates:", len(C), ") ===")
S = R.groupby("db").agg(candidates_with_hits=("n", lambda s: int((s > 0).sum())), median_n_structures=("n", "median"), assigned_in_db=("in_db", "sum"),
                        rank1=("rank", lambda s: int((s == 1).sum())), top3=("rank", lambda s: int((s <= 3).sum())), mean_expl_frac_assigned=("expl_frac", "mean"), mean_expl_frac_top1=("top1_expl_frac", "mean"))
print(S.to_string())
if len(P):
    e = P.groupby("cand_id").agg(mf=("explained_by_metfrag", "mean"), bi=("explained_by_builtin", "mean"), n=("mz", "size"))
    print(f"\nproduct ions: share of MS2 peaks (>=1% base) explained for the assigned structure: MetFrag {e.mf.mean():.0%}, built-in predictor {e.bi.mean():.0%} (mean over {len(e)} candidates)")
print("->", out, "|", os.path.join(W, "metfrag_product_ions.csv"))
