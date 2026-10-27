"""MetFrag second opinion on the screening candidates.

For each candidate (water-identified chemical with a consensus MS2) every structure in SusDat + PubChemLite within 5 ppm of the neutral mass is scored by MetFrag
(FragmenterScore + offline MassBank spectral score).  Reported: rank/score of OUR assignment, MetFrag's top structure, and the same rank when the candidate list is
restricted to a drug-rich NORMAN SusDat collection (default: source code S1, which contains every one of 14 reference drugs we tested).

    python tools/metfrag_verify.py work/pec_screen_israel_pos/candidates.csv --names "NPN,Isobutyl*" | --ids file.txt | --all  [--out work/metfrag_verify.csv]"""
import argparse, fnmatch, os, re, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hrms_nts.common import load_config, resolve, adduct_table
from hrms_nts.database import build_database
from hrms_nts.ms2 import string_to_peaks
from hrms_nts import insilico as IS

ap = argparse.ArgumentParser(); ap.add_argument("candidates", help="candidates.csv of a pec-screen run (its _ident folder is used)"); ap.add_argument("--names"); ap.add_argument("--all", action="store_true"); ap.add_argument("--ids")
ap.add_argument("--pharma-source", default="S1"); ap.add_argument("--out"); ap.add_argument("--ppm", type=float, default=5.0); ap.add_argument("--max-cands", type=int, default=2500)
ap.add_argument("--min-level", default="2b"); a = ap.parse_args()
cfg = load_config()
import glob, pickle
from hrms_nts import pecscreen as PS
d = os.path.dirname(os.path.abspath(a.candidates))                      # candidates are rebuilt from <run>/_ident/*.pkl so that the consensus MS2 is available
ident = {os.path.basename(f)[:-4]: pickle.load(open(f, "rb")) for f in glob.glob(os.path.join(d, "_ident", "*.pkl"))}
C = PS.build_candidates(ident, cfg, a.min_level).fillna("")
C = C[C.ms2_peaks.astype(str) != ""]
if a.names:
    pats = [p.strip().lower() for p in a.names.split(",")]
    C = C[[any(fnmatch.fnmatch(str(n).lower(), p if "*" in p else f"*{p}*") for p in pats) for n in C.name]]
elif a.ids:
    C = C[C.cand_id.isin([x.strip() for x in open(a.ids) if x.strip()])]
elif not a.all:
    sys.exit("choose --names, --ids or --all")
print(f"{len(C)} candidates")
db = build_database(cfg).df.drop_duplicates("ik1"); db = db[db.smiles.astype(str).str.len() > 1]
raw = pd.read_csv(resolve(cfg, cfg["paths"]["susdat"]), sep="\t", engine="python", dtype=str, usecols=["StdInChIKey", "Source"])
phar = set(raw.loc[raw.Source.fillna("").map(lambda s: a.pharma_source in [t.strip() for t in str(s).split(";")]), "StdInChIKey"].dropna().str.split("-").str[0])
print(f"drug-rich NORMAN subset ({a.pharma_source}): {len(phar)} structures")
ok, jar, java = IS.metfrag_available(cfg)
if not ok: sys.exit("MetFrag / Java not available (paths.metfrag_jar, java)")
mbf = resolve(cfg, cfg["paths"].get("metfrag_mb_folder")); work = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "work", "metfrag_verify")
rows = []
for i, c in enumerate(C.itertuples(), 1):
    ad = adduct_table(c.polarity)[c.adduct] if c.adduct in adduct_table(c.polarity) else (1.007276 if c.polarity == "positive" else -1.007276)
    M = c.mz - ad
    cands = db[(db.exact_mass - M).abs() / M * 1e6 <= a.ppm].copy()
    if c.ik1 not in set(cands.ik1):
        cands = pd.concat([cands, db[db.ik1 == c.ik1]])
    if len(cands) > a.max_cands: cands = cands.sample(a.max_cands, random_state=0) if c.ik1 not in set(cands.head(0).ik1) else cands
    mzp, itp = string_to_peaks(c.ms2_peaks); o = np.argsort(mzp); mzp, itp = mzp[o], itp[o]
    rec = dict(cand_id=c.cand_id, name=c.name, level=c.level, n_candidates=len(cands))
    for tag, cs in (("all", cands), ("drug", cands[cands.ik1.isin(phar) | (cands.ik1 == c.ik1)] if c.ik1 in phar else cands[cands.ik1.isin(phar)])):
        if cs.empty: rec.update({f"{tag}_n": 0}); continue
        t0 = time.time()
        r = IS.run_metfrag(cs.rename(columns={"db_id": "db_id"}), mzp, itp, c.mz, c.polarity, c.adduct or "[M+H]+", work, jar, java, cfg["insilico"]["frag_tol_ppm"],
                           cfg["insilico"]["frag_tol_da"], cfg["insilico"]["max_depth"], mbf, tag=re.sub(r"\W", "_", f"{c.cand_id}_{tag}"))
        rec[f"{tag}_n"] = len(cs)
        if r is None or r.empty: rec[f"{tag}_note"] = "no MetFrag result"; continue
        r = r.merge(cs[["db_id", "ik1", "name"]], on="db_id", how="left").sort_values("mf_score", ascending=False).reset_index(drop=True)
        mine = r.index[r.ik1 == c.ik1]
        rec[f"{tag}_rank"] = int(mine[0]) + 1 if len(mine) else np.nan
        rec[f"{tag}_score"] = float(r.loc[mine[0], "mf_score"]) if len(mine) else np.nan
        rec[f"{tag}_top1"] = str(r.loc[0, "name"])[:50]; rec[f"{tag}_top1_score"] = float(r.loc[0, "mf_score"])
        rec[f"{tag}_sec"] = float(r.loc[1, "mf_score"]) if len(r) > 1 else np.nan
        rec[f"{tag}_s"] = round(time.time() - t0, 1)
    rows.append(rec); print(f"[{i}/{len(C)}] {str(c.name)[:40]:40s} all: n={rec.get('all_n')} rank={rec.get('all_rank')} top1={rec.get('all_top1')} | drug: n={rec.get('drug_n')} rank={rec.get('drug_rank')}", flush=True)
R = pd.DataFrame(rows)
out = a.out or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "work", "metfrag_verify.csv")
R.to_csv(out, index=False); print("->", out)
