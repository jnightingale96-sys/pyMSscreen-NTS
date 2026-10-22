"""Benchmark the built-in product-ion scorer on real library spectra (target-decoy).

For N random MassBank spectra of [M+H]+ / [M-H]- we score the TRUE structure and K decoy structures (other library compounds whose
precursor m/z lies within +-window Da, different planar structure) against the same experimental peaks.  Reports:
  * top-1 accuracy among 1+K candidates, AUC (true vs decoy)
  * 2b-style pass rate for true structures vs decoys (= empirical false-positive rate of the evidence criterion)
Usage: python benchmark_insilico.py [N=300] [K=7] [polarity=positive]"""
import sys, os, sqlite3, random, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hrms_nts.common import load_config, resolve
from hrms_nts import insilico as IS
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
K = int(sys.argv[2]) if len(sys.argv) > 2 else 7
POL = sys.argv[3] if len(sys.argv) > 3 else "positive"
ADD = "[M+H]+" if POL == "positive" else "[M-H]-"
cfg = load_config(None)
con = sqlite3.connect(resolve(cfg, cfg["paths"]["library_db"]))
rows = con.execute("SELECT id,name,ik1,smiles,precursor_mz,peaks FROM spectra WHERE polarity=? AND adduct=? AND smiles!='' AND precursor_mz BETWEEN 120 AND 600",
                   (POL, ADD)).fetchall()
random.seed(1)
random.shuffle(rows)
by_mass = sorted(rows, key=lambda r: r[4])
pm = np.array([r[4] for r in by_mass])
ic = cfg["insilico"]; lv = cfg["levels"]
res = []; t0 = time.time(); used = 0
for r in rows:
    if used >= N: break
    lid, name, ik1, smi, prec, blob = r
    a = np.frombuffer(blob, dtype=np.float32).reshape(-1, 2).astype(float)
    k = a[:, 0] < prec - 0.5
    if k.sum() < 4: continue
    lo, hi = np.searchsorted(pm, [prec - 30, prec + 30])
    pool = [x for x in by_mass[lo:hi] if x[2] != ik1]
    if len(pool) < K: continue
    decoys = random.sample(pool, K)
    def sc(s):
        try:
            return IS.score_candidate(s, ADD, a[:, 0], a[:, 1], prec, ic["frag_tol_da"], ic["frag_tol_ppm"], ic["rel_min"], ic["max_depth"])
        except Exception:
            return None
    t = sc(smi)
    if t is None: continue
    ds = [sc(d[3]) for d in decoys]; ds = [d for d in ds if d is not None]
    if len(ds) < 3: continue
    used += 1
    res.append(dict(name=name, true_score=t["bi_score"], true_adj=t["bi_expl_adj_n"], dec_scores=[d["bi_score"] for d in ds], dec_adj=[d["bi_expl_adj_n"] for d in ds]))
    if used % 50 == 0: print(f"  {used}/{N}  ({time.time()-t0:.0f}s)", flush=True)

def passes(score, adj):
    return adj >= lv["min_explained_peaks_2b"] and score >= lv["insilico_min_score"]
top1 = np.mean([r["true_score"] >= max(r["dec_scores"]) and r["true_score"] > 0 for r in res])
auc = np.mean([np.mean([(r["true_score"] > d) + 0.5 * (r["true_score"] == d) for d in r["dec_scores"]]) for r in res])
tp = np.mean([passes(r["true_score"], r["true_adj"]) for r in res])
fp = np.mean([passes(s, a) for r in res for s, a in zip(r["dec_scores"], r["dec_adj"])])
print(f"\n=== built-in scorer benchmark ({POL}, {ADD}); {len(res)} spectra, {K} decoys each ===")
print(f"top-1 accuracy among {K+1}: {top1:.2f}   (chance {1/(K+1):.2f})")
print(f"AUC true vs decoy: {auc:.2f}   (0.5 = random)")
print(f"2b-style evidence pass: true structures {tp:.2f}  vs  decoys {fp:.2f}  (empirical false-positive rate of the criterion)")
pd.DataFrame([dict(polarity=POL, n=len(res), k=K, top1=top1, auc=auc, pass_true=tp, pass_decoy=fp,
                   min_expl=lv["min_explained_peaks_2b"], min_score=lv["insilico_min_score"])]).to_csv(os.path.join(cfg["_root"], "work", f"benchmark_insilico_{POL}.csv"), index=False)
