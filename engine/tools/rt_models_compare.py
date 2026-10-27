"""Which RT model suits THIS run?  Candidate improvements taken from the three papers, evaluated by leave-one-out on the robust Level 1/2a anchors:
  A  RT ~ RTI                       (Aalizadeh 2021: unified RTI; the plain linear RTI->tR map, here Huber-robust)
  B  sqrt(RT) ~ RTI                 (Celma 2022: model sqrt(RT) so the prediction-interval variance is constant)
  C  RT ~ RTI + TPSA                (descriptor correction; TPSA from Updated_NORMAN_with_TPSA)
  D  RT ~ RTI + TPSA + Crippen logP
  E  A restricted to anchors inside the RTI model's applicability domain ('Covered by Model')
and report how many Level 2b candidates each model puts inside its own 95 % band.   (MultiConditionRT, Souihi 2022, needs training data for other columns/eluents and PaDEL
descriptors - not reproducible from the paper alone; its gradient-mapping step is what A/B already do with the anchors.)

    python tools/rt_models_compare.py work/pec_screen_israel_pos --norman Downloads/Updated_NORMAN_with_TPSA.xlsx"""
import argparse, glob, os, pickle, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from sklearn.linear_model import HuberRegressor
from hrms_nts.common import load_config
from hrms_nts import pecscreen as PS
from hrms_nts.idfilter import SERIES

ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--norman", required=True); ap.add_argument("--min-cos", type=float, default=0.8); a = ap.parse_args()
cfg = load_config(); ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); W = os.path.join(ROOT, "work")
ident = {os.path.basename(f)[:-4]: pickle.load(open(f, "rb")) for f in glob.glob(os.path.join(a.run, "_ident", "*.pkl"))}
C = PS.build_candidates(ident, cfg, "2b").fillna(""); C = C[C.polarity == "positive"].copy(); C["rt"] = C.rt.astype(float)
nw = pd.read_pickle(os.path.join(W, "norman_tpsa.pkl")) if os.path.exists(os.path.join(W, "norman_tpsa.pkl")) else pd.read_excel(a.norman, dtype=str)
nw["ik1"] = nw.StdInChIKey.str.split("-").str[0]
N = pd.DataFrame({"ik1": nw.ik1, "rti": pd.to_numeric(nw.Pred_RTI_Positive_ESI, errors="coerce"), "ad": nw.Uncertainty_RTI_pos.fillna("").astype(str), "tpsa": pd.to_numeric(nw.TPSA_calculated, errors="coerce")}).drop_duplicates("ik1")
# SusDat full for the compounds the workbook does not contain
sus = pd.read_csv(cfg["paths"]["susdat"], sep="\t", engine="python", dtype=str, usecols=["StdInChIKey", "Pred_RTI_Positive_ESI", "Uncertainty_RTI_pos"])
S = pd.DataFrame({"ik1": sus.StdInChIKey.str.split("-").str[0], "rti": pd.to_numeric(sus.Pred_RTI_Positive_ESI, errors="coerce"), "ad": sus.Uncertainty_RTI_pos.fillna("").astype(str)}).drop_duplicates("ik1")
R = pd.concat([N, S[~S.ik1.isin(N.ik1)]]).drop_duplicates("ik1")
C = C.merge(R, on="ik1", how="left"); C = C[C.rti.notna()]
from rdkit import Chem
from rdkit.Chem import Crippen, rdMolDescriptors
def desc(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None: return (np.nan, np.nan)
    m = max(Chem.GetMolFrags(m, asMols=True), key=lambda x: x.GetNumHeavyAtoms()); return (rdMolDescriptors.CalcTPSA(m), Crippen.MolLogP(m))
d = np.array([desc(s) for s in C.smiles]); C["tpsa_c"] = np.where(np.isnan(C.tpsa), d[:, 0], C.tpsa); C["logp"] = d[:, 1]
C["series"] = C.name.map(lambda n: bool(SERIES.search(str(n))))
anc = C[C.level.isin(["1", "2a"]) & (C.lib_score.astype(float) >= a.min_cos) & ~C.series].copy(); b2 = C[C.level == "2b"].copy()
print(f"anchors {len(anc)} (in applicability domain: {int(anc.ad.str.contains('Covered').sum())}) | 2b with RTI {len(b2)} (in AD: {int(b2.ad.str.contains('Covered').sum())})")

def feats(d, cols): return d[cols].values
def loo(cols, sqrt=False, data=None):
    d = anc if data is None else data; X = feats(d, cols); y = d.rt.values; y = np.sqrt(y) if sqrt else y; e = []
    for i in range(len(d)):
        m = np.ones(len(d), bool); m[i] = False
        mu, sd = X[m].mean(0), X[m].std(0) + 1e-9
        h = HuberRegressor(max_iter=500).fit((X[m] - mu) / sd, y[m]); p = float(h.predict(((X[[i]] - mu) / sd))[0])
        e.append((p ** 2 if sqrt else p) - d.rt.values[i])
    return np.array(e)
models = {"A RT~RTI": (["rti"], False, None), "B sqrt(RT)~RTI": (["rti"], True, None), "C RT~RTI+TPSA": (["rti", "tpsa_c"], False, None),
          "D RT~RTI+TPSA+logP": (["rti", "tpsa_c", "logp"], False, None), "E A, AD-covered anchors only": (["rti"], False, anc[anc.ad.str.contains("Covered")])}
print(f"\n{'model':32s} {'n':>3s} {'LOO RMSE':>9s} {'median|e|':>9s} {'<=1 min':>8s} {'95% band(min)':>13s} | 2b inside band")
res = {}
for k, (cols, sq, dd) in models.items():
    e = loo(cols, sq, dd); n = len(e); band = float(np.quantile(np.abs(e), 0.95))
    d = anc if dd is None else dd; X = d[cols].values; y = np.sqrt(d.rt.values) if sq else d.rt.values; mu, sd = X.mean(0), X.std(0) + 1e-9
    h = HuberRegressor(max_iter=500).fit((X - mu) / sd, y); pr = h.predict((b2[cols].values - mu) / sd); pr = pr ** 2 if sq else pr
    inside = float((np.abs(pr - b2.rt.values) <= band).mean())
    print(f"{k:32s} {n:3d} {np.sqrt((e ** 2).mean()):9.2f} {np.median(np.abs(e)):9.2f} {(np.abs(e) <= 1).mean():8.0%} {band:13.2f} | {inside:.0%}"); res[k] = (e, band)
b2in = b2[b2.ad.str.contains("Covered")]; h = HuberRegressor().fit(anc[["rti"]], anc.rt); err = np.abs(h.predict(b2in[["rti"]]) - b2in.rt.values)
print(f"\n2b candidates inside the RTI applicability domain: {len(b2in)}; median |error| {np.median(err):.2f} min; within 2.1 min: {(err <= 2.1).mean():.0%}  (outside AD: {len(b2) - len(b2in)} - RT should not veto these)")
ad_in = anc[anc.ad.str.contains("Covered")]; ad_out = anc[~anc.ad.str.contains("Covered")]
print("anchor error inside AD vs outside AD (model A, all-anchor fit):", end=" ")
for nm, dd in (("in", ad_in), ("out", ad_out)):
    if len(dd): print(f"{nm}: n={len(dd)} median|e| {np.median(np.abs(h.predict(dd[['rti']]) - dd.rt.values)):.2f}", end="  ")
print()
