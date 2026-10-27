"""Retention-time model WITHOUT authentic standards: use the confident library-matched identifications (Level 2a, high cosine) of the run itself as RT anchors,
plus the few true standards (carbamazepine, alprazolam on 0_01.mzML), and fit the toolkit's RT models (NORMAN Pred_RTI -> RT map + RDKit QSAR).
Reports leave-one-out accuracy and writes user_data/standards_anchor.csv (same format as standards.csv) for `paths.standards`.

    python tools/rt_anchor_fit.py work/pec_screen_israel_pos/candidates.csv [--min-cos 0.8] [--polarity positive]"""
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hrms_nts.common import load_config, resolve
from hrms_nts.database import build_database
from hrms_nts.rt_prediction import RTModel

ap = argparse.ArgumentParser(); ap.add_argument("candidates"); ap.add_argument("--min-cos", type=float, default=0.8); ap.add_argument("--polarity", default="positive")
ap.add_argument("--out", default=None); a = ap.parse_args()
cfg = load_config()
C = pd.read_csv(a.candidates).fillna("")
anc = C[(C.level.isin(["1", "2a"])) & (C.lib_score.astype(float) >= a.min_cos) & (C.polarity == a.polarity) & (C.smiles != "")].copy()
# known standards on this method (0_01.mzML mix): RT from the EIC apex
std = pd.DataFrame([dict(name="Carbamazepine", smiles="NC(=O)N1c2ccccc2C=Cc2ccccc12", ik1="FFGPTBGBLSHEPO", rt=6.27, kind="standard"),
                    dict(name="Alprazolam", smiles="Cc1nnc2n1-c1ccc(Cl)cc1C(=NC2)c1ccccc1", ik1="VREFGVBLTWBCJP", rt=12.86, kind="standard")])
A = pd.concat([std, anc.assign(rt=anc.rt.astype(float), kind="Level " + anc.level)[["name", "smiles", "ik1", "rt", "kind"]]], ignore_index=True).drop_duplicates("ik1", keep="first")
print(f"anchors: {len(A)} ({(A.kind == 'standard').sum()} standards + {(A.kind != 'standard').sum()} library-matched IDs, cosine >= {a.min_cos})")
db = build_database(cfg).df
rti = db.drop_duplicates("ik1")[["ik1", "name", "pred_rti_pos", "pred_rti_neg"]]
A = A.merge(rti[["ik1", "pred_rti_pos"]], on="ik1", how="left"); print("anchors with a NORMAN Pred_RTI:", int(A.pred_rti_pos.notna().sum()))

# leave-one-out for the RTI map (robust linear) and for the full model
from sklearn.linear_model import HuberRegressor
R = A.dropna(subset=["pred_rti_pos"]).copy()
loo = []
for i in R.index:
    tr = R.drop(i); h = HuberRegressor().fit(tr[["pred_rti_pos"]], tr.rt); loo.append(float(h.predict(R.loc[[i], ["pred_rti_pos"]])[0]))
R["rt_loo_rti"] = loo; R["err"] = R.rt_loo_rti - R.rt
print(f"RTI->RT (Huber) leave-one-out: RMSE {np.sqrt((R.err ** 2).mean()):.2f} min, median |err| {R.err.abs().median():.2f}, within +-0.5: {(R.err.abs() <= 0.5).mean():.0%}, within +-1.0: {(R.err.abs() <= 1.0).mean():.0%}")
print(R.sort_values("err", key=abs, ascending=False)[["name", "kind", "rt", "pred_rti_pos", "rt_loo_rti", "err"]].head(8).round(2).to_string(index=False))
m = RTModel(a.polarity, tol=cfg["rt"]["tolerance_min"]).fit(A[["name", "smiles", "rt"]], db)
print(f"RTModel: n_train {m.n_train}  QSAR CV-RMSE {m.qsar_rmse:.2f}  RTI CV-RMSE {m.rti_rmse:.2f}")
out = a.out or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "user_data", "standards_anchor.csv")
S = A[["name", "smiles", "ik1", "rt"]].rename(columns={"ik1": "inchikey", "rt": "rt_pos" if a.polarity == "positive" else "rt_neg"})
S["rt_neg" if a.polarity == "positive" else "rt_pos"] = np.nan
S.to_csv(out, index=False); print("->", out)
