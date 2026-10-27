"""How much does retention time tell us about the Level 2b identifications?
Fit NORMAN Pred_RTI -> RT on the robust Level 1/2a anchors (oligomer series and surfactant homologues excluded - one library entry matches every chain length), then compare
the prediction error of the held-out anchors (leave-one-out) with the error of every Level 2b candidate.  If 2b errors are much wider than the anchor errors, a share of the
2b identifications is probably wrong; the share beyond +-2 x anchor-RMSE is reported.

    python tools/rt_check_levels.py work/pec_screen_israel_pos/candidates.csv"""
import os, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hrms_nts.common import load_config
from hrms_nts.database import build_database
from sklearn.linear_model import HuberRegressor

C = pd.read_csv(sys.argv[1]).fillna("")
C = C[C.polarity == "positive"].copy(); C["rt"] = C.rt.astype(float)
db = build_database(load_config()).df.drop_duplicates("ik1")[["ik1", "pred_rti_pos"]]
C = C.merge(db, on="ik1", how="left"); C = C[C.pred_rti_pos.notna()]
series = re.compile(r"glycol|PEG|PPG|ethoxy|polyoxy|EO\b|oxa|ether|alkyl|amide of|dodecyl|hexadecyl|tetradecyl", re.I)
C["series"] = C.name.map(lambda n: bool(series.search(str(n))))
anc = C[C.level.isin(["1", "2a"]) & (C.lib_score.astype(float) >= 0.8) & ~C.series & (C.name != "Alprazolam")]
b = C[C.level == "2b"]
print(f"anchors (2a, cos>=0.8, non-series): {len(anc)}   2b candidates with RTI: {len(b)} ({int(b.series.sum())} series/surfactant-like)")
loo = []
for i in anc.index:
    tr = anc.drop(i); h = HuberRegressor().fit(tr[["pred_rti_pos"]], tr.rt); loo.append(float(h.predict(anc.loc[[i], ["pred_rti_pos"]])[0]) - anc.loc[i, "rt"])
e = np.array(loo); rm = float(np.sqrt((e ** 2).mean()))
print(f"anchors LOO: RMSE {rm:.2f} min | median |err| {np.median(np.abs(e)):.2f} | within 0.5: {(np.abs(e) <= .5).mean():.0%} | within 1.0: {(np.abs(e) <= 1).mean():.0%} | within 2.0: {(np.abs(e) <= 2).mean():.0%}")
h = HuberRegressor().fit(anc[["pred_rti_pos"]], anc.rt); print(f"map: RT = {h.intercept_:.2f} + {h.coef_[0]:.4f} x RTI")
for nm, d in (("2b (all)", b), ("2b non-series", b[~b.series]), ("2b series-like", b[b.series])):
    if len(d):
        err = h.predict(d[["pred_rti_pos"]]) - d.rt.values
        print(f"{nm:15s} n={len(d):3d} median |err| {np.median(np.abs(err)):.2f} | within 1.0: {(np.abs(err) <= 1).mean():.0%} | beyond +-{2*rm:.1f} min: {(np.abs(err) > 2 * rm).mean():.0%}")
bb = b.assign(err=h.predict(b[["pred_rti_pos"]]) - b.rt.values)
bb["abs"] = bb.err.abs()
out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "work", "rt_check_2b.csv")
bb.sort_values("abs", ascending=False)[["name", "level", "rt", "pred_rti_pos", "err"]].to_csv(out, index=False); print("->", out)
