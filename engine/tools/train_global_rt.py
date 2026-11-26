"""Train the global (method-independent) retention model on the METLIN SMRT dataset (80 038 compounds, reversed-phase C18,
water/acetonitrile + formic acid; Domingo-Almenara et al. 2019, Nat. Commun. 10:5811, CC BY 4.0 / Apache 2.0).

    python tools/train_global_rt.py <SMRT_dataset.csv> <out_dir>

The model predicts the SMRT retention time (min) from RDKit descriptors. It is NOT used directly: for each lab method the
engine fits a monotonic map from this 'global retention index' to the lab's own RTs, using anchor compounds (the lab's
standards and the internal standards found in the targeted screen). Non-retained compounds (RT < 5 min, ~2 % of SMRT)
are dropped, as in the original paper. Writes global_rt_model.pkl + global_rt_model_cv.txt.
"""
from __future__ import annotations
import os
import pickle
import sys
import numpy as np
import pandas as pd


def main(src, out_dir):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from hrms_nts.rt_prediction import global_descriptors, GLOBAL_DESC_NAMES
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.model_selection import KFold, cross_val_predict
    df = pd.read_csv(src, sep=";")
    df["rt_min"] = df.rt / 60.0
    df = df[df.rt_min >= 5.0].reset_index(drop=True)
    os.makedirs(out_dir, exist_ok=True)
    cache = os.path.join(out_dir, "_smrt_descriptors.npz")
    if os.path.exists(cache):
        z = np.load(cache); X, y = z["X"], z["y"]
    else:
        X, y = [], []
        for inchi, rt in zip(df.inchi, df.rt_min):
            m = Chem.MolFromInchi(inchi)
            if m is None:
                continue
            d = global_descriptors(Chem.MolToSmiles(m))
            if d is not None:
                X.append(d); y.append(rt)
        X, y = np.array(X, float), np.array(y, float)
        np.savez_compressed(cache, X=X, y=y)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    print(f"{len(y)} retained SMRT compounds, {X.shape[1]} descriptors", flush=True)
    mdl = HistGradientBoostingRegressor(max_iter=800, learning_rate=0.06, max_leaf_nodes=63, min_samples_leaf=20,
                                        l2_regularization=0.5, random_state=0)
    pred = cross_val_predict(mdl, X, y, cv=KFold(5, shuffle=True, random_state=0))
    err = np.abs(pred - y)
    rep = (f"5-fold CV on {len(y)} compounds: MAE {err.mean():.2f} min, median AE {np.median(err):.2f} min, "
           f"RMSE {np.sqrt(np.mean(err ** 2)):.2f} min, R2 {1 - np.sum((pred - y) ** 2) / np.sum((y - y.mean()) ** 2):.3f} "
           f"(SMRT RT range {y.min():.1f}-{y.max():.1f} min)")
    print(rep, flush=True)
    mdl.fit(X, y)
    with open(os.path.join(out_dir, "global_rt_model.pkl"), "wb") as f:
        pickle.dump(dict(model=mdl, descriptors=GLOBAL_DESC_NAMES, source="METLIN SMRT (retained, RT >= 5 min)",
                         n=len(y), cv=rep), f)
    open(os.path.join(out_dir, "global_rt_model_cv.txt"), "w").write(rep + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
