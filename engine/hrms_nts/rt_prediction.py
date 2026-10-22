"""Retention-time plausibility for candidates (supports the +-0.5 min criterion used in the project).

Two independent predictors, both calibrated on YOUR method's standards (standards.csv: name, smiles, rt_pos, rt_neg):
  * RTI  - linear/robust map from NORMAN SusDat 'Pred_RTI' (retention-time index) to run RT
  * QSAR - RDKit-descriptor regression (GradientBoosting / Ridge) trained on the same standards
Candidates are 'rt_agree' when ANY available model is within tolerance and 'rt_conflict' when ALL available models
are off by > conflict_factor x tolerance."""
from __future__ import annotations
import logging, pickle, os
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")

_SMARTS = {
    "acid": "[CX3](=O)[OX2H1,OX1-]", "sulfonic": "S(=O)(=O)[OH,O-]", "tetrazole": "c1nnn[nH]1",
    "phenol": "c[OH]", "amine_basic": "[NX3;H2,H1,H0;!$(NC=O);!$(N-a);!$(N-S(=O)=O);!$(N=*);!$(N#*)]",
    "amidine": "[NX3][CX3]=[NX2]", "quat": "[N+;!$(N=O);!$(N-[O-])]", "halogen": "[F,Cl,Br,I]",
}


def descriptors(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, Lipinski
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return None
    # largest fragment (drop counter-ions)
    frags = Chem.GetMolFrags(m, asMols=True)
    m = max(frags, key=lambda x: x.GetNumHeavyAtoms())
    d = [Crippen.MolLogP(m), Crippen.MolMR(m), rdMolDescriptors.CalcTPSA(m), Descriptors.MolWt(m),
         Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m), rdMolDescriptors.CalcNumRotatableBonds(m),
         rdMolDescriptors.CalcNumAromaticRings(m), rdMolDescriptors.CalcNumRings(m), rdMolDescriptors.CalcFractionCSP3(m),
         m.GetNumHeavyAtoms()]
    for k, s in _SMARTS.items():
        p = Chem.MolFromSmarts(s)
        d.append(len(m.GetSubstructMatches(p)))
    return d


DESC_NAMES = ["logP", "MR", "TPSA", "MW", "HBD", "HBA", "RotB", "AroRings", "Rings", "Fsp3", "HeavyAtoms"] + list(_SMARTS)


class RTModel:
    def __init__(self, polarity: str, tol=0.5, conflict_factor=2.0, adaptive_k=0.0):
        self.polarity = polarity; self.tol = tol; self.conflict_factor = conflict_factor; self.adaptive_k = adaptive_k
        self.qsar = None; self.qsar_rmse = np.nan
        self.rti = None; self.rti_rmse = np.nan; self.n_train = 0; self.rt_range = (0, 99)

    # ---- training -------------------------------------------------------------------------------------------
    def fit(self, standards: pd.DataFrame, rti_lookup: pd.DataFrame | None = None):
        """standards: columns name, smiles, rt (min, for this polarity); rti_lookup: db with ik1/name + pred_rti_{pos,neg}."""
        from sklearn.ensemble import GradientBoostingRegressor
        from sklearn.linear_model import Ridge, HuberRegressor
        from sklearn.model_selection import cross_val_predict, KFold
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        st = standards.dropna(subset=["rt"]).copy()
        X, y, keep = [], [], []
        for i, r in st.iterrows():
            d = descriptors(r.get("smiles"))
            if d is not None:
                X.append(d); y.append(float(r["rt"])); keep.append(i)
        X = np.array(X, float); y = np.array(y)
        self.n_train = len(y); self.rt_range = (float(y.min()), float(y.max())) if len(y) else (0, 99)
        if len(y) >= 8:
            if len(y) >= 25:
                mdl = GradientBoostingRegressor(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8, random_state=0)
            else:
                mdl = make_pipeline(StandardScaler(), Ridge(alpha=3.0))
            k = min(5, len(y))
            pred = cross_val_predict(mdl, X, y, cv=KFold(k, shuffle=True, random_state=0))
            self.qsar_rmse = float(np.sqrt(np.mean((pred - y) ** 2)))
            mdl.fit(X, y); self.qsar = mdl
            log.info("RT QSAR (%s): n=%d, CV RMSE=%.2f min", self.polarity, len(y), self.qsar_rmse)
        else:
            log.warning("RT QSAR needs >= 8 standards with SMILES (have %d)", len(y))
        # RTI calibration
        if rti_lookup is not None and len(st):
            col = "pred_rti_pos" if self.polarity == "positive" else "pred_rti_neg"
            ik = st["smiles"].map(_smiles_ik1)
            tab = rti_lookup[rti_lookup[col].notna()].drop_duplicates("ik1").set_index("ik1")[col]
            rti = ik.map(tab)
            ok = rti.notna() & st["rt"].notna()
            if ok.sum() >= 6:
                xr, yr = rti[ok].values.astype(float), st.loc[ok, "rt"].values.astype(float)
                hub = HuberRegressor().fit(xr.reshape(-1, 1), yr)
                pr = cross_val_predict(HuberRegressor(), xr.reshape(-1, 1), yr, cv=KFold(min(5, len(yr)), shuffle=True, random_state=0))
                self.rti_rmse = float(np.sqrt(np.mean((pr - yr) ** 2))); self.rti = hub
                log.info("RT RTI-map (%s): n=%d, CV RMSE=%.2f min", self.polarity, int(ok.sum()), self.rti_rmse)
        return self

    # ---- prediction -----------------------------------------------------------------------------------------
    def predict_qsar(self, smiles_list):
        out = np.full(len(smiles_list), np.nan)
        if self.qsar is None:
            return out
        cache = {}
        for i, s in enumerate(smiles_list):
            if s not in cache:
                d = descriptors(s); cache[s] = None if d is None else float(self.qsar.predict(np.array([d]))[0])
            out[i] = cache[s] if cache[s] is not None else np.nan
        return out

    def predict_rti(self, rti_values):
        v = np.asarray(rti_values, float)
        out = np.full(len(v), np.nan)
        if self.rti is None:
            return out
        ok = np.isfinite(v)
        if ok.any():
            out[ok] = self.rti.predict(v[ok].reshape(-1, 1))
        return out

    def tolerance(self, rmse):
        if self.adaptive_k and np.isfinite(rmse):
            return max(self.tol, self.adaptive_k * rmse)
        return self.tol

    def annotate(self, hits: pd.DataFrame) -> pd.DataFrame:
        """Adds rt_pred_qsar, rt_pred_rti, rt_diff_qsar/rti, rt_agree, rt_conflict, rt_available."""
        h = hits.copy()
        n = len(h)
        h["rt_pred_qsar"] = self.predict_qsar(h["smiles"].tolist()) if n else np.array([])
        col = "pred_rti_pos" if self.polarity == "positive" else "pred_rti_neg"
        h["rt_pred_rti"] = self.predict_rti(h[col].values) if (col in h.columns and n) else np.nan
        h["rt_diff_qsar"] = h["rt_apex"] - h["rt_pred_qsar"]
        h["rt_diff_rti"] = h["rt_apex"] - h["rt_pred_rti"]
        tq, tr = self.tolerance(self.qsar_rmse), self.tolerance(self.rti_rmse)
        aq = h["rt_diff_qsar"].abs() <= tq; ar = h["rt_diff_rti"].abs() <= tr
        have_q = h["rt_pred_qsar"].notna(); have_r = h["rt_pred_rti"].notna()
        h["rt_available"] = have_q | have_r
        h["rt_agree"] = (have_q & aq) | (have_r & ar)
        cq = have_q & (h["rt_diff_qsar"].abs() > self.conflict_factor * tq)
        cr = have_r & (h["rt_diff_rti"].abs() > self.conflict_factor * tr)
        h["rt_conflict"] = h["rt_available"] & ((cq | ~have_q) & (cr | ~have_r))
        return h

    def save(self, path):
        with open(path, "wb") as f: pickle.dump(self, f)

    @staticmethod
    def load(path):
        with open(path, "rb") as f: return pickle.load(f)


def _smiles_ik1(smi):
    try:
        from rdkit import Chem
        from rdkit.Chem.inchi import MolToInchiKey
        m = Chem.MolFromSmiles(smi)
        return MolToInchiKey(m).split("-")[0] if m is not None else None
    except Exception:
        return None


def known_standard_rt(hits: pd.DataFrame, standards: pd.DataFrame, polarity: str) -> pd.DataFrame:
    """Join authentic-standard RTs (same method) by InChIKey first block / name.  Adds std_rt, std_rt_diff."""
    h = hits.copy()
    col = "rt_pos" if polarity == "positive" else "rt_neg"
    if standards is None or col not in standards.columns:
        h["std_rt"] = np.nan; h["std_rt_diff"] = np.nan; return h
    s = standards.dropna(subset=[col]).copy()
    s["ik1"] = s["inchikey"].astype(str).str.split("-").str[0] if "inchikey" in s.columns else ""
    if "smiles" in s.columns:
        miss = s["ik1"].str.len() < 5
        s.loc[miss, "ik1"] = s.loc[miss, "smiles"].map(_smiles_ik1)
    m_ik = s.drop_duplicates("ik1").set_index("ik1")[col]
    h["std_rt"] = h["ik1"].map(m_ik)
    nm = s.assign(_n=s["name"].str.lower().str.strip()).drop_duplicates("_n").set_index("_n")[col]
    miss = h["std_rt"].isna()
    h.loc[miss, "std_rt"] = h.loc[miss, "name"].astype(str).str.lower().str.strip().map(nm)
    h["std_rt_diff"] = h["rt_apex"] - h["std_rt"]
    return h
