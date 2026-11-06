"""IS-surrogate semi-quantification (engine 0.6.0).

Evidence (semiquant_eval, SCIEX X500R ESI+, 75 labelled drug IS, 2026-11): in a neat IS mix the molar response of
Good-peak IS varied only ~1.7-fold (SD 0.23 log); the median of all IS predicted a held-out IS within 1.49x (median)
and 3.2x (worst). Single-reference calibration and single-nearest-neighbour choices were the WORST options; descriptor
regressions did not beat the median. Matrix effects (up to 35x, RT/hydrophobicity dependent) dominated.

Method
1. is_response_factors : molar RF of every IS in the IS-mix injection(s) = area / (conc / MW_labelled).
   Only Good peaks, resolved (non-isobaric), outside the void (apex >= t0 + void_margin) are used.
2. physchem_table      : RDKit structure descriptors (MW, McGowan volume, TPSA, HBD, HBA, Crippen logP) plus,
   when a curated table is given (semiquant_eval 08_comptox_physchem_curated.csv format or any csv with
   smiles/inchikey + logKow + pKa_a/pKa_b), curated logKow and pKa -> cationic fraction at the mobile-phase pH.
3. semiquantify_is     : target molar RF = inverse-distance-weighted mean of the k nearest IS in z-scored descriptor
   space (k = 5 default; method='median' uses all IS). conc = area / RF * MW_target.
   Matrix correction: divide by the IS recovery in that sample (sample IS area / IS-mix area), median of the IS
   eluting within +-rt_window of the target (fallback: all IS in the sample).
   Uncertainty: x/÷ uf by distance to the nearest IS (applicability domain), x/÷ 10 for void-eluting targets and
   for ionisation classes absent from the IS set (e.g. acids / neutrals when the IS are bases).
4. rt_property_model   : column-specific log k = a + b*logP + c*TPSA + d*charged fitted on the IS; gives a 95 %
   prediction interval used as WEAK identification evidence (in the evaluation it kept ~90 % of true structures and
   rejected ~28 % of same-formula decoys).
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

MCG = {"C": 16.35, "H": 8.71, "N": 14.39, "O": 12.43, "F": 10.48, "Cl": 20.95, "Br": 26.21, "I": 34.53,
       "S": 22.91, "P": 24.87}
DESCRIPTORS = ["logKow", "TPSA", "McGowan_V", "MW", "HBD", "HBA", "f_cation"]
_PKA_DEFAULT = {"aliphatic_amine": (np.nan, 9.5), "zwitterion": (4.5, 9.5), "weak_base": (np.nan, 3.5),
                "acid": (4.5, np.nan), "neutral": (np.nan, np.nan)}


# ------------------------------------------------------------------------------------------------ descriptors
def _rdkit_desc(smiles: str) -> dict | None:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
    from rdkit.Chem.inchi import MolToInchiKey
    m = Chem.MolFromSmiles(str(smiles)) if smiles else None
    if m is None:
        return None
    mh = Chem.AddHs(m)
    v = (sum(MCG.get(a.GetSymbol(), 20.0) for a in mh.GetAtoms()) - 6.56 * mh.GetNumBonds()) / 100.0
    pat = lambda s: m.HasSubstructMatch(Chem.MolFromSmarts(s))
    amine = pat("[NX3;!$(N-a);!$(N-[C,S,P]=[O,N,S]);!$(N-S(=O));!$(N-[N,O]);!$(N-C#N)]") or \
        pat("[NX3;!$(N-C=O)]-C=[NX2;!a]")
    cooh = pat("C(=O)[OH1,O-]")
    acid = cooh or pat("[NH1]S(=O)(=O)")
    weak = pat("[nX2;r6]") or pat("[CX3;!$(C-[O,N]);!a]=[NX2;!a]") or pat("[NX3;H2,H1;!$(N-C=O);!$(N-S)]-a")
    cls = ("zwitterion" if amine and cooh else "aliphatic_amine" if amine else "acid" if acid else
           "weak_base" if weak else "neutral")
    return dict(MW=Descriptors.MolWt(m), logP=Crippen.MolLogP(m), TPSA=rdMolDescriptors.CalcTPSA(m), McGowan_V=v,
                HBD=Lipinski.NumHDonors(m), HBA=Lipinski.NumHAcceptors(m), ionisable_class=cls,
                ik1=MolToInchiKey(m).split("-")[0])


def _fractions(pka_a, pka_b, ph):
    a = 10 ** (ph - pka_a) if np.isfinite(pka_a) else 0.0
    b = 10 ** (pka_b - ph) if np.isfinite(pka_b) else 0.0
    den = 1 + a + b + a * b
    return b / den, 1 / den          # cationic, neutral


def physchem_table(compounds: pd.DataFrame, curated: str | pd.DataFrame | None = None, ph: float = 2.7) -> pd.DataFrame:
    """compounds: name + smiles. curated: optional csv/DataFrame with smiles or inchikey/ik_comptox + logKow and
    pKa_a/pKa_b (or pKa_a_final/pKa_b_final). Returns one row per name with the descriptor columns used here."""
    cur = None
    if curated is not None:
        cur = pd.read_csv(curated) if isinstance(curated, str) else curated.copy()
        ikc = next((c for c in ("ik_comptox", "inchikey", "InChIKey", "ik_list") if c in cur.columns), None)
        if ikc is not None:
            cur["ik1"] = cur[ikc].astype(str).str[:14]
        elif "smiles" in cur.columns:
            cur["ik1"] = cur.smiles.map(lambda s: (_rdkit_desc(s) or {}).get("ik1", ""))
        pa = "pKa_a_final" if "pKa_a_final" in cur.columns else "pKa_a"
        pb = "pKa_b_final" if "pKa_b_final" in cur.columns else "pKa_b"
        cur = cur.rename(columns={pa: "_pka_a", pb: "_pka_b"}).drop_duplicates("ik1").set_index("ik1")
    rows = []
    for r in compounds.itertuples():
        d = _rdkit_desc(getattr(r, "smiles", ""))
        if d is None:
            continue
        pa, pb = _PKA_DEFAULT[d["ionisable_class"]]
        logkow, src, psrc = d["logP"], "RDKit Crippen", "class default"
        if cur is not None and d["ik1"] in cur.index:
            c = cur.loc[d["ik1"]]
            if pd.notna(c.get("logKow", np.nan)):
                logkow, src = float(c["logKow"]), "curated"
            if "_pka_a" in c.index or "_pka_b" in c.index:
                pa, pb = float(c.get("_pka_a", np.nan)), float(c.get("_pka_b", np.nan))
                psrc = "curated"
        fc, fn = _fractions(pa, pb, ph)
        rows.append(dict(name=r.name, smiles=r.smiles, **d, logKow=logkow, logKow_source=src, pKa_a=pa, pKa_b=pb,
                         pKa_source=psrc, f_cation=fc, logD=logkow + np.log10(max(fn, 1e-12))))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ IS response factors
def is_response_factors(assigned: pd.DataFrame, comp: pd.DataFrame, is_mix_files, conc_ug_ml=None,
                        t0: float = 1.0, void_margin: float = 0.5, ratings=("Good",)) -> pd.DataFrame:
    """Molar RF of each IS in the IS-mix injection(s). conc_ug_ml: scalar (same for all IS), dict name->conc, or None
    to read comp['conc'] (list column 'conc', ug/mL). Returns name, smiles, rt, mw_labelled, area, logRF, use, reason."""
    a = assigned[assigned.file.isin(list(is_mix_files))].copy()
    if a.empty:
        raise ValueError(f"no targeted results for the IS-mix file(s) {list(is_mix_files)}")
    C = comp.set_index("name")
    rows = []
    for name, g in a.groupby("name"):
        if name not in C.index:
            continue
        c = C.loc[name]
        if isinstance(conc_ug_ml, dict):
            conc = conc_ug_ml.get(name, np.nan)
        elif conc_ug_ml is not None:
            conc = float(conc_ug_ml)
        else:
            conc = float(c.get("conc", np.nan)) if "conc" in C.columns else np.nan
        area = float(g.int_area.mean())
        rating = "/".join(g.int_rating.astype(str))
        rt = float(g.int_apex_rt.mean())
        level = "/".join(g.level.astype(str))
        why = []
        if not all(str(x) in ratings for x in g.int_rating):
            why.append(f"peak rating {rating}")
        if "unresolved" in level:
            why.append("isobaric, unresolved")
        if rt < t0 + void_margin:
            why.append(f"void (apex {rt:.2f} < {t0 + void_margin:.2f} min)")
        if not (np.isfinite(conc) and conc > 0 and area > 0):
            why.append("no concentration / area")
        mw = float(c.avg_mw)
        rows.append(dict(name=name, smiles=c.smiles, rt=rt, mw_labelled=mw, conc_ug_ml=conc, area=area, rating=rating,
                         level=level, logRF=np.log10(area / (conc / mw)) if (np.isfinite(conc) and conc > 0 and area > 0)
                         else np.nan,
                         use=not why, reason="; ".join(why)))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ matrix correction
def is_recovery(assigned: pd.DataFrame, is_rf: pd.DataFrame, is_mix_files) -> pd.DataFrame:
    """Per sample and IS: recovery = area in sample / area in the IS mix (only IS used for RFs)."""
    ok = is_rf[is_rf.use].set_index("name")
    a = assigned[assigned.name.isin(ok.index) & ~assigned.file.isin(list(is_mix_files))]
    out = a.assign(recovery=a.int_area / a.name.map(ok.area), rt=a.int_apex_rt)[["file", "name", "rt", "recovery"]]
    return out[np.isfinite(out.recovery) & (out.recovery > 0)]


def _local_recovery(rec: pd.DataFrame, sample: str, rt: float, window: float) -> tuple[float, int, str]:
    s = rec[rec.file == sample]
    if s.empty:
        return 1.0, 0, "no IS in sample - not corrected"
    loc = s[np.abs(s.rt - rt) <= window]
    if len(loc) >= 2:
        return float(10 ** np.median(np.log10(loc.recovery))), len(loc), f"IS within +-{window} min"
    return float(10 ** np.median(np.log10(s.recovery))), len(s), "all IS in sample (none near this RT)"


# ------------------------------------------------------------------------------------------------ prediction
def semiquantify_is(targets: pd.DataFrame, is_rf: pd.DataFrame, props: pd.DataFrame, method: str = "knn", k: int = 5,
                    recovery: pd.DataFrame | None = None, rt_window: float = 1.0, t0: float = 1.0,
                    void_margin: float = 0.5, uf: tuple = (3.0, 5.0, 10.0), ad_cuts: tuple = (1.0, 2.0)) -> pd.DataFrame:
    """targets: name, smiles, rt, area, sample. Returns targets + sq_conc_ug_ml (in the injected solution, before any
    matrix/enrichment factor), sq_low, sq_high, sq_rf_logpred, sq_nn_distance, sq_nearest_IS, sq_recovery, sq_flags."""
    cal = is_rf[is_rf.use].merge(props.drop(columns=["smiles"], errors="ignore"), on="name", how="inner")
    out = targets.copy()
    for c in ("sq_conc_ug_ml", "sq_low", "sq_high", "sq_rf_logpred", "sq_nn_distance", "sq_recovery"):
        out[c] = np.nan
    for c in ("sq_nearest_IS", "sq_method", "sq_flags", "sq_recovery_basis"):
        out[c] = ""
    if len(cal) < 3:
        log.warning("IS-surrogate semi-quant: fewer than 3 usable IS")
        return out
    P = props.set_index("name")
    X = cal[DESCRIPTORS].to_numpy(float)
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1
    Z = (X - mu) / sd
    y = cal.logRF.to_numpy()
    is_classes = set(cal.ionisable_class)
    for i, r in out.iterrows():
        if r.get("name") not in P.index or not np.isfinite(r.get("area", np.nan)) or r.get("area", 0) <= 0:
            continue
        p = P.loc[r["name"]]
        if isinstance(p, pd.DataFrame):
            p = p.iloc[0]
        z = (p[DESCRIPTORS].to_numpy(float) - mu) / sd
        dist = np.sqrt(((Z - z) ** 2).sum(1))
        j = np.argsort(dist)
        if method == "median":
            lr, meth = float(np.median(y)), "median of all IS"
        else:
            kk = j[:k]
            w = 1 / np.maximum(dist[kk], 1e-6)
            lr, meth = float(np.sum(w * y[kk]) / w.sum()), f"kNN k={k} (IDW)"
        flags = []
        dmin = float(dist[j[0]])
        f = uf[0] if dmin <= ad_cuts[0] else uf[1] if dmin <= ad_cuts[1] else uf[2]
        if dmin > ad_cuts[1]:
            flags.append("outside IS domain")
        rt = float(r.get("rt", np.nan))
        if np.isfinite(rt) and rt < t0 + void_margin:
            flags.append("void-eluting"); f = max(f, uf[2])
        if p.ionisable_class not in is_classes:
            flags.append(f"class '{p.ionisable_class}' not represented by the IS"); f = max(f, uf[2])
        rec, nrec, basis = (1.0, 0, "not corrected")
        if recovery is not None and "sample" in out.columns:
            rec, nrec, basis = _local_recovery(recovery, r["sample"], rt, rt_window)
        conc = r["area"] / 10 ** lr / rec * float(p.MW) / 1.0     # area / (area per umol/mL) -> umol/mL -> ug/mL
        out.at[i, "sq_conc_ug_ml"] = conc
        out.at[i, "sq_low"], out.at[i, "sq_high"] = conc / f, conc * f
        out.at[i, "sq_rf_logpred"] = round(lr, 3)
        out.at[i, "sq_nn_distance"] = round(dmin, 2)
        out.at[i, "sq_nearest_IS"] = cal.name.iloc[j[0]]
        out.at[i, "sq_recovery"] = round(rec, 3)
        out.at[i, "sq_recovery_basis"] = f"{basis} (n={nrec})"
        out.at[i, "sq_method"] = meth
        out.at[i, "sq_flags"] = "; ".join(flags)
    return out


# ------------------------------------------------------------------------------------------------ RT <-> properties
class RTPropertyModel:
    """Column/method-specific retention regression on the IS: log k = b0 + b.x, x = (logP, TPSA, charged).
    k = (tR - t0)/t0 (gradient: an empirical retention index for this method). predict_interval() gives the 95 %
    prediction interval; plausible() tests a candidate structure against an observed RT (weak evidence only)."""

    COLS = ["logP", "TPSA", "charged"]

    def __init__(self, t0: float = 1.0):
        self.t0 = t0
        self.coef = None

    def fit(self, rt, props: pd.DataFrame):
        X = self._x(props)
        lk = np.log10((np.asarray(rt, float) - self.t0) / self.t0)
        ok = np.isfinite(lk) & np.isfinite(X).all(1)
        X, lk = X[ok], lk[ok]
        A = np.column_stack([np.ones(len(X)), X])
        self.coef, *_ = np.linalg.lstsq(A, lk, rcond=None)
        res = lk - A @ self.coef
        self.dof = len(lk) - A.shape[1]
        self.s = float(np.sqrt(np.sum(res ** 2) / max(self.dof, 1)))
        self.XtX_inv = np.linalg.pinv(A.T @ A)
        self.n = len(lk)
        return self

    def _x(self, props):
        cl = props.get("ionisable_class", pd.Series([""] * len(props)))
        charged = np.isin(np.asarray(cl), ["aliphatic_amine", "zwitterion"]).astype(float)
        return np.column_stack([props["logP"].to_numpy(float), props["TPSA"].to_numpy(float), charged])

    def predict_interval(self, props: pd.DataFrame, level: float = 0.95):
        from scipy.stats import t as tdist
        A = np.column_stack([np.ones(len(props)), self._x(props)])
        pred = A @ self.coef
        se = self.s * np.sqrt(1 + np.einsum("ij,jk,ik->i", A, self.XtX_inv, A))
        h = tdist.ppf(0.5 + level / 2, max(self.dof, 1)) * se
        to_rt = lambda lk: self.t0 + self.t0 * 10 ** lk
        return pd.DataFrame(dict(rt_pred=to_rt(pred), rt_lo=to_rt(pred - h), rt_hi=to_rt(pred + h)), index=props.index)

    def plausible(self, rt_obs, props: pd.DataFrame, level: float = 0.95) -> pd.Series:
        iv = self.predict_interval(props, level)
        return (np.asarray(rt_obs) >= iv.rt_lo) & (np.asarray(rt_obs) <= iv.rt_hi)
