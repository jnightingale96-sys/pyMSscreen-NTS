"""Suspect / candidate databases: NORMAN SusDat, PubChemLite, user list (pyMSscreen-style compounds.csv)."""
from __future__ import annotations
import logging, os
import numpy as np
import pandas as pd
from .common import formula_mass

log = logging.getLogger("hrms_nts")

STD_COLS = ["db_id", "name", "formula", "exact_mass", "smiles", "inchikey", "ik1", "cas", "source",
            "pred_rti_pos", "pred_rti_neg", "direct_precursor_mz", "db_level"]


def _finish(df: pd.DataFrame) -> pd.DataFrame:
    for c in STD_COLS:
        if c not in df.columns:
            df[c] = np.nan if c in ("exact_mass", "pred_rti_pos", "pred_rti_neg", "direct_precursor_mz") else ""
    df["exact_mass"] = pd.to_numeric(df["exact_mass"], errors="coerce")
    df["inchikey"] = df["inchikey"].fillna("").astype(str)
    df["ik1"] = df["inchikey"].str.split("-").str[0]
    return df[STD_COLS]


def load_susdat(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", engine="python", dtype=str)
    ren = {"Norman_SusDat_ID": "db_id", "Name": "name", "Molecular_Formula": "formula", "Monoiso_Mass": "exact_mass",
           "SMILES": "smiles", "StdInChIKey": "inchikey", "CAS_RN_Dashboard": "cas", "Validation_Level": "db_level",
           "Pred_RTI_Positive_ESI": "pred_rti_pos", "Pred_RTI_Negative_ESI": "pred_rti_neg"}
    df = df.rename(columns=ren)
    for c in ("pred_rti_pos", "pred_rti_neg"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["source"] = "NORMAN-SusDat"
    df = df.dropna(subset=["exact_mass"])
    return _finish(df)


def load_pubchemlite(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    ren = {"Identifier": "db_id", "CompoundName": "name", "MolecularFormula": "formula", "MonoisotopicMass": "exact_mass",
           "SMILES": "smiles", "InChIKey": "inchikey", "Synonym": "cas"}
    df = df.rename(columns=ren)
    df["cas"] = df["cas"].astype(str).str.split(r"\s|\|").str[0]
    df["db_id"] = "PCL" + df["db_id"].astype(str)
    df["source"] = "PubChemLite"
    return _finish(df.dropna(subset=["exact_mass"]))


def load_sle(path: str) -> pd.DataFrame:
    """NORMAN Suspect List Exchange structures (tools/build_reference_index.py -> sle_candidates.csv.gz)."""
    df = pd.read_csv(path, dtype=str)
    df["db_id"] = "SLE:" + df["lists"].fillna("").str.split(";").str[0]
    df["source"] = "NORMAN-SLE"
    df["cas"] = ""
    return _finish(df.dropna(subset=["exact_mass"]))


def load_user_compounds(path: str) -> pd.DataFrame:
    """pyMSscreen-style compounds.csv: ID, Name, SMILES, mz (neutral mass), precursor_mz (any of the last three)."""
    df = pd.read_csv(path, dtype=str)
    low = {c.lower().strip(): c for c in df.columns}
    g = lambda k: df[low[k]] if k in low else pd.Series([""] * len(df))
    out = pd.DataFrame({"db_id": g("id"), "name": g("name"), "smiles": g("smiles"),
                        "inchikey": g("inchikey") if "inchikey" in low else ""})
    out["formula"] = g("formula") if "formula" in low else ""
    em = pd.to_numeric(g("mz") if "mz" in low else g("exact_mass"), errors="coerce")
    out["direct_precursor_mz"] = pd.to_numeric(g("precursor_mz") if "precursor_mz" in low else g("precursormz"), errors="coerce")
    # SMILES takes priority over mz (pyMSscreen rule)
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        from rdkit.Chem.inchi import MolToInchiKey
        masses, forms, iks = [], [], []
        for smi, m0, f0, k0 in zip(out["smiles"], em, out["formula"], out["inchikey"]):
            mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) and smi.strip() else None
            if mol is not None:
                masses.append(Descriptors.ExactMolWt(mol)); forms.append(rdMolDescriptors.CalcMolFormula(mol))
                iks.append(k0 if isinstance(k0, str) and k0 else MolToInchiKey(mol))
            else:
                masses.append(m0); forms.append(f0); iks.append(k0)
        out["exact_mass"] = masses; out["formula"] = forms; out["inchikey"] = iks
    except Exception as e:  # rdkit missing -> use supplied masses only
        log.warning("RDKit unavailable for user compounds (%s); using supplied masses", e)
        out["exact_mass"] = em
    out["source"] = "user"
    out["cas"] = ""
    out = out[(out["exact_mass"].notna()) | (out["direct_precursor_mz"].notna())]
    return _finish(out)


class CompoundDB:
    """Merged, mass-sorted candidate table with fast ppm range search."""

    def __init__(self, df: pd.DataFrame):
        df = df.copy()
        df["has_ik"] = df["inchikey"].str.len() > 10
        df = df.sort_values(["has_ik", "source"], ascending=[False, True])
        # merge duplicates by full InChIKey (keeps first record, concatenates source names)
        dup_key = np.where(df["has_ik"], df["inchikey"], df["db_id"].astype(str) + "|" + df["name"].astype(str))
        df["_k"] = dup_key
        src = df.groupby("_k")["source"].agg(lambda s: "+".join(sorted(set(s))))
        df = df.drop_duplicates("_k", keep="first").copy()
        df["source"] = df["_k"].map(src)
        df = df.drop(columns=["_k", "has_ik"])
        df = df[df["exact_mass"].notna() | df["direct_precursor_mz"].notna()].sort_values("exact_mass", na_position="last").reset_index(drop=True)
        self.df = df
        self.masses = df["exact_mass"].values
        log.info("candidate database: %d unique compounds (%s)", len(df), ", ".join(sorted(set("+".join(df['source']).split("+")))))

    def search(self, neutral_mass: np.ndarray, ppm: float):
        """For an array of neutral masses -> (query_idx, db_idx) pairs within ppm."""
        m = np.asarray(neutral_mass, dtype=float)
        valid = self.masses[~np.isnan(self.masses)]
        n_valid = len(valid)  # NaNs are sorted last
        lo = np.searchsorted(self.masses[:n_valid], m * (1 - ppm * 1e-6), side="left")
        hi = np.searchsorted(self.masses[:n_valid], m * (1 + ppm * 1e-6), side="right")
        cnt = hi - lo
        qi = np.repeat(np.arange(len(m)), cnt)
        di = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi) if b > a]) if cnt.sum() else np.array([], int)
        return qi, di


def build_database(cfg: dict) -> CompoundDB:
    from .common import resolve
    frames = []
    p = cfg["paths"]
    for key, fn in (("susdat", load_susdat), ("pubchemlite", load_pubchemlite), ("sle_candidates", load_sle)):
        path = resolve(cfg, p.get(key))
        if path and os.path.exists(path):
            frames.append(fn(path)); log.info("loaded %s: %d", key, len(frames[-1]))
        elif path:
            log.warning("%s not found: %s", key, path)
    upath = resolve(cfg, p.get("user_compounds"))
    if upath and os.path.exists(upath):
        frames.append(load_user_compounds(upath)); log.info("loaded user compounds: %d", len(frames[-1]))
    if not frames:
        raise FileNotFoundError("No candidate database found - set paths.susdat / paths.pubchemlite / paths.user_compounds in config.yaml")
    return CompoundDB(pd.concat(frames, ignore_index=True))
