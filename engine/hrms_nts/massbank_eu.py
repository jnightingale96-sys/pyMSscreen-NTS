"""NORMAN MassBank-EU compound table (MassBankEU_Cmpds_*_wMS_DTXSIDs_*.csv): compounds with MEASURED MS/MS records and their fragment m/z lists.

File quirks handled here
  * the NAME column holds several synonyms joined by '|' ("Lidocaine|Lignocaine|2-(diethylamino)-N-...") and some rows start with a record-code token such as
    'LID_235.1805_10.1' or 'VEN_278.2114_12.5' (code_mz_RT of the contributing lab, NOT a name) - code tokens are dropped, the first real synonym becomes the name;
  * the same chemical appears in several rows (one per contributing institute) - merged on the first InChIKey block, fragments united;
  * 'Simplified_fragments_Positive/Negative' are ';'-separated m/z lists (no intensities) - used for fragment support, not for cosine scoring.
Use: independent, measured-fragment evidence for an identification (does the feature's MS2 contain the fragments MassBank recorded for that structure?)."""
from __future__ import annotations
import logging, re
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")
CODE_RE = re.compile(r"^[A-Za-z0-9]{2,8}_\d+(\.\d+)?(_\d+(\.\d+)?)?$")


def split_names(s) -> list[str]:
    """'LID_235.1805_10.1|Lidocaine| 2-(diethylamino)-N-(2,6-dimethylphenyl)acetamide' -> ['Lidocaine', '2-(diethylamino)-...']"""
    if not isinstance(s, str): return []
    out = []
    for t in s.split("|"):
        t = t.strip()
        if t and not CODE_RE.match(t) and t.lower() != "nan" and t not in out: out.append(t)
    return out


def _frags(series) -> np.ndarray:
    v = []
    for s in series.dropna().astype(str):
        for t in s.replace("|", ";").split(";"):
            try: v.append(float(t))
            except ValueError: pass
    return np.unique(np.round(np.array(v, float), 4)) if v else np.array([])


def load_massbank_eu(path: str) -> pd.DataFrame:
    d = pd.read_csv(path, dtype=str, encoding="utf-8-sig")
    d["ik1"] = d["InChIKey"].fillna("").str.split("-").str[0]
    d = d[d.ik1 != ""]
    rows = []
    for ik1, g in d.groupby("ik1", sort=False):
        names = []
        for n in g["NAME"]:
            for t in split_names(n):
                if t not in names: names.append(t)
        first = g.iloc[0]
        rows.append(dict(ik1=ik1, name=names[0] if names else "", aliases=" | ".join(names[1:12]), formula=first.get("FORMULA", ""), smiles=first.get("SMILES", ""),
                         exact_mass=float(first["NEUTRAL_EXACT_MASS"]) if pd.notna(first.get("NEUTRAL_EXACT_MASS")) else np.nan,
                         dtxsid=first.get("DTXSID", ""), cas=str(first.get("CAS_RN", "")).replace("CAS_RN:", "").strip(),
                         pos_frags=_frags(g["Simplified_fragments_Positive"]), neg_frags=_frags(g["Simplified_fragments_Negative"]),
                         n_rows=len(g), has_pos=g["MB_Records Positive Ionization"].notna().any(), has_neg=g["MB_Records Negative Ionization"].notna().any(),
                         institutes=";".join(sorted({x for s in list(g["Institute Positive Ionization"].dropna()) + list(g["Institute Negative Ionization"].dropna()) for x in str(s).split(";")}))))
    return pd.DataFrame(rows).set_index("ik1", drop=False)


def fragment_support(obs_mz, ref_mz, precursor_mz=None, tol_da=0.005, tol_ppm=10.0, precursor_excl=1.5) -> tuple[int, int]:
    """(n reference fragments found among the observed MS2 peaks, n reference fragments) - precursor and its isotopes excluded from the reference."""
    obs = np.sort(np.asarray(obs_mz, float)); ref = np.asarray(ref_mz, float)
    if precursor_mz is not None and len(ref): ref = ref[ref < precursor_mz - precursor_excl]
    if len(ref) == 0 or len(obs) == 0: return 0, len(ref)
    hit = 0
    for r in ref:
        tol = max(tol_da, r * tol_ppm * 1e-6); i = np.searchsorted(obs, [r - tol, r + tol])
        if i[1] > i[0]: hit += 1
    return hit, len(ref)
