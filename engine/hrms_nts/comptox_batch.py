"""CompTox Chemicals Dashboard batch search for the project (no API key needed).

1. identifier_list(compounds)   names + InChIKeys of the compounds to look up (identified compounds, IS, calibrants)
2. fetch_export(...)             runs the Dashboard 'Batch Search' export automatically (the same request the web page sends:
                                 identifiers, OPERA/TEST predictions, 'Physicochemical Property Values' sheet), or
   the user runs it by hand at https://comptox.epa.gov/dashboard/batch-search and uploads the Excel file
3. curate(exports, compounds)   structure check (InChIKey skeleton must match the compound's SMILES - name searches can
                                 return a different structure), then per compound: logKow and pKa (experimental median >
                                 OPERA 2.8 > TEST), water solubility, logD5.5/7.4, with sources -> table read by
                                 surrogate.physchem_table (semi-quantification).
"""
from __future__ import annotations
import json
import logging
import time
import urllib.request
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")
BASE = "https://comptox.epa.gov/dashboard-api/batchsearch/export"
BATCH_URL = "https://comptox.epa.gov/dashboard/batch-search"
ITEMS = ["CASRN", "INCHIKEY", "SMILES", "QSAR_READY_SMILES", "MOLECULAR_FORMULA", "MONOISOTOPIC_MASS", "QC_LEVEL",
         "CHEMICAL_PROPERTIES_DETAILS", "OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED", "OPERA_PKAA_OPERA_PRED", "OPERA_PKAB_OPERA_PRED",
         "WATER_SOLUBILITY_MOL/L_OPERA_PRED", "LOGD5.5", "LOGD7.4"]
MANUAL_STEPS = (f"Manual alternative: open {BATCH_URL} -> tick 'Chemical Name' and 'InChIKey' -> paste the list -> "
                "'Choose export options' -> Excel -> tick 'CAS-RN', 'InChIKey', 'SMILES', 'OPERA Model Predictions' and "
                "'Physicochemical Property Values' -> Download export file -> upload that .xlsx in the notebook.")


def _ik(smiles):
    from rdkit import Chem
    m = Chem.MolFromSmiles(str(smiles)) if isinstance(smiles, str) and smiles else None
    return Chem.MolToInchiKey(m) if m is not None else ""


def identifier_list(compounds: pd.DataFrame) -> pd.DataFrame:
    """compounds: name, smiles -> name, smiles, inchikey (one row per unique structure)."""
    c = compounds[["name", "smiles"]].dropna().drop_duplicates("name").copy()
    c["inchikey"] = c.smiles.map(_ik)
    return c[c.inchikey != ""].drop_duplicates("inchikey").reset_index(drop=True)


def fetch_export(identifiers, id_type, out_path, timeout_s=300):
    """id_type: 'chemical_name' or 'INCHIKEY'. Writes the Dashboard Excel export to out_path."""
    body = {"identifierTypes": [id_type], "massError": 0, "downloadItems": ITEMS, "searchItems": "\n".join(identifiers),
            "inputType": "IDENTIFIER", "downloadType": "EXCEL"}
    hdr = {"Content-Type": "application/json", "Accept": "application/json", "User-Agent": "hrms-nts/0.6"}
    req = urllib.request.Request(BASE + "/", data=json.dumps(body).encode(), headers=hdr, method="POST")
    job = urllib.request.urlopen(req, timeout=120).read().decode().strip().strip('"')
    t = time.time()
    while time.time() - t < timeout_s:
        st = urllib.request.urlopen(urllib.request.Request(f"{BASE}/status/{job}", headers=hdr), timeout=60).read().decode()
        if "true" in st.lower():
            break
        time.sleep(4)
    data = urllib.request.urlopen(urllib.request.Request(f"{BASE}/content/{job}", headers={"User-Agent": hdr["User-Agent"]}),
                                  timeout=300).read()
    open(out_path, "wb").write(data)
    return out_path


def _read_export(path):
    x = pd.read_excel(path, sheet_name=None)
    main = next((v for k, v in x.items() if k.lower().startswith("main")), pd.DataFrame())
    props = next((v for k, v in x.items() if "propert" in k.lower()), pd.DataFrame())
    return main, props


def curate(export_paths, compounds: pd.DataFrame) -> pd.DataFrame:
    ids = identifier_list(compounds)
    ids["ik1"] = ids.inchikey.str[:14]
    mains, props = [], []
    for p in export_paths:
        m, pr = _read_export(p)
        mains.append(m); props.append(pr)
    main = pd.concat(mains, ignore_index=True) if mains else pd.DataFrame()
    prop = pd.concat(props, ignore_index=True).drop_duplicates() if props else pd.DataFrame()
    main = main[main.get("DTXSID").notna()] if len(main) else main
    main["ik1"] = main["INCHIKEY"].astype(str).str[:14] if len(main) else ""
    rows = []
    for _, c in ids.iterrows():
        hits = main[main.ik1 == c.ik1] if len(main) else main
        rec = dict(name=c["name"], smiles=c.smiles, inchikey=c.inchikey)
        if hits.empty:
            rec["comptox_status"] = "not found (or only a different structure under this name)"
            rows.append(rec); continue
        dtx = sorted(hits.DTXSID.unique())
        P = prop[prop.DTXSID.isin(dtx)] if len(prop) else prop
        nexp = P[P.TYPE == "experimental"].groupby("DTXSID").size() if len(P) else pd.Series(dtype=int)
        best = max(dtx, key=lambda d: nexp.get(d, 0))
        h = hits[hits.DTXSID == best].iloc[0]
        rec.update(DTXSID=best, preferred_name=h.get("PREFERRED_NAME"), CASRN=h.get("CASRN"), ik_comptox=h.get("INCHIKEY"),
                   comptox_status="OK (structure-checked)")
        for pname, short, opera_col in (("LogKow: Octanol-Water", "logKow", "OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED"),
                                         ("pKa Acidic Apparent", "pKa_a", "OPERA_PKAA_OPERA_PRED"),
                                         ("pKa Basic Apparent", "pKa_b", "OPERA_PKAB_OPERA_PRED")):
            e = P[(P.NAME == pname) & (P.TYPE == "experimental")].VALUE.astype(float) if len(P) else pd.Series(dtype=float)
            op = pd.to_numeric(h.get(opera_col), errors="coerce")
            if len(e):
                rec[short], rec[f"{short}_source"], rec[f"{short}_n_exp"] = float(e.median()), "experimental (median)", len(e)
            elif np.isfinite(op):
                rec[short], rec[f"{short}_source"] = float(op), "OPERA 2.8 prediction"
        ws = pd.to_numeric(h.get("WATER_SOLUBILITY_MOL/L_OPERA_PRED"), errors="coerce")
        rec["logS_molL"] = float(np.log10(ws)) if np.isfinite(ws) and ws > 0 else np.nan
        rec["logD5.5"] = pd.to_numeric(h.get("LOGD5.5"), errors="coerce")
        rec["logD7.4"] = pd.to_numeric(h.get("LOGD7.4"), errors="coerce")
        rows.append(rec)
    return pd.DataFrame(rows)
