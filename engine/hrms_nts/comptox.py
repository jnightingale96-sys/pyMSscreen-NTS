"""Physico-chemical properties from the US EPA CompTox Chemicals Dashboard (CTX API) for semi-quantification.

API: https://comptox.epa.gov/ctx-api (free key from the CTX API team; set CTX_API_KEY or pass api_key). Endpoints used
(from the published OpenAPI spec): GET /chemical/search/equal/{word} (InChIKey / CAS / name -> DTXSID) and
GET /chemical/property/summary/search/by-dtxsid/{dtxsid} (experimental and predicted (OPERA etc.) medians per property).
Experimental medians are preferred over predicted ones; every value carries its source. Responses are cached on disk.
Without a key: a Dashboard batch-search export (CSV/XLSX, read by hrms_nts.properties) or RDKit estimates are used."""
from __future__ import annotations
import json, logging, os, time
from urllib.parse import quote
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")
BASE_URL = "https://comptox.epa.gov/ctx-api"
# property-name fragment (lower case) -> output column
PROPS = [("octanol-water", "logkow"), ("logkow", "logkow"), ("pka acidic", "pka_acid"), ("pka basic", "pka_base"),
         ("water solubility", "solubility_mol_l"), ("molar refractivity", "molar_refractivity"), ("polarizability", "polarizability"),
         ("henry", "henry_atm_m3_mol"), ("vapor pressure", "vp_mmhg"), ("boiling point", "bp_c"), ("logd", "logd")]


class CompTox:
    def __init__(self, api_key=None, cache_path=None, delay=0.15):
        self.key = api_key or os.environ.get("CTX_API_KEY") or ""
        self.cache_path = cache_path
        self.cache = {}
        if cache_path and os.path.exists(cache_path):
            try:
                self.cache = json.load(open(cache_path, encoding="utf-8"))
            except Exception:
                self.cache = {}
        self.delay = delay

    @property
    def available(self):
        return bool(self.key)

    def _get(self, path):
        import requests
        r = requests.get(BASE_URL + path, headers={"x-api-key": self.key, "accept": "application/json"}, timeout=30)
        time.sleep(self.delay)
        if r.status_code == 404:
            return None
        r.raise_for_status()
        return r.json()

    def save(self):
        if self.cache_path:
            os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
            json.dump(self.cache, open(self.cache_path, "w", encoding="utf-8"))

    def properties(self, identifier: str) -> dict:
        """identifier: InChIKey (preferred), CAS or name -> dict(dtxsid, preferred_name, logkow, ..., <col>_source)."""
        if not identifier:
            return {}
        if identifier in self.cache:
            return self.cache[identifier]
        out = {}
        try:
            hits = self._get(f"/chemical/search/equal/{quote(identifier, safe='')}") or []
            hits = hits if isinstance(hits, list) else [hits]
            if hits:
                h = hits[0]
                out["dtxsid"] = h.get("dtxsid"); out["preferred_name"] = h.get("preferredName")
                summ = self._get(f"/chemical/property/summary/search/by-dtxsid/{out['dtxsid']}") or []
                for p in summ if isinstance(summ, list) else [summ]:
                    nm = str(p.get("propName", "")).lower()
                    col = next((c for frag, c in PROPS if frag in nm), None)
                    if col is None or col in out:
                        continue
                    if (p.get("experimentalCount") or 0) > 0 and p.get("experimentalMedian") is not None:
                        out[col] = float(p["experimentalMedian"]); out[col + "_source"] = f"CompTox experimental median (n={p['experimentalCount']})"
                    elif p.get("predictedMedian") is not None:
                        out[col] = float(p["predictedMedian"]); out[col + "_source"] = f"CompTox predicted median (n={p.get('predictedCount')})"
        except Exception as e:
            log.warning("CompTox lookup failed for %s: %s", identifier, e)
            return out                      # not cached: retry next time
        self.cache[identifier] = out
        return out


def property_table(compounds: pd.DataFrame, api_key=None, cache_path=None, batch_export=None) -> pd.DataFrame:
    """compounds: name, smiles, inchikey (or ik1). -> one row per compound with logkow/pka/... and their sources.
    Order: CompTox API (if a key) -> Dashboard batch export (if given) -> RDKit (logP = Crippen, flagged 'RDKit estimate')."""
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    from rdkit.Chem.inchi import MolToInchiKey
    ct = CompTox(api_key, cache_path)
    exp = None
    if batch_export:
        from .properties import read_property_table
        exp = read_property_table(batch_export)
    rows = []
    for r in compounds.drop_duplicates(subset=["smiles"]).itertuples():
        m = Chem.MolFromSmiles(r.smiles) if isinstance(r.smiles, str) else None
        ik = getattr(r, "inchikey", "") or (MolToInchiKey(m) if m else "")
        row = dict(name=r.name, smiles=r.smiles, inchikey=ik, ik1=str(ik).split("-")[0])
        if ct.available and ik:
            row.update(ct.properties(ik))
        if exp is not None and not exp.empty and "logkow" not in row:
            e = exp[exp.get("ik1", pd.Series(dtype=str)) == row["ik1"]] if "ik1" in exp.columns else pd.DataFrame()
            if len(e) and pd.notna(e.iloc[0].get("logkow")):
                row["logkow"] = float(e.iloc[0]["logkow"]); row["logkow_source"] = "CompTox batch export"
                for c in ("pka_acid", "pka_base"):
                    if pd.notna(e.iloc[0].get(c)): row[c] = float(e.iloc[0][c])
        if m is not None:
            row.setdefault("mw", Descriptors.MolWt(m))
            row["tpsa"] = rdMolDescriptors.CalcTPSA(m)
            row["hba"] = rdMolDescriptors.CalcNumHBA(m); row["hbd"] = rdMolDescriptors.CalcNumHBD(m)
            if "logkow" not in row:
                row["logkow"] = Crippen.MolLogP(m); row["logkow_source"] = "RDKit Crippen estimate"
        rows.append(row)
    ct.save()
    df = pd.DataFrame(rows)
    if len(df):
        log.info("properties: %d compounds; logKow from %s", len(df), df["logkow_source"].value_counts().to_dict() if "logkow_source" in df else {})
    return df
