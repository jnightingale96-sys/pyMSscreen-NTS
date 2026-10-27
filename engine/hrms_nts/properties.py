"""Chemical properties for the fate model (PEC) of screening candidates.

Sources, in priority order
  1. a properties table (CSV/XLSX) supplied by the user - own columns (see PROP_COLS) and/or a CompTox Dashboard batch export
     (OPERA predictions: AVERAGE_MASS, VAPOR_PRESSURE_MMHG_*, WATER_SOLUBILITY_MOL/L_*, HENRYS_LAW_ATM-M3/MOLE_*, OCTANOL_WATER_*,
     OPERA_PKAA/PKAB_*, BIODEGRADATION_HALF_LIFE_DAYS_*), matched by InChIKey first block, then CAS, then name;
  2. structure-based estimates from the SMILES (RDKit Crippen logP, ionisable-group SMARTS) - always flagged 'estimated'.
A chemical whose table row is incomplete (no Kaw / DT50 / logKow) is run over a property GRID and reported as a PEC range,
never as a single invented value."""
from __future__ import annotations
import logging, os, re
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")

PROP_COLS = ["name", "ik1", "cas", "smiles", "mw", "logkow", "pka_acid", "pka_base", "type", "kaw", "vp_pa", "solubility_mg_L",
             "henry_atm_m3_mol", "biowin4", "dt50_soil_d", "k_fruit", "source"]
R_ATM_T = 8.2057e-5 * 298.15            # m3 atm / mol at 25 C
R_PA_T = 8.314 * 298.15                  # J / mol
MMHG_TO_PA = 133.322
SOIL_PH = 7.2

# CompTox/OPERA export column -> internal column (matched on the column-name prefix, case-insensitive)
_COMPTOX = [("INCHIKEY", "inchikey"), ("CASRN", "cas"), ("PREFERRED_NAME", "name"), ("INPUT", "input"), ("SMILES", "smiles"),
            ("AVERAGE_MASS", "mw"), ("OCTANOL_WATER_PARTITION_COEFF_LOGKOW_OPERA", "logkow"), ("OCTANOL_WATER_PARTITION", "logkow"),
            ("OPERA_PKAA", "pka_acid"), ("OPERA_PKAB", "pka_base"), ("HENRYS_LAW_ATM-M3/MOLE_OPERA", "henry_atm_m3_mol"),
            ("HENRYS_LAW_ATM-M3/MOLE", "henry_atm_m3_mol"), ("VAPOR_PRESSURE_MMHG_OPERA", "vp_mmhg"), ("VAPOR_PRESSURE_MMHG", "vp_mmhg"),
            ("WATER_SOLUBILITY_MOL/L_OPERA", "sol_mol_l"), ("WATER_SOLUBILITY_MOL/L", "sol_mol_l"),
            ("BIODEGRADATION_HALF_LIFE_DAYS", "opera_biodeg_d")]


def _num(x):
    if x is None: return np.nan
    s = str(x).strip().replace(",", ".")
    if s in ("", "-", "nan", "NaN", "N/A", "#N/A"): return np.nan
    m = re.fullmatch(r"([0-9.]+)E10-([0-9]+)", s, re.I)               # '3.17E10-5' typed for 3.17E-5 in the project workbook
    if m: return float(m.group(1)) * 10 ** (-int(m.group(2)))
    try: return float(s)
    except ValueError: return np.nan


def read_property_tables(paths, default_sheet: str | None = None) -> pd.DataFrame | None:
    """Several tables, FIRST match wins (put curated project values before predicted ones).  'file.xlsx::Sheet' selects a sheet."""
    if not paths: return None
    paths = [paths] if isinstance(paths, str) else list(paths)
    parts = []
    for p in paths:
        f, _, sh = str(p).partition("::")
        parts.append(read_property_table(f, sh or default_sheet))
    return pd.concat(parts, ignore_index=True) if parts else None


def read_property_table(path: str, sheet: str | None = None) -> pd.DataFrame:
    """Read a properties CSV/XLSX in the toolkit format or a CompTox batch export; returns the internal column set."""
    if not path or not os.path.exists(path):
        return pd.DataFrame(columns=PROP_COLS)
    if path.lower().endswith((".xlsx", ".xlsm", ".xls")):
        raw = pd.read_excel(path, sheet_name=sheet or 0, dtype=object)
    else:
        from .common import read_table
        raw = read_table(path).astype(object)
    cols = {}
    for c in raw.columns:
        cu = str(c).strip().upper().replace(" ", "_")
        for pref, tgt in _COMPTOX:
            if cu.startswith(pref) and tgt not in cols.values():
                cols[c] = tgt; break
        else:
            if str(c).strip().lower() in PROP_COLS:
                cols[c] = str(c).strip().lower()
    df = raw.rename(columns=cols)
    df = df.loc[:, ~df.columns.duplicated()]
    if "ik1" not in df.columns:
        df["ik1"] = df["inchikey"].astype(str).str.split("-").str[0] if "inchikey" in df.columns else ""
    for c in ("mw", "logkow", "pka_acid", "pka_base", "kaw", "vp_pa", "solubility_mg_L", "henry_atm_m3_mol", "biowin4", "dt50_soil_d",
              "k_fruit", "vp_mmhg", "sol_mol_l", "opera_biodeg_d"):
        df[c] = df[c].map(_num) if c in df.columns else np.nan
    if "vp_pa" in df.columns:
        df["vp_pa"] = df["vp_pa"].fillna(df["vp_mmhg"] * MMHG_TO_PA)
    df["solubility_mg_L"] = df["solubility_mg_L"].fillna(df["sol_mol_l"] * df["mw"] * 1000.0)
    for c in ("name", "cas", "smiles", "type", "source"):
        if c not in df.columns: df[c] = ""
        df[c] = df[c].fillna("").astype(str)
    df["ik1"] = df["ik1"].fillna("").astype(str).str.strip()
    df["source"] = df["source"].where(df["source"] != "", os.path.basename(path))
    return df


def _lookup(table: pd.DataFrame, ik1: str, cas: str, name: str):
    if table is None or table.empty: return None
    for col, key in (("ik1", ik1), ("cas", cas)):
        if key and col in table.columns:
            m = table[table[col].astype(str).str.lower() == str(key).lower()]
            if len(m): return m.iloc[0]
    if name:
        for col in ("name", "input"):
            if col in table.columns:
                m = table[table[col].astype(str).str.lower() == name.lower()]
                if len(m): return m.iloc[0]
    return None


# ----------------------------------------------------------------------------------------------------------------------------------
# structure-based estimates (always flagged)
# ----------------------------------------------------------------------------------------------------------------------------------
_ACID = ["[CX3](=O)[OX2H1,OX1-]", "[SX4](=O)(=O)[OX2H1,OX1-]", "[PX4](=O)[OX2H1,OX1-]"]
_BASE = ["[NX3;H2,H1,H0;!$(NC=O);!$(N-a);!$(N-S(=O)=O);!$(N-C=[N,O,S]);!$(N#*);!$(N=*)]([#6])"]
_QUAT = ["[NX4+;!$([N+][O-])]", "[n+;!$([n+][O-])]", "[P+]", "[S+;!$([S+][O-])]"]
# zwitterions: an acid (carboxylic / sulfonic / phosphonic / tetrazole) together with a base that is protonated at pH ~7 (aliphatic amine, amidine, guanidine),
# or a charge-separated N-oxide / betaine.  Table 'type' entries only know Acid/Base/Neutral and would hide these.
_ACID_Z = _ACID + ["[#6]c1nn[nH]n1", "[#6]c1n[nH]nn1"]
_BASE_Z = _BASE + ["[NX3;!$(NC=O);!$(N-a)]C(=[NX2;!$(N-a);!$(N-C=O)])", "[NX3]C(=[NX2])[NX3]", "[NX3;!$(NC=O)]C(=[NX2;H1,H0;!$(N-a)])[#6,#7]", "[NX3;H2]C(=[NH])"]
_NOXIDE = ["[NX4+;!R]([#6])([#6])([#6])[O-]", "[N+;!$([N+]=O);!$([N+]([O-])=O)][O-]"]


def structure_flags(smiles: str) -> dict:
    out = dict(has_acid=False, has_base=False, quaternary=False, zwitterion=False, mw=np.nan, logkow_est=np.nan, parsed=False)
    if not isinstance(smiles, str) or not smiles.strip(): return out
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors
        m = Chem.MolFromSmiles(smiles)
        if m is None: return out
        frags = sorted(Chem.GetMolFrags(m, asMols=True), key=lambda x: x.GetNumHeavyAtoms(), reverse=True)
        m = frags[0]
        has = lambda pats: any(m.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in pats)
        acid_z, base_z = has(_ACID_Z), has(_BASE_Z)
        betaine = has(_QUAT) and acid_z
        n_oxide = any(m.HasSubstructMatch(Chem.MolFromSmarts(q)) and not m.HasSubstructMatch(Chem.MolFromSmarts("[N+](=O)[O-]")) for q in _NOXIDE[:1])
        out["zwitterion"] = bool((acid_z and base_z) or betaine or n_oxide)
        out.update(has_acid=has(_ACID), has_base=has(_BASE), quaternary=has(_QUAT), mw=float(Descriptors.MolWt(m)),
                   logkow_est=float(Crippen.MolLogP(m)), parsed=True)
    except Exception as e:
        log.debug("structure_flags failed for %s: %s", smiles, e)
    return out


def species_type(pka_acid, pka_base, ph=SOIL_PH):
    """'Acid' | 'Base' | 'Neutral' | 'Zwitterion' from % ionised at soil pH (pKa of the conjugate acid for bases)."""
    a = 1 / (1 + 10 ** (pka_acid - ph)) if np.isfinite(pka_acid) else 0.0
    b = 1 / (1 + 10 ** (ph - pka_base)) if np.isfinite(pka_base) else 0.0
    if a > 0.5 and b > 0.5: return "Zwitterion"
    if a > 0.5: return "Acid"
    if b > 0.5: return "Base"
    return "Neutral"


def resolve_properties(cands: pd.DataFrame, table: pd.DataFrame | None, cfg_pec: dict) -> pd.DataFrame:
    """One row per candidate: the property set the fate model needs + 'prop_status' (table | table-partial | estimated | not-applicable)
    and the property GRID to use when something is missing.  cands needs: cand_id, name, smiles, ik1 (cas optional)."""
    grid = cfg_pec.get("property_grid", {})
    d_logkow = grid.get("logkow_delta", [-0.5, 0.0, 0.5]); kaw_g = grid.get("kaw", [1e-9, 1e-6, 1e-4]); dt50_g = grid.get("dt50_soil_d", [10, 60, 365, 1000])
    rows = []
    for c in cands.itertuples():
        t = _lookup(table, getattr(c, "ik1", ""), getattr(c, "cas", ""), getattr(c, "name", ""))
        sf = structure_flags(getattr(c, "smiles", ""))
        g = (lambda k: float(t[k]) if (t is not None and k in t.index and pd.notna(t[k])) else np.nan)
        mw = g("mw"); mw = sf["mw"] if not np.isfinite(mw) else mw
        logkow = g("logkow"); lk_src = "table"
        if not np.isfinite(logkow): logkow, lk_src = sf["logkow_est"], "estimated (Crippen)"
        pka_a, pka_b = g("pka_acid"), g("pka_base")
        ttype = str(t["type"]) if (t is not None and str(t["type"]) in ("Acid", "Base", "Neutral")) else ""
        t_src = "table" if ttype else ""
        if not ttype:
            if np.isfinite(pka_a) or np.isfinite(pka_b):
                ttype, t_src = species_type(pka_a, pka_b), "table pKa"
            else:
                a, b = sf["has_acid"], sf["has_base"]
                ttype = "Zwitterion" if (a and b) else "Acid" if a else "Base" if b else "Neutral"; t_src = "estimated (SMARTS)"
                pka_a = 4.5 if a else np.nan; pka_b = 9.5 if b else np.nan
        if sf["zwitterion"] and ttype != "Zwitterion": ttype, t_src = "Zwitterion", (t_src + " + structure (acid and base / N-oxide)").strip(" +")
        # Kaw
        kaw = g("kaw"); kaw_src = "table"
        if not np.isfinite(kaw):
            h = g("henry_atm_m3_mol")
            if np.isfinite(h): kaw, kaw_src = h / R_ATM_T, "Henry (table)"
            else:
                vp, sol = g("vp_pa"), g("solubility_mg_L")
                if np.isfinite(vp) and np.isfinite(sol) and sol > 0 and np.isfinite(mw): kaw, kaw_src = vp * mw / (sol * R_PA_T), "VP/solubility (table)"
                else: kaw_src = ""
        # degradation
        bio, dt50 = g("biowin4"), g("dt50_soil_d"); dsrc = "BIOWIN4" if np.isfinite(bio) else ("table DT50" if np.isfinite(dt50) else "")
        if not dsrc and t is not None and "opera_biodeg_d" in t.index and pd.notna(t["opera_biodeg_d"]):
            dt50, dsrc = float(t["opera_biodeg_d"]), "OPERA biodegradation half-life used as soil DT50"
        status = "table" if (t is not None and np.isfinite(logkow) and np.isfinite(kaw) and dsrc) else ("table-partial" if t is not None else "estimated")
        na = ""
        if sf["quaternary"]: na = "permanent charge (quaternary) - fate model not applicable"
        elif ttype == "Zwitterion": na = "zwitterion - fate model not applicable"
        elif not np.isfinite(mw) or not np.isfinite(logkow): na = "no structure/properties"
        if na: status = "not-applicable"
        rows.append(dict(cand_id=c.cand_id, mw=mw, logkow=logkow, logkow_src=lk_src, pka_acid=pka_a, pka_base=pka_b, species=ttype, species_src=t_src,
                         kaw=kaw, kaw_src=kaw_src, biowin4=bio, dt50_soil_d=dt50, dt50_src=dsrc, k_fruit_user=g("k_fruit"), prop_status=status, not_applicable=na,
                         grid_logkow=";".join(f"{logkow + dl:.2f}" for dl in d_logkow) if np.isfinite(logkow) else "",
                         grid_kaw=";".join(f"{k:g}" for k in kaw_g), grid_dt50=";".join(f"{x:g}" for x in dt50_g),
                         prop_source=(str(t["source"]) if t is not None else "SMILES-based estimate")))
    return pd.DataFrame(rows)


def parse_grid(s) -> list[float]:
    return [float(x) for x in str(s).split(";") if x != ""]
