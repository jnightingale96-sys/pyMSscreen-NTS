"""Build the reference index from all free structure databases (one row per InChIKey first block).

    python tools/build_reference_index.py <databases_dir> <out_dir>

<databases_dir> is the folder written by the database installer (layout: norman_sle/ coconut/ lotus/ npatlas/ chebi/
lipidmaps/ pubchemlite/ ...). Writes to <out_dir>:
    reference_index.csv.gz    ik1, name, norman_lists, n_norman_lists, known_tp, np_sources, lotus_n_organisms,
                              npatlas_origin, chebi_roles, lipid_category, lipid_class
    tp_pairs.csv.gz           ik1_tp, tp_name, parent_name, parent_ik1, transformation, source  (literature TP -> parent)
    sle_candidates.csv.gz     NORMAN-SLE structures (name, formula, exact_mass, smiles, inchikey, lists) for the candidate
                              database (most are also in SusDat; the loader de-duplicates by InChIKey)
Sources and licences: see DATA_SOURCES.md. Everything here is free to use; NP Atlas is CC BY-NC (academic use).
"""
from __future__ import annotations

import glob
import gzip
import os
import re
import sys
import zipfile
from collections import defaultdict

import pandas as pd

TP_LISTS = {"S038", "S066", "S068", "S074", "S078", "S081", "S121"}       # transformation-product collections
CHEBI_ROLES = {"CHEBI:76924": "plant metabolite", "CHEBI:77746": "human metabolite", "CHEBI:76946": "fungal metabolite",
               "CHEBI:76969": "bacterial metabolite", "CHEBI:75771": "mouse metabolite", "CHEBI:76507": "marine metabolite"}


def _rdkit():
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    return Chem


def ik_from_smiles(smiles, cache={}):
    Chem = _rdkit()
    out = []
    for s in smiles:
        if not isinstance(s, str) or not s.strip():
            out.append(""); continue
        if s not in cache:
            m = Chem.MolFromSmiles(s)
            cache[s] = Chem.MolToInchiKey(m) if m is not None else ""
        out.append(cache[s])
    return out


def _col(cols, *names):
    low = {c.lower().replace(" ", "").replace("_", ""): c for c in cols}
    for n in names:
        if n in low:
            return low[n]
    return None


def read_csv_any(p, **kw):
    for enc in ("utf-8", "latin-1"):
        try:
            return pd.read_csv(p, encoding=enc, low_memory=False, on_bad_lines="skip", **kw)
        except UnicodeDecodeError:
            continue
        except Exception:
            return None
    return None


# ------------------------------------------------------------------------------------------------ NORMAN SLE
def norman_sle(d):
    members = defaultdict(set)          # ik1 -> {S066, ...}
    names, rows, pairs = {}, [], []
    cid_ik = {}
    for folder in sorted(glob.glob(os.path.join(d, "norman_sle", "S*"))):
        code = os.path.basename(folder).split("_")[0]
        for p in sorted(glob.glob(os.path.join(folder, "*.csv"))):
            df = read_csv_any(p)
            if df is None or df.empty:
                continue
            ikc = _col(df.columns, "inchikey", "stdinchikey", "standardinchikey")
            smc = _col(df.columns, "smiles", "isomericsmiles", "canonicalsmiles")
            nmc = _col(df.columns, "name", "compoundname", "chemicalname", "cpd")
            cidc = _col(df.columns, "pubchemcid", "cid")
            if ikc is None and smc is None:
                continue
            if "transformation" in os.path.basename(p).lower():
                continue                                                   # handled below
            iks = df[ikc].astype(str) if ikc else pd.Series(ik_from_smiles(df[smc]), index=df.index)
            iks = iks.where(iks.str.len() >= 14, "")
            for i, ik in iks.items():
                if not ik:
                    continue
                k1 = ik[:14]
                members[k1].add(code)
                if nmc and k1 not in names and isinstance(df.at[i, nmc], str):
                    names[k1] = df.at[i, nmc]
                if cidc is not None and pd.notna(df.at[i, cidc]):
                    try:
                        cid_ik[int(float(df.at[i, cidc]))] = k1
                    except Exception:
                        pass
                if smc is not None:
                    rows.append((ik, df.at[i, smc], names.get(k1, ""), code))
        # transformation tables: predecessor -> successor
        for p in glob.glob(os.path.join(folder, "*ransformation*.csv")):
            t = read_csv_any(p)
            if t is None or t.empty:
                continue
            ps, ss = _col(t.columns, "predecessorsmiles"), _col(t.columns, "successorsmiles")
            pc, sc = _col(t.columns, "predecessorcid"), _col(t.columns, "successorcid")
            pn, sn = _col(t.columns, "predecessorname"), _col(t.columns, "successorname")
            tr = _col(t.columns, "transformation")
            pik = ik_from_smiles(t[ps]) if ps else [cid_ik.get(int(float(x)), "") if pd.notna(x) else "" for x in t[pc]] if pc else []
            sik = ik_from_smiles(t[ss]) if ss else [cid_ik.get(int(float(x)), "") if pd.notna(x) else "" for x in t[sc]] if sc else []
            for j, (a, b) in enumerate(zip(pik, sik)):
                if b:
                    pairs.append(dict(ik1_tp=b[:14], tp_name=t[sn].iloc[j] if sn else "", parent_name=t[pn].iloc[j] if pn else "",
                                      parent_ik1=a[:14] if a else "", transformation=t[tr].iloc[j] if tr else "", source=code))
                    members[b[:14]].add(code)
    # Eawag S66: ParentOrTP with Name_of_pair
    p = glob.glob(os.path.join(d, "norman_sle", "S066_*", "*.csv"))
    if p:
        e = read_csv_any(p[0])
        if e is not None and "ParentOrTP" in e:
            tp = e[e.ParentOrTP.astype(str).str.upper().str.startswith("TP")]
            for _, r in tp.iterrows():
                if isinstance(r.get("InChIKey"), str):
                    pairs.append(dict(ik1_tp=r.InChIKey[:14], tp_name=r.get("Name", ""), parent_name=r.get("Name_of_pair", ""),
                                      parent_ik1="", transformation=r.get("Transformation", ""), source="S066"))
    cand = pd.DataFrame(rows, columns=["inchikey", "smiles", "name", "list"]).drop_duplicates("inchikey")
    return members, names, cand, pd.DataFrame(pairs).drop_duplicates()


# ------------------------------------------------------------------------------------------------ natural products, lipids, ChEBI
def coconut(d):
    z = glob.glob(os.path.join(d, "coconut", "*.zip"))
    if not z:
        return set(), {}
    zf = zipfile.ZipFile(z[0])
    ks, nm = set(), {}
    for ch in pd.read_csv(zf.open(zf.namelist()[0]), usecols=["standard_inchi_key", "name"], chunksize=200_000):
        for k, n in zip(ch.standard_inchi_key.astype(str), ch["name"]):
            k1 = k[:14]
            ks.add(k1)
            if isinstance(n, str) and k1 not in nm:
                nm[k1] = n
    return ks, nm


def lotus(d):
    p = glob.glob(os.path.join(d, "lotus", "*frozen.csv.gz"))
    if not p:
        return {}
    df = pd.read_csv(p[0], usecols=["structure_inchikey", "organism_name"])
    df["ik1"] = df.structure_inchikey.astype(str).str[:14]
    return df.groupby("ik1").organism_name.nunique().to_dict()


def npatlas(d):
    p = glob.glob(os.path.join(d, "npatlas", "*.tsv"))
    if not p:
        return {}, {}
    df = pd.read_csv(p[0], sep="\t", usecols=["compound_inchikey", "compound_name", "origin_type"])
    df["ik1"] = df.compound_inchikey.astype(str).str[:14]
    g = df.groupby("ik1")
    return g.origin_type.agg(lambda s: ";".join(sorted(set(map(str, s))))).to_dict(), g.compound_name.first().to_dict()


def _sdf_records(fh, tags):
    rec = {}
    tag = None
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("$$$$"):
            yield rec; rec = {}; tag = None; continue
        m = re.match(r">\s*<([^>]+)>", line)
        if m:
            tag = m.group(1) if m.group(1) in tags else None
            continue
        if tag and line and tag not in rec:
            rec[tag] = line


def lipidmaps(d):
    z = glob.glob(os.path.join(d, "lipidmaps", "*.zip"))
    if not z:
        return {}
    zf = zipfile.ZipFile(z[0])
    out = {}
    with zf.open(zf.namelist()[0]) as fb:
        fh = (l.decode("utf-8", "ignore") for l in fb)
        for r in _sdf_records(fh, {"INCHI_KEY", "CATEGORY", "MAIN_CLASS", "NAME"}):
            k = r.get("INCHI_KEY", "")[:14]
            if k and k not in out:
                out[k] = (r.get("CATEGORY", ""), r.get("MAIN_CLASS", ""), r.get("NAME", ""))
    return out


def chebi(d):
    roles = defaultdict(set)
    obo = glob.glob(os.path.join(d, "chebi", "*.obo.gz"))
    if obo:
        cur = None
        with gzip.open(obo[0], "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("id: CHEBI:"):
                    cur = line[4:].strip()
                elif line.startswith("relationship: has_role ") and cur:
                    r = line.split()[2]
                    if r in CHEBI_ROLES:
                        roles[cur].add(CHEBI_ROLES[r])
    out, names = {}, {}
    sdf = glob.glob(os.path.join(d, "chebi", "*.sdf.gz"))
    if sdf:
        with gzip.open(sdf[0], "rt", encoding="utf-8", errors="ignore") as f:
            for r in _sdf_records(f, {"ChEBI ID", "INCHIKEY", "ChEBI NAME"}):
                k = r.get("INCHIKEY", "")[:14]
                if not k:
                    continue
                cid = r.get("ChEBI ID", "")
                out.setdefault(k, set()).update(roles.get(cid, set()) or {"in ChEBI"})
                names.setdefault(k, r.get("ChEBI NAME", ""))
    return out, names


# ------------------------------------------------------------------------------------------------ main
def main(dbdir, outdir):
    os.makedirs(outdir, exist_ok=True)
    print("NORMAN SLE ...", flush=True)
    members, names, cand, pairs = norman_sle(dbdir)
    print(f"  {len(members)} structures in {len({c for s in members.values() for c in s})} lists; {len(pairs)} TP pairs", flush=True)
    print("COCONUT ...", flush=True); coco, coco_nm = coconut(dbdir); print(" ", len(coco))
    print("LOTUS ...", flush=True); lot = lotus(dbdir); print(" ", len(lot))
    print("NP Atlas ...", flush=True); npa, npa_nm = npatlas(dbdir); print(" ", len(npa))
    print("LIPID MAPS ...", flush=True); lmsd = lipidmaps(dbdir); print(" ", len(lmsd))
    print("ChEBI ...", flush=True); che, che_nm = chebi(dbdir); print(" ", len(che))
    keys = set(members) | coco | set(lot) | set(npa) | set(lmsd) | set(che)
    rows = []
    for k in keys:
        if len(k) != 14:
            continue
        lists = sorted(members.get(k, ()))
        src = [s for s, ok in (("COCONUT", k in coco), ("LOTUS", k in lot), ("NPAtlas", k in npa), ("LIPIDMAPS", k in lmsd),
                               ("ChEBI", k in che)) if ok]
        lm = lmsd.get(k, ("", "", ""))
        rows.append(dict(ik1=k, name=names.get(k) or lm[2] or che_nm.get(k) or npa_nm.get(k) or coco_nm.get(k, ""),
                         norman_lists=";".join(lists), n_norman_lists=len(lists), known_tp=bool(set(lists) & TP_LISTS),
                         np_sources=";".join(src), lotus_n_organisms=lot.get(k, 0), npatlas_origin=npa.get(k, ""),
                         chebi_roles=";".join(sorted(r for r in che.get(k, ()) if r != "in ChEBI")),
                         lipid_category=lm[0], lipid_class=lm[1]))
    idx = pd.DataFrame(rows).sort_values("ik1")
    idx.to_csv(os.path.join(outdir, "reference_index.csv.gz"), index=False, compression="gzip")
    pairs.to_csv(os.path.join(outdir, "tp_pairs.csv.gz"), index=False, compression="gzip")
    # SLE candidates with mass / formula
    Chem = _rdkit()
    from rdkit.Chem import Descriptors, rdMolDescriptors
    out = []
    for r in cand.itertuples():
        m = Chem.MolFromSmiles(r.smiles) if isinstance(r.smiles, str) else None
        if m is None:
            continue
        out.append(dict(name=r.name, formula=rdMolDescriptors.CalcMolFormula(m), exact_mass=Descriptors.ExactMolWt(m),
                        smiles=r.smiles, inchikey=r.inchikey, lists=";".join(sorted(members.get(r.inchikey[:14], ())))))
    pd.DataFrame(out).to_csv(os.path.join(outdir, "sle_candidates.csv.gz"), index=False, compression="gzip")
    print(f"reference index: {len(idx)} structures -> {outdir}")
    print(idx[["n_norman_lists"]].gt(0).sum().to_dict(), "known TPs:", int(idx.known_tp.sum()),
          "natural products:", int((idx.np_sources != "").sum()), "lipids:", int((idx.lipid_category != "").sum()))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
