"""Targeted screen of a compound list (internal standards, targets, calibrants) in any mzML (DDA or DIA/SWATH, centroid or profile).

Per compound and file:  Q1 (precursor m/z error, wide +-40 ppm search so a calibration offset cannot hide a compound) ->
candidate peaks -> MS2 (DIA: correlation pseudo-MS2; DDA: precursor-matched spectra) -> product ions annotated as
sub-formulae of the (isotope-labelled) formula (Q4 error) -> label-aware library match against the installed spectral library
(a library fragment may match with k deuteriums / 13C retained) -> Schymanski level -> +-10 ppm re-integration with bounds,
baseline and S/N -> blank check, calibration curves, quantification, RT anchors for the RT models.

Validated 2026-10 on SCIEX X500R SWATH data with 75 deuterated drug internal standards (see engine/VALIDATION.md)."""
from __future__ import annotations
import os, re, logging
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

from . import mzml_io
from .common import config_for_file, resolve, guess_polarity, read_table, PROTON, ELECTRON

log = logging.getLogger("hrms_nts")

# ---------------------------------------------------------------------------------------------------- masses / labels
MONO = {"C": 12.0, "13C": 13.0033548378, "H": 1.00782503207, "D": 2.01410177812, "N": 14.0030740048, "15N": 15.0001088982,
        "O": 15.99491461956, "S": 31.97207100, "P": 30.97376163, "F": 18.99840322, "Cl": 34.96885268, "Br": 78.9183371,
        "I": 126.904473, "Na": 22.9897692809, "K": 38.96370668, "Si": 27.9769265325}
AVG = {"C": 12.0106, "13C": 13.0033548, "H": 1.00794, "D": 2.0141018, "N": 14.0067, "15N": 15.0001089, "O": 15.9994, "S": 32.065,
       "P": 30.973762, "F": 18.9984032, "Cl": 35.453, "Br": 79.904, "I": 126.90447, "Na": 22.98977, "K": 39.0983, "Si": 28.0855}
DD = MONO["D"] - MONO["H"]; D13 = MONO["13C"] - MONO["C"]; D15 = MONO["15N"] - MONO["N"]
ORDER = ["C", "13C", "H", "D", "Br", "Cl", "F", "I", "N", "15N", "Na", "K", "O", "P", "S", "Si"]
ADDUCTS = {"positive": {"[M+H]+": PROTON, "[M+Na]+": MONO["Na"] - ELECTRON, "[M+NH4]+": MONO["N"] + 3 * MONO["H"] + PROTON,
                        "[M+K]+": MONO["K"] - ELECTRON},
           "negative": {"[M-H]-": -PROTON, "[M+Cl]-": MONO["Cl"] + ELECTRON, "[M+HCOO]-": MONO["C"] + MONO["H"] + 2 * MONO["O"] + ELECTRON}}
MAIN = {"positive": "[M+H]+", "negative": "[M-H]-"}
SALT_WORDS = r"\b(hcl|hydrochloride|hbr|oxalate|perchlorate|sodium salt|potassium salt|sulfate salt|tartrate|citrate|maleate|fumarate|acetate salt|solution)\b"


def parse_formula(f: str) -> dict:
    """'C4H[2H]4NO4S', 'C5[13C]3H10N4O2', 'C10H6D9N', 'C9H7(2)H6N' -> element counts with D / 13C / 15N as own elements."""
    if not isinstance(f, str):
        return {}
    f = f.replace(" ", "").replace("(2)H", "[2H]").replace("(13)C", "[13C]").replace("(15)N", "[15N]")
    d = {}
    for el, n in re.findall(r"(\[13C\]|\[2H\]|\[15N\]|Cl|Br|Na|Si|[A-Z][a-z]?)(\d*)", f):
        el = {"[13C]": "13C", "[2H]": "D", "[15N]": "15N"}.get(el, el)
        d[el] = d.get(el, 0) + (int(n) if n else 1)
    return {k: v for k, v in d.items() if v}


def fstr(d: dict) -> str:
    return "".join({"13C": "[13C]", "D": "[2H]", "15N": "[15N]"}.get(e, e) + (str(d[e]) if d[e] != 1 else "") for e in ORDER if d.get(e))


def mono(d): return float(sum(MONO[k] * v for k, v in d.items()))
def avg_mw(d): return float(sum(AVG[k] * v for k, v in d.items()))


def unlabel(d: dict) -> dict:
    u = dict(d)
    u["H"] = u.get("H", 0) + u.pop("D", 0); u["C"] = u.get("C", 0) + u.pop("13C", 0); u["N"] = u.get("N", 0) + u.pop("15N", 0)
    return {k: v for k, v in u.items() if v}


def label_from_text(text: str) -> tuple[int, int, int]:
    """nD, n13C, n15N from a label string or compound name: 'D4', '-d10', '13C3', '¹³C,D2', '13C6,15N2', 'd3 (MAM-D3)'."""
    t = str(text).replace("¹³", "13").replace("¹⁵", "15")
    nD = n13 = n15 = 0
    m = re.findall(r"(?:^|[-_ ,(])[dD](\d{1,2})\b", t)
    if m: nD = int(m[0])
    m = re.search(r"13C(\d*)", t)
    if m: n13 = int(m.group(1) or 1)
    m = re.search(r"15N(\d*)", t)
    if m: n15 = int(m.group(1) or 1)
    return nD, n13, n15


def parent_name(name: str) -> str:
    """'Fentanyl-D5 oxalate' -> 'Fentanyl';  '6-Acetylmorphine-D3 (MAM-D3)' -> '6-Acetylmorphine'."""
    n = re.sub(r"\(.*?[dD]\d+.*?\)", "", str(name))
    n = re.sub(r"[*†‡]", "", n)
    n = re.sub(SALT_WORDS, "", n, flags=re.I)
    n = re.sub(r"[-_ ]?(?:¹³C|13C)\d*(?:,\s*)?", "", n)
    n = re.sub(r"[-_ ]?(?:2,4,5,6-)?[dD]\d{1,2}\b", "", n)
    n = re.sub(r"[-_ ,]+$", "", n.strip())
    return re.sub(r"\s{2,}", " ", n).strip()


def _name_variants(n):
    n = n.replace("′", "'").replace("’", "'").strip()
    v = [n.lower()]
    s = re.sub(r"^(\(\s*[±+\-/]+\s*\)-|rac-|\(rs\)-)", "", n, flags=re.I)
    s2 = re.sub(r"^(cis|trans|threo|erythro|dl|d|l|r|s)-", "", s, flags=re.I)
    s3 = re.sub(r"^(cis|trans|threo|erythro)-", "", s2, flags=re.I)
    for x in (s, s2, s3, n.replace("-", " ")):
        if x.lower() not in v: v.append(x.lower())
    return v


# ---------------------------------------------------------------------------------------------------- compound list
_COLS = {"name": ("name", "compound", "compound name", "analyte"), "formula": ("formula", "chemical formula", "molecular formula"),
         "smiles": ("smiles", "parent_smiles", "parent smiles"), "label": ("label", "isotope label", "labelling"),
         "role": ("role", "type", "kind"), "q1": ("q1", "precursor", "precursor (q1) mass (da)", "precursor m/z", "mz", "m/z"),
         "rt": ("rt", "rt_expected", "expected rt", "retention time"), "polarity": ("polarity", "mode", "esi"),
         "is_name": ("is_name", "internal standard", "is"), "adduct": ("adduct", "adduct/charge")}


def load_compounds(path: str, db=None) -> pd.DataFrame:
    """CSV/XLSX compound list -> standard table. Accepts the vendor-style IS sheet (Name, Chemical formula, Adduct, Q1) or
    a list with SMILES (+ label / nD / n13C). Missing SMILES are looked up by parent name in the installed compound database."""
    if str(path).lower().endswith((".xlsx", ".xls")):
        raw = pd.read_excel(path)
    else:
        raw = read_table(path)
    low = {str(c).strip().lower(): c for c in raw.columns}
    pick = {k: next((low[a] for a in al if a in low), None) for k, al in _COLS.items()}
    if pick["name"] is None:
        raise ValueError(f"{path}: needs a name column (one of {_COLS['name']})")
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    from rdkit.Chem.inchi import MolToInchiKey
    rows = []
    dbn = None
    if db is not None:
        dbn = db.df.assign(_n=db.df["name"].astype(str).str.lower().str.strip())
    for _, r in raw.iterrows():
        name = r[pick["name"]]
        if not isinstance(name, str) or not name.strip():
            continue
        name = name.strip()
        get = lambda k: (r[pick[k]] if pick[k] is not None and pd.notna(r[pick[k]]) else None)
        smi = get("smiles"); fml = get("formula"); lab = get("label")
        nD = int(r[low["nd"]]) if "nd" in low and pd.notna(r[low["nd"]]) else None
        n13 = int(r[low["n13c"]]) if "n13c" in low and pd.notna(r[low["n13c"]]) else None
        d = parse_formula(fml) if fml else {}
        labelled_formula = any(k in d for k in ("D", "13C", "15N"))
        if not labelled_formula:
            lD, l13, l15 = label_from_text(lab if lab else name)
            nD = nD if nD is not None else lD; n13 = n13 if n13 is not None else l13
            if not d and smi:
                m = Chem.MolFromSmiles(str(smi))
                d = parse_formula(CalcMolFormula(m).rstrip("+-")) if m else {}
            if d and (nD or n13 or l15):
                d = dict(d)
                d["H"] -= nD or 0; d["C"] -= n13 or 0; d["N"] = d.get("N", 0) - l15
                if nD: d["D"] = nD
                if n13: d["13C"] = n13
                if l15: d["15N"] = l15
                d = {k: v for k, v in d.items() if v}
        pname = parent_name(name)
        pf = unlabel(d) if d else {}
        src = "list"
        if not smi and dbn is not None:
            for v in _name_variants(pname):
                hit = dbn[dbn["_n"] == v]
                if pf and len(hit):
                    hit = hit[hit["formula"].astype(str).map(lambda x: parse_formula(x) == pf)]
                if len(hit):
                    smi = hit.iloc[0]["smiles"]; src = "database (name)"
                    break
        if not d and smi:
            m = Chem.MolFromSmiles(str(smi))
            if m: d = parse_formula(CalcMolFormula(m).rstrip("+-")); pf = d
        ik1 = ""
        if smi:
            m = Chem.MolFromSmiles(str(smi))
            if m is not None:
                ik1 = MolToInchiKey(m).split("-")[0]
                if pf and parse_formula(CalcMolFormula(m).rstrip("+-")) != pf:
                    log.warning("%s: SMILES formula %s does not match the (unlabelled) list formula %s - SMILES ignored", name,
                                CalcMolFormula(m), fstr(pf))
                    smi, ik1, src = None, "", "list"
        if not d:
            log.warning("%s: no formula and no SMILES - skipped", name); continue
        role = str(get("role") or ("IS" if any(k in d for k in ("D", "13C", "15N")) else "target")).strip()
        rows.append(dict(name=name, parent=pname, formula=fstr(d), formula_dict=d, parent_formula=fstr(pf), nD=d.get("D", 0),
                         n13C=d.get("13C", 0), n15N=d.get("15N", 0), mono=mono(d), avg_mw=avg_mw(d), parent_mono=mono(pf),
                         smiles=smi or "", ik1=ik1, smiles_source=src if smi else "", role=role,
                         q1_list=float(get("q1")) if get("q1") is not None and str(get("q1")).replace(".", "", 1).isdigit() else np.nan,
                         rt_expected=float(get("rt")) if get("rt") is not None else np.nan,
                         polarity=str(get("polarity") or "").lower(), is_name=str(get("is_name") or ""),
                         adduct_list=str(get("adduct") or "")))
    df = pd.DataFrame(rows)
    # isobaric groups (identical labelled formula)
    df["isobaric_group"] = ""
    for k, (f, g) in enumerate(df.groupby("formula")):
        if len(g) > 1:
            df.loc[g.index, "isobaric_group"] = " / ".join(g["name"])
    log.info("compound list: %d compounds (%d labelled), %d SMILES from the database, %d isobaric groups", len(df),
             int((df.nD + df.n13C + df.n15N > 0).sum()), int((df.smiles_source == "database (name)").sum()),
             df.loc[df.isobaric_group != "", "isobaric_group"].nunique())
    return df


def ion_mz(c, adduct, pol):
    return c["mono"] + ADDUCTS[pol][adduct]


def likely_negative(smiles: str) -> bool:
    """acidic group and no basic nitrogen -> ionises mainly in ESI-."""
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        has = lambda s: m.HasSubstructMatch(Chem.MolFromSmarts(s))
        no_n = not any(a.GetSymbol() == "N" for a in m.GetAtoms())
        acid = (any(has(s) for s in ("C(=O)[OX2H1]", "S(=O)(=O)[OX2H1]", "OS(=O)(=O)[OX2H1,OX1-]", "C(=O)NS(=O)(=O)O"))
                or (no_n and has("c[OX2H1]"))                                                      # phenols without any N (triclosan)
                or len(m.GetSubstructMatches(Chem.MolFromSmarts("[NX3]S(=O)(=O)"))) >= 2)            # bis-sulfonamides (thiazides)
        basic = has("[NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(N-a);!$(N=*);!$(N#*)]")
        return acid and not basic
    except Exception:
        return False


# ---------------------------------------------------------------------------------------------------- peaks
def find_candidate_peaks(rt, I, MZ, mz0, min_apex, max_peaks=4):
    sm = gaussian_filter1d(I, 1.0)
    pk, _ = find_peaks(sm, height=min_apex, prominence=min_apex * 0.5, distance=5)
    out = []
    for p in pk:
        h = sm[p]
        a = p
        while a > 0 and sm[a - 1] > 0.05 * h and sm[a - 1] <= sm[a] * 1.05: a -= 1
        b = p
        while b < len(sm) - 1 and sm[b + 1] > 0.05 * h and sm[b + 1] <= sm[b] * 1.05: b += 1
        top = np.arange(a, b + 1)
        top = top[(I[top] >= 0.3 * I[a:b + 1].max()) & np.isfinite(MZ[top])]
        if len(top) == 0:
            continue
        w = I[top]; mzw = float(np.sum(MZ[top] * w) / w.sum())
        out.append(dict(rt=float(rt[p]), height=float(I[a:b + 1].max()), n_scans=int(b - a + 1), mz_meas=mzw,
                        ppm=(mzw - mz0) / mz0 * 1e6, ppm_sd=float(np.std((MZ[top] - mz0) / mz0 * 1e6)),
                        rt_start=float(rt[a]), rt_end=float(rt[b])))
    return sorted(out, key=lambda d: -d["height"])[:max_peaks]


def integrate(rt, I, rt0, floor, search=0.1):
    """+-search min apex search, bounds at baseline + 5 % of apex (or a valley), trapezoid area above baseline (counts x s)."""
    sm = gaussian_filter1d(I, 1.0)
    w = np.flatnonzero(np.abs(rt - rt0) <= search)
    if len(w) == 0 or sm[w].max() <= 0:
        return dict(apex_rt=rt0, start_rt=np.nan, end_rt=np.nan, height=0.0, area=0.0, baseline=0.0, noise=0.0, sn=0.0, fwhm_s=np.nan,
                    n_points=0, valley_limited=False)
    p = w[np.argmax(sm[w])]; h = sm[p]
    out = (np.abs(rt - rt[p]) <= 1.5) & (np.abs(rt - rt[p]) > 0.3)
    base = float(np.median(sm[out])) if out.any() else 0.0
    noise = float(np.median(np.abs(sm[out] - base)) * 1.4826) if out.any() else 0.0
    thr = base + 0.05 * (h - base)
    a = p
    while a > 0 and sm[a - 1] > thr and sm[a - 1] <= sm[a] * 1.02: a -= 1
    b = p
    while b < len(sm) - 1 and sm[b + 1] > thr and sm[b + 1] <= sm[b] * 1.02: b += 1
    seg = slice(a, b + 1)
    area = float(np.trapezoid(np.clip(I[seg] - base, 0, None), rt[seg]) * 60)
    half = np.flatnonzero(sm[seg] >= base + 0.5 * (h - base))
    return dict(apex_rt=float(rt[p]), start_rt=float(rt[a]), end_rt=float(rt[b]), height=float(I[p]), area=area, baseline=base, noise=noise,
                sn=float((h - base) / max(noise, floor)), fwhm_s=float((rt[a + half[-1]] - rt[a + half[0]]) * 60) if len(half) > 1 else np.nan,
                n_points=int(b - a + 1), valley_limited=bool(sm[a] > base + 0.2 * (h - base) or sm[b] > base + 0.2 * (h - base)))


def integration_rating(r):
    if r["sn"] >= 100 and r["n_points"] >= 8 and not r["valley_limited"]: return "Good"
    if r["sn"] >= 10 and r["n_points"] >= 5: return "Fair" + (" (valley-limited)" if r["valley_limited"] else "")
    return "Poor"


def iso_theory(d):
    c = d.get("C", 0); h = d.get("H", 0); n = d.get("N", 0); o = d.get("O", 0); s = d.get("S", 0); cl = d.get("Cl", 0); br = d.get("Br", 0)
    m1 = c * 0.010816 + h * 0.000115 + n * 0.003654 + o * 0.000381 + s * 0.00800
    m2 = (c * 0.010816) ** 2 / 2 + o * 0.002005 + s * 0.04517 + cl * 0.31978 + (cl * (cl - 1) / 2) * 0.31978 ** 2 + br * 0.9728
    return m1, m2


# ---------------------------------------------------------------------------------------------------- product ions
def subformula_table(d):
    els = [e for e in ["C", "13C", "H", "D", "N", "15N", "O", "S", "P", "Cl", "Br", "F", "I", "Si"] if d.get(e)]
    grids = np.array(np.meshgrid(*[np.arange(d[e] + 1) for e in els], indexing="ij")).reshape(len(els), -1).T
    masses = grids @ np.array([MONO[e] for e in els])
    g = lambda e: grids[:, els.index(e)] if e in els else 0
    nC = g("C") + g("13C"); nH = g("H") + g("D"); nN = g("N") + g("15N"); nX = g("Cl") + g("F") + g("I") + g("Br")
    rdbe = nC - (nH + nX) / 2 + nN / 2 + 1
    k = (nC >= 1) & (rdbe >= -0.5)
    o = np.argsort(masses[k])
    return els, grids[k][o], masses[k][o]


def annotate_ion(mz, sub, pol, ppm=15.0):
    """Best sub-formula for a product ion: even-electron [F+H]+ / [F-H]- preferred, radical F+. / F-. allowed. err = obs - theo."""
    els, g, m = sub
    best = None
    shifts = ((PROTON, "[F+H]+"), (-ELECTRON, "F+.")) if pol == "positive" else ((-PROTON, "[F-H]-"), (ELECTRON, "F-."))
    for sh, tag in shifts:
        neutral = mz - sh
        lo, hi = np.searchsorted(m, neutral * (1 - ppm * 1e-6)), np.searchsorted(m, neutral * (1 + ppm * 1e-6))
        for i in range(lo, hi):
            theo = m[i] + sh
            err = (mz - theo) / theo * 1e6
            even = tag.startswith("[")
            if best is None or abs(err) < abs(best[1]) - 1 or (even and not best[2].startswith("[") and abs(err) < abs(best[1]) + 2):
                best = ({e: int(v) for e, v in zip(els, g[i]) if v}, err, tag, theo)
    if best is None:
        return None
    f, err, tag, theo = best
    f = dict(f)
    if tag == "[F+H]+": f["H"] = f.get("H", 0) + 1
    if tag == "[F-H]-": f["H"] = f.get("H", 0) - 1
    return dict(formula=fstr({k: v for k, v in f.items() if v}) + ("+" if pol == "positive" else "-"), theo=float(theo), ppm=float(err),
                ion_type=tag, frag_nD=f.get("D", 0), frag_n13C=f.get("13C", 0))


def ms2_for_peak(run, mz0, pk, cfg):
    """(window/scan description, [dict(mz, inten, r)]) for a candidate peak: DIA correlation or DDA precursor match."""
    if run.acquisition == "DIA" and run.dia is not None:
        from .dia import pseudo_ms2
        dc = cfg.get("dia", {})
        rt, I, _ = run.xic_max(mz0, 10.0)
        w, fr = pseudo_ms2(run, mz0, pk["rt_start"], pk["rt_end"], rt, I, dc.get("corr_min", 0.7), dc.get("apex_tol_scans", 3),
                           dc.get("frag_ppm", 20.0), dc.get("max_candidates", 150), dc.get("min_frag_int", 30.0))
        return (f"DIA {w[1]:.1f}-{w[2]:.1f} ({w[3]:g} eV)" if w else "no DIA window"), fr
    from .ms2 import MS2Index
    idx = MS2Index(run)
    m2 = cfg["ms2"]
    sp = idx.best_spectrum(mz0, pk["rt"], pk["rt_start"], pk["rt_end"], m2["precursor_tol_da"], m2["rt_pad_min"], m2["merge_scans"],
                           m2["rel_min"], m2["top_n"])
    if sp is None:
        return "no DDA MS2", []
    fr = [dict(mz=float(a), inten=float(b), r=np.nan) for a, b in zip(sp["mz"], sp["inten"])]
    return f"DDA {sp['n_scans']} scan(s), CE {sp['ce']}", sorted(fr, key=lambda f: -f["inten"])


def library_spectra(lib, ik1, pmz_label, pol, parent_mz):
    """Parent spectra (InChIKey block 1) + labelled-standard spectra at the labelled precursor m/z."""
    if lib is None:
        return []
    con = lib.con
    rows = []
    if ik1:
        rows += [("parent",) + r for r in con.execute(
            "SELECT id,name,adduct,ce,instrument,source,precursor_mz,peaks FROM spectra WHERE ik1=? AND polarity=?", (ik1, pol))]
    tol = pmz_label * 10e-6
    for r in con.execute("SELECT id,name,adduct,ce,instrument,source,precursor_mz,peaks FROM spectra WHERE polarity=? AND precursor_mz BETWEEN ? AND ?",
                         (pol, pmz_label - tol, pmz_label + tol)):
        if re.search(r"(-d\d|d\d+\b|13c|15n|deuter|_d\d)", (r[1] or "").lower()):
            rows.append(("labelled",) + r)
    out = []
    for kind, i, nm, add, ce, ins, src, pmz, blob in rows:
        a = np.frombuffer(blob, dtype=np.float32).reshape(-1, 2).astype(float)
        out.append(dict(kind=kind, id=i, name=nm, adduct=add, ce=ce, instrument=ins, source=src, pmz=pmz or parent_mz, mz=a[:, 0], it=a[:, 1]))
    return out


def label_cosine(frags, lib, nD, n13, prec_obs, tol_da=0.01, tol_ppm=20.0, shift=True):
    """Square-root-intensity cosine, greedy one-to-one matching; with shift=True a library product ion of the UNLABELLED parent
    may match an observed ion at +k*dD + j*d13C (k <= nD, j <= n13C), i.e. the fragment kept k deuteriums / j 13C."""
    qm = np.array([f["mz"] for f in frags]); qi = np.sqrt(np.array([f["inten"] for f in frags]))
    k = qm < prec_obs - 1.5; qm, qi = qm[k], qi[k]
    lm, li = lib["mz"], np.sqrt(np.clip(lib["it"], 0, None))
    k = lm < (lib["pmz"] or prec_obs) - 1.5; lm, li = lm[k], li[k]
    if len(qm) == 0 or len(lm) == 0:
        return 0.0, 0, ""
    sh = [a * DD + b * D13 for a in range(nD + 1) for b in range(n13 + 1)] if shift else [0.0]
    pairs = []
    for a, m in enumerate(lm):
        for s in sh:
            t = m + s; tol = max(tol_da, t * tol_ppm * 1e-6)
            for b in np.flatnonzero(np.abs(qm - t) <= tol):
                pairs.append((li[a] * qi[b], a, b, s))
    pairs.sort(key=lambda p: -p[0])
    ua, ub, sc, used = set(), set(), 0.0, []
    for w, a, b, s in pairs:
        if a in ua or b in ub: continue
        ua.add(a); ub.add(b); sc += w; used.append(f"{lm[a]:.3f}->{qm[b]:.4f}(+{round(s / DD) if s else 0}D)")
    den = np.sqrt((qi ** 2).sum() * (li ** 2).sum())
    return (float(sc / den) if den else 0.0), len(used), "; ".join(used[:10])


# ---------------------------------------------------------------------------------------------------- per file screen
def screen_file(path, comp: pd.DataFrame, cfg_base: dict, lib=None, rt_models=None, standards=None, polarity="auto"):
    cfg = config_for_file(cfg_base, path)
    tc = cfg["targeted"]
    run = mzml_io.read_for_config(path, cfg)
    pol = polarity if polarity in ("positive", "negative") else (guess_polarity(path, cfg.get("polarity_patterns")) or run.polarity or "positive")
    run.polarity = pol
    floor = tc.get("noise_floor") or cfg.get("centroid", {}).get("min_height_ms1", 150.0)
    info = dict(file=run.name, path=path, instrument=cfg.get("_instrument", "generic"), acquisition=run.acquisition, polarity=pol,
                start_time=run.start_time, n_ms1=run.n_ms1, rt_end=float(run.ms1_rt[-1]) if run.n_ms1 else np.nan,
                profile_centroided=run.centroided_here, truncated=run.truncated,
                gaps="; ".join(f"{a:.2f}-{b:.2f} min" for a, b in run.gaps),
                dia_windows=len(run.dia["windows"]) if run.dia else 0)
    log.info("targeted: %s | %s | %s | %s | %d compounds", run.name, info["instrument"], run.acquisition, pol, len(comp))
    peaks, frags = [], []
    for _, c in comp.iterrows():
        if c.polarity in ("positive", "negative") and c.polarity != pol:
            continue
        d = c.formula_dict
        sub = None
        libs = None
        m1t, m2t = iso_theory(d)
        for add in ADDUCTS[pol]:
            mz0 = ion_mz(c, add, pol)
            rt, I, MZ = run.xic_max(mz0, tc["wide_ppm"])
            for rank, pk in enumerate(find_candidate_peaks(rt, I, MZ, mz0, tc["min_apex"], tc["max_peaks"])):
                row = dict(name=c["name"], file=run.name, adduct=add, mz_theo=mz0, peak_rank=rank + 1, **pk)
                if add == MAIN[pol]:
                    sel = (rt >= pk["rt_start"]) & (rt <= pk["rt_end"])
                    base = I[sel].sum() or 1.0
                    _, I1, _ = run.xic_max(mz0 + 1.003355, 15)
                    _, I2, _ = run.xic_max(mz0 + (1.997050 if d.get("Cl") or d.get("Br") else 2.006710), 15)
                    row.update(m1_obs=I1[sel].sum() / base, m1_theo=m1t, m2_obs=I2[sel].sum() / base, m2_theo=m2t)
                    if c.nD:
                        _, Im, _ = run.xic_max(mz0 - DD, 15); row["d_minus1"] = Im[sel].sum() / base
                    if c.nD or c.n13C or c.n15N:
                        _, In, _ = run.xic_max(c.parent_mono + ADDUCTS[pol][add], 15); row["unlabelled_same_rt"] = In[sel].sum() / base
                    if rank < 3:
                        if sub is None:
                            sub = subformula_table(d)
                            libs = library_spectra(lib, c.ik1, mz0, pol, c.parent_mono + ADDUCTS[pol][add])
                        desc, fr = ms2_for_peak(run, mz0, pk, cfg)
                        ann = []
                        for f in fr:
                            a = annotate_ion(f["mz"], sub, pol, tc["annotate_ppm"])
                            f = dict(f, **(a or {})); ann.append(f)
                            frags.append(dict(name=c["name"], file=run.name, rt=pk["rt"], peak_rank=rank + 1, frag_mz=f["mz"], inten=f["inten"],
                                              corr_r=f.get("r"), formula=f.get("formula"), theo=f.get("theo"), ppm=f.get("ppm"),
                                              ion_type=f.get("ion_type"), frag_nD=f.get("frag_nD"), frag_n13C=f.get("frag_n13C")))
                        expl = [x for x in ann if x.get("formula") and 60 <= x["mz"] < mz0 - 1.5]
                        q4 = [x for x in expl if abs(x["ppm"]) <= tc["q4_ppm"]]
                        best = (0.0, 0, "", None)
                        for L in libs:
                            s = label_cosine(ann, L, int(c.nD), int(c.n13C), mz0, cfg["library"]["frag_tol_da"], cfg["library"]["frag_tol_ppm"],
                                             shift=(L["kind"] == "parent"))
                            if (s[0], s[1]) > (best[0], best[1]):
                                best = (*s, L)
                        L = best[3]
                        row.update(ms2_source=desc, n_frags=len(ann), n_expl=len(expl), n_q4=len(q4),
                                   q4_top3_all=bool(expl) and all(abs(x["ppm"]) <= tc["q4_ppm"] for x in expl[:3]),
                                   q4_median_abs_ppm=float(np.median([abs(x["ppm"]) for x in q4[:5]])) if q4 else np.nan,
                                   q4_ions="; ".join(f"{x['theo']:.4f} {x['formula']} (meas {x['mz']:.4f}, {x['ppm']:+.1f} ppm)" for x in expl[:3]),
                                   top_frags="; ".join(f"{x['mz']:.4f}({x.get('formula', '?')})" for x in ann[:8]),
                                   n_lib=len(libs), lib_score=best[0], lib_matched=best[1], lib_matches=best[2],
                                   lib_kind=L["kind"] if L else "", lib_name=L["name"] if L else "", lib_source=L["source"] if L else "",
                                   lib_ce=L["ce"] if L else "")
                peaks.append(row)
    P = pd.DataFrame(peaks); Fr = pd.DataFrame(frags)
    return dict(info=info, peaks=P, frags=Fr, run=run, cfg=cfg)


# ---------------------------------------------------------------------------------------------------- assignment / levels
def assign_peaks(P: pd.DataFrame, comp: pd.DataFrame, cfg: dict, standards=None, pol="positive", batch_rt=None, batch_tol=0.25):
    """Best peak per compound per file. RT is reported; it only steers the choice when an expected RT exists (list column rt or
    standards.csv) or to split isobaric compounds: by their MS2 library match in this file, or (batch_rt) by the RT at which the
    member was library-identified in another file of the same batch. batch_rt never creates Level 1."""
    tc = cfg["targeted"]
    if P.empty:
        return P
    M = P[P.adduct.isin(MAIN.values()) & (P.ppm.abs() <= 10)].copy()
    if M.empty:
        return M
    M["batch_rt_claim"] = False
    for c in ("lib_score", "lib_matched", "n_expl", "n_q4"):
        M[c] = M[c].fillna(0) if c in M.columns else 0
    M["score"] =(np.where(M.lib_matched >= 2, M.lib_score, 0) + 0.15 * M.n_expl.clip(upper=8) + 0.1 * np.log10(M.height.clip(lower=1))
                  - 0.03 * M.ppm.abs())
    exp_rt = comp.set_index("name")["rt_expected"].to_dict()
    if standards is not None:
        col = "rt_pos" if pol == "positive" else "rt_neg"
        if col in standards.columns:
            srt = standards.assign(_n=standards["name"].astype(str).str.lower()).dropna(subset=[col]).set_index("_n")[col].to_dict()
            for n in comp["name"]:
                if not np.isfinite(exp_rt.get(n, np.nan)) and n.lower() in srt:
                    exp_rt[n] = srt[n.lower()]
    M["rt_expected"] = M.name.map(exp_rt)
    tol = cfg["rt"]["tolerance_min"]
    M["rt_diff_expected"] = M.rt - M.rt_expected
    M.loc[M.rt_diff_expected.abs() <= tol, "score"] += 1.0
    lib_ok = (M.lib_score >= tc["lib_cosine"]) & (M.lib_matched >= tc["lib_min_matched"])
    # isobaric groups: a peak whose library match identifies one member is removed from the others
    groups = comp.loc[comp.isobaric_group != "", ["name", "isobaric_group"]]
    for grp, g in groups.groupby("isobaric_group"):
        members = list(g["name"])
        for f in M.file.unique():
            sub = M[(M.file == f) & M.name.isin(members)]
            for rtv, pg in sub.groupby(sub.rt.round(2)):
                ok = pg[lib_ok.loc[pg.index]]
                rest = pg.lib_score.drop(ok.lib_score.idxmax()) if len(ok) else pg.lib_score
                if len(ok) and ok.lib_score.max() - (rest.max() if len(rest) else 0.0) >= 0.1:
                    winner = ok.loc[ok.lib_score.idxmax(), "name"]
                    M = M.drop(pg[pg.name != winner].index)
                    continue
                # batch RT transfer: exactly one member was library-identified at this RT in another file of the batch
                if batch_rt:
                    near = [n for n in pg.name.unique() if n in batch_rt and abs(float(rtv) - batch_rt[n]) <= batch_tol]
                    if len(near) == 1:
                        M = M.drop(pg[pg.name != near[0]].index)
                        M.loc[pg[pg.name == near[0]].index, "batch_rt_claim"] = True
            # a member identified by library MS2 in this file keeps only that peak; its other peaks are left to the other members
            sub = M[(M.file == f) & M.name.isin(members)]
            sub_ok = sub[lib_ok.reindex(sub.index).fillna(False).astype(bool) | sub.batch_rt_claim.astype(bool)]
            for nm in sub_ok.name.unique():
                M = M.drop(sub[(sub.name == nm) & ~sub.index.isin(sub_ok.index)].index)
    # height competition AFTER the isobaric split (a peak claimed by another isomer must not suppress this compound's own peak);
    # library-confirmed peaks are always eligible
    M["rel_h"] = M.height / M.groupby(["name", "file"]).height.transform("max")
    M = M[(M.rel_h >= tc["rel_height_compete"]) | lib_ok.reindex(M.index).fillna(False).astype(bool) | M.batch_rt_claim.astype(bool)]
    M.loc[M.batch_rt_claim.astype(bool), "score"] += 2.0
    ch = M.sort_values("score", ascending=False).groupby(["name", "file"]).head(1).reset_index(drop=True)
    return ch


def level_for(r, tc, level1_tol):
    q1ok = abs(r.ppm) <= tc["q1_ppm"]
    lib = r.lib_score >= tc["lib_cosine"] and r.lib_matched >= tc["lib_min_matched"]
    if not q1ok:
        return "Q1 > %g ppm" % tc["q1_ppm"]
    if lib and np.isfinite(r.rt_expected) and abs(r.rt - r.rt_expected) <= level1_tol:
        return "1"
    if lib:
        return "2a"
    if r.n_q4 >= 3:
        return "2b"
    if r.n_q4 >= 1 or (r.lib_matched >= 2 and r.lib_score >= 0.5):
        return "3"
    return "4"


# ---------------------------------------------------------------------------------------------------- calibration
def fit_calibration(x, y, weight="1/x"):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x, y = x[ok], y[ok]
    if len(x) < 2:
        return None
    w = {"1/x": 1 / x, "1/x2": 1 / x ** 2}.get(weight, np.ones_like(x))
    W = np.sqrt(w)
    A = np.c_[x, np.ones_like(x)] * W[:, None]
    (k, b), *_ = np.linalg.lstsq(A, y * W, rcond=None)
    pred = k * x + b
    ss = ((y - pred) ** 2).sum(); r2 = 1 - ss / ((y - y.mean()) ** 2).sum() if len(x) > 2 else np.nan
    sres = np.sqrt(ss / max(len(x) - 2, 1))
    back = (y - b) / k
    return dict(slope=k, intercept=b, r2=r2, n=len(x), lod=3.3 * sres / k if k else np.nan, loq=10 * sres / k if k else np.nan,
                accuracy_pct=";".join(f"{v:.0f}" for v in back / x * 100), range=(x.min(), x.max()))


# ---------------------------------------------------------------------------------------------------- orchestration
def run_targeted(files, list_path, out, cfg, sheet=None, blanks=None, polarity="auto", use_library=True, use_db=True,
                 update_standards=False, plots=True, workers=1):
    from .provenance import write_provenance
    os.makedirs(out, exist_ok=True)
    db = lib = None; rt_models = {}; standards = None
    if use_db:
        try:
            from .pipeline import Context
            ctx = Context(cfg, need_library=use_library)
            db, lib, standards, rt_models = ctx.db, ctx.lib, ctx.standards, ctx.rt_models
        except Exception as e:
            log.warning("compound database / RT models unavailable (%s) - continuing with the list only", e)
    if lib is None and use_library:
        lp = resolve(cfg, cfg["paths"].get("library_db"))
        if lp and os.path.exists(lp):
            from .library import SpectralLibrary
            lib = SpectralLibrary(lp)
    comp = load_compounds(list_path, db)
    roles = _file_roles(files, sheet, blanks)
    res = [screen_file(f, comp, cfg, lib, rt_models, standards, polarity) for f in files]
    info = pd.DataFrame([r["info"] for r in res])
    P = pd.concat([r["peaks"] for r in res if not r["peaks"].empty], ignore_index=True) if any(not r["peaks"].empty for r in res) else pd.DataFrame()
    Fr = pd.concat([r["frags"] for r in res if not r["frags"].empty], ignore_index=True) if any(not r["frags"].empty for r in res) else pd.DataFrame()
    info["role"] = info.file.map(lambda f: roles.get(f, {}).get("type", "sample"))
    info["conc"] = info.file.map(lambda f: roles.get(f, {}).get("conc", np.nan))
    tc = res[0]["cfg"]["targeted"] if res else cfg["targeted"]
    cfgF = res[0]["cfg"] if res else cfg
    blank_files = set(info.loc[info.role == "blank", "file"])
    # pass 1: isobaric members library-identified anywhere in the batch -> their RT (median) steers the other files (pass 2)
    batch_rt = {}
    iso_names = set(comp.loc[comp.isobaric_group != "", "name"])
    if iso_names and not P.empty:
        for r in res:
            a1 = assign_peaks(P[P.file == r["info"]["file"]], comp, r["cfg"], standards, r["info"]["polarity"])
            if len(a1):
                okl = a1[a1.name.isin(iso_names) & (a1.lib_score >= tc["lib_cosine"]) & (a1.lib_matched >= tc["lib_min_matched"])]
                for n, g in okl.groupby("name"):
                    batch_rt.setdefault(n, []).extend(g.rt.tolist())
        batch_rt = {n: float(np.median(v)) for n, v in batch_rt.items()}
        if batch_rt:
            log.info("isobaric compounds library-identified in the batch (RT used to split them in the other files): %s",
                     ", ".join(f"{n} {v:.2f} min" for n, v in batch_rt.items()))
    ch = []
    for r in res:
        pol = r["info"]["polarity"]
        sub = P[P.file == r["info"]["file"]] if not P.empty else P
        a = assign_peaks(sub, comp, r["cfg"], standards, pol, batch_rt)
        if len(a):
            a["level"] = a.apply(lambda x: level_for(x, r["cfg"]["targeted"], r["cfg"]["levels"]["level1_rt_tol"]), axis=1)
            iso = comp.set_index("name")["isobaric_group"]
            amb = a.name.map(iso).fillna("") != ""
            # isobaric members still sharing one peak in this file -> unresolved
            dup = a[amb].groupby(a.loc[amb, "rt"].round(2)).name.transform("size") > 1 if amb.any() else pd.Series(dtype=bool)
            a.loc[dup[dup].index, "level"] = a.loc[dup[dup].index, "level"] + " (isobaric, unresolved)"
            # integration
            run = r["run"]; fl = r["cfg"]["targeted"].get("noise_floor") or r["cfg"].get("centroid", {}).get("min_height_ms1", 150.0)
            integ = []
            for x in a.itertuples():
                rt, I, _ = run.xic_max(x.mz_theo, r["cfg"]["targeted"]["integrate_ppm"])
                g = integrate(rt, I, x.rt, fl); g["rating"] = integration_rating(g); integ.append(g)
            a = pd.concat([a.reset_index(drop=True), pd.DataFrame(integ).add_prefix("int_")], axis=1)
            ch.append(a)
    ch = pd.concat(ch, ignore_index=True) if ch else pd.DataFrame()
    # RT model plausibility for the unlabelled parent (QSAR/RTI models trained on standards.csv)
    if not ch.empty:
        ch = _rt_plausibility(ch, comp, rt_models, db)
    blank = _blank_check(ch, res, comp, blank_files, tc) if blank_files and not ch.empty else pd.DataFrame()
    cal, quant, cal_export = _calibration(ch, comp, info, tc) if (info.role == "cal").any() and not ch.empty else (pd.DataFrame(),) * 3
    natives = _natives(ch, res, comp, lib) if not ch.empty else pd.DataFrame()
    anchors = _rt_anchors(ch, comp, info)
    if update_standards and not anchors.empty:
        _merge_standards(anchors, resolve(cfg, cfg["paths"]["standards"]))
    from .targeted_report import write_report
    xlsx = write_report(out, comp, info, P, Fr, ch, blank, cal, quant, natives, anchors, tc, cfgF)
    if cal_export is not None and len(cal_export):
        cal_export.to_csv(os.path.join(out, "calibrants_for_semiquant.csv"), index=False)
    if not anchors.empty:
        anchors.to_csv(os.path.join(out, "rt_anchors_for_standards.csv"), index=False)
    if plots and not ch.empty:
        from .targeted_report import eic_pdf
        eic_pdf(os.path.join(out, "EIC_integration_QC.pdf"), ch, res, comp, blank_files, tc)
    write_provenance(out, cfg, files, None, extra=dict(mode="targeted", list=list_path, sample_sheet=sheet, blanks=sorted(blank_files)))
    return dict(xlsx=xlsx, assigned=ch, info=info, blank=blank, calibration=cal, anchors=anchors)


def _file_roles(files, sheet, blanks):
    roles = {}
    if sheet:
        s = read_table(sheet) if not str(sheet).lower().endswith(("xlsx", "xls")) else pd.read_excel(sheet)
        s.columns = [str(c).strip().lower() for c in s.columns]
        for _, r in s.iterrows():
            stem = os.path.splitext(os.path.basename(str(r["file"])))[0]
            roles[stem] = dict(type=str(r.get("type", "sample")).strip().lower(), conc=float(r["conc"]) if "conc" in s.columns and pd.notna(r.get("conc")) else np.nan)
    import fnmatch
    for f in files:
        stem = os.path.splitext(os.path.basename(f))[0]
        if blanks and any(fnmatch.fnmatch(stem.lower(), b.lower()) or fnmatch.fnmatch(os.path.basename(f).lower(), b.lower()) for b in blanks):
            roles.setdefault(stem, {})["type"] = "blank"
    for k, v in roles.items():
        t = v.get("type", "sample")
        v["type"] = {"calibration": "cal", "calibrant": "cal", "standard": "cal", "std": "cal", "solvent blank": "blank", "procedural blank": "blank"}.get(t, t)
    return roles


def _rt_plausibility(ch, comp, rt_models, db):
    ch = ch.copy()
    for c in ("rt_pred_qsar", "rt_pred_rti", "rt_model_agree", "rt_model_conflict"):
        ch[c] = np.nan
    if not rt_models:
        return ch
    cm = comp.set_index("name")
    for pol, mdl in rt_models.items():
        if mdl is None:
            continue
        sel = ch.adduct.map(lambda a: ("positive" if a.endswith("+") else "negative") == pol)
        if not sel.any():
            continue
        h = pd.DataFrame({"smiles": ch.loc[sel, "name"].map(cm["smiles"]).values, "rt_apex": ch.loc[sel, "rt"].values})
        if db is not None:
            col = "pred_rti_pos" if pol == "positive" else "pred_rti_neg"
            rti = db.df.dropna(subset=[col]).drop_duplicates("ik1").set_index("ik1")[col]
            h[col] = ch.loc[sel, "name"].map(cm["ik1"]).map(rti).values
        a = mdl.annotate(h)
        ch.loc[sel, "rt_pred_qsar"] = a["rt_pred_qsar"].values; ch.loc[sel, "rt_pred_rti"] = a["rt_pred_rti"].values
        ch.loc[sel, "rt_model_agree"] = a["rt_agree"].values; ch.loc[sel, "rt_model_conflict"] = a["rt_conflict"].values
    return ch


def _blank_check(ch, res, comp, blank_files, tc):
    runs = {r["info"]["file"]: r["run"] for r in res}
    rows = []
    samp = ch[~ch.file.isin(blank_files)]
    for n, g in samp.groupby("name"):
        ref = g.sort_values("int_area", ascending=False).iloc[0]
        for bf in blank_files:
            run = runs[bf]
            rt, I, MZ = run.xic_max(ref.mz_theo, tc["integrate_ppm"])
            lo, hi = ref.int_start_rt - 0.1, ref.int_end_rt + 0.1
            m = (rt >= lo) & (rt <= hi)
            area = float(np.trapezoid(I[m], rt[m]) * 60) if m.any() else 0.0
            k = np.flatnonzero(m)[np.argmax(I[m])] if m.any() and I[m].max() > 0 else None
            pts = int((I[m] > 0).sum())
            pct = area / g.int_area.mean() * 100 if g.int_area.mean() else np.nan
            if area <= 0: v = "Clean - no signal"
            elif pts < 3: v = f"Clean - single-scan noise spike ({int(I[k])} counts)"
            elif pct < tc["blank_trace_pct"]: v = f"Minimal trace (< {tc['blank_trace_pct']:g} %)"
            elif pct < tc["blank_low_pct"]: v = f"Low trace ({tc['blank_trace_pct']:g}-{tc['blank_low_pct']:g} %)"
            else: v = f"Significant (> {tc['blank_low_pct']:g} %) - carry-over or contamination"
            rows.append(dict(name=n, blank_file=bf, rt_window=f"{lo:.2f}-{hi:.2f} (from {ref.file})", mean_sample_area=g.int_area.mean(),
                             min_sample_area=g.int_area.min(), blank_area=area, blank_height=float(I[k]) if k is not None else 0.0,
                             blank_points=pts, blank_apex_rt=float(rt[k]) if k is not None else np.nan,
                             blank_ppm=float((MZ[k] - ref.mz_theo) / ref.mz_theo * 1e6) if k is not None and np.isfinite(MZ[k]) else np.nan,
                             blank_pct=pct, verdict=v))
    return pd.DataFrame(rows)


def _calibration(ch, comp, info, tc):
    cal_files = info[info.role == "cal"].set_index("file")["conc"]
    cm = comp.set_index("name")
    is_area = ch.set_index(["name", "file"])["int_area"]
    rows, qrows, export = [], [], []
    for n, g in ch.groupby("name"):
        isn = cm.at[n, "is_name"] if n in cm.index else ""
        def resp(f, a):
            if isn and (isn, f) in is_area.index and is_area[(isn, f)] > 0:
                return a / is_area[(isn, f)]
            return np.nan if isn else a
        gc = g[g.file.isin(cal_files.index)]
        x = gc.file.map(cal_files).values; y = np.array([resp(f, a) for f, a in zip(gc.file, gc.int_area)])
        usable = np.isfinite(x.astype(float)) & np.isfinite(y) & (x.astype(float) > 0)
        if len(np.unique(x[usable])) >= tc["min_cal_levels"]:
            fit = fit_calibration(x, y, tc["calibration_weight"])
            if fit:
                warn = "; ".join(w for w in (("slope <= 0 - response does not rise with concentration" if fit["slope"] <= 0 else ""),
                                             (f"R2 {fit['r2']:.3f} < 0.99" if np.isfinite(fit["r2"]) and fit["r2"] < 0.99 else "")) if w)
                if fit["slope"] <= 0:
                    fit["lod"] = fit["loq"] = np.nan
                rows.append(dict(name=n, response="area/IS area (" + isn + ")" if isn else "area", weighting=tc["calibration_weight"], **{k: v for k, v in fit.items() if k != "range"},
                                 range_low=fit["range"][0], range_high=fit["range"][1], warning=warn))
                if fit["slope"] <= 0:
                    continue
                for x_ in g[~g.file.isin(cal_files.index)].itertuples():
                    yv = resp(x_.file, x_.int_area); conc = (yv - fit["intercept"]) / fit["slope"] if fit["slope"] else np.nan
                    flag = "below LOQ" if conc < fit["loq"] else ("above range" if conc > fit["range"][1] else ("below range" if conc < fit["range"][0] else ""))
                    qrows.append(dict(name=n, file=x_.file, response=yv, conc=conc, flag=flag, level=x_.level))
        for x_ in gc.itertuples():
            export.append(dict(name=cm.at[n, "parent"] if n in cm.index else n, smiles=cm.at[n, "smiles"] if n in cm.index else "",
                               rt=round(x_.int_apex_rt, 3), conc_ug_ml=cal_files[x_.file], area=x_.int_area,
                               polarity="positive" if x_.adduct.endswith("+") else "negative", file=x_.file))
    return pd.DataFrame(rows), pd.DataFrame(qrows), pd.DataFrame(export)


def _natives(ch, res, comp, lib):
    """Unlabelled analogue of each labelled IS at the IS RT (co-elution with its labelled standard + exact mass + library MS2)."""
    runs = {r["info"]["file"]: (r["run"], r["cfg"]) for r in res}
    cm = comp.set_index("name")
    rows = []
    for x in ch.itertuples():
        c = cm.loc[x.name]
        if not (c.nD or c.n13C or c.n15N):
            continue
        run, cfg = runs[x.file]
        pol = "positive" if x.adduct.endswith("+") else "negative"
        mz0 = c.parent_mono + ADDUCTS[pol][x.adduct]
        rt, I, MZ = run.xic_max(mz0, cfg["targeted"]["integrate_ppm"])
        g = integrate(rt, I, x.int_apex_rt, cfg["targeted"].get("noise_floor") or 150.0, search=0.08)
        if g["height"] < cfg["targeted"]["min_apex"] / 2 or g["n_points"] < 3:
            continue
        k = np.argmin(np.abs(rt - g["apex_rt"]))
        ppm = (MZ[k] - mz0) / mz0 * 1e6 if np.isfinite(MZ[k]) else np.nan
        pk = dict(rt=g["apex_rt"], rt_start=g["start_rt"], rt_end=g["end_rt"])
        desc, fr = ms2_for_peak(run, mz0, pk, cfg)
        best = (0.0, 0)
        if lib is not None and c.ik1:
            for L in library_spectra(lib, c.ik1, mz0, pol, mz0):
                if L["kind"] == "parent":
                    s = label_cosine(fr, L, 0, 0, mz0, shift=False)
                    best = max(best, s[:2])
        rows.append(dict(IS=x.name, native=c.parent, file=x.file, native_mz=mz0, IS_rt=x.int_apex_rt, native_rt=g["apex_rt"],
                         drt=g["apex_rt"] - x.int_apex_rt, ppm=ppm, native_area=g["area"], IS_area=x.int_area,
                         native_to_IS=g["area"] / x.int_area if x.int_area else np.nan, lib_score=best[0], lib_matched=best[1],
                         evidence=("co-elutes with labelled IS + library MS2" if best[0] >= 0.7 and best[1] >= 3 else "co-elutes with labelled IS, exact mass")))
    return pd.DataFrame(rows)


def _rt_anchors(ch, comp, info):
    """RT of each labelled IS / target detected at Level 1-2b -> standards.csv rows for the RT models (deuterated IS co-elute
    with their native analogue within a few seconds on RP-LC)."""
    if ch.empty:
        return pd.DataFrame()
    ok = ch[ch.level.str.match(r"^(1|2a|2b)$") & ~ch.file.isin(info.loc[info.role == "blank", "file"])]
    cm = comp.set_index("name")
    rows = []
    for (n, ad), g in ok.groupby(["name", "adduct"]):
        c = cm.loc[n]
        if not c.smiles:
            continue
        from rdkit import Chem
        from rdkit.Chem.inchi import MolToInchiKey
        m = Chem.MolFromSmiles(c.smiles)
        ik = MolToInchiKey(m) if m else ""
        pol = "positive" if ad.endswith("+") else "negative"
        rows.append(dict(name=c.parent, smiles=c.smiles, inchikey=ik, rt_pos=round(g.int_apex_rt.median(), 3) if pol == "positive" else np.nan,
                         rt_neg=round(g.int_apex_rt.median(), 3) if pol == "negative" else np.nan, from_compound=n, n_files=len(g),
                         rt_spread_min=round(g.int_apex_rt.max() - g.int_apex_rt.min(), 3)))
    a = pd.DataFrame(rows)
    if a.empty:
        return a
    return a.groupby(["name", "smiles", "inchikey"], as_index=False).agg(rt_pos=("rt_pos", "max"), rt_neg=("rt_neg", "max"),
                                                                        from_compound=("from_compound", "first"), n_files=("n_files", "max"),
                                                                        rt_spread_min=("rt_spread_min", "max"))


def _merge_standards(anchors, path):
    """Append IS-derived RTs to user_data/standards.csv (existing rows / RTs are never overwritten)."""
    cols = ["name", "smiles", "inchikey", "rt_pos", "rt_neg"]
    if path and os.path.exists(path):
        with open(path, "r", encoding="utf-8-sig") as f:
            comments = [ln for ln in f if ln.lstrip().startswith("#")]
        old = read_table(path)
    else:
        comments, old = [], pd.DataFrame(columns=cols)
    have = set(old["inchikey"].astype(str).str[:14]) | set(old["name"].astype(str).str.lower())
    new = anchors[~(anchors.inchikey.str[:14].isin(have) | anchors.name.str.lower().isin(have))][cols]
    out = pd.concat([old, new], ignore_index=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(",".join(cols) + "\n")
        for c in comments: f.write(c)
        out[cols].to_csv(f, index=False, header=False)
    log.info("standards.csv: %d RT anchors added from the targeted screen (%d already present) -> %s", len(new), len(anchors) - len(new), path)
