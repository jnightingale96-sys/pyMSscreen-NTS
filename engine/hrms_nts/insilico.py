"""In-silico MS/MS: (a) built-in rule + bond-cleavage product-ion predictor (pure RDKit, no external tools),
(b) optional MetFrag command-line wrapper (Java + MetFragCommandLine jar), (c) candidate scoring against experimental MS2.

The built-in predictor is a transparent, screening-grade approach (neutral-loss rules + 1-2 bond cleavages with
decoy-corrected matching).  It is NOT a substitute for CFM-ID/MetFrag-quality fragmentation and is reported as such;
use engine 'both' to add MetFrag scores for the final Level 2b assignment."""
from __future__ import annotations
import os, re, subprocess, shutil, logging, itertools
from functools import lru_cache
import numpy as np
import pandas as pd
from .common import ELEMENT_MASS, NEUTRAL_LOSSES, PROTON, normalize_adduct

log = logging.getLogger("hrms_nts")
MH = ELEMENT_MASS["H"]
CARRIER = {"[M+H]+": 1.007276, "[M+Na]+": 22.989218, "[M+K]+": 38.963158, "[M+NH4]+": 18.033823}
NEG_ADDUCT_SHIFT = {"[M-H]-": -1.007276, "[M+Cl]-": 34.969402, "[M+HCOO]-": 44.998201, "[M+CH3COO]-": 59.013851}

_NL_SMARTS = {
    "H2O": "[OX2H]", "NH3": "[NX3;H2,H1;!$(NC=O)]", "CO": "[CX3]=O", "CO2": "C(=O)[OX2H1,OX1-]", "HCOOH": "C(=O)[OX2H]",
    "CH3OH": "[OX2][CH3]", "CH2O": "[OX2][CH3]", "SO2": "S(=O)=O", "SO3": "S(=O)(=O)[OX2,OX1-]", "HF": "F", "HCl": "Cl",
    "H2S": "[SX2H]", "HCN": "C#N", "NO2": "[N+](=O)[O-]", "NO": "[N+](=O)[O-]", "C2H2O": "[CH3]C(=O)", "C2H4": "[CH2][CH3]",
    "C3H6": "[CH3][CH]([CH3])", "CH3": "[CH3]", "N2": "N=N",
}


def _prep(smiles):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return None
    m = max(Chem.GetMolFrags(m, asMols=True), key=lambda x: x.GetNumHeavyAtoms())
    return m


def _atom_data(m):
    n = m.GetNumAtoms()
    sym = [a.GetSymbol() for a in m.GetAtoms()]
    base = np.array([ELEMENT_MASS.get(s, 0.0) for s in sym])
    hcount = np.array([a.GetTotalNumHs() for a in m.GetAtoms()], float)
    adj = [[] for _ in range(n)]
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        adj[i].append((j, b.GetIdx())); adj[j].append((i, b.GetIdx()))
    return sym, base, hcount, adj


def _components(adj, removed):
    n = len(adj); seen = [-1] * n; comps = []
    for s in range(n):
        if seen[s] >= 0:
            continue
        stack = [s]; seen[s] = len(comps); cur = [s]
        while stack:
            u = stack.pop()
            for v, bi in adj[u]:
                if bi in removed or seen[v] >= 0:
                    continue
                seen[v] = len(comps); stack.append(v); cur.append(v)
        comps.append(cur)
    return comps


def fragment_masses(smiles, max_depth=2, max_pairs=250):
    """Closed-shell (R-H) neutral fragment masses from cleaving 1 (or 2) bonds.  Returns dict mass -> (n_bonds, label)."""
    m = _prep(smiles)
    if m is None:
        return {}, None
    sym, base, hcount, adj = _atom_data(m)
    full = float(base.sum() + hcount.sum() * MH)
    cl = [b.GetIdx() for b in m.GetBonds() if str(b.GetBondType()) in ("SINGLE",) and not (b.GetIsAromatic())]
    out = {}
    def add(comp, removed_bonds, nb):
        mass = base[comp].sum() + (hcount[comp].sum()) * MH
        # one H per cleaved bond endpoint inside this component (closed-shell R-H)
        extra = 0
        for bi in removed_bonds:
            b = m.GetBondWithIdx(bi)
            if b.GetBeginAtomIdx() in comp: extra += 1
            if b.GetEndAtomIdx() in comp: extra += 1
        mass += extra * MH
        key = round(float(mass), 4)
        if key not in out or out[key][0] > nb:
            out[key] = (nb, f"cleave{nb}")
    combos = [(b,) for b in cl]
    if max_depth >= 2:
        pairs = list(itertools.combinations(cl, 2))
        if len(pairs) > max_pairs:
            step = len(pairs) / max_pairs
            pairs = [pairs[int(i * step)] for i in range(max_pairs)]
        combos += pairs
    # ring bonds: cleave pairs inside the same ring system (ring opening + fragmentation)
    ring_b = [b.GetIdx() for b in m.GetBonds() if b.IsInRing() and str(b.GetBondType()) == "SINGLE"]
    if max_depth >= 2 and 1 < len(ring_b) <= 30:
        combos += [p for p in itertools.combinations(ring_b, 2) if p not in set(combos)]
    for rem in combos:
        rs = set(rem)
        for comp in _components(adj, rs):
            if len(comp) < len(sym):
                add(comp, rs, len(rem))
    return out, full


def neutral_loss_masses(smiles, max_losses=2):
    """Allowed (sequential) neutral losses -> dict(mass -> label), only for functional groups present."""
    from rdkit import Chem
    m = _prep(smiles)
    if m is None:
        return {}
    ok = [k for k, s in _NL_SMARTS.items() if m.HasSubstructMatch(Chem.MolFromSmarts(s))]
    out = {}
    for r in range(1, max_losses + 1):
        for combo in itertools.combinations_with_replacement(ok, r):
            if len(set(combo)) < len(combo) and combo[0] not in ("H2O", "CO", "CH3"):
                continue
            mass = sum(NEUTRAL_LOSSES[c] for c in combo)
            out["+".join(combo)] = mass
    return out


def predict_product_ions(smiles, adduct="[M+H]+", max_depth=2, min_mz=40.0, top=40):
    return _predict_cached(smiles, normalize_adduct(adduct), max_depth, min_mz, top).copy()


@lru_cache(maxsize=60000)
def _predict_cached(smiles, adduct, max_depth, min_mz, top):
    """Rank-ordered predicted product ions for a candidate.  Returns DataFrame(mz, label, kind, rank_score)."""
    adduct = normalize_adduct(adduct)
    m = _prep(smiles)
    if m is None:
        return pd.DataFrame(columns=["mz", "label", "kind", "rank_score"])
    if any(a.GetFormalCharge() != 0 for a in m.GetAtoms()):
        log.debug("charged species not supported by built-in fragmenter: %s", smiles)
    frags, full = fragment_masses(smiles, max_depth)
    pos = adduct.endswith("+")
    if adduct in CARRIER:
        prec = full + CARRIER[adduct]
    elif adduct in NEG_ADDUCT_SHIFT:
        prec = full + NEG_ADDUCT_SHIFT[adduct]
    else:
        return pd.DataFrame(columns=["mz", "label", "kind", "rank_score"])
    rows = []
    for f0, (nb, lab) in frags.items():
        if pos:
            c = CARRIER.get(adduct, PROTON)
            rows.append((f0 + c, f"{lab} [F{adduct[2:-2]}]+", "bond", 1.0 / nb))
            if adduct in ("[M+H]+", "[M+NH4]+"):
                rows.append((f0 - MH - 0.000549, f"{lab} [F-H]+", "bond", 0.8 / nb))
        else:
            rows.append((f0 - PROTON, f"{lab} [F-H]-", "bond", 1.0 / nb))
    for lab, ms in neutral_loss_masses(smiles).items():
        n = lab.count("+") + 1
        radical = any(x in ("CH3", "NO") for x in lab.split("+"))      # odd-electron losses are rare from even-electron ESI ions
        rows.append((prec - ms, "-" + lab.replace("+", "-"), "neutral_loss", (0.5 if radical else 1.6) / n))
        if adduct == "[M+NH4]+":
            rows.append((prec - 17.026549 - ms, f"-NH3-{lab}", "neutral_loss", 1.2 / n))
    if adduct == "[M+NH4]+":
        rows.append((prec - 17.026549, "-NH3", "neutral_loss", 2.0))
    df = pd.DataFrame(rows, columns=["mz", "label", "kind", "rank_score"])
    df = df[(df["mz"] >= min_mz) & (df["mz"] < prec - 0.5)]
    prec_ = prec
    df["rank_score"] = df["rank_score"] * (0.4 + 0.6 * df["mz"] / prec_)   # favour larger, more diagnostic ions
    df["_r"] = df["mz"].round(3)
    df = df.sort_values("rank_score", ascending=False).drop_duplicates("_r").drop(columns="_r").head(top * 3)
    df["precursor_mz"] = prec
    return df.sort_values("rank_score", ascending=False).reset_index(drop=True).head(top)


def suggest_transitions(smiles, adduct="[M+H]+", n=3):
    """Top-n predicted precursor->product transitions for MRM/PRM method building (prediction only - confirm with a standard)."""
    df = predict_product_ions(smiles, adduct)
    if df.empty:
        return df
    df = df[df["mz"] > 50].copy()
    return df.head(n)[["precursor_mz", "mz", "label", "kind"]]


def _explained(ion_masses, peaks_mz, tol_da, tol_ppm, shift=0.0):
    tol = np.maximum(tol_da, peaks_mz * tol_ppm * 1e-6)
    idx = np.searchsorted(ion_masses, peaks_mz + shift)
    hit = np.zeros(len(peaks_mz), bool)
    for off in (-1, 0):
        j = np.clip(idx + off, 0, len(ion_masses) - 1)
        hit |= np.abs(ion_masses[j] - (peaks_mz + shift)) <= tol
    return hit


def score_candidate(smiles, adduct, peaks_mz, peaks_int, precursor_mz, tol_da=0.005, tol_ppm=10.0, rel_min=0.02, max_depth=2,
                    min_diag_mz=60.0):
    """Decoy-corrected explained-intensity score of an experimental MS2 against the built-in predicted ions."""
    pred = predict_product_ions(smiles, adduct, max_depth=max_depth, top=2000)
    peaks_mz = np.asarray(peaks_mz, float); peaks_int = np.asarray(peaks_int, float)
    k = (peaks_mz < precursor_mz - 0.5) & (peaks_int >= rel_min * (peaks_int.max() if len(peaks_int) else 1))
    if pred.empty or k.sum() == 0:
        return dict(bi_expl_n=0, bi_expl_adj_n=0.0, bi_expl_frac=0.0, bi_decoy_frac=0.0, bi_score=0.0, bi_n_considered=int(k.sum()), bi_expl_peaks="")
    ions = np.sort(pred["mz"].values)
    mz, it = peaks_mz[k], peaks_int[k]
    hit = _explained(ions, mz, tol_da, tol_ppm)
    frac = it[hit].sum() / it.sum()
    dec = []
    for sh in (0.0237, 0.0431, -0.0317, -0.0529, 0.0713, -0.0841):
        dec.append(it[_explained(ions, mz, tol_da, tol_ppm, sh)].sum() / it.sum())
    dfrac = float(np.mean(dec))
    diag = mz >= min_diag_mz            # tiny fragments are explained by chance far too easily
    dec_n = float(np.mean([(_explained(ions, mz, tol_da, tol_ppm, sh) & diag).sum() for sh in (0.0237, -0.0317, 0.0431, -0.0529)]))
    adj_n = max(0.0, float((hit & diag).sum()) - dec_n)
    lab = pred.assign(r=pred["mz"].round(3)).drop_duplicates("r").set_index("r")["label"]
    expl = [f"{a:.4f}={lab.get(round(ions[np.argmin(np.abs(ions - a))], 3), '?')}" for a in mz[hit]][:15]
    return dict(bi_expl_n=int(hit.sum()), bi_expl_adj_n=float(adj_n), bi_expl_frac=float(frac), bi_decoy_frac=dfrac,
                bi_score=float(max(0.0, frac - dfrac)), bi_n_considered=int(k.sum()), bi_expl_peaks=";".join(expl))


# ------------------------------------------------------------------------------------------------ MetFrag wrapper
_METFRAG_MODE = {"[M+H]+": 1, "[M+NH4]+": 18, "[M+Na]+": 23, "[M+K]+": 39, "[M-H]-": -1, "[M+Cl]-": 35, "[M+HCOO]-": 45,
                 "[M+CH3COO]-": 59}


def metfrag_available(cfg):
    jar = (cfg.get("insilico", {}) or {}).get("metfrag_jar") or cfg.get("paths", {}).get("metfrag_jar")
    java = cfg.get("paths", {}).get("java", "java")
    return bool(jar and os.path.exists(jar) and shutil.which(java)), jar, java


def run_metfrag(cands: pd.DataFrame, peaks_mz, peaks_int, ion_mz, polarity, adduct, workdir, jar, java="java",
                frag_ppm=10.0, frag_abs=0.001, depth=2, mb_folder=None, timeout=300, tag="x"):
    """cands: columns db_id,name,formula,exact_mass,smiles,inchikey.  Returns DataFrame(db_id, mf_score, mf_frag_score, mf_expl, mf_used)."""
    os.makedirs(workdir, exist_ok=True)
    c = cands.drop_duplicates("db_id")
    cpath = os.path.join(workdir, f"{tag}_cands.csv"); ppath = os.path.join(workdir, f"{tag}_peaks.txt")
    def _inchi(smi):
        try:
            from rdkit import Chem
            from rdkit.Chem.inchi import MolToInchi
            return MolToInchi(Chem.MolFromSmiles(smi))
        except Exception:
            return ""
    c_out = pd.DataFrame({"Identifier": c["db_id"], "MonoisotopicMass": c["exact_mass"], "MolecularFormula": c["formula"],
                          "InChI": [_inchi(x) for x in c["smiles"]], "SMILES": c["smiles"], "InChIKey": c["inchikey"],
                          "CompoundName": c["name"]})
    c_out = c_out[c_out["InChI"].astype(str).str.startswith("InChI=")]      # MetFrag cannot process structures without InChI
    if c_out.empty:
        return pd.DataFrame()
    c_out.to_csv(cpath, index=False)
    with open(ppath, "w") as f:
        for a, b in zip(peaks_mz, peaks_int):
            f.write(f"{a:.5f}\t{b:.1f}\n")
    pos = polarity == "positive"
    P = [f"PeakListPath = {ppath}", "MetFragDatabaseType = LocalCSV", f"LocalDatabasePath = {cpath}",
         f"IonizedPrecursorMass = {ion_mz}", "DatabaseSearchRelativeMassDeviation = 1000000",
         f"PrecursorIonMode = {_METFRAG_MODE.get(normalize_adduct(adduct), 1 if pos else -1)}",
         f"IsPositiveIonMode = {pos}", f"FragmentPeakMatchAbsoluteMassDeviation = {frag_abs}",
         f"FragmentPeakMatchRelativeMassDeviation = {frag_ppm}",
         "MetFragCandidateWriter = CSV", f"SampleName = {tag}", f"ResultsPath = {workdir}",
         f"MaximumTreeDepth = {depth}", "MetFragPreProcessingCandidateFilter = UnconnectedCompoundFilter,IsotopeFilter",
         "MetFragPostProcessingCandidateFilter = InChIKeyFilter", "NumberThreads = 1", "UseSmiles = true"]
    if mb_folder and os.path.exists(mb_folder):
        P += ["MetFragScoreTypes = FragmenterScore,OfflineIndividualMoNAScore", "MetFragScoreWeights = 1.0,1.0",
              f"OfflineSpectralDatabaseFile = {mb_folder}"]
    else:
        P += ["MetFragScoreTypes = FragmenterScore", "MetFragScoreWeights = 1.0"]
    prm = os.path.join(workdir, f"{tag}_params.txt")
    with open(prm, "w") as f:
        f.write("\n".join(P) + "\n")
    try:
        # MetFrag mis-parses an absolute params path containing spaces -> run inside workdir with the bare file name
        rp = subprocess.run([java, "-jar", jar, os.path.basename(prm)], cwd=workdir, capture_output=True, timeout=timeout, check=False, text=True)
        if rp.returncode != 0:
            log.debug("MetFrag stderr: %s", (rp.stderr or "")[-400:])
    except Exception as e:
        log.warning("MetFrag failed (%s): %s", tag, e)
        return pd.DataFrame()
    res = os.path.join(workdir, f"{tag}.csv")
    if not os.path.exists(res):
        return pd.DataFrame()
    r = pd.read_csv(res)
    out = pd.DataFrame({"db_id": r["Identifier"].astype(str), "mf_score": r["Score"], "mf_frag_score": r.get("FragmenterScore"),
                        "mf_mona": r.get("OfflineIndividualMoNAScore", 0.0), "mf_expl": r.get("NoExplPeaks"),
                        "mf_used": r.get("NumberPeaksUsed")})
    return out
