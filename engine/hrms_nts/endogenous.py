"""Endogenous / biogenic metabolites: chemicals a maize plant (or any plant tissue) makes itself, so a grain concentration cannot be read as uptake from irrigation
water.  Used by pec-screen / cd-screen to take them out of the PEC-vs-MEC ranking (they stay in the tables, labelled).

Three independent routes (any one flags the chemical):
  1. curated names (primary metabolism, hormones, phenolics, plant lipids ...), matched on the normalised name,
  2. structure rules (SMARTS): free amino acids / dipeptides, long-chain free fatty acids and oxylipins, glycosides, glyco/phospholipids, tetrapyrroles, flavonoid cores,
  3. keywords in the Compound Discoverer 'Class' text ('metabolic intermediate', 'amino acid', 'peptide', 'nucleoside', 'vitamin', ...).
`user_data/endogenous_extra.csv` (name or ik1, category) adds / overrides entries.  Confidence: 'high' (definite plant metabolite) or 'probable' (could also be industrial,
e.g. small dicarboxylic acids, or an unusual plant natural product that is more likely biogenic than irrigation-derived in a plant matrix)."""
from __future__ import annotations
import os, re
import numpy as np
import pandas as pd

HIGH = {
    # central metabolism: organic acids
    "butanedioic acid": "organic acid (TCA)", "(r)-lactic acid": "organic acid", "(s)-lactic acid": "organic acid", "lactic acid": "organic acid", "aconitate": "organic acid (TCA)", "aconitic acid": "organic acid (TCA)",
    "threonate": "sugar acid", "d-glucitol": "sugar alcohol", "glyceraldehyde, d-": "sugar", "2-methylcitric acid": "organic acid (TCA-like)", "n-carbamoylaspartic acid": "pyrimidine biosynthesis",
    "arabinofuranosyluracil": "nucleoside", "laevosan": "anhydro sugar (levoglucosan)", "d-erythro-l-talo-octono-1,4-lactone": "sugar lactone", "nonanedioic acid": "dicarboxylic acid (plant lipid oxidation / cutin)",
    "undecanedioic acid": "dicarboxylic acid (cutin)", "tridecanedioic acid": "dicarboxylic acid (cutin)", "2,6-resorcylic acid": "phenolic acid", "4-hydroxybenzoic acid": "phenolic acid", "syringic acid": "phenolic acid",
    "uric acid": "purine", "tp2/syringic acid": "phenolic acid", "jasmonic acid": "phytohormone", "chlorogenic acid": "phenolic acid (phenylpropanoid)",
    "citric acid": "organic acid (TCA)", "l-(-)-malic acid": "organic acid (TCA)", "malic acid": "organic acid (TCA)", "succinic acid": "organic acid (TCA)",
    "fumaric acid": "organic acid (TCA)", "maleic acid": "organic acid (TCA)", "trans-aconitic acid": "organic acid (TCA)", "itaconic acid": "organic acid", "citraconic acid": "organic acid",
    "glycolic acid": "organic acid", "l-(+)-lactic acid": "organic acid", "mevalonic acid": "isoprenoid precursor", "2-isopropylmalic acid": "amino-acid biosynthesis (leucine)",
    "3-hydroxy-3-methylglutaric acid": "isoprenoid precursor", "2-methylglutaric acid": "organic acid", "phosphonosuccinic acid": "organic acid (phosphonate)",
    # amino acids, peptides
    "dl-tryptophan": "amino acid", "l-tyrosine": "amino acid", "l-aspartic acid": "amino acid", "n-acetyl-l-methionine": "amino acid derivative", "n-isovalerylglycine": "amino acid derivative",
    "leucylproline": "dipeptide", "valylproline": "dipeptide", "norophthalmic acid": "glutathione-type peptide", "s-acetylglutathione": "glutathione", "l-glutathione (reduced)": "glutathione",
    "methyl l-tryptophanate": "amino acid derivative", "goralatide": "peptide",
    # sugars and derivatives
    "sucrose": "sugar", "cellobiose": "sugar", "maltose": "sugar", "lactose": "sugar", "trehalose": "sugar", "raffinose": "sugar", "glucose": "sugar", "d-glucose": "sugar", "fructose": "sugar", "d-fructose": "sugar", "galactose": "sugar", "mannose": "sugar", "xylose": "sugar", "melibiose": "sugar", "turanose": "sugar", "isomaltose": "sugar", "gentiobiose": "sugar", "sorbitol": "sugar alcohol", "mannitol": "sugar alcohol", "myo-inositol": "sugar alcohol", "ascorbic acid": "vitamin C", "muramic acid": "amino sugar", "3,4,5-trihydroxycyclohex-1-ene-1-carboxylic acid": "shikimate pathway",
    # nucleosides / bases / cytokinins
    "adenosine": "nucleoside", "2'-deoxyadenosine": "nucleoside", "xanthine": "purine", "trans-zeatin": "phytohormone (cytokinin)", "cordycepin": "nucleoside",
    # vitamins / cofactors
    "4-pyridoxic acid": "vitamin B6 catabolite",
    # phenylpropanoids / phenolic acids / flavonoids
    "ferulic acid": "phenolic acid", "sinapinic acid": "phenolic acid", "vanillic acid": "phenolic acid", "gallic acid": "phenolic acid", "gentisic acid": "phenolic acid",
    "2,4-dihydroxybenzoic acid": "phenolic acid", "salicylic acid": "phytohormone (salicylate)", "pyrogallol": "phenolic", "prunin": "flavonoid glycoside",
    "3',5,7-trihydroxy-4'-methoxyflavanone": "flavonoid (hesperetin)",
    # oxylipins / plant fatty acids
    "12-oxo phytodienoic acid": "phytohormone (jasmonate precursor)", "9s,13r-12-oxophytodienoic acid": "phytohormone (jasmonate precursor)",
    "gamolenic acid": "plant fatty acid", "stearidonic acid": "plant fatty acid", "alpha-eleostearic acid": "plant fatty acid", "9z,12z,15e-octadecatrienoic acid": "plant fatty acid",
    "1-linoleoyl glycerol": "plant glycerolipid", "azelaic acid": "dicarboxylic acid (plant defence/cutin)", "dodecanedioic acid": "dicarboxylic acid (cutin)",
    "tetradecanedioic acid": "dicarboxylic acid (cutin)",
}
PROBABLE = {   # could also be industrial / not specific to maize, but unlikely to be irrigation-derived in a plant matrix
    "crotonic acid": "small organic acid", "acrylic acid": "small organic acid", "2-furoic acid": "small organic acid", "adipic acid": "dicarboxylic acid", "suberic acid": "dicarboxylic acid",
    "dimethyl maleate": "small organic ester", "3-tert-butyladipic acid": "dicarboxylic acid",
    "andrographolide": "plant natural product", "cafestol": "plant natural product", "cicutoxin": "plant natural product", "xanthumin": "plant natural product",
    "acetyllythosenine": "plant natural product", "abietic acid": "plant resin acid", "5-(4-carboxy-3-methylbutyl)-5,6,8a-trimethyl-3,4,4a,5,6,7,8,8a-octahydronaphthalene-1-carboxylic acid": "plant terpenoid",
}
_NAME_KEYWORDS = [(r"porphyrin|pheophor|phorbide|chlorin\b|bacteriochlor", "tetrapyrrole (chlorophyll catabolite)", "high"),
                  (r"tetrahydroxycyclohex|trihydroxycyclohex", "quinic / shikimic acid", "high"),
                  (r"phosphoglycer|galactose-6-phosphate|glucose-6-phosphate|fructose-6-phosphate|ribose-5-phosphate|\bphosphoenolpyruvate", "sugar phosphate", "high"),
                  (r"^(nona|deca|undeca|dodeca|trideca|tetradeca|hepta|octa)nedioic acid$", "dicarboxylic acid", "probable"),
                  (r"\bglucoside\b|\bglucuronide\b|\bgalactoside\b", "glycoside", "probable"),
                  (r"phosphatidyl|lysophosph|glyceroph", "phospholipid", "high")]
_CD_CLASS = [(r"metabolic (intermediate|pathway)|energy metabolism|krebs", "metabolic intermediate"), (r"amino acid", "amino acid"), (r"peptide", "peptide"),
             (r"nucleoside|nucleotide", "nucleoside"), (r"vitamin", "vitamin"), (r"plant (metabolite|hormone|growth)|phytohormone", "plant metabolite")]
_SMARTS = {
    "dipeptide / peptide (free C-terminus)": ("[NX3][CX4][CX3](=O)[NX3][CX4][CX3](=O)[OX2H1,OX1-]", "high"),
    "alpha-amino acid": ("[NX3;H2,H1;!$(NC=O)][CX4][CX3](=O)[OX2H1,OX1-]", "high"),
    "flavonoid core (flavanone/flavone/isoflavone)": ("[#6]1(=O)[#6]~[#6](-c2ccccc2)[#8]c2ccccc12", "probable"),
    "flavonoid core (flavone/isoflavone)": ("O=c1cc(-c2ccccc2)oc2ccccc12", "probable"),
    "plant glycoside (>= 2 pyranose rings)": ("[OX2]1[CX4][CX4]([OX2H,OX2])[CX4]([OX2H,OX2])[CX4]([OX2H,OX2])[CX4]1", "probable"),
}
_ACID = "[CX3](=O)[OX2H1,OX1-]"
_CHAIN12 = "[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]"


def _norm(n: str) -> str:
    n = str(n or "").split(" (=")[0].strip().lower()
    return re.sub(r"\s+", " ", n)


def _mol(smiles):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) and smiles else None
    if m is None: return None
    return max(Chem.GetMolFrags(m, asMols=True), key=lambda x: x.GetNumHeavyAtoms())


def load_extra(path: str | None) -> dict:
    if not path or not os.path.exists(path): return {}
    from .common import read_table
    t = read_table(path); t.columns = [c.lower().strip() for c in t.columns]
    out = {}
    for r in t.itertuples():
        for k in ("ik1", "name"):
            v = str(getattr(r, k, "") or "").strip()
            if v and v != "nan": out[v.lower() if k == "name" else v] = (str(getattr(r, "category", "user-listed") or "user-listed"), "high")
    return out


def classify(name: str, smiles: str = "", cls_text: str = "", ik1: str = "", extra: dict | None = None, names_all: str = "") -> tuple[bool, str, str]:
    """-> (is_endogenous, category, confidence)."""
    extra = extra or {}
    if ik1 and ik1 in extra: return True, extra[ik1][0], extra[ik1][1]
    for n in {_norm(name), *[_norm(x) for x in str(names_all).split("|") if x.strip()]}:
        if n in extra: return True, extra[n][0], extra[n][1]
        if n in HIGH: return True, HIGH[n], "high"
        if n in PROBABLE: return True, PROBABLE[n], "probable"
        for pat, cat, conf in _NAME_KEYWORDS:
            if re.search(pat, n): return True, cat, conf
    for pat, cat in _CD_CLASS:
        if re.search(pat, str(cls_text or "").lower()): return True, f"CD class: {cat}", "probable"
    m = _mol(smiles)
    if m is None: return False, "", ""
    if not any(a.GetSymbol() == "C" for a in m.GetAtoms()): return True, "inorganic (no carbon)", "high"
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    q = lambda s: m.HasSubstructMatch(Chem.MolFromSmarts(s))
    cnt = lambda s: len(m.GetSubstructMatches(Chem.MolFromSmarts(s)))
    heavy = {a.GetSymbol() for a in m.GetAtoms()}
    nO = sum(a.GetSymbol() == "O" for a in m.GetAtoms()); nN = sum(a.GetSymbol() == "N" for a in m.GetAtoms())
    if q(_SMARTS["dipeptide / peptide (free C-terminus)"][0]) and Descriptors.MolWt(m) < 800 and not q("[a]-[OX2]-[CX4]CN"): return True, "dipeptide / peptide", "high"
    if q(_SMARTS["alpha-amino acid"][0]) and Descriptors.MolWt(m) < 400 and rdMolDescriptors.CalcNumRings(m) <= 2 and not q("[CX3](=O)[NX3]c"): return True, "alpha-amino acid", "high"
    if q(_ACID) and q(_CHAIN12) and not q("a") and nN == 0 and nO <= 5 and heavy <= {"C", "O", "H"}: return True, "free fatty acid / oxylipin", "high"
    if q(_SMARTS["flavonoid core (flavanone/flavone/isoflavone)"][0]) or q(_SMARTS["flavonoid core (flavone/isoflavone)"][0]): return True, "flavonoid", "probable"
    if cnt(_SMARTS["plant glycoside (>= 2 pyranose rings)"][0]) >= 2 and Descriptors.MolWt(m) > 450: return True, "plant glycoside", "probable"
    if cnt(_SMARTS["plant glycoside (>= 2 pyranose rings)"][0]) >= 1 and q("[CX3](=O)O[CX4]") and q(_CHAIN12) and nN == 0: return True, "glycolipid", "high"
    ri = m.GetRingInfo()
    pyrrole_like = sum(1 for r in ri.AtomRings() if len(r) == 5 and sum(m.GetAtomWithIdx(i).GetSymbol() == "N" for i in r) == 1 and sum(m.GetAtomWithIdx(i).GetSymbol() == "C" for i in r) == 4)
    if cnt("[CX4][OX2H]") >= 3 and nN == 0 and not q("a") and Descriptors.MolWt(m) < 420 and heavy <= {"C", "O", "H", "P"}: return True, "sugar / polyol", "high"
    if nN == 4 and pyrrole_like >= 4 and 400 < Descriptors.MolWt(m) < 800: return True, "tetrapyrrole (chlorophyll catabolite)", "high"
    return False, "", ""


def flag_table(df: pd.DataFrame, extra_path: str | None = None, cls_col="cls") -> pd.DataFrame:
    """df needs name, smiles, ik1 (and optionally 'names' with ' | '-joined synonyms and a CD class column).  Adds endogenous, endogenous_category, endogenous_confidence."""
    extra = load_extra(extra_path)
    res = [classify(r.get("name", ""), r.get("smiles", ""), r.get(cls_col, ""), r.get("ik1", ""), extra, r.get("names", "")) for r in df.to_dict("records")]
    out = df.copy()
    out["endogenous"] = [r[0] for r in res]; out["endogenous_category"] = [r[1] for r in res]; out["endogenous_confidence"] = [r[2] for r in res]
    return out
