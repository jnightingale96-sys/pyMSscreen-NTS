"""Polyether oligomer series (PEG / PPG) - ubiquitous plastic / lab / surfactant contaminants that appear in every ESI(+) run as evenly spaced series
(44.026 Da ethylene oxide, 58.042 Da propylene oxide).  They are not environmental micropollutants of interest, so they are flagged and left out of the PEC-vs-MEC ranking.

Rule: name says PEG / PPG / polyethylene|polypropylene glycol, or the structure is a pure polyether (only C and O, acyclic, no carbonyl, >= 3 consecutive -O-C-C- units
and O/C >= 0.3, i.e. glycols and their alkyl ethers; esters and ring-bearing ethoxylates are not flagged)."""
from __future__ import annotations
import re
import pandas as pd

_NAME = re.compile(r"(^|[\s(\[])(peg|ppg|ppo|peo)(\b|\d|-)|poly\s*(ethylene|propylene)\s*(glycol|oxide)|polyoxy(ethylene|propylene)|polyethylene glycol|polypropylene glycol|"
                   r"\b(tri|tetra|penta|hexa|hepta|octa|nona|deca|undeca|dodeca)(ethylene|propylene) glycol\b(?! (diacrylate|dimethacrylate))", re.I)
_CHAIN = None


def _chain_patt():
    global _CHAIN
    if _CHAIN is None:
        from rdkit import Chem
        _CHAIN = Chem.MolFromSmarts("[OX2;!R]-[CX4;!R]-[CX4;!R]-[OX2;!R]-[CX4;!R]-[CX4;!R]-[OX2;!R]-[CX4;!R]-[CX4;!R]-[OX2;!R]")      # >= 3 repeat units incl. both end oxygens
    return _CHAIN


def _pure_polyether(smiles: str) -> bool:
    if not smiles or str(smiles) in ("nan", "None"): return False
    from rdkit import Chem
    m = Chem.MolFromSmiles(str(smiles))
    if m is None: return False
    if any(a.GetSymbol() not in ("C", "O") for a in m.GetAtoms()): return False
    hits = m.GetSubstructMatches(_chain_patt())
    if not hits: return False
    nC = sum(a.GetSymbol() == "C" for a in m.GetAtoms()); nO = sum(a.GetSymbol() == "O" for a in m.GetAtoms())
    carbonyl = m.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=[OX1]"))
    return not any(a.GetIsAromatic() or a.IsInRing() for a in m.GetAtoms()) and not carbonyl and nO / max(nC, 1) >= 0.3     # PEG O/C = 0.5, PPG 0.33; alkyl-capped / branched variants stay above 0.3


def classify_polymer(name: str, smiles: str = "", names: str = "") -> tuple[bool, str]:
    txt = " | ".join(str(x) for x in (name, names) if x and str(x) != "nan")
    if any(_NAME.search(t) for t in txt.split(" | ")): return True, "PEG/PPG oligomer series"
    if _pure_polyether(smiles): return True, "PEG/PPG-type polyether"
    return False, ""


def flag_polymers(df: pd.DataFrame) -> pd.DataFrame:
    res = [classify_polymer(r.get("name", ""), r.get("smiles", ""), r.get("names", "")) for r in df.to_dict("records")]
    out = df[["cand_id"]].copy(); out["polymer"] = [r[0] for r in res]; out["polymer_category"] = [r[1] for r in res]
    return out
