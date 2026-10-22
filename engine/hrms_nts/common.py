"""Shared constants, adduct tables, formula/isotope helpers, config loading."""
from __future__ import annotations
import os, re, logging
import numpy as np
import yaml

log = logging.getLogger("hrms_nts")

PROTON = 1.007276466
ELECTRON = 0.000548580

# neutral-mass -> ion shifts (add to M to obtain m/z for z=1). Electron mass already included.
ADDUCTS_POS = {"[M+H]+": 1.007276, "[M+Na]+": 22.989218, "[M+K]+": 38.963158, "[M+NH4]+": 18.033823}
ADDUCTS_NEG = {"[M-H]-": -1.007276, "[M+Cl]-": 34.969402, "[M+HCOO]-": 44.998201, "[M+CH3COO]-": 59.013851}

ISO_SPACING = 1.0033548  # 13C-12C

# common neutral losses (formula -> exact mass) used by diagnostic checks / built-in fragmenter
NEUTRAL_LOSSES = {
    "H2O": 18.010565, "NH3": 17.026549, "CO": 27.994915, "CO2": 43.989829, "HCOOH": 46.005479,
    "CH3OH": 32.026215, "CH2O": 30.010565, "SO2": 63.961900, "SO3": 79.956815, "HF": 20.006228,
    "HCl": 35.976678, "C2H4": 28.031300, "C3H6": 42.046950, "CH3": 15.023475, "HCN": 27.010899,
    "H2S": 33.987721, "N2": 28.006148, "NO": 29.997989, "NO2": 45.992904, "C2H2O": 42.010565,
}

ELEMENT_MASS = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956, "S": 31.972071,
                "P": 30.97376163, "F": 18.99840322, "Cl": 34.96885268, "Br": 78.9183371, "I": 126.904473,
                "Na": 22.9897692809, "K": 38.96370668, "Si": 27.9769265325, "B": 11.0093054, "Se": 79.9165218}

# isotopologue abundances for +0, +1, +2 nominal mass offsets
_ISO = {
    "C": (0.9893, 0.0107, 0.0), "H": (0.999885, 0.000115, 0.0), "N": (0.99636, 0.00364, 0.0),
    "O": (0.99757, 0.00038, 0.00205), "S": (0.9493, 0.0076, 0.0429), "Cl": (0.7576, 0.0, 0.2424),
    "Br": (0.5069, 0.0, 0.4931), "P": (1.0, 0, 0), "F": (1.0, 0, 0), "I": (1.0, 0, 0),
    "Si": (0.9223, 0.0467, 0.0310), "Na": (1.0, 0, 0), "K": (0.932581, 0.000117, 0.067302), "B": (0.199, 0.0, 0.0),
    "Se": (0.0089 + 0.0937 + 0.0763, 0.0, 0.2377),
}

_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> dict:
    """'C15H12N2O' -> {'C':15,'H':12,'N':2,'O':1}; tolerant of spaces and charges, ignores '.' salts (takes largest)."""
    if not isinstance(formula, str) or not formula.strip():
        return {}
    parts = [p for p in re.split(r"[.\s]", formula) if p]
    best = {}
    for p in parts:
        d = {}
        for el, n in _FORMULA_RE.findall(re.sub(r"[+\-]\d*$", "", p)):
            if el:
                d[el] = d.get(el, 0) + (int(n) if n else 1)
        if sum(d.values()) > sum(best.values()):
            best = d
    return best


def formula_mass(formula: str | dict) -> float:
    d = parse_formula(formula) if isinstance(formula, str) else formula
    return float(sum(ELEMENT_MASS.get(el, 0.0) * n for el, n in d.items()))


def isotope_pattern(formula, n: int = 3) -> np.ndarray:
    if isinstance(formula, str):
        return _iso_cached(formula, n).copy()
    return _iso_calc(formula, n)


def _iso_calc(d, n=3):
    """Relative abundance (M0=1) of the first n nominal isotopologues (M, M+1, M+2)."""
    pat = np.array([1.0, 0.0, 0.0])
    for el, cnt in d.items():
        a = np.array(_ISO.get(el, (1.0, 0.0, 0.0)))
        for _ in range(min(int(cnt), 400)):
            pat = np.convolve(pat, a)[:3]
    pat = pat / pat[0]
    return pat[:n]


from functools import lru_cache


@lru_cache(maxsize=500000)
def _iso_cached(formula, n):
    return _iso_calc(parse_formula(formula), n)


def ppm_err(obs, theo):
    return (np.asarray(obs) - np.asarray(theo)) / np.asarray(theo) * 1e6


def normalize_adduct(a: str) -> str:
    if not isinstance(a, str):
        return ""
    a = a.strip().replace(" ", "")
    if not a:
        return ""
    if not a.startswith("["):
        a = "[" + a
    if not (a.endswith("]+") or a.endswith("]-") or a.endswith("]2+") or a.endswith("]2-")):
        # e.g. "M+H" -> "[M+H]+" ; "M-H" -> "[M-H]-"
        a = a.rstrip("]")
        a += "]-" if re.search(r"(-H|\+Cl|\+HCOO|\+CH3COO|\+Br|-H2O-H)$", a) else "]+"
    return a


def guess_polarity(name: str, patterns: dict | None = None) -> str | None:
    n = os.path.basename(name).upper()
    patterns = patterns or {"positive": r"(^|[_\-\. ])(POS|POSITIVE|ESI\+)([_\-\. 0-9]|$)",
                            "negative": r"(^|[_\-\. ])(NEG|NEGATIVE|ESI\-)([_\-\. 0-9]|$)"}
    for pol, pat in patterns.items():
        if re.search(pat, n):
            return pol
    return None


def adduct_table(polarity: str, selected: list[str] | None = None) -> dict:
    base = ADDUCTS_POS if polarity == "positive" else ADDUCTS_NEG
    if selected:
        sel = {normalize_adduct(s) for s in selected}
        base = {k: v for k, v in base.items() if k in sel}
    return dict(base)


def _deep_merge(a: dict, b: dict) -> dict:
    out = dict(a)
    for k, v in (b or {}).items():
        out[k] = _deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def load_config(path: str | None = None, overrides: dict | None = None) -> dict:
    """config.yaml (portable defaults) <- config.local.yaml (this machine, git-ignored) <- $HRMS_NTS_CONFIG <- explicit overrides."""
    here = os.path.dirname(os.path.abspath(__file__))
    default = os.path.join(os.path.dirname(here), "config.yaml")
    path = os.path.abspath(path or default)
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    layers = [os.path.join(os.path.dirname(path), "config.local.yaml"), os.environ.get("HRMS_NTS_CONFIG")]
    used = [path]
    for lp in layers:
        if lp and os.path.exists(lp) and os.path.abspath(lp) != path:
            with open(lp, "r", encoding="utf-8") as f:
                cfg = _deep_merge(cfg, yaml.safe_load(f) or {}); used.append(os.path.abspath(lp))
    if overrides:
        cfg = _deep_merge(cfg, overrides)
    cfg["_config_path"] = path
    cfg["_config_layers"] = used
    cfg["_root"] = os.path.dirname(path)
    return cfg


def resolve(cfg: dict, p: str | None) -> str | None:
    if not p:
        return None
    p = os.path.expandvars(os.path.expanduser(str(p)))
    return p if os.path.isabs(p) else os.path.join(cfg["_root"], p)


def setup_logging(level=logging.INFO, logfile: str | None = None):
    fmt = "%(asctime)s %(levelname)-7s %(message)s"
    handlers = [logging.StreamHandler()]
    if logfile:
        handlers.append(logging.FileHandler(logfile, encoding="utf-8"))
    logging.basicConfig(level=level, format=fmt, datefmt="%H:%M:%S", handlers=handlers, force=True)


def read_table(path):
    """CSV reader that ignores full-line '#' comments (SMILES may contain '#', so pandas comment= is unsafe)."""
    import io, pandas as pd
    with open(path, "r", encoding="utf-8-sig") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#") and ln.strip()]
    return pd.read_csv(io.StringIO("".join(lines)))
