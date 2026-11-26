"""Reference index annotations from the free structure databases (built by tools/build_reference_index.py).

For every identified compound (matched on the InChIKey first block):
  norman_lists     NORMAN Suspect List Exchange collections it appears in (e.g. S010;S066), n_norman_lists
  known_tp         listed in a transformation-product collection (S38, S66, S68, S74, S78, S81, S121)
  tp_of            parent compound(s) of a literature transformation product (+ transformation type)
  np_sources       natural-product / lipid / ontology databases it appears in (COCONUT, LOTUS, NPAtlas, LIPIDMAPS, ChEBI)
  lotus_n_organisms, npatlas_origin, chebi_roles (plant / human / fungal / bacterial metabolite), lipid_category / lipid_class
These are context labels, not identification evidence: they help to separate anthropogenic contaminants from
natural / endogenous compounds and to recognise known transformation products. Levels are not changed.
"""
from __future__ import annotations
import logging
import os
import pandas as pd

log = logging.getLogger("hrms_nts")
COLS = ["norman_lists", "n_norman_lists", "known_tp", "tp_of", "np_sources", "lotus_n_organisms", "npatlas_origin",
        "chebi_roles", "lipid_category", "lipid_class"]
_CACHE = {}


def _paths(cfg):
    from .common import resolve
    p = cfg.get("paths", {})
    return resolve(cfg, p.get("reference_index")), resolve(cfg, p.get("tp_pairs"))


def load(cfg):
    idx_p, tp_p = _paths(cfg)
    key = (idx_p, tp_p)
    if key in _CACHE:
        return _CACHE[key]
    idx = pd.read_csv(idx_p, dtype={"ik1": str}).set_index("ik1") if idx_p and os.path.exists(idx_p) else None
    tp = None
    if tp_p and os.path.exists(tp_p):
        t = pd.read_csv(tp_p, dtype=str).fillna("")
        t["txt"] = t.parent_name.where(t.parent_name != "", "parent " + t.parent_ik1) + \
            t.transformation.map(lambda s: f" ({s})" if s else "") + " [" + t.source + "]"
        tp = t.groupby("ik1_tp").txt.agg(lambda s: "; ".join(sorted(set(s))[:5]))
    if idx is None:
        log.info("reference index not configured (paths.reference_index) - NORMAN list / natural-product labels skipped")
    _CACHE[key] = (idx, tp)
    return idx, tp


def annotate(df: pd.DataFrame, cfg, ik1_col="ik1", inchikey_col="inchikey") -> pd.DataFrame:
    """Adds the COLS to a table of identifications (in place on a copy). Rows without a structure stay empty."""
    idx, tp = load(cfg)
    out = df.copy()
    if idx is None or out.empty:
        return out
    if ik1_col in out.columns:
        k = out[ik1_col].astype(str).str[:14]
    elif inchikey_col in out.columns:
        k = out[inchikey_col].astype(str).str[:14]
    else:
        return out
    for c in COLS:
        if c == "tp_of":
            out[c] = k.map(tp) if tp is not None else ""
        elif c in idx.columns:
            out[c] = k.map(idx[c])
    for c in ("norman_lists", "tp_of", "np_sources", "npatlas_origin", "chebi_roles", "lipid_category", "lipid_class"):
        out[c] = out[c].fillna("")
    # a known TP is a compound the literature parent -> TP tables give as the PRODUCT (TP collections also list parents)
    out["known_tp"] = out["tp_of"] != ""
    out["n_norman_lists"] = out["n_norman_lists"].fillna(0).astype(int)
    out["lotus_n_organisms"] = out["lotus_n_organisms"].fillna(0).astype(int)
    out["origin_hint"] = [_hint(r) for _, r in out.iterrows()]
    return out


def is_natural(r) -> bool:
    """Natural / endogenous evidence: reported from >= 2 organisms (LOTUS) or from NP Atlas, a LIPID MAPS lipid, or a ChEBI metabolite
    role. COCONUT alone and plain ChEBI membership are NOT enough - both also contain synthetic drugs and contaminants."""
    src = str(r.get("np_sources", ""))
    lotus = int(r.get("lotus_n_organisms", 0) or 0) >= 2          # one organism is often a contamination record
    return lotus or any(s in src for s in ("NPAtlas", "LIPIDMAPS")) or "metabolite" in str(r.get("chebi_roles", ""))


def _hint(r):
    """Origin reading of the labels (transformation-product status is the separate known_tp / tp_of columns)."""
    natural = is_natural(r)
    anthro = int(r.get("n_norman_lists", 0) or 0) >= 1
    if natural and anthro:
        return "natural product, also on contaminant lists"
    if natural:
        return "natural / endogenous"
    if anthro:
        return "anthropogenic (NORMAN lists)"
    return ""
