"""QA/QC across injections: replicate presence rule (e.g. 2/3), QC stability, blank rule, internal-standard QC.

Sample sheet columns: file (stem or name), type (sample | blank | qc | cal), group (biological sample / site / treatment),
replicate (optional label), conc (cal level).  Replicate injections of one sample share a group."""
from __future__ import annotations
import math, os, re
import numpy as np
import pandas as pd

TYPE_ALIASES = {"calibration": "cal", "calibrant": "cal", "standard": "cal", "std": "cal", "solvent blank": "blank",
                "procedural blank": "blank", "pooled qc": "qc", "pool": "qc", "qc pool": "qc", "spiked blank": "qc"}


def read_sample_sheet(path_or_df) -> pd.DataFrame:
    s = path_or_df.copy() if isinstance(path_or_df, pd.DataFrame) else (
        pd.read_excel(path_or_df) if str(path_or_df).lower().endswith(("xlsx", "xls")) else pd.read_csv(path_or_df))
    s.columns = [str(c).strip().lower() for c in s.columns]
    s["file"] = s["file"].astype(str).map(lambda f: os.path.splitext(os.path.basename(f))[0])
    s["type"] = s.get("type", "sample").astype(str).str.strip().str.lower().map(lambda t: TYPE_ALIASES.get(t, t))
    if "group" not in s.columns: s["group"] = s["file"]
    s["group"] = s["group"].fillna(s["file"]).astype(str)
    if "conc" not in s.columns: s["conc"] = np.nan
    return s


def guess_sample_sheet(files) -> pd.DataFrame:
    """Starting point from file names: blank / QC / calibration patterns, replicate suffixes (_1, _2, -R1 ...) -> group."""
    rows = []
    for f in files:
        stem = os.path.splitext(os.path.basename(f))[0]
        low = stem.lower()
        if re.search(r"(blank|blk|solvent|eq\d+|^mb\d*|procedural)", low): t = "blank"
        elif re.search(r"(^|[_\-\s])(qc|pool(ed)?qc|qcpool)([_\-\s\d]|$)", low): t = "qc"
        elif re.search(r"(cal|std|standard|level)[_\-\s]?\d", low): t = "cal"
        else: t = "sample"
        grp = re.sub(r"([_\-\s](r|rep|inj)?\d{1,2})$", "", stem, flags=re.I)
        m = re.search(r"(\d+(?:\.\d+)?)\s*(ng|ug|µg|pg)", low)
        rows.append(dict(file=stem, type=t, group=grp, replicate=stem[len(grp):].strip("_- ") or "1",
                         conc=float(m.group(1)) if (t == "cal" and m) else np.nan))
    return pd.DataFrame(rows)


def parse_rule(rule) -> tuple[float, str]:
    """'2/3' -> (0.667, 'fraction'); '2' or '1' -> at least that many injections (count); '0.5' / '50%' -> fraction."""
    s = str(rule).strip()
    if "/" in s:
        a, b = s.split("/"); return float(a) / float(b), "fraction"
    if s.endswith("%"):
        return float(s[:-1]) / 100, "fraction"
    if s.isdigit():
        return float(s), "count"
    v = float(s)
    return (v, "fraction") if v <= 1 else (v, "count")


def required(n, rule):
    v, kind = parse_rule(rule)
    return max(1, math.ceil(v * n - 1e-9)) if kind == "fraction" else min(int(v), n)


def replicate_presence(aligned: pd.DataFrame, sheet: pd.DataFrame, rule="2/3", value="area"):
    """aligned: one row per feature per sample with consensus_id (aggregate.align_features) and blank_flag.
    -> long table consensus_id x group: n_detected, n_injections, required, present, mean/RSD of the value."""
    s = sheet[sheet.type == "sample"]
    if s.empty or aligned.empty:
        return pd.DataFrame()
    a = aligned[aligned["sample"].isin(s.file)].copy()
    if "blank_flag" in a.columns:
        a = a[~a["blank_flag"].fillna(False).astype(bool)]
    a = a.merge(s[["file", "group"]], left_on="sample", right_on="file", how="left")
    size = s.groupby("group").file.nunique()
    g = a.groupby(["consensus_id", "group"]).agg(n_detected=("sample", "nunique"), mean_value=(value, "mean"), sd_value=(value, "std")).reset_index()
    g["n_injections"] = g.group.map(size)
    g["required"] = [required(n, rule) for n in g.n_injections]
    g["present"] = g.n_detected >= g.required
    g["rsd_pct"] = 100 * g.sd_value / g.mean_value
    g["rule"] = str(rule)
    return g.drop(columns="sd_value")


def qc_stability(aligned: pd.DataFrame, sheet: pd.DataFrame, max_rsd=30.0, min_frac=0.5, value="area"):
    """Pooled-QC injections: detection rate and area RSD per consensus feature -> qc_pass (RSD <= max_rsd, detected in >= min_frac)."""
    q = sheet[sheet.type == "qc"].file
    if len(q) == 0 or aligned.empty:
        return pd.DataFrame()
    a = aligned[aligned["sample"].isin(q)]
    g = a.groupby("consensus_id").agg(qc_detected=("sample", "nunique"), qc_mean=(value, "mean"), qc_sd=(value, "std")).reset_index()
    g["qc_injections"] = len(q)
    g["qc_rsd_pct"] = 100 * g.qc_sd / g.qc_mean
    g["qc_pass"] = (g.qc_detected / len(q) >= min_frac) & ((g.qc_rsd_pct <= max_rsd) | (len(q) < 2))
    return g.drop(columns="qc_sd")


def is_qc(assigned: pd.DataFrame, sheet: pd.DataFrame | None = None, max_rsd=30.0):
    """Internal-standard QC from the targeted screen: per IS area RSD, RT range and mass-error spread across injections
    (blanks excluded); per injection the median IS response relative to all injections (injection / matrix effect)."""
    if assigned is None or assigned.empty:
        return pd.DataFrame(), pd.DataFrame()
    a = assigned
    if sheet is not None:
        a = a[~a.file.isin(sheet.loc[sheet.type == "blank", "file"])]
    per_is = a.groupby("name").agg(n_injections=("file", "nunique"), area_mean=("int_area", "mean"), area_sd=("int_area", "std"),
                                   rt_min=("int_apex_rt", "min"), rt_max=("int_apex_rt", "max"), ppm_median=("ppm", "median"),
                                   ppm_max_abs=("ppm", lambda v: v.abs().max())).reset_index()
    per_is["area_rsd_pct"] = 100 * per_is.area_sd / per_is.area_mean
    per_is["rt_range_min"] = per_is.rt_max - per_is.rt_min
    per_is["qc_pass"] = (per_is.area_rsd_pct <= max_rsd) | (per_is.n_injections < 2)
    rel = a.assign(rel=a.int_area / a.groupby("name").int_area.transform("median"))
    per_inj = rel.groupby("file").agg(n_IS=("name", "nunique"), median_IS_response_rel=("rel", "median"),
                                      median_ppm=("ppm", "median")).reset_index()
    per_inj["flag"] = np.where(per_inj.median_IS_response_rel < 0.5, "IS response < 50 % of the batch median (suppression / injection problem)", "")
    return per_is.drop(columns=["area_sd", "rt_min", "rt_max"]), per_inj
