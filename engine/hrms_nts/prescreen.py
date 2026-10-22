"""Targeted prescreen (pyMSscreen-style): given a compound list, extract EIC + MS2 from each mzML and report QA flags.

Concept and column names follow pyMSscreen (Univ. Luxembourg LCSB-ECI, Artistic License 2.0; independent re-implementation).
Differences: noise estimate uses the EIC median/MAD (not mean), apex is a smoothed peak pick, results also carry ppm error."""
from __future__ import annotations
import logging
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from .common import adduct_table, ppm_err
from .ms2 import MS2Index, peaks_to_string

log = logging.getLogger("hrms_nts")


def _apex(rt, y):
    if len(y) < 5 or y.max() <= 0:
        return None
    sm = np.convolve(y, np.ones(5) / 5, mode="same")
    pk, _ = find_peaks(sm, width=2)
    if len(pk) == 0:
        return int(np.argmax(y))
    prom = [i for i in pk if sm[i] >= 0.05 * sm.max()]
    return int(max(prom or list(pk), key=lambda i: y[i]))


def prescreen(run, compounds: pd.DataFrame, adducts: list[str] | None, ppm=5.0, ms1_thresh=1e5, ms2_thresh=2500.0, s2n=3.0,
              rt_shift_tol=0.5, rt_expected_col=None, ms2_tol_da=0.02) -> pd.DataFrame:
    pol = run.polarity or "positive"
    ad = adduct_table(pol, adducts or (["[M+H]+"] if pol == "positive" else ["[M-H]-"]))
    idx = MS2Index(run)
    rows = []
    for c in compounds.itertuples():
        targets = []
        if pd.notna(getattr(c, "direct_precursor_mz", np.nan)):
            targets.append(("direct", float(c.direct_precursor_mz)))
        if pd.notna(getattr(c, "exact_mass", np.nan)):
            targets += [(a, float(c.exact_mass) + s) for a, s in ad.items()]
        for a, mz in targets:
            rt, y = run.eic(mz, ppm)
            r = dict(ID=c.db_id, Name=c.name, adduct=a, tag=run.name, mz=mz, qa_ms1_exists=False, qa_ms1_good_int=False,
                     qa_ms1_above_noise=False, qa_ms2_exists=False, qa_ms2_good_int=False, qa_ms2_near=False, qa_pass=False,
                     ms1_rt=np.nan, ms1_apex_int=0.0, ms1_snr=np.nan, ms2_rt=np.nan, peak_list="", rt_expected=np.nan)
            ai = _apex(rt, y)
            if ai is not None and y[ai] > 0:
                r.update(qa_ms1_exists=True, ms1_rt=float(rt[ai]), ms1_apex_int=float(y[ai]))
                nz = y[y > 0]; med = np.median(nz); mad = np.median(np.abs(nz - med)) * 1.4826 or max(med * 0.1, 1.0)
                r["ms1_snr"] = float((y[ai] - med) / mad)
                r["qa_ms1_good_int"] = bool(y[ai] > ms1_thresh)
                r["qa_ms1_above_noise"] = bool(r["ms1_snr"] > s2n)
                # peak boundaries for MS2 search
                lo = ai
                while lo > 0 and y[lo - 1] < y[lo] and y[lo - 1] > 0.1 * y[ai]: lo -= 1
                hi = ai
                while hi < len(y) - 1 and y[hi + 1] < y[hi] and y[hi + 1] > 0.1 * y[ai]: hi += 1
                sp = idx.best_spectrum(mz, rt[ai], rt[lo], rt[hi], tol_da=ms2_tol_da, rt_pad=rt_shift_tol)
                if sp is not None:
                    r["qa_ms2_exists"] = True; r["ms2_rt"] = sp["rt"]; r["peak_list"] = peaks_to_string(sp["mz"], sp["inten"])
                    r["qa_ms2_good_int"] = bool(len(sp["inten"]) and sp["inten"].max() > ms2_thresh)
                    r["qa_ms2_near"] = abs(sp["rt"] - rt[ai]) < rt_shift_tol
            if rt_expected_col and hasattr(c, rt_expected_col) and pd.notna(getattr(c, rt_expected_col)):
                r["rt_expected"] = float(getattr(c, rt_expected_col))
            rt_ok = (not np.isfinite(r["rt_expected"])) or (np.isfinite(r["ms1_rt"]) and abs(r["ms1_rt"] - r["rt_expected"]) <= rt_shift_tol)
            r["qa_pass"] = bool(r["qa_ms1_exists"] and r["qa_ms1_good_int"] and r["qa_ms1_above_noise"] and r["qa_ms2_exists"]
                                and r["qa_ms2_good_int"] and r["qa_ms2_near"] and rt_ok)
            rows.append(r)
    return pd.DataFrame(rows)
