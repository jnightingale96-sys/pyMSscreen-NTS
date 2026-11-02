"""Peak quality rating (0-10, in the spirit of Compound Discoverer's Peak Rating) and presence confidence (0-100 %).

Peak rating = 10 x weighted mean of five shape/strength scores, each 0-1 (all reported, so a low rating can be explained):
  points      >= 8 data points across the peak = 1 (linear below, 0 at 3)
  fwhm_base   FWHM / base width close to a Gaussian (~0.48): 1 inside 0.25-0.75, falling outside (spikes / plateaus)
  jaggedness  1 - fraction of direction changes of the raw trace inside the peak beyond the single rise and fall
  modality    1 for one maximum; 0.5 for two maxima above 20 % of the apex; 0 for more
  sn          log-scaled signal-to-noise: S/N 3 -> 0, 10 -> 0.5, >= 100 -> 1
Presence confidence = weighted mean of the available evidence (missing terms are left out and the weights renormalised):
  peak rating (0.25), S/N (0.15), Q1 mass accuracy (0.15), isotope pattern (0.10), MS2 evidence (0.20: library cosine or
  explained product ions), replicate frequency (0.10), blank ratio (0.05).  Default pass threshold 70 %.
These are transparent ranking scores, not calibrated probabilities - validate the threshold with spiked/blank samples."""
from __future__ import annotations
import math
import numpy as np
from scipy.ndimage import gaussian_filter1d


def peak_rating(rt, I, start_rt, end_rt, sn=None, noise_floor=150.0):
    rt = np.asarray(rt, float); I = np.asarray(I, float)
    m = (rt >= start_rt) & (rt <= end_rt)
    y = I[m]; n = int(m.sum())
    if n < 3 or y.max() <= 0:
        return dict(peak_rating=0.0, pr_points=0.0, pr_fwhm_base=0.0, pr_jaggedness=0.0, pr_modality=0.0, pr_sn=0.0)
    s_points = float(np.clip((n - 3) / 5, 0, 1))
    sm = gaussian_filter1d(y, 0.8)
    h = sm.max()
    half = np.flatnonzero(sm >= 0.5 * h)
    ratio = (half[-1] - half[0] + 1) / n if len(half) else 0.0
    s_fwhm = 1.0 if 0.25 <= ratio <= 0.75 else float(max(0.0, 1 - abs(ratio - (0.25 if ratio < 0.25 else 0.75)) / 0.25))
    d = np.sign(np.diff(y)); d = d[d != 0]
    changes = int((d[1:] != d[:-1]).sum()) if len(d) > 1 else 0
    s_jag = float(np.clip(1 - max(changes - 1, 0) / max(n - 2, 1), 0, 1))
    pk = [i for i in range(1, len(sm) - 1) if sm[i] >= sm[i - 1] and sm[i] > sm[i + 1] and sm[i] >= 0.2 * h]
    s_mod = 1.0 if len(pk) <= 1 else (0.5 if len(pk) == 2 else 0.0)
    if sn is None:
        sn = h / noise_floor
    s_sn = float(np.clip((math.log10(max(sn, 1e-9)) - math.log10(3)) / (2 - math.log10(3)), 0, 1))
    w = dict(points=0.2, fwhm=0.2, jag=0.2, mod=0.15, sn=0.25)
    rating = 10 * (w["points"] * s_points + w["fwhm"] * s_fwhm + w["jag"] * s_jag + w["mod"] * s_mod + w["sn"] * s_sn)
    return dict(peak_rating=round(float(rating), 1), pr_points=round(s_points, 2), pr_fwhm_base=round(s_fwhm, 2),
                pr_jaggedness=round(s_jag, 2), pr_modality=round(s_mod, 2), pr_sn=round(s_sn, 2))


WEIGHTS = dict(rating=0.25, sn=0.15, ppm=0.15, iso=0.10, ms2=0.20, replicates=0.10, blank=0.05)


def presence_confidence(peak_rating=None, sn=None, ppm=None, ppm_tol=5.0, iso_score=None, lib_score=None, lib_matched=0,
                        n_explained=None, replicate_frac=None, blank_ratio=None, weights=None):
    """-> (confidence %, dict of the component scores used)."""
    c = {}
    if peak_rating is not None and np.isfinite(peak_rating): c["rating"] = float(np.clip(peak_rating / 10, 0, 1))
    if sn is not None and np.isfinite(sn): c["sn"] = float(np.clip((math.log10(max(sn, 1e-9)) - math.log10(3)) / (2 - math.log10(3)), 0, 1))
    if ppm is not None and np.isfinite(ppm): c["ppm"] = float(np.clip(1 - max(abs(ppm) - ppm_tol / 2, 0) / (ppm_tol * 1.5), 0, 1))
    if iso_score is not None and np.isfinite(iso_score): c["iso"] = float(np.clip(iso_score, 0, 1))
    ms2 = []
    if lib_score is not None and np.isfinite(lib_score) and (lib_matched or 0) >= 2: ms2.append(float(np.clip(lib_score, 0, 1)))
    if n_explained is not None and np.isfinite(n_explained): ms2.append(float(np.clip(n_explained / 3, 0, 1)))
    if ms2: c["ms2"] = max(ms2)
    if replicate_frac is not None and np.isfinite(replicate_frac): c["replicates"] = float(np.clip(replicate_frac, 0, 1))
    if blank_ratio is not None:
        br = float(blank_ratio)
        c["blank"] = 1.0 if not np.isfinite(br) or br >= 10 else float(np.clip((br - 1) / 9, 0, 1))
    W = {**WEIGHTS, **(weights or {})}
    c = {k: v for k, v in c.items() if W.get(k, 0) > 0}
    if not c:
        return 0.0, c
    wsum = sum(W[k] for k in c)
    return round(100 * sum(W[k] * v for k, v in c.items()) / wsum, 1), c
