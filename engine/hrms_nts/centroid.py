"""Profile -> centroid peak picking for TOF/Orbitrap profile spectra (validated on SCIEX X500R SWATH data, 2026-10).

Vendor-neutral and vectorised.  Profile exports usually omit zero-intensity points, so the spectrum is split into
contiguous segments (gap > 3x the local point spacing, which on a TOF scales with sqrt(m/z)).  Maxima are found on a
5-point binomial-smoothed trace (raw low-count profiles are spiky and otherwise split one peak into 3-4 centroids);
the centroid m/z is the intensity-weighted mean of the raw points >= 50 % of the apex within +-3 points, and the
reported intensity is the summed intensity of the +-3 points (area proxy)."""
from __future__ import annotations
import numpy as np


def centroid_profile(mz, it, min_height: float = 150.0, half_width: int = 3):
    mz = np.asarray(mz, dtype=np.float64); it = np.asarray(it, dtype=np.float64)
    if len(mz) < 3:
        return np.empty(0), np.empty(0, np.float32)
    o = np.argsort(mz, kind="mergesort")
    mz, it = mz[o], it[o]
    d = np.diff(mz)
    pos = d > 0
    if not pos.any():
        return np.empty(0), np.empty(0, np.float32)
    loc = np.median(d[pos] / np.sqrt(mz[1:][pos])) * np.sqrt(mz[1:])
    n = len(mz)
    segid = np.cumsum(np.r_[0, d > 3 * loc])

    def nb(arr, k, fill=0.0):
        idx = np.arange(n) + k
        ok = (idx >= 0) & (idx < n)
        idx = np.clip(idx, 0, n - 1)
        ok &= segid[idx] == segid
        return np.where(ok, arr[idx], fill), ok

    sm = 6 * it
    for k, w in ((-2, 1), (-1, 4), (1, 4), (2, 1)):
        sm = sm + w * nb(it, k)[0]
    sm /= 16
    sl, _ = nb(sm, -1); sr, _ = nb(sm, 1); sl2, _ = nb(sm, -2); sr2, _ = nb(sm, 2)
    pk = (sm >= sl) & (sm > sr) & (sm >= sl2) & (sm > sr2) & (sm >= min_height * 0.6)
    num = np.zeros(n); den = np.zeros(n); tot = np.zeros(n)
    for k in range(-half_width, half_width + 1):
        y, ok = nb(it, k); x, _ = nb(mz, k)
        use = ok & (y >= 0.5 * sm)
        num += np.where(use, y * x, 0); den += np.where(use, y, 0); tot += np.where(ok, y, 0)
    keep = pk & (den > 0) & (tot >= min_height)
    return num[keep] / den[keep], tot[keep].astype(np.float32)
