"""mzML reading (pyteomics) into compact numpy structures, with an on-disk cache.

MS1: flattened peak arrays + per-scan offsets (fast vectorised feature picking / EIC).
MS2: list of dicts (precursor, RT, peaks) for DDA data.
Assumes centroided data (MSConvert 'peak picking' on, or ThermoRawFileParser default)."""
from __future__ import annotations
import os, pickle, logging
from dataclasses import dataclass, field
import numpy as np

log = logging.getLogger("hrms_nts")


@dataclass
class Run:
    path: str
    name: str
    polarity: str | None
    ms1_rt: np.ndarray              # (n_scans,) minutes
    ms1_off: np.ndarray             # (n_scans+1,) offsets into flat arrays
    mz: np.ndarray                  # flat float64
    inten: np.ndarray               # flat float32
    ms2: list = field(default_factory=list)
    profile_warning: bool = False
    truncated: bool = False         # file ended mid-spectrum (interrupted conversion/copy); only the readable part was used

    @property
    def n_ms1(self): return len(self.ms1_rt)

    def scan_index_of_peaks(self):
        return np.repeat(np.arange(self.n_ms1), np.diff(self.ms1_off))

    def ms1_scan(self, i):
        a, b = self.ms1_off[i], self.ms1_off[i + 1]
        return self.mz[a:b], self.inten[a:b]

    def nearest_scan(self, rt):
        return int(np.argmin(np.abs(self.ms1_rt - rt)))

    def eic(self, mz0, ppm=5.0, rt_min=None, rt_max=None):
        """Summed-intensity EIC (rt, intensity) in a ppm window. Vectorised over the flat arrays."""
        lo, hi = mz0 * (1 - ppm * 1e-6), mz0 * (1 + ppm * 1e-6)
        sel = (self.mz >= lo) & (self.mz <= hi)
        idx = np.flatnonzero(sel)
        sc = np.searchsorted(self.ms1_off, idx, side="right") - 1
        out = np.zeros(self.n_ms1, dtype=float)
        np.add.at(out, sc, self.inten[idx])
        rt = self.ms1_rt
        if rt_min is not None or rt_max is not None:
            m = np.ones(len(rt), bool)
            if rt_min is not None: m &= rt >= rt_min
            if rt_max is not None: m &= rt <= rt_max
            return rt[m], out[m]
        return rt, out


def _rt_minutes(scan):
    t = scan.get("scan start time", scan.get("scan time"))
    if t is None:
        return None
    unit = getattr(t, "unit_info", None)
    t = float(t)
    return t / 60.0 if unit == "second" else t


def _spectra_tolerant(reader, path):
    """Yield spectra; on an XML error after at least one spectrum (truncated file) yield None once and stop."""
    n = 0
    it = iter(reader)
    while True:
        try:
            sp = next(it)
        except StopIteration:
            return
        except Exception as e:
            if n == 0:
                raise
            log.debug("parse stopped in %s after %d spectra: %s", path, n, e)
            yield None
            return
        n += 1
        yield sp


def read_mzml(path: str, polarity: str | None = None, store_min_intensity: float = 1e4,
              ms2_min_intensity: float = 0.0, use_cache: bool = True, cache_dir: str | None = None) -> Run:
    from pyteomics import mzml
    st = os.stat(path)
    cdir = cache_dir or os.path.join(os.path.dirname(path), ".nts_cache")
    cpath = os.path.join(cdir, f"{os.path.basename(path)}.{st.st_size}.{int(st.st_mtime)}.{int(store_min_intensity)}.pkl")
    if use_cache and os.path.exists(cpath):
        try:
            with open(cpath, "rb") as f:
                run = pickle.load(f)
            run.polarity = polarity or run.polarity
            log.info("cache hit: %s", os.path.basename(path))
            return run
        except Exception:
            pass

    rts, offs, mzs, ints, ms2 = [], [0], [], [], []
    pol_seen = set(); profile = False; truncated = False
    with mzml.read(path, read_schema=False) as rd:
        for sp in _spectra_tolerant(rd, path):
            if sp is None:
                truncated = True
                break
            lvl = sp.get("ms level")
            mz = sp.get("m/z array"); it = sp.get("intensity array")
            if mz is None or it is None or len(mz) == 0:
                continue
            rt = _rt_minutes(sp.get("scanList", {}).get("scan", [{}])[0])
            if rt is None:
                continue
            if "negative scan" in sp: pol_seen.add("negative")
            elif "positive scan" in sp: pol_seen.add("positive")
            if "profile spectrum" in sp: profile = True
            if lvl == 1:
                keep = it >= store_min_intensity
                m = np.asarray(mz[keep], dtype=np.float64); i = np.asarray(it[keep], dtype=np.float32)
                order = np.argsort(m)
                rts.append(rt); mzs.append(m[order]); ints.append(i[order]); offs.append(offs[-1] + len(m))
            elif lvl == 2:
                pre = sp.get("precursorList", {}).get("precursor", [{}])[0]
                si = pre.get("selectedIonList", {}).get("selectedIon", [{}])[0]
                pmz = si.get("selected ion m/z")
                if pmz is None:
                    continue
                iso = pre.get("isolationWindow", {})
                act = pre.get("activation", {})
                ce = act.get("collision energy", act.get("normalized collision energy"))
                keep = it >= ms2_min_intensity
                ms2.append(dict(scan_id=str(sp.get("id", "")), rt=rt, precursor_mz=float(pmz),
                                precursor_int=float(si.get("peak intensity", 0.0) or 0.0),
                                iso_width=float(iso.get("isolation window lower offset", 0.5) or 0.5) +
                                          float(iso.get("isolation window upper offset", 0.5) or 0.5),
                                ce=float(ce) if ce is not None else None,
                                mz=np.asarray(mz[keep], dtype=np.float64), inten=np.asarray(it[keep], dtype=np.float32),
                                tic=float(np.sum(it))))
    run = Run(path=path, name=os.path.splitext(os.path.basename(path))[0],
              polarity=polarity or (list(pol_seen)[0] if len(pol_seen) == 1 else None),
              ms1_rt=np.asarray(rts), ms1_off=np.asarray(offs, dtype=np.int64),
              mz=np.concatenate(mzs) if mzs else np.array([]), inten=np.concatenate(ints) if ints else np.array([], np.float32),
              ms2=ms2, profile_warning=profile, truncated=truncated)
    if truncated:
        log.warning("%s is TRUNCATED (file ends mid-spectrum - interrupted conversion or copy). Using the %d MS1 scans up to "
                    "RT %.2f min; re-convert the .raw file for complete results.", run.name, run.n_ms1, rts[-1] if rts else 0.0)
    if profile:
        log.warning("%s contains profile-mode spectra; centroid in MSConvert/ThermoRawFileParser for reliable results.", run.name)
    log.info("read %s: %d MS1 scans, %d MS1 peaks kept (>=%g), %d MS2 scans, polarity=%s",
             run.name, run.n_ms1, len(run.mz), store_min_intensity, len(ms2), run.polarity)
    if use_cache:
        try:
            os.makedirs(cdir, exist_ok=True)
            with open(cpath, "wb") as f:
                pickle.dump(run, f, protocol=4)
        except Exception as e:
            log.debug("cache write failed: %s", e)
    return run
