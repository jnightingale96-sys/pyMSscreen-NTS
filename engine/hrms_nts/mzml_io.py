"""mzML reading (pyteomics) into compact numpy structures, with an on-disk cache.

MS1: flattened peak arrays + per-scan offsets (fast vectorised feature picking / EIC).
MS2: list of dicts (precursor, RT, peaks) for DDA data.
Assumes centroided data (MSConvert 'peak picking' on, or ThermoRawFileParser default)."""
from __future__ import annotations
import os, pickle, logging
from dataclasses import dataclass, field
import numpy as np

log = logging.getLogger("hrms_nts")


def _offline_psims():
    """pyteomics/psims download the PSI-MS and unit ontologies from purl.obolibrary.org every time a new mzML is parsed; a dropped connection then kills the run.
    psims ships vendored copies - register them as the resolvers so parsing never needs the network."""
    try:
        from psims.controlled_vocabulary import controlled_vocabulary as C
        from psims.controlled_vocabulary.vendor import _use_vendored_psims_obo, _use_vendored_unit_obo
        memo = {}

        def make(vendored):
            def resolver(cache):
                if vendored not in memo:
                    fh = vendored()
                    with fh:
                        memo[vendored] = C.ControlledVocabulary.from_obo(fh)
                return memo[vendored]
            return resolver
        for uri in ("http://purl.obolibrary.org/obo/ms/psi-ms.obo", "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"):
            C.obo_cache.set_resolver(uri, make(_use_vendored_psims_obo))
        for uri in ("http://purl.obolibrary.org/obo/uo.obo", "http://ontologies.berkeleybop.org/uo.obo"):
            C.obo_cache.set_resolver(uri, make(_use_vendored_unit_obo))
    except Exception as e:                              # older psims: fall back to the library default
        log.debug("offline PSI-MS vocabulary not available: %s", e)


_offline_psims()


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
    instrument: str = ""            # detect_instrument(): 'sciex_qtof' | 'orbitrap' | 'generic'
    acquisition: str = "DDA"        # 'DDA' | 'DIA' (SWATH / AIF / MSE: isolation windows wider than dia_min_window)
    centroided_here: bool = False   # profile spectra were centroided by hrms_nts.centroid on reading
    dia: dict | None = None         # DIA MS2: rt, win (window key), off, mz, inten (flat) + windows table (key, lo, hi, ce)
    start_time: str = ""            # acquisition start (mzML run startTimeStamp)
    gaps: list = field(default_factory=list)   # (from_rt, to_rt) MS1 gaps > 1 min inside the run (incomplete export)

    @property
    def n_ms1(self): return len(self.ms1_rt)

    def xic_max(self, mz0, ppm=10.0):
        """Per MS1 scan the most intense centroid inside +-ppm -> (rt, intensity, m/z). Unlike eic() it does not sum, so a
        neighbouring ion inside a wide window cannot inflate the trace; the m/z array gives the measured mass per scan."""
        k = self.__dict__.get("_ms1_key")
        if k is None or len(k) != len(self.mz):
            k = self.__dict__["_ms1_key"] = scan_key(self.ms1_off, self.mz)
        return _xic_max(self.ms1_rt, self.ms1_off, self.mz, self.inten, np.arange(self.n_ms1), mz0, ppm, k)

    # ---- DIA helpers -----------------------------------------------------------------------------------------------
    def dia_window(self, mz0):
        """Key of the DIA window that best contains mz0 (centre-most when windows overlap), or None."""
        if not self.dia:
            return None
        W = self.dia["windows"]
        ok = [w for w in W if w[1] + 0.5 <= mz0 <= w[2] - 0.5] or [w for w in W if w[1] <= mz0 <= w[2]]
        return min(ok, key=lambda w: abs((w[1] + w[2]) / 2 - mz0)) if ok else None

    def dia_xic(self, win_key, mz0, ppm=20.0, scans=None):
        d = self.dia
        sc = np.flatnonzero(d["win"] == win_key) if scans is None else scans
        if "_key" not in d or len(d["_key"]) != len(d["mz"]):
            d["_key"] = scan_key(d["off"], d["mz"])
        return _xic_max(d["rt"], d["off"], d["mz"], d["inten"], sc, mz0, ppm, d["_key"])


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


_KEY_STRIDE = 1e5          # scan index * stride + m/z: one sorted key array -> vectorised per-scan searchsorted (m/z < 1e5)


def scan_key(off, mz):
    """Monotonic key over a flat centroid array (each scan's m/z sorted): scan_index * stride + m/z."""
    return np.repeat(np.arange(len(off) - 1, dtype=np.float64), np.diff(off)) * _KEY_STRIDE + mz


def _xic_max(rt_all, off, mz, inten, scans, mz0, ppm, key=None):
    """Most intense centroid within +-ppm in each of `scans` -> (rt, intensity, m/z). Vectorised (no per-scan Python loop)."""
    scans = np.asarray(scans, dtype=np.int64)
    I = np.zeros(len(scans)); MZ = np.full(len(scans), np.nan)
    if len(scans) == 0 or len(mz) == 0:
        return rt_all[scans], I, MZ
    if key is None:
        key = scan_key(off, mz)
    lo, hi = mz0 * (1 - ppm * 1e-6), mz0 * (1 + ppm * 1e-6)
    base = scans.astype(np.float64) * _KEY_STRIDE
    i0 = np.searchsorted(key, base + lo, side="left"); i1 = np.searchsorted(key, base + hi, side="right")
    n = i1 - i0
    one = np.flatnonzero(n == 1)
    I[one] = inten[i0[one]]; MZ[one] = mz[i0[one]]
    for j in np.flatnonzero(n > 1):                 # rare: several centroids inside the window
        k = i0[j] + int(np.argmax(inten[i0[j]:i1[j]]))
        I[j] = inten[k]; MZ[j] = mz[k]
    return rt_all[scans], I, MZ


_INSTRUMENT_PATTERNS = (("sciex_qtof", ("SCIEX", "WIFF", "X500", "TripleTOF", "ZenoTOF", "Analyst", "QSTAR")),
                        ("orbitrap", ("Orbitrap", "Exploris", "Q Exactive", "Thermo", "Astral", "LTQ")),
                        ("bruker_qtof", ("Bruker", "timsTOF", "maXis", "impact", "compact")),
                        ("agilent_qtof", ("Agilent", "6545", "6546", "6550", "6560")),
                        ("waters_qtof", ("Waters", "Xevo", "Synapt", "MassLynx")))


def detect_instrument(path: str) -> tuple[str, str]:
    """(instrument family, acquisition start time) from the mzML header (first ~40 kB)."""
    try:
        with open(path, "rb") as f:
            head = f.read(40000).decode("utf-8", "ignore")
    except OSError:
        return "generic", ""
    import re
    ts = re.search(r'startTimeStamp="([^"]+)"', head)
    hl = head.lower()
    for fam, keys in _INSTRUMENT_PATTERNS:
        if any(k.lower() in hl for k in keys):
            return fam, ts.group(1) if ts else ""
    return "generic", ts.group(1) if ts else ""

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
              ms2_min_intensity: float = 0.0, use_cache: bool = True, cache_dir: str | None = None,
              centroid_min_ms1: float = 150.0, centroid_min_ms2: float = 80.0, dia_min_window: float = 5.0) -> Run:
    """Read one mzML.  Profile spectra are centroided on the fly (hrms_nts.centroid); MS2 scans whose isolation window is
    wider than dia_min_window (SWATH / AIF / MSE) are stored as DIA in run.dia instead of run.ms2."""
    from pyteomics import mzml
    from .centroid import centroid_profile
    st = os.stat(path)
    cdir = cache_dir or os.path.join(os.path.dirname(path), ".nts_cache")
    cpath = os.path.join(cdir, f"{os.path.basename(path)}.{st.st_size}.{int(st.st_mtime)}.{int(store_min_intensity)}"
                               f".c{int(centroid_min_ms1)}-{int(centroid_min_ms2)}.dia{dia_min_window:g}.v2.pkl")
    if use_cache and os.path.exists(cpath):
        try:
            with open(cpath, "rb") as f:
                run = pickle.load(f)
            run.polarity = polarity or run.polarity
            log.info("cache hit: %s", os.path.basename(path))
            return run
        except Exception:
            pass

    instrument, start_time = detect_instrument(path)
    rts, offs, mzs, ints, ms2 = [], [0], [], [], []
    d_rt, d_win, d_off, d_mz, d_int, windows = [], [], [0], [], [], {}
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
            is_profile = "profile spectrum" in sp
            if is_profile: profile = True
            if lvl == 1:
                if is_profile:
                    m, i = centroid_profile(mz, it, centroid_min_ms1)
                    keep = i >= store_min_intensity
                    m, i = m[keep], i[keep]
                else:
                    keep = it >= store_min_intensity
                    m = np.asarray(mz[keep], dtype=np.float64); i = np.asarray(it[keep], dtype=np.float32)
                    order = np.argsort(m); m, i = m[order], i[order]
                rts.append(rt); mzs.append(m); ints.append(i); offs.append(offs[-1] + len(m))
            elif lvl == 2:
                pre = sp.get("precursorList", {}).get("precursor", [{}])[0]
                si = pre.get("selectedIonList", {}).get("selectedIon", [{}])[0]
                iso = pre.get("isolationWindow", {})
                act = pre.get("activation", {})
                ce = act.get("collision energy", act.get("normalized collision energy"))
                tgt = iso.get("isolation window target m/z", si.get("selected ion m/z"))
                lo_off = float(iso.get("isolation window lower offset", 0.5) or 0.5)
                hi_off = float(iso.get("isolation window upper offset", 0.5) or 0.5)
                if is_profile:
                    mz, it = centroid_profile(mz, it, centroid_min_ms2)
                if tgt is not None and lo_off + hi_off > dia_min_window:          # DIA / SWATH window
                    key = round(float(tgt), 1)
                    windows.setdefault(key, (key, float(tgt) - lo_off, float(tgt) + hi_off, float(ce) if ce is not None else np.nan))
                    o = np.argsort(mz)
                    d_rt.append(rt); d_win.append(key); d_mz.append(np.asarray(mz, np.float64)[o]); d_int.append(np.asarray(it, np.float32)[o])
                    d_off.append(d_off[-1] + len(o))
                    continue
                pmz = si.get("selected ion m/z")
                if pmz is None:
                    continue
                keep = it >= ms2_min_intensity
                ms2.append(dict(scan_id=str(sp.get("id", "")), rt=rt, precursor_mz=float(pmz),
                                precursor_int=float(si.get("peak intensity", 0.0) or 0.0),
                                iso_width=lo_off + hi_off,
                                ce=float(ce) if ce is not None else None,
                                mz=np.asarray(mz[keep], dtype=np.float64), inten=np.asarray(it[keep], dtype=np.float32),
                                tic=float(np.sum(it))))
    dia = None
    if d_rt:
        dia = dict(rt=np.asarray(d_rt), win=np.asarray(d_win), off=np.asarray(d_off, dtype=np.int64),
                   mz=np.concatenate(d_mz), inten=np.concatenate(d_int), windows=sorted(windows.values()))
    rts_a = np.asarray(rts)
    gaps = [(float(a), float(b)) for a, b in zip(rts_a[:-1], rts_a[1:]) if b - a > 1.0] if len(rts_a) > 1 else []
    run = Run(path=path, name=os.path.splitext(os.path.basename(path))[0],
              polarity=polarity or (list(pol_seen)[0] if len(pol_seen) == 1 else None),
              ms1_rt=rts_a, ms1_off=np.asarray(offs, dtype=np.int64),
              mz=np.concatenate(mzs) if mzs else np.array([]), inten=np.concatenate(ints) if ints else np.array([], np.float32),
              ms2=ms2, profile_warning=profile, truncated=truncated, instrument=instrument,
              acquisition="DIA" if dia is not None and len(d_rt) > len(ms2) else "DDA", centroided_here=profile, dia=dia,
              start_time=start_time, gaps=gaps)
    if truncated:
        log.warning("%s is TRUNCATED (file ends mid-spectrum - interrupted conversion or copy). Using the %d MS1 scans up to "
                    "RT %.2f min; re-convert the .raw file for complete results.", run.name, run.n_ms1, rts[-1] if rts else 0.0)
    for a, b in gaps:
        log.warning("%s has NO MS1 scans between %.2f and %.2f min (incomplete export?) - compounds eluting there cannot be found; "
                    "re-export the raw file.", run.name, a, b)
    if profile:
        log.info("%s: profile-mode spectra centroided on reading (hrms_nts.centroid, MS1 floor %g, MS2 floor %g)",
                 run.name, centroid_min_ms1, centroid_min_ms2)
    log.info("read %s: %s, %s, %d MS1 scans, %d MS1 peaks kept (>=%g), %d DDA MS2, %d DIA MS2 in %d windows, polarity=%s",
             run.name, instrument, run.acquisition, run.n_ms1, len(run.mz), store_min_intensity, len(ms2), len(d_rt),
             len(windows), run.polarity)
    if use_cache:
        try:
            os.makedirs(cdir, exist_ok=True)
            with open(cpath, "wb") as f:
                pickle.dump(run, f, protocol=4)
        except Exception as e:
            log.debug("cache write failed: %s", e)
    return run


def read_for_config(path: str, cfg: dict, polarity: str | None = None) -> Run:
    """read_mzml with the thresholds of the instrument profile that applies to this file (see common.config_for_file)."""
    from .common import resolve
    c = cfg
    return read_mzml(path, polarity, store_min_intensity=c["ms1"]["store_min_intensity"], cache_dir=resolve(c, c["paths"].get("cache_dir")),
                     centroid_min_ms1=c.get("centroid", {}).get("min_height_ms1", 150.0),
                     centroid_min_ms2=c.get("centroid", {}).get("min_height_ms2", 80.0),
                     dia_min_window=c.get("dia", {}).get("min_window_da", 5.0))
