"""Experimental MS/MS spectral library (MSP -> SQLite) and cosine search.  Supports MassBank / GNPS / MoNA NIST-format MSP
and your own in-house standard spectra (source='inhouse' -> enables Level 1 together with a standard's RT)."""
from __future__ import annotations
import os, re, sqlite3, logging, time
import numpy as np
import pandas as pd
from .common import normalize_adduct

log = logging.getLogger("hrms_nts")

SCHEMA = """
CREATE TABLE IF NOT EXISTS spectra(
  id INTEGER PRIMARY KEY, name TEXT, inchikey TEXT, ik1 TEXT, smiles TEXT, formula TEXT, adduct TEXT, polarity TEXT,
  precursor_mz REAL, ce TEXT, instrument TEXT, source TEXT, n_peaks INTEGER, peaks BLOB);
CREATE INDEX IF NOT EXISTS idx_pm ON spectra(polarity, precursor_mz);
CREATE INDEX IF NOT EXISTS idx_ik ON spectra(ik1);
"""


def _iter_msp(path):
    meta, peaks, npk = {}, [], None
    with open(path, "r", encoding="utf-8", errors="replace", buffering=1 << 20) as f:
        for line in f:
            s = line.strip()
            if npk is not None and (npk > 0 or not s):
                if not s:
                    if peaks or npk == 0:
                        yield meta, peaks
                    meta, peaks, npk = {}, [], None
                    continue
                parts = s.replace(",", " ").split()
                if len(parts) >= 2:
                    try:
                        peaks.append((float(parts[0]), float(parts[1])))
                    except ValueError:
                        pass
                if len(peaks) >= npk:
                    yield meta, peaks
                    meta, peaks, npk = {}, [], None
                continue
            if not s:
                continue
            if ":" in s:
                k, v = s.split(":", 1)
                k = k.strip().lower().replace(" ", "_"); v = v.strip()
                if k in ("num_peaks", "numpeaks"):
                    try:
                        npk = int(float(v))
                    except ValueError:
                        npk = 0
                else:
                    if k not in meta:
                        meta[k] = v
    if peaks:
        yield meta, peaks


def _pol(meta):
    m = (meta.get("ion_mode") or meta.get("ionmode") or "").lower()
    if m.startswith("pos"): return "positive"
    if m.startswith("neg"): return "negative"
    a = (meta.get("precursor_type") or meta.get("precursortype") or "")
    if a.endswith("+"): return "positive"
    if a.endswith("-"): return "negative"
    return ""


def build_library(msp_paths, db_path, source_names=None, max_peaks=100, min_peaks=3, limit=None, append=False,
                  only_esi=True):
    if os.path.exists(db_path) and not append:
        os.remove(db_path)
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
    con = sqlite3.connect(db_path)
    con.executescript(SCHEMA)
    ins = ("INSERT INTO spectra(name,inchikey,ik1,smiles,formula,adduct,polarity,precursor_mz,ce,instrument,source,n_peaks,peaks) "
           "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)")
    total = 0
    for pi, p in enumerate(msp_paths):
        src = (source_names[pi] if source_names else os.path.splitext(os.path.basename(p))[0])
        t0 = time.time(); batch = []; kept = skipped = 0
        for meta, peaks in _iter_msp(p):
            try:
                pmz = float(meta.get("precursormz") or meta.get("precursor_mz") or meta.get("parent") or 0)
            except ValueError:
                pmz = 0
            pol = _pol(meta)
            if pmz <= 0 or len(peaks) < min_peaks or not pol:
                skipped += 1; continue
            ins_type = (meta.get("instrument_type") or meta.get("instrumenttype") or "")
            if only_esi and ins_type and "ESI" not in ins_type.upper() and "LC" not in ins_type.upper():
                skipped += 1; continue
            if (meta.get("spectrum_type") or "MS2").upper() not in ("MS2", "MS/MS", "MSN"):
                skipped += 1; continue
            a = np.array(peaks, dtype=np.float64)
            a = a[a[:, 1] > 0]
            if len(a) > max_peaks:
                a = a[np.argsort(a[:, 1])[::-1][:max_peaks]]
            a = a[np.argsort(a[:, 0])]
            a[:, 1] = a[:, 1] / a[:, 1].max() * 1000.0
            ik = (meta.get("inchikey") or "").strip()
            add = normalize_adduct(meta.get("precursor_type") or meta.get("precursortype") or "")
            batch.append((meta.get("name", "")[:200], ik, ik.split("-")[0], meta.get("smiles", ""), meta.get("formula", ""), add, pol,
                          pmz, meta.get("collision_energy") or meta.get("collisionenergy") or "", ins_type[:60], src, len(a),
                          a.astype(np.float32).tobytes()))
            kept += 1
            if len(batch) >= 5000:
                con.executemany(ins, batch); con.commit(); batch = []
            if limit and kept >= limit:
                break
        if batch:
            con.executemany(ins, batch); con.commit()
        total += kept
        log.info("library %-28s kept %7d skipped %7d (%.0fs)", src, kept, skipped, time.time() - t0)
    con.close()
    return total


def cosine(q_mz, q_int, l_mz, l_int, tol_da=0.01, tol_ppm=20.0, int_pow=0.5, mz_pow=0.0):
    if len(q_mz) == 0 or len(l_mz) == 0:
        return 0.0, 0
    wq = np.power(q_int, int_pow) * np.power(q_mz, mz_pow); wl = np.power(l_int, int_pow) * np.power(l_mz, mz_pow)
    nq, nl = np.linalg.norm(wq), np.linalg.norm(wl)
    if nq == 0 or nl == 0:
        return 0.0, 0
    pairs = []
    for i, m in enumerate(q_mz):
        tol = max(tol_da, m * tol_ppm * 1e-6)
        lo, hi = np.searchsorted(l_mz, [m - tol, m + tol])
        for j in range(lo, hi):
            pairs.append((wq[i] * wl[j], i, j))
    if not pairs:
        return 0.0, 0
    pairs.sort(reverse=True)
    ui, uj, tot, n = set(), set(), 0.0, 0
    for p, i, j in pairs:
        if i in ui or j in uj:
            continue
        ui.add(i); uj.add(j); tot += p; n += 1
    return float(tot / (nq * nl)), n


class SpectralLibrary:
    def __init__(self, db_path):
        if not os.path.exists(db_path):
            raise FileNotFoundError(f"spectral library not built: {db_path}  (run: build-library)")
        self.con = sqlite3.connect(db_path)
        self.n = self.con.execute("SELECT COUNT(*) FROM spectra").fetchone()[0]

    def candidates(self, precursor_mz, polarity, ppm=20.0, tol_da=0.005):
        tol = max(tol_da, precursor_mz * ppm * 1e-6)
        cur = self.con.execute("SELECT id,name,inchikey,ik1,adduct,ce,instrument,source,n_peaks,peaks,smiles,formula,precursor_mz FROM spectra "
                               "WHERE polarity=? AND precursor_mz BETWEEN ? AND ?", (polarity, precursor_mz - tol, precursor_mz + tol))
        return cur.fetchall()

    def search(self, precursor_mz, polarity, q_mz, q_int, ppm=20.0, frag_tol_da=0.01, frag_tol_ppm=20.0,
               min_score=0.3, remove_precursor=True, max_return=15):
        q_mz = np.asarray(q_mz, float); q_int = np.asarray(q_int, float)
        if remove_precursor:
            k = q_mz < precursor_mz - 1.5
            if k.sum() >= 1:
                q_mz, q_int = q_mz[k], q_int[k]
        rows = []
        for (lid, name, ik, ik1, add, ce, ins, src, npk, blob, smi, frm, lpm) in self.candidates(precursor_mz, polarity, ppm):
            a = np.frombuffer(blob, dtype=np.float32).reshape(-1, 2).astype(float)
            lm, li = a[:, 0], a[:, 1]
            if remove_precursor:
                k = lm < precursor_mz - 1.5
                if k.sum() >= 1:
                    lm, li = lm[k], li[k]
            sc, n = cosine(q_mz, q_int, lm, li, frag_tol_da, frag_tol_ppm)
            if sc >= min_score:
                rows.append(dict(lib_id=lid, lib_name=name, lib_inchikey=ik, ik1=ik1, lib_adduct=add, lib_ce=ce,
                                 lib_instrument=ins, lib_source=src, lib_n_peaks=npk, score=sc, matched=n,
                                 lib_smiles=smi, lib_formula=frm, lib_precursor_mz=lpm))
        if not rows:
            return pd.DataFrame(columns=["lib_id", "lib_name", "lib_inchikey", "ik1", "lib_adduct", "lib_ce", "lib_instrument",
                                         "lib_source", "lib_n_peaks", "score", "matched", "lib_smiles", "lib_formula", "lib_precursor_mz"])
        return pd.DataFrame(rows).sort_values(["score", "matched"], ascending=False).head(max_return).reset_index(drop=True)

    def spectrum(self, lib_id):
        r = self.con.execute("SELECT peaks FROM spectra WHERE id=?", (int(lib_id),)).fetchone()
        a = np.frombuffer(r[0], dtype=np.float32).reshape(-1, 2).astype(float)
        return a[:, 0], a[:, 1]
