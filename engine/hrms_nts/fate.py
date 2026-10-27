"""Fate model access for PEC: drives the project's Excel Trapp/SimpleTreat workbook (Model_Worst_Case_PEC_PLANT .xlsm) through COM,
overriding the effluent/irrigation-water concentration and the chemical inputs, and reading the fruit/grain uptake factor.

    k_fruit  [mg/kg fw per mg/L irrigation water]   ->   PEC grain (ng/g dw) = k_fruit x C_water(ng/L) x 0.005      (0.005 = 1e-6 x 1000 / 0.2 fw->dw)

Every k_fruit is cached in a CSV keyed by the full input set, so reruns (and colleagues without Excel) reuse it.  A `backend` callable can
replace Excel (tests, other fate models).  Chemicals with missing properties are evaluated over a GRID and returned as (low, high)."""
from __future__ import annotations
import hashlib, itertools, logging, os
import numpy as np
import pandas as pd
from .properties import parse_grid

log = logging.getLogger("hrms_nts")


class ExcelFate:
    """Same cell map as the project's fate_driver.py (Input!B6/B7/B13/B17-B24/B29-B31/B64/F11, Koc prediction!C4/C9, PEC_Soil_PEC_Plant!D20)."""

    def __init__(self, model_path, depth_m=0.1, irrigation_l_m2_d=1.5, khsa=10000.0):
        import win32com.client, pythoncom
        pythoncom.CoInitialize()
        self.xl = win32com.client.DispatchEx("Excel.Application"); self.xl.Visible = False; self.xl.DisplayAlerts = False
        self.wb = self.xl.Workbooks.Open(model_path, 0, True)
        self.inp = self.wb.Worksheets("Input"); self.koc = self.wb.Worksheets("Koc prediction"); self.pec = self.wb.Worksheets("PEC_Soil_PEC_Plant")
        self.inp.Range("F11").Value = depth_m; self.inp.Range("B64").Value = irrigation_l_m2_d; self.khsa = khsa

    def close(self):
        try: self.wb.Close(False); self.xl.Quit()
        except Exception: pass

    def __call__(self, mw, logkow, pka, base, kaw, dt50_soil_d=None, biowin4=None, quantity="fruit"):
        i = self.inp
        i.Range("B17").Value = mw; i.Range("B18").Value = logkow; self.koc.Range("C4").Value = logkow
        i.Range("B23").Value = pka; i.Range("B24").Value = pka
        self.koc.Range("C9").Value = 1 if base else -1
        i.Range("B22").Value = "Base" if base else ("Acid" if pka < 12 else "Neutral")
        i.Range("B30").Value = kaw; i.Range("B29").Value = self.khsa; i.Range("B13").Value = 0; i.Range("B31").Value = 0
        if biowin4 is not None and np.isfinite(biowin4):
            i.Range("B6").Value = biowin4; i.Range("B7").Formula = "=(10^(6-B6))"            # aquatic DT50 h; soil DT50 = 0.5 x aquatic
        elif dt50_soil_d is not None:
            i.Range("B7").Value = dt50_soil_d * 24 / 0.5
        self.pec.Range("D20").Value = 1.0
        self.xl.Calculate()
        return float(i.Range("G39" if quantity == "soil" else "G45").Value)       # G45 fruit, G39 steady-state soil (both per 1 mg/L irrigation water)


class FateRunner:
    def __init__(self, cfg: dict, backend=None, cache_path: str | None = None):
        from .common import resolve
        self.cfg = cfg; self.pc = cfg.get("pec_screen", {}); fc = self.pc.get("fate", {})
        self.backend = backend; self._excel = None; self.depth = fc.get("depth_m", 0.1); self.irr = fc.get("irrigation_L_m2_d", 1.5)
        self.khsa = fc.get("khsa", 10000.0); self.model = resolve(cfg, cfg["paths"].get("fate_model"))
        self.cache_path = cache_path or resolve(cfg, fc.get("cache", "work/kfruit_cache.csv"))
        self.cache = {}; self.quantity = self.pc.get("_quantity", "fruit")
        if self.cache_path and os.path.exists(self.cache_path):
            c = pd.read_csv(self.cache_path); self.cache = dict(zip(c["key"], c["k_fruit"]))
        self.n_calls = 0

    def _key(self, *vals):
        s = "|".join(f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v) for v in vals + ((self.quantity,) if self.quantity != "fruit" else ()) + (self.depth, self.irr, self.khsa))
        return hashlib.md5(s.encode()).hexdigest()

    def _backend(self):
        if self.backend is not None: return self.backend
        if self._excel is None:
            if not self.model or not os.path.exists(self.model):
                raise RuntimeError("fate model workbook not found (set paths.fate_model in config.local.yaml) and k_fruit not in cache - "
                                   "supply k_fruit in the properties table or run on a machine with Excel + the model workbook")
            self._excel = ExcelFate(self.model, self.depth, self.irr, self.khsa)
        return self._excel

    def save(self):
        if self.cache_path and self.cache:
            os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
            if os.path.exists(self.cache_path):                      # another run (e.g. grain and soil in parallel) may have added keys since we loaded
                try:
                    c = pd.read_csv(self.cache_path); self.cache = {**dict(zip(c["key"], c["k_fruit"])), **self.cache}
                except Exception: pass
            pd.DataFrame({"key": list(self.cache), "k_fruit": list(self.cache.values())}).to_csv(self.cache_path, index=False)

    def close(self):
        if self._excel is not None: self._excel.close(); self._excel = None
        self.save()

    def k_point(self, mw, logkow, species, pka_acid, pka_base, kaw, dt50=None, biowin4=None):
        base = species == "Base"
        pka = pka_base if base else (pka_acid if species == "Acid" else 14.0)
        if not np.isfinite(pka): pka = 14.0
        key = self._key(float(mw), float(logkow), species, float(pka), float(kaw), "bio" if (biowin4 is not None and np.isfinite(biowin4)) else "dt50",
                        float(biowin4) if (biowin4 is not None and np.isfinite(biowin4)) else float(dt50))
        if key not in self.cache:
            self.cache[key] = float(self._backend()(mw, logkow, pka, base, max(float(kaw), 1e-12),
                                                    dt50_soil_d=None if (biowin4 is not None and np.isfinite(biowin4)) else float(dt50),
                                                    biowin4=biowin4 if (biowin4 is not None and np.isfinite(biowin4)) else None,
                                                    **({"quantity": self.quantity} if self.quantity != "fruit" else {})))
            self.n_calls += 1
        return self.cache[key]

    def k_range(self, p) -> dict:
        """p: one row of resolve_properties().  Returns k_fruit low/high/central, n grid points and a note."""
        if p["not_applicable"]:
            return dict(k_low=np.nan, k_high=np.nan, k_central=np.nan, n_grid=0, k_mode="not applicable")
        if np.isfinite(p.get("k_fruit_user", np.nan)):
            k = float(p["k_fruit_user"]); return dict(k_low=k, k_high=k, k_central=k, n_grid=0, k_mode="user k_fruit")
        have_bio = np.isfinite(p["biowin4"]); have_dt = np.isfinite(p["dt50_soil_d"])
        lk = [p["logkow"]] if str(p["logkow_src"]) == "table" else parse_grid(p["grid_logkow"])
        kaw = [p["kaw"]] if np.isfinite(p["kaw"]) else parse_grid(p["grid_kaw"])
        deg = [(None, p["biowin4"])] if have_bio else ([(p["dt50_soil_d"], None)] if have_dt else [(d, None) for d in parse_grid(p["grid_dt50"])])
        ks = []
        for l, k_, (dt, bw) in itertools.product(lk, kaw, deg):
            ks.append(self.k_point(p["mw"], l, p["species"], p["pka_acid"], p["pka_base"], k_, dt50=dt, biowin4=bw))
        ks = np.array(ks)
        return dict(k_low=float(ks.min()), k_high=float(ks.max()), k_central=float(np.exp(np.mean(np.log(np.clip(ks, 1e-30, None))))) if len(ks) > 1 else float(ks[0]),
                    n_grid=len(ks), k_mode="point" if len(ks) == 1 else "grid")
