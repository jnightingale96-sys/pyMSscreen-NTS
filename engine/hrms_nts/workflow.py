"""Step-by-step NTS project API (used by the Google Colab notebook and usable from any Python session).

    P = Project(data_dir, out_dir)              # finds .raw / .wiff / .mzML
    P.convert_raw()                             # Thermo .raw -> centroided indexed mzML (ThermoRawFileParser)
    P.set_sample_sheet(df)                      # file, type (sample|blank|qc|cal), group, replicate, conc
    P.calibrate(cal_list)                       # parents + IS + concentrations -> linear / quadratic curves (+ plots)
    P.targeted(is_list)                         # IS / targets in every injection: Q1/Q3, levels, integration, blanks
    P.screen()                                  # NTS Level 5 -> 1 for samples, QCs and blanks
    P.qa(rule="2/3")                            # blank rule, replicate presence, QC stability, peak rating, presence confidence
    P.semiquant(api_key=...)                    # elucidation-distance semi-quant with CompTox properties
    P.report()                                  # NTS_report.xlsx
    P.chromatogram("Carbamazepine", "S1_rep1")  # Q1/Q3 traces, integration, rating (plotly figure)
Every step stores its results in out_dir (pickles + CSV), so a notebook can be re-opened and continued."""
from __future__ import annotations
import glob, json, logging, os, pickle, shutil, subprocess, sys
import numpy as np
import pandas as pd

from .common import load_config, config_for_file, resolve, setup_logging
from . import qa as QA

log = logging.getLogger("hrms_nts")
LEVEL_ORDER = ["1", "2a", "2b", "3", "4", "5"]


class Project:
    def __init__(self, data_dir, out_dir=None, config=None, overrides=None):
        self.data_dir = os.path.abspath(data_dir)
        self.out_dir = os.path.abspath(out_dir or os.path.join(self.data_dir, "NTS_project"))
        os.makedirs(self.out_dir, exist_ok=True)
        self.cfg = load_config(config, overrides)
        setup_logging(logfile=os.path.join(self.out_dir, "project.log"))
        try:
            from rdkit import RDLogger
            RDLogger.DisableLog("rdApp.*")            # library SMILES that RDKit rejects are skipped; avoid flooding the notebook
        except Exception:
            pass
        self.sheet = None
        self._state = {}
        p = os.path.join(self.out_dir, "sample_sheet.csv")
        if os.path.exists(p):
            self.sheet = QA.read_sample_sheet(p)

    # ------------------------------------------------------------------------------------------------ files
    def inputs(self):
        def f(*ext):                       # case variants of an extension match the same file on Windows -> de-duplicate
            seen = {}
            for e in ext:
                for x in glob.glob(os.path.join(self.data_dir, "**", "*" + e), recursive=True):
                    if "NTS_project" not in x:
                        seen.setdefault(os.path.normcase(os.path.abspath(x)), x)
            return sorted(seen.values())
        return dict(raw=f(".raw", ".RAW"), sciex=f(".wiff", ".wiff2", ".WIFF", ".WIFF2"), mzml=f(".mzML", ".mzml"))

    def mzml_files(self):
        out = {}
        for p in self.inputs()["mzml"]:
            out.setdefault(os.path.splitext(os.path.basename(p))[0], p)
        return out

    def convert_raw(self, trfp=None, out=None):
        """Thermo .raw -> centroided indexed mzML with ThermoRawFileParser (only files without an mzML yet).
        SCIEX .wiff/.wiff2 need ProteoWizard msconvert on Windows (vendor libraries) - convert those first and upload the mzML."""
        ins = self.inputs(); have = set(self.mzml_files())
        todo = [r for r in ins["raw"] if os.path.splitext(os.path.basename(r))[0] not in have]
        if ins["sciex"]:
            log.warning("%d SCIEX .wiff/.wiff2 file(s): convert them with ProteoWizard msconvert on Windows "
                        "(Peak picking: Vendor, MS levels 1-2) and put the mzML files in %s", len(ins["sciex"]), self.data_dir)
        if not todo:
            return []
        cmd = _trfp_command(trfp or resolve(self.cfg, self.cfg["paths"].get("thermorawfileparser")))
        out = out or os.path.join(self.data_dir, "mzML"); os.makedirs(out, exist_ok=True)
        done = []
        for r in todo:
            log.info("converting %s", os.path.basename(r))
            subprocess.run(cmd + ["-i", r, "-o", out, "-f", "2"], check=True)
            done.append(os.path.join(out, os.path.splitext(os.path.basename(r))[0] + ".mzML"))
        return done

    # ------------------------------------------------------------------------------------------------ sample sheet
    def guess_sample_sheet(self):
        return QA.guess_sample_sheet(list(self.mzml_files()))

    def set_sample_sheet(self, df_or_path):
        s = QA.read_sample_sheet(df_or_path)
        missing = sorted(set(s.file) - set(self.mzml_files()))
        if missing:
            raise ValueError(f"sample sheet lists files without an mzML: {missing}")
        self.sheet = s
        s.to_csv(os.path.join(self.out_dir, "sample_sheet.csv"), index=False)
        log.info("sample sheet: %s", s.type.value_counts().to_dict())
        return s

    def files_of(self, *types):
        f = self.mzml_files()
        if self.sheet is None:
            return list(f.values())
        return [f[x] for x in self.sheet.loc[self.sheet.type.isin(types), "file"] if x in f]

    def _save(self, key, obj):
        self._state[key] = obj
        with open(os.path.join(self.out_dir, f"_{key}.pkl"), "wb") as fh:
            pickle.dump(obj, fh, protocol=4)

    def _load(self, key):
        if key in self._state:
            return self._state[key]
        p = os.path.join(self.out_dir, f"_{key}.pkl")
        if os.path.exists(p):
            with open(p, "rb") as fh:
                self._state[key] = pickle.load(fh)
            return self._state[key]
        return None

    # ------------------------------------------------------------------------------------------------ calibration
    def calibrate(self, cal_list, weight="1/x", min_levels_quadratic=5, units="ng/mL"):
        """cal_list (DataFrame/CSV): name, smiles or formula (+ label), role (target|IS), is_name (IS of each target).
        Concentrations come from the sample sheet (type 'cal', column conc). Response = area / IS area when is_name is set."""
        from .targeted import run_targeted
        from . import calibration as CAL
        if self.sheet is None or not (self.sheet.type == "cal").any():
            raise ValueError("no 'cal' files in the sample sheet")
        lst = cal_list if isinstance(cal_list, str) else self._write_list(cal_list, "calibration_list.csv")
        out = os.path.join(self.out_dir, "calibration")
        files = self.files_of("cal", "blank")
        sheet = os.path.join(self.out_dir, "sample_sheet.csv")
        r = run_targeted(files, lst, out, self.cfg, sheet=sheet, plots=True)
        ch = r["assigned"]
        conc = self.sheet.set_index("file")["conc"]
        comp = pd.read_csv(lst) if isinstance(lst, str) else lst
        # per-compound level concentrations (template columns L1..L8 + sample-sheet cal_level) take precedence over the
        # single 'conc' per calibration injection
        from .templates import level_concentrations
        lvl = level_concentrations(comp)
        lvl_conc = {(n, l): c for n, l, c in lvl[["name", "cal_level", "conc"]].itertuples(index=False)}
        file_level = self.sheet.set_index("file")["cal_level"].astype(str) if "cal_level" in self.sheet.columns else pd.Series(dtype=str)
        comp.columns = [c.lower() for c in comp.columns]
        isn = dict(zip(comp["name"], comp["is_name"].fillna("").astype(str))) if "is_name" in comp.columns else {}
        role = dict(zip(comp["name"], comp["role"].fillna("target").astype(str).str.lower())) if "role" in comp.columns else {}
        area = ch.set_index(["name", "file"])["int_area"]
        rows, abs_rows = [], []
        for x in ch[ch.file.isin(self.sheet.loc[self.sheet.type == "cal", "file"])].itertuples():
            c = lvl_conc.get((x.name, file_level.get(x.file, "")), conc.get(x.file, np.nan))
            if not np.isfinite(c):
                continue
            isname = isn.get(x.name, "")
            if isname:
                if (isname, x.file) not in area.index or area[(isname, x.file)] <= 0:
                    continue
                rows.append(dict(name=x.name, file=x.file, conc=c, response=x.int_area / area[(isname, x.file)]))
            elif role.get(x.name, "target") != "is":
                rows.append(dict(name=x.name, file=x.file, conc=c, response=x.int_area))
            abs_rows.append(dict(name=x.name, file=x.file, conc=c, response=x.int_area, rt=x.int_apex_rt))
        cal = pd.DataFrame(rows)
        summ, fits = CAL.calibrate_table(cal, weight, min_levels_quadratic) if len(cal) else (pd.DataFrame(), {})
        pdf = CAL.plot_curves(fits, os.path.join(out, "calibration_curves.pdf"), units) if fits else None
        # absolute-area response factors (slope of area vs conc) for the semi-quantification
        A = pd.DataFrame(abs_rows); rf = []
        for n, g in A.groupby("name") if len(A) else []:
            f = CAL.fit(g.conc.values, g.response.values, weight, "linear")
            if f and f["coef"][1] > 0:
                smi = comp.set_index("name").get("smiles", pd.Series(dtype=str)).get(n, "")
                rf.append(dict(name=n, smiles=smi, rt=float(g.rt.median()), rf=f["coef"][1], r2=f["r2"], n_levels=f["n"]))
        res = dict(summary=summ, fits=fits, points=cal, rf=pd.DataFrame(rf), targeted=r["xlsx"], pdf=pdf, units=units)
        self._save("calibration", res)
        summ.to_csv(os.path.join(out, "calibration_summary.csv"), index=False)
        return res

    def _write_list(self, df, name):
        p = os.path.join(self.out_dir, name); pd.DataFrame(df).to_csv(p, index=False); return p

    # ------------------------------------------------------------------------------------------------ targeted IS check
    def targeted(self, is_list, q1_ppm=None, q3_ppm=None):
        from .targeted import run_targeted
        if q1_ppm: self.cfg["targeted"]["q1_ppm"] = q1_ppm; self.cfg["_overrides"].setdefault("targeted", {})["q1_ppm"] = q1_ppm
        if q3_ppm: self.cfg["targeted"]["q4_ppm"] = q3_ppm; self.cfg["_overrides"].setdefault("targeted", {})["q4_ppm"] = q3_ppm
        lst = is_list if isinstance(is_list, str) else self._write_list(is_list, "targeted_list.csv")
        files = list(self.mzml_files().values()) if self.sheet is None else self.files_of("sample", "qc", "blank", "cal", "is_mix", "solvent_std", "matrix_spike")
        sheet = os.path.join(self.out_dir, "sample_sheet.csv") if self.sheet is not None else None
        r = run_targeted(files, lst, os.path.join(self.out_dir, "targeted"), self.cfg, sheet=sheet)
        per_is, per_inj = QA.is_qc(r["assigned"], self.sheet)
        res = dict(assigned=r["assigned"], blank=r["blank"], info=r["info"], xlsx=r["xlsx"], is_qc=per_is, injection_qc=per_inj,
                   list_path=lst)
        self._save("targeted", res)
        return res

    # ------------------------------------------------------------------------------------------------ NTS screening
    def screen(self, matrix="water", insilico="builtin", polarity="auto", types=("sample", "qc", "blank")):
        from .pipeline import Context, process_sample
        from .report import write_all
        from .aggregate import blank_filter, align_features
        self.cfg["insilico"]["engine"] = insilico
        files = self.files_of(*types) if self.sheet is not None else list(self.mzml_files().values())
        out = os.path.join(self.out_dir, "nts"); os.makedirs(out, exist_ok=True)
        ctx = Context(self.cfg)
        results = []
        for f in files:
            try:
                r = process_sample(f, ctx, matrix, polarity, out, None, False)
                r.pop("run", None)                      # keep memory low: runs are re-read from the cache when needed
                results.append(r)
            except Exception:
                log.exception("sample failed: %s", f)
        allf = pd.concat([r["features"] for r in results if not r["features"].empty], ignore_index=True)
        blanks = list(self.sheet.loc[self.sheet.type == "blank", "file"]) if self.sheet is not None else []
        allf = align_features(blank_filter(allf, blanks, self.cfg["blank"]["ratio"]))
        for r in results:
            sub = allf[allf["sample"] == r["sample"]][["feature_id", "consensus_id", "blank_flag", "blank_ratio"]]
            for k in ("features", "top", "hits"):
                if not r[k].empty:
                    r[k] = r[k].drop(columns=[c for c in ("consensus_id", "blank_flag", "blank_ratio") if c in r[k].columns]).merge(sub, on="feature_id", how="left")
        res = write_all(results, out, self.cfg, ctx, params=dict(matrix=matrix, insilico=insilico, files=";".join(os.path.basename(f) for f in files)))
        res["aligned"] = allf
        res["hits_all"] = pd.concat([r["hits"] for r in results if not r["hits"].empty], ignore_index=True, sort=False)
        self._save("screen", {k: v for k, v in res.items() if k in ("top", "compounds", "level_summary", "aligned", "hits_all", "xlsx")})
        return res

    # ------------------------------------------------------------------------------------------------ QA / QC
    def qa(self, rule="2/3", max_qc_rsd=30.0, min_level="3", ppm_tol=None, presence_threshold=70.0):
        """Replicate presence (rule, e.g. '2/3'), QC stability, blank rule, and per identification: peak rating (0-10) and
        presence confidence (%) from the re-integrated Q1 trace + MS2 evidence + replicate frequency + blank ratio."""
        from . import mzml_io
        from .targeted import integrate
        from .quality import peak_rating, presence_confidence
        S = self._load("screen")
        if S is None:
            raise RuntimeError("run screen() first")
        al, top = S["aligned"], S["top"]
        pres = QA.replicate_presence(al, self.sheet, rule) if self.sheet is not None else pd.DataFrame()
        qc = QA.qc_stability(al, self.sheet, max_qc_rsd) if self.sheet is not None else pd.DataFrame()
        ok_lv = LEVEL_ORDER[:LEVEL_ORDER.index(min_level) + 1]
        T = top[top.level.isin(ok_lv)].copy()
        rep_frac = {}
        if len(pres):
            for cid, g in pres.groupby("consensus_id"):
                rep_frac[cid] = float((g.n_detected / g.n_injections).max())
        ppm_tol = ppm_tol or self.cfg["level4"]["ppm"]
        files = self.mzml_files()
        rows = []
        for smp, g in T.groupby("sample"):
            if smp not in files:
                continue
            cfg = config_for_file(self.cfg, files[smp]); run = mzml_io.read_for_config(files[smp], cfg)
            floor = cfg.get("centroid", {}).get("min_height_ms1", 150.0) if run.centroided_here else max(cfg["ms1"]["store_min_intensity"], 1.0)
            for r in g.itertuples():
                rt, I, MZ = run.xic_max(float(r.mz_apex), 10.0)
                it = integrate(rt, I, float(r.rt_apex), floor)
                pr = peak_rating(rt, I, it["start_rt"], it["end_rt"], it["sn"], floor) if np.isfinite(it["start_rt"]) else dict(peak_rating=0.0)
                conf, comp = presence_confidence(pr["peak_rating"], it["sn"], getattr(r, "ppm_error", np.nan), ppm_tol, getattr(r, "iso_score", np.nan),
                                                 getattr(r, "lib_score", np.nan), getattr(r, "lib_matched", 0), getattr(r, "bi_expl_adj_n", np.nan),
                                                 rep_frac.get(getattr(r, "consensus_id", None), np.nan), getattr(r, "blank_ratio", np.nan),
                                                 cfg.get("quality", {}).get("weights"))
                rows.append(dict(sample=smp, feature_id=r.feature_id, int_area=it["area"], int_sn=it["sn"], int_start=it["start_rt"],
                                 int_end=it["end_rt"], **pr, presence_confidence_pct=conf,
                                 presence_pass=conf >= presence_threshold, confidence_components=json.dumps(comp)))
            del run
        Q = pd.DataFrame(rows)
        T = T.merge(Q, on=["sample", "feature_id"], how="left") if len(Q) else T
        res = dict(presence=pres, qc=qc, identifications=T, rule=rule, presence_threshold=presence_threshold)
        self._save("qa", res)
        return res

    # ------------------------------------------------------------------------------------------------ semi-quantification
    def semiquant(self, api_key=None, batch_export=None, ed_max=0.5, min_level="2b", matrix_factor=1.0, unit="ng/mL",
                  method="auto", is_mix=None, is_conc_ug_ml=None, physchem=None, k=None, matrix_correction=True):
        """method: 'is' = IS-surrogate (kNN of the labelled IS response factors, engine 0.6.0, recommended when an
        IS mix was injected and targeted() was run); 'ed' = elucidation distance to calibration standards (0.5.0);
        'auto' = 'is' when targeted results + an IS-mix injection exist, else 'ed'.
        is_mix: file name(s) of the IS-mix injection(s) (default: sample-sheet type 'qc'); is_conc_ug_ml: IS
        concentration in that mix (scalar or dict; default list column 'conc'); physchem: curated property csv
        (logKow + pKa) e.g. from semiquant_eval/08_comptox_physchem_curated.csv."""
        if method == "auto":
            T0 = self._load("targeted")
            has_mix = is_mix is not None or (self.sheet is not None and self.sheet.type.isin(["is_mix", "solvent_std", "qc"]).any())
            method = "is" if (T0 is not None and has_mix) else "ed"
        if method == "is":
            return self._semiquant_is(min_level, matrix_factor, unit, is_mix, is_conc_ug_ml, physchem, k, matrix_correction)
        from .comptox import property_table
        from .semiquant import semiquantify_ed
        Q = self._load("qa"); C = self._load("calibration")
        if Q is None:
            raise RuntimeError("run qa() first")
        if C is None or C["rf"].empty:
            raise RuntimeError("semi-quantification needs calibrants: run calibrate() with calibration injections first")
        ok_lv = LEVEL_ORDER[:LEVEL_ORDER.index(min_level) + 1]
        T = Q["identifications"]
        T = T[T.level.isin(ok_lv) & T.smiles.notna() & (T.smiles.astype(str) != "")].copy()
        T["area_sq"] = T["int_area"].fillna(T.get("area"))
        allc = pd.concat([T[["name", "smiles", "inchikey"]], C["rf"].assign(inchikey="")[["name", "smiles", "inchikey"]]], ignore_index=True)
        props = property_table(allc[allc.smiles.astype(str) != ""], api_key, os.path.join(self.out_dir, "comptox_cache.json"), batch_export)
        out = semiquantify_ed(T.rename(columns={"rt_apex": "rt"}).assign(area=T["area_sq"]), C["rf"], props, ed_max,
                              self.cfg["semiquant"].get("uncertainty_factor", 10))
        for c in ("sq_conc", "sq_low", "sq_high"):
            out[c] = out[c] * matrix_factor
        out["sq_unit"] = unit
        res = dict(table=out, properties=props, method="ed")
        self._save("semiquant", res)
        return res

    def _semiquant_is(self, min_level, matrix_factor, unit, is_mix, is_conc, physchem, k, matrix_correction):
        from . import surrogate as SG
        from .targeted import load_compounds
        Q = self._load("qa"); T = self._load("targeted")
        if Q is None or T is None:
            raise RuntimeError("IS-surrogate semi-quant needs targeted() (IS list incl. the IS-mix injection) and qa()")
        scfg = self.cfg.get("semiquant", {}).get("is_surrogate", {})
        t0 = scfg.get("t0_min", 1.0); vm = scfg.get("void_margin_min", 0.5)
        if is_mix is None:
            is_mix = []
            for t in ("is_mix", "solvent_std", "qc"):          # first sample-sheet type present
                if self.sheet is not None and (self.sheet.type == t).any():
                    is_mix = list(self.sheet.loc[self.sheet.type == t, "file"])
                    break
        is_mix = [os.path.splitext(os.path.basename(f))[0] for f in ([is_mix] if isinstance(is_mix, str) else is_mix)]
        comp = load_compounds(T["list_path"])
        if is_conc is None and "conc" not in comp.columns:
            raw = pd.read_csv(T["list_path"]) if str(T["list_path"]).endswith(".csv") else pd.read_excel(T["list_path"])
            cc = next((c for c in raw.columns if str(c).lower() in ("conc", "conc_ug_ml", "is_conc")), None)
            if cc is None:
                raise ValueError("give is_conc_ug_ml= or a 'conc' (ug/mL) column in the IS list")
            comp = comp.merge(raw[[raw.columns[0], cc]].rename(columns={raw.columns[0]: "name", cc: "conc"}), on="name", how="left")
        rf = SG.is_response_factors(T["assigned"], comp, is_mix, is_conc, t0, vm)
        ok_lv = LEVEL_ORDER[:LEVEL_ORDER.index(min_level) + 1]
        I = Q["identifications"]
        I = I[I.level.isin(ok_lv) & I.smiles.notna() & (I.smiles.astype(str) != "")].copy()
        I["area"] = I["int_area"].fillna(I.get("area")) if "int_area" in I.columns else I.get("area")
        I = I.rename(columns={"rt_apex": "rt"})
        allc = pd.concat([rf[["name", "smiles"]], I[["name", "smiles"]]]).drop_duplicates("name")
        if physchem is None and os.path.exists(os.path.join(self.out_dir, "comptox_properties.csv")):
            physchem = os.path.join(self.out_dir, "comptox_properties.csv")      # from comptox()
        props = SG.physchem_table(allc, physchem, scfg.get("mobile_phase_ph", 2.7))
        # matrix correction: True / "is_recovery" = IS recovery in each sample near the target RT;
        # "me_map" = matrix-effect map from matrix_effect() (spiked control extract vs solvent standard); False/"none" = off
        mode = {True: "is_recovery", False: "none", None: "none"}.get(matrix_correction, matrix_correction)
        rec = SG.is_recovery(T["assigned"], rf, is_mix) if mode == "is_recovery" else None
        out = SG.semiquantify_is(I, rf, props, scfg.get("method", "knn"), k or scfg.get("k", 5), rec,
                                 scfg.get("rt_window_min", 1.0), t0, vm, tuple(scfg.get("uncertainty", (3, 5, 10))),
                                 tuple(scfg.get("ad_cuts", (1.0, 2.0))))
        out["sq_matrix_correction"] = mode
        if mode == "me_map":
            from .matrix import me_at_rt
            ME = self._load("matrix_effect")
            if ME is None:
                raise RuntimeError("matrix_correction='me_map' needs matrix_effect() first")
            me = np.array([me_at_rt(ME["map"], float(t)) for t in out["rt"]])
            out["sq_ME_pct_at_rt"] = me
            f = np.where(np.isfinite(me), 1 + me / 100, 1.0)
            for c in ("sq_conc_ug_ml", "sq_low", "sq_high"):
                out[c] = out[c] / f
        # RT <-> property plausibility (column-specific regression fitted on the IS; weak identification evidence)
        used = rf[rf.use].merge(props, on="name", suffixes=("", "_p"))
        out["rt_pred_from_properties"] = np.nan
        out["rt_property_plausible"] = None
        rtm = None
        if len(used) >= 6:
            rtm = SG.RTPropertyModel(t0).fit(used.rt.to_numpy(), used)
            pp = props.set_index("name").reindex(out["name"]).reset_index()
            okp = pp.logP.notna().to_numpy() & np.isfinite(out["rt"].to_numpy(float))
            if okp.any():
                iv = rtm.predict_interval(pp[okp])
                idx = out.index[okp]
                out.loc[idx, "rt_pred_from_properties"] = iv.rt_pred.to_numpy()
                rr = out.loc[idx, "rt"].to_numpy(float)
                out.loc[idx, "rt_property_plausible"] = (rr >= iv.rt_lo.to_numpy()) & (rr <= iv.rt_hi.to_numpy())
        for c in ("sq_conc_ug_ml", "sq_low", "sq_high"):
            out[c.replace("_ug_ml", "")] = out[c] * matrix_factor
        out["sq_unit"] = unit
        res = dict(table=out, properties=props, is_rf=rf, recovery=rec, method="is",
                   rt_model=None if rtm is None else dict(coef=rtm.coef.tolist(), s=rtm.s, n=rtm.n))
        self._save("semiquant", res)
        return res

    # ------------------------------------------------------------------------------------------------ chromatogram viewer
    def compounds(self, min_level="3"):
        Q = self._load("qa"); S = self._load("screen")
        T = (Q or {}).get("identifications") if Q else (S or {}).get("top")
        if T is None:
            return pd.DataFrame()
        ok = LEVEL_ORDER[:LEVEL_ORDER.index(min_level) + 1]
        return T[T.level.isin(ok)].sort_values(["name", "sample"])

    def chromatogram(self, name, sample, n_q3=6, window=1.0, show=True):
        from . import mzml_io, chromatogram as CH
        from .ms2 import string_to_peaks
        T = self.compounds("5")
        row = T[(T.name == name) & (T["sample"] == sample)]
        files = self.mzml_files()
        if row.empty:
            # compound not identified in this sample: trace it at the RT / m/z where it was found elsewhere
            other = T[T.name == name]
            if other.empty:
                raise KeyError(f"{name} not in the identifications")
            row = other.sort_values("level").head(1)
        r = row.iloc[0]
        cfg = config_for_file(self.cfg, files[sample]); run = mzml_io.read_for_config(files[sample], cfg)
        ions = None
        mzp, itp = string_to_peaks(r.get("ms2_peaks", ""))
        if len(mzp):
            expl = set(str(r.get("bi_expl_peaks", "")).split(";"))
            o = np.argsort(-itp)[:n_q3]
            ions = [dict(mz=float(mzp[k]), label=f"{mzp[k]:.4f}",
                         source=("observed, explained in-silico" if any(e.startswith(f"{mzp[k]:.3f}") for e in expl) else "observed")) for k in o
                    if mzp[k] < float(r.mz_apex) - 1.5]
        tr = CH.extract(run, float(r.mz_apex), float(r.rt_apex), ions or None, cfg, smiles=r.get("smiles"),
                        adduct=r.get("adduct", "[M+H]+") if r.get("adduct") != "direct" else "[M+H]+", n_q3=n_q3, window=window,
                        noise_floor=cfg.get("centroid", {}).get("min_height_ms1", 150.0))
        summ = CH.summary(tr, r.get("lib_score"), r.get("lib_matched", 0), r.get("bi_expl_adj_n"), None, r.get("blank_ratio"), r.get("iso_score"),
                          self.cfg["level4"]["ppm"], cfg.get("quality", {}).get("weights"))
        fig = CH.plot(tr, f"{name} in {sample} - Level {r.level}, {r.get('adduct', '')}")
        if show:
            fig.show()
        return dict(trace=tr, summary=summ, figure=fig)

    # ------------------------------------------------------------------------------------------------ report
    # ------------------------------------------------------------------------------------------------ CompTox batch search
    def comptox(self, auto=True, exports=None, min_level="2b"):
        """Physico-chemical properties from the CompTox Dashboard batch search for every identified compound (>= min_level),
        IS and calibrant. auto=True runs the Dashboard export automatically (no key); exports = uploaded .xlsx file(s)
        from a manual batch search. Result (structure-checked) -> comptox_properties.csv, used by semiquant()."""
        from . import comptox_batch as CB
        comp = []
        Q = self._load("qa"); T = self._load("targeted"); C = self._load("calibration")
        if Q is not None:
            ok = LEVEL_ORDER[:LEVEL_ORDER.index(min_level) + 1]
            comp.append(Q["identifications"].loc[Q["identifications"].level.isin(ok), ["name", "smiles"]])
        if T is not None:
            from .targeted import load_compounds
            comp.append(load_compounds(T["list_path"])[["name", "smiles"]])
        if C is not None and len(C["rf"]):
            comp.append(C["rf"][["name", "smiles"]])
        comp = pd.concat(comp, ignore_index=True) if comp else pd.DataFrame(columns=["name", "smiles"])
        comp = comp[comp.smiles.astype(str).str.len() > 0]
        ids = CB.identifier_list(comp)
        d = os.path.join(self.out_dir, "comptox"); os.makedirs(d, exist_ok=True)
        lst = os.path.join(d, "comptox_batch_identifiers.txt")
        open(lst, "w", encoding="utf-8").write("\n".join(list(ids["name"]) + list(ids.inchikey)))
        paths = list(exports or [])
        if auto and not paths:
            try:
                paths = [CB.fetch_export(list(ids["name"]), "chemical_name", os.path.join(d, "CCD_batch_by_name.xlsx")),
                         CB.fetch_export(list(ids.inchikey), "INCHIKEY", os.path.join(d, "CCD_batch_by_inchikey.xlsx"))]
            except Exception as e:
                raise RuntimeError(f"automatic CompTox export failed ({e}). {CB.MANUAL_STEPS} Identifier list: {lst}")
        tab = CB.curate(paths, comp)
        out = os.path.join(self.out_dir, "comptox_properties.csv")
        tab.to_csv(out, index=False)
        res = dict(table=tab, path=out, identifiers=lst, exports=paths)
        self._save("comptox", res)
        return res

    # ------------------------------------------------------------------------------------------------ templates / matrix effect
    def templates(self, out=None):
        """Write fill-in templates (sample sheet, calibration list with levels L1..L8, IS list) + README."""
        from .templates import write_templates
        return write_templates(out or os.path.join(self.out_dir, "templates"))

    def check_list(self, path_or_df, kind="is"):
        """Check a filled-in IS ('is') or calibration ('cal') list: formula vs monoisotopic mass (5 dp), roles, levels."""
        from .templates import check_compound_list
        return check_compound_list(path_or_df, kind)

    def matrix_effect(self, window_min=1.0):
        """Matrix effect after the run: 'matrix_spike' (mix in control extract) vs 'solvent_std' (same mix in solvent),
        from the targeted results. Saves a table, an RT map and a figure; semiquant() can use the map."""
        from . import matrix as MX
        T = self._load("targeted")
        if T is None:
            raise RuntimeError("run targeted() with an IS / standard list first")
        me = MX.matrix_effect(T["assigned"], self.sheet)
        mp = MX.me_map(me, window_min)
        d = os.path.join(self.out_dir, "matrix_effect"); os.makedirs(d, exist_ok=True)
        me.to_csv(os.path.join(d, "matrix_effect_per_compound.csv"), index=False)
        mp.to_csv(os.path.join(d, "matrix_effect_map.csv"), index=False)
        png = MX.plot_me(me, mp, os.path.join(d, "matrix_effect_map.png"))
        res = dict(table=me, map=mp, png=png)
        self._save("matrix_effect", res)
        return res

    def report(self, path=None):
        """NTS_report.xlsx (all tables, cover sheet) + NTS_report.html (one-file summary with figures). Returns both paths."""
        from .report_project import write_project_report
        from .report_html import write_html_report
        xlsx = write_project_report(self, path or os.path.join(self.out_dir, "NTS_report.xlsx"))
        htm = write_html_report(self, os.path.splitext(xlsx)[0] + ".html", getattr(self, "_report_notes", None))
        return dict(xlsx=xlsx, html=htm)


def _trfp_command(trfp=None):
    """Command prefix for ThermoRawFileParser: Windows .exe, a native binary on PATH, 'dotnet <dll>' or 'mono <exe>'."""
    if trfp and os.path.exists(trfp):
        if trfp.lower().endswith(".dll"):
            return ["dotnet", trfp]
        if trfp.lower().endswith(".exe") and sys.platform != "win32":
            return ["mono", trfp]
        return [trfp]
    for name in ("ThermoRawFileParser", "thermorawfileparser", "ThermoRawFileParser.sh"):
        p = shutil.which(name)
        if p:
            return [p]
    raise FileNotFoundError("ThermoRawFileParser not found: run the Colab setup cell (install_thermorawfileparser) or set "
                            "paths.thermorawfileparser; or convert .raw to mzML elsewhere and upload the mzML files")


def install_thermorawfileparser(dest="/content/trfp"):
    """Colab / Linux: download the latest ThermoRawFileParser release from GitHub and return its path (needs mono or dotnet:
    the setup cell installs mono-runtime via apt when the release is a .NET Framework build)."""
    import io, urllib.request, zipfile
    rel = json.load(urllib.request.urlopen("https://api.github.com/repos/compomics/ThermoRawFileParser/releases/latest"))
    asset = next(a for a in rel["assets"] if a["name"].lower().endswith(".zip") and "osx" not in a["name"].lower() and "win" not in a["name"].lower()) \
        if any("linux" in a["name"].lower() for a in rel["assets"]) else next(a for a in rel["assets"] if a["name"].lower().endswith(".zip"))
    os.makedirs(dest, exist_ok=True)
    zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(asset["browser_download_url"]).read())).extractall(dest)
    for cand in ("ThermoRawFileParser", "ThermoRawFileParser.dll", "ThermoRawFileParser.exe"):
        hits = glob.glob(os.path.join(dest, "**", cand), recursive=True)
        if hits:
            p = hits[0]
            if not p.endswith((".dll", ".exe")):
                os.chmod(p, 0o755)
            log.info("ThermoRawFileParser %s -> %s", rel.get("tag_name"), p)
            return p
    raise FileNotFoundError("ThermoRawFileParser binary not found in the release archive")
