"""End-to-end sample processing: mzML -> features -> Level 4 -> RT/MS2/library/in-silico -> levels -> semi-quant."""
from __future__ import annotations
import os, glob, re, logging, time
from functools import lru_cache
import numpy as np
import pandas as pd

from . import common, mzml_io, features as F, level4 as L4, ms2 as M2, insilico as IS, levels as LV, semiquant as SQ
from .common import resolve, adduct_table, guess_polarity, ppm_err
from .database import build_database
from .rt_prediction import RTModel, known_standard_rt

log = logging.getLogger("hrms_nts")


class Context:
    """Heavy shared objects (database, library, RT models, standards) loaded once per run."""

    def __init__(self, cfg, need_library=True):
        self.cfg = cfg
        self.db = build_database(cfg)
        self.lib = None
        lp = resolve(cfg, cfg["paths"].get("library_db"))
        if need_library and lp and os.path.exists(lp):
            from .library import SpectralLibrary
            self.lib = SpectralLibrary(lp); log.info("spectral library: %d spectra", self.lib.n)
        elif need_library:
            log.warning("no spectral library at %s - Level 2a/1 unavailable until you run: build-library", lp)
        sp = resolve(cfg, cfg["paths"].get("standards"))
        self.standards = common.read_table(sp) if sp and os.path.exists(sp) else None
        cp = resolve(cfg, cfg["paths"].get("calibrants"))
        self.calibrants = SQ.calibrants_table(cp) if cp and os.path.exists(cp) else None
        self.rt_models = {}
        for pol in ("positive", "negative"):
            self.rt_models[pol] = self._rt_model(pol)

    def _rt_model(self, pol):
        cfg = self.cfg; rc = cfg["rt"]
        if not rc.get("enabled", True) or self.standards is None:
            return None
        col = "rt_pos" if pol == "positive" else "rt_neg"
        if col not in self.standards.columns:
            return None
        tr = self.standards.rename(columns={col: "rt"}).dropna(subset=["rt"])
        if len(tr) < 6:
            log.warning("RT model (%s): only %d standards - skipped", pol, len(tr)); return None
        m = RTModel(pol, tol=rc["tolerance_min"], conflict_factor=rc["conflict_factor"], adaptive_k=rc.get("adaptive_k", 0.0))
        return m.fit(tr, self.db.df)


@lru_cache(maxsize=200000)
def _bi_cached(smiles, adduct, mz_key):
    return None  # placeholder to keep API stable (spectra differ per feature; fragment sets are cached inside insilico)


def sample_polarity(path, cfg, run, forced=None):
    if forced and forced != "auto":
        return forced
    return guess_polarity(path, cfg.get("polarity_patterns")) or run.polarity or "positive"


def process_sample(path, ctx: Context, matrix="water", polarity="auto", out_dir=None, engine=None, semiquant=True) -> dict:
    cfg = ctx.cfg; t0 = time.time()
    name0 = os.path.splitext(os.path.basename(path))[0]
    run = mzml_io.read_mzml(path, None, store_min_intensity=cfg["ms1"]["store_min_intensity"],
                            cache_dir=resolve(cfg, cfg["paths"].get("cache_dir")))
    run.polarity = sample_polarity(path, cfg, run, polarity)
    pol = run.polarity
    m1 = cfg["ms1"]
    feat = F.detect_features(run, m1["feature_ppm"], m1["rt_gap_min"], m1["min_scans"], m1["min_intensity"],
                             tuple(m1["rt_range"]) if m1.get("rt_range") else None,
                             m1["split_broad_min"], m1["split_prominence"], m1["flag_width_min"],
                             m1.get("bg_occupancy", 0.30), m1.get("bg_peak_to_median", 15.0))
    if feat.empty:
        return dict(sample=name0, polarity=pol, features=feat, hits=pd.DataFrame(), top=pd.DataFrame())
    feat = F.annotate_relations(feat, m1["feature_ppm"] * 0.7, cfg["ms1"]["relation_rt_tol"])
    m2c = cfg["ms2"]
    feat = M2.assign_ms2_to_features(feat, run, m2c["precursor_tol_da"], m2c["rt_pad_min"], m2c["merge_scans"], m2c["rel_min"],
                                     m2c["top_n"], m2c["good_int"])
    ad = adduct_table(pol, cfg["level4"].get("adducts_" + pol))
    hits = L4.match_features(feat, ctx.db, ad, cfg["level4"]["ppm"], run, cfg["ms1"]["store_min_intensity"],
                             cfg["level4"]["drop_isotope_peaks"], cfg["level4"]["min_iso_score"], cfg["level4"]["use_isotopes"])
    if not hits.empty:
        hits = hits.merge(feat[[c for c in feat.columns if c.startswith("ms2_")] + ["feature_id"]], on="feature_id", how="left",
                          suffixes=("", "_dup"))
        hits = hits.loc[:, ~hits.columns.str.endswith("_dup")]
        # ---- RT --------------------------------------------------------------------------------------------
        rm = ctx.rt_models.get(pol)
        if rm is not None:
            hits = rm.annotate(hits)
        else:
            for c in ("rt_pred_qsar", "rt_pred_rti", "rt_diff_qsar", "rt_diff_rti"): hits[c] = np.nan
            hits["rt_available"] = False; hits["rt_agree"] = False; hits["rt_conflict"] = False
        hits = known_standard_rt(hits, ctx.standards, pol)
        # ---- limit candidate load for MS2-heavy steps -------------------------------------------------------
        # candidate priority for the expensive steps: RT-consistent, curated sources first, then smallest mass error
        hits["_ppm_abs"] = hits["ppm_error"].abs()
        hits["_prio"] = (hits["rt_agree"].astype(int) * 3 + hits["source"].astype(str).str.contains("user|NORMAN", regex=True).astype(int) * 2
                         + (hits["std_rt_diff"].abs() <= cfg["levels"]["level1_rt_tol"]).astype(int) * 5 - hits["rt_conflict"].astype(int) * 3)
        hits = hits.sort_values(["feature_id", "_prio", "_ppm_abs"], ascending=[True, False, True]).reset_index(drop=True)
        hits["cand_rank"] = hits.groupby("feature_id").cumcount()
        # ---- library (per feature) ---------------------------------------------------------------------------
        hits = _attach_library(hits, feat, ctx, pol)
        # ---- in-silico ----------------------------------------------------------------------------------------
        hits = _attach_insilico(hits, feat, ctx, pol, run, engine or cfg["insilico"]["engine"], out_dir)
        hits = LV.mark_clear_top(hits, cfg["levels"]["clear_top_margin"], cfg["levels"]["insilico_min_score"])
        hits = LV.assign_levels(hits, cfg)
    # ---- features with no candidate: library-only + Level 5 -----------------------------------------------------
    matched_ids = set(hits["feature_id"]) if not hits.empty else set()
    unk = feat[~feat["feature_id"].isin(matched_ids) & (feat["iso_of"].astype(str) == "")].copy()
    unk["level"] = "5"; unk["level_basis"] = "exact mass only"
    hits_all = pd.concat([hits, unk[[c for c in unk.columns if c not in hits.columns or True]]], ignore_index=True, sort=False) if not unk.empty else hits
    top = hits_all[hits_all.get("is_top_candidate", pd.Series(False, index=hits_all.index)).fillna(False).astype(bool) | (hits_all["level"] == "5")].copy()
    if semiquant and cfg["semiquant"].get("mode", "single") != "none" and not top.empty:
        try:
            top = SQ.semiquantify(top, cfg, matrix, ctx.calibrants, pol, cfg["semiquant"].get("area_column", "area"))
        except Exception as e:
            log.error("semi-quantification failed: %s", e)
    log.info("%s: %d features, %d candidate rows; top levels: %s (%.0fs)", name0, len(feat), len(hits),
             top["level"].value_counts().to_dict() if not top.empty else {}, time.time() - t0)
    return dict(sample=name0, polarity=pol, run=run, features=feat, hits=hits, top=top)


def _attach_library(hits, feat, ctx, pol):
    cfg = ctx.cfg; lc = cfg["library"]
    for c, d in (("lib_score", 0.0), ("lib_matched", 0), ("lib_name", ""), ("lib_source", ""), ("lib_adduct", ""), ("lib_inhouse", False)):
        hits[c] = d
    if ctx.lib is None or hits.empty:
        return hits
    fmap = feat.set_index("feature_id")
    sel = hits[(hits["ms2_n_peaks"] >= cfg["levels"]["min_ms2_peaks"])]["feature_id"].unique()
    best = {}
    for fid in sel:
        f = fmap.loc[fid]
        mzp, itp = M2.string_to_peaks(f["ms2_peaks"])
        res = ctx.lib.search(float(f["mz_apex"]), pol, mzp, itp, lc["precursor_ppm"], lc["frag_tol_da"], lc["frag_tol_ppm"],
                             lc["min_report_score"], lc["remove_precursor"])
        if not res.empty:
            best[fid] = res
    if best:
        L = pd.concat([r.assign(feature_id=fid) for fid, r in best.items()], ignore_index=True)
        L = L[["feature_id", "ik1", "score", "matched", "lib_name", "lib_source", "lib_adduct"]].rename(
            columns={"score": "_ls", "matched": "_lm", "lib_name": "_ln", "lib_source": "_lsrc", "lib_adduct": "_la"})
        hits = hits.reset_index(drop=True); hits["_hid"] = np.arange(len(hits))
        J = hits[["_hid", "feature_id", "ik1", "adduct"]].merge(L, on=["feature_id", "ik1"], how="inner")
        # library precursor type must agree with the candidate adduct (when the library states one)
        J = J[(J["_la"].astype(str) == "") | (J["_la"].astype(str) == J["adduct"].astype(str)) | (J["adduct"] == "direct")]
        J = J.sort_values(["_ls", "_lm"], ascending=False).drop_duplicates("_hid")
        hits = hits.merge(J[["_hid", "_ls", "_lm", "_ln", "_lsrc", "_la"]], on="_hid", how="left").drop(columns="_hid")
        m = hits["_ls"].notna()
        hits.loc[m, "lib_score"] = hits.loc[m, "_ls"]; hits.loc[m, "lib_matched"] = hits.loc[m, "_lm"]
        hits.loc[m, "lib_name"] = hits.loc[m, "_ln"]; hits.loc[m, "lib_source"] = hits.loc[m, "_lsrc"]; hits.loc[m, "lib_adduct"] = hits.loc[m, "_la"]
        hits["lib_inhouse"] = hits["lib_source"].astype(str).str.lower().str.startswith("inhouse")
        hits = hits.drop(columns=["_ls", "_lm", "_ln", "_lsrc", "_la"])
    # library-only discoveries (structure not in suspect lists)
    if lc.get("discovery", True):
        extra = []
        for fid, res in best.items():
            f = fmap.loc[fid]
            known = set(hits.loc[hits["feature_id"] == fid, "ik1"])
            for _, b in res.iterrows():
                if (b["ik1"] and b["ik1"] not in known and b["score"] >= cfg["levels"]["lib_cosine"] and b["matched"] >= cfg["levels"]["lib_min_matched"]
                    and abs(float(ppm_err(f["mz_apex"], b["lib_precursor_mz"]))) <= cfg["level4"]["ppm"]):
                    row = {c: f[c] for c in feat.columns if c in f.index}
                    row.update(db_id=f"LIB{int(b['lib_id'])}", name=b["lib_name"], formula=b["lib_formula"], smiles=b["lib_smiles"],
                               inchikey=b["lib_inchikey"], ik1=b["ik1"], exact_mass=np.nan, source="library-only", adduct=b["lib_adduct"],
                               calc_mz=b["lib_precursor_mz"], ppm_error=float(ppm_err(f["mz_apex"], b["lib_precursor_mz"])),
                               lib_score=b["score"], lib_matched=b["matched"], lib_name=b["lib_name"], lib_source=b["lib_source"],
                               lib_adduct=b["lib_adduct"], lib_inhouse=str(b["lib_source"]).lower().startswith("inhouse"),
                               n_adducts_seen=1, n_candidates_for_feature=1, cand_rank=99)
                    extra.append(row); known.add(b["ik1"])
        if extra:
            ex = pd.DataFrame(extra)
            for c in hits.columns:
                if c not in ex.columns: ex[c] = np.nan
            hits = pd.concat([hits, ex[hits.columns]], ignore_index=True)
            if ctx.standards is not None:
                hits = known_standard_rt(hits.drop(columns=["std_rt", "std_rt_diff"], errors="ignore"), ctx.standards, pol)
            log.info("library-only identifications added: %d", len(extra))
    return hits


def _attach_insilico(hits, feat, ctx, pol, run, engine, out_dir):
    cfg = ctx.cfg; ic = cfg["insilico"]
    for c, d in (("bi_expl_n", 0), ("bi_expl_adj_n", 0.0), ("bi_expl_frac", 0.0), ("bi_decoy_frac", 0.0), ("bi_score", 0.0), ("bi_expl_peaks", ""),
                 ("mf_score", np.nan), ("mf_frag_score", np.nan), ("mf_expl", np.nan), ("mf_used", np.nan)):
        hits[c] = d
    if engine == "none" or hits.empty:
        return hits
    fmap = feat.set_index("feature_id")
    ms2_feats = hits[(hits["ms2_n_peaks"] >= cfg["levels"]["min_ms2_peaks"])]
    done = 0
    for fid, g in ms2_feats.groupby("feature_id"):
        f = fmap.loc[fid]
        mzp, itp = M2.string_to_peaks(f["ms2_peaks"])
        g = g[g["cand_rank"] < ic["max_candidates_per_feature"]]
        if engine in ("builtin", "both"):
            for i, r in g.iterrows():
                if not isinstance(r["smiles"], str) or not r["smiles"]:
                    continue
                sc = IS.score_candidate(r["smiles"], r["adduct"] if r["adduct"] != "direct" else ("[M+H]+" if pol == "positive" else "[M-H]-"),
                                        mzp, itp, float(f["mz_apex"]), ic["frag_tol_da"], ic["frag_tol_ppm"], ic["rel_min"], ic["max_depth"])
                for k, v in sc.items():
                    hits.at[i, k] = v
        if engine in ("metfrag", "both"):
            ok, jar, java = IS.metfrag_available(cfg)
            if ok and len(g) >= 1 and f["max_int"] >= ic.get("metfrag_min_intensity", 0):
                work = os.path.join(out_dir or ".", "metfrag_work")
                adduct = g["adduct"].iloc[0]
                res = IS.run_metfrag(g, mzp, itp, float(f["mz_apex"]), pol, adduct, work, jar, java, ic["frag_tol_ppm"], ic["frag_tol_da"],
                                     ic["max_depth"], resolve(cfg, cfg["paths"].get("metfrag_mb_folder")), tag=re.sub(r"\W", "_", fid))
                if not res.empty:
                    rr = res.set_index("db_id")
                    for i, r in g.iterrows():
                        if str(r["db_id"]) in rr.index:
                            hits.at[i, "mf_score"] = rr.at[str(r["db_id"]), "mf_score"]
                            hits.at[i, "mf_frag_score"] = rr.at[str(r["db_id"]), "mf_frag_score"]
                            hits.at[i, "mf_expl"] = rr.at[str(r["db_id"]), "mf_expl"]
                            hits.at[i, "mf_used"] = rr.at[str(r["db_id"]), "mf_used"]
            elif not ok and done == 0:
                log.warning("MetFrag/Java not available - using built-in predictor only"); done = 1
    return hits
