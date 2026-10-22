"""Outputs: CSV tables, a multi-sheet Excel workbook, and MS2 mirror plots for identified candidates."""
from __future__ import annotations
import os, json, logging, datetime
import numpy as np
import pandas as pd
from . import ms2 as M2

log = logging.getLogger("hrms_nts")

LEVEL_ORDER = ["1", "2a", "2b", "3", "4", "5"]


def _clean(df: pd.DataFrame, drop=("apex_scan",)):
    d = df.drop(columns=[c for c in drop if c in df.columns]).copy()
    for c in d.columns:
        if d[c].dtype == object:
            d[c] = d[c].astype(str).str.slice(0, 32000)
    return d


def level_summary(top: pd.DataFrame) -> pd.DataFrame:
    if top.empty:
        return pd.DataFrame(columns=["sample", "level", "n_features"])
    return (top.groupby(["sample", "level"]).size().unstack("level", fill_value=0)
              .reindex(columns=LEVEL_ORDER, fill_value=0).reset_index())


def mirror_plot(path, title, q_mz, q_int, l_mz=None, l_int=None, pred_mz=None, precursor=None):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.8))
    if len(q_mz):
        ax.vlines(q_mz, 0, q_int / q_int.max() * 100, color="#1b6ca8", lw=1.2, label="experimental MS2")
    if l_mz is not None and len(l_mz):
        ax.vlines(l_mz, 0, -l_int / l_int.max() * 100, color="#c0392b", lw=1.2, label="library")
    if pred_mz is not None and len(pred_mz):
        ax.plot(pred_mz, np.full(len(pred_mz), -105), "|", color="#7f8c8d", ms=8, label="predicted ions")
    if precursor:
        ax.axvline(precursor, color="k", ls=":", lw=0.8)
    ax.axhline(0, color="k", lw=0.6); ax.set_xlabel("m/z"); ax.set_ylabel("relative intensity (%)")
    ax.set_title(title, fontsize=9); ax.legend(fontsize=7, loc="upper left"); ax.set_ylim(-115 if (l_mz is not None or pred_mz is not None) else 0, 110)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def write_all(results: list[dict], out_dir: str, cfg: dict, ctx=None, params: dict | None = None, plots_max=60):
    os.makedirs(out_dir, exist_ok=True)
    feats = pd.concat([r["features"] for r in results if not r["features"].empty], ignore_index=True) if results else pd.DataFrame()
    hits = pd.concat([r["hits"] for r in results if not r["hits"].empty], ignore_index=True, sort=False) if any(not r["hits"].empty for r in results) else pd.DataFrame()
    tops = [r["top"] for r in results if not r["top"].empty]
    top = pd.concat(tops, ignore_index=True, sort=False) if tops else pd.DataFrame()
    if not top.empty:
        top["level_rank"] = top["level"].map({k: i for i, k in enumerate(LEVEL_ORDER)})
        top = top.sort_values(["sample", "level_rank", "max_int"], ascending=[True, True, False])
    feats.to_csv(os.path.join(out_dir, "features.csv"), index=False)
    if not hits.empty:
        hits.drop(columns=["_ppm_abs"], errors="ignore").to_csv(os.path.join(out_dir, "all_candidates.csv"), index=False)
    top.to_csv(os.path.join(out_dir, "top_candidates_per_feature.csv"), index=False)
    ident = top[top["level"].isin(["1", "2a", "2b"])] if not top.empty else top
    ident.to_csv(os.path.join(out_dir, "identified_level1_2b.csv"), index=False)
    lsum = level_summary(top)
    lsum.to_csv(os.path.join(out_dir, "level_summary.csv"), index=False)

    # compound-level table across samples (best level per compound), the table you usually want for a paper
    comp = pd.DataFrame()
    if not top.empty:
        t = top[top["level"] != "5"].copy()
        conc_cols = [c for c in t.columns if c.startswith("conc_") and not c.endswith(("_low", "_high")) and c != "conc_ug_mL_extract"]
        agg = {"level": "min", "confidence_score": "max", "area": "max", "max_int": "max", "rt_apex": "median", "ppm_error": "median", "sample": "nunique"}
        for c in conc_cols: agg[c] = "median"
        t["level_num"] = t["level"].map({k: i for i, k in enumerate(LEVEL_ORDER)})
        g = t.sort_values("level_num").groupby(["name", "ik1"], as_index=False)
        comp = g.agg({**{k: v for k, v in agg.items() if k != "level"}, "level_num": "min", "formula": "first", "smiles": "first"})
        comp["best_level"] = comp["level_num"].map(dict(enumerate(LEVEL_ORDER)))
        comp = comp.rename(columns={"sample": "n_samples_detected"}).drop(columns="level_num").sort_values(["best_level", "max_int"], ascending=[True, False])
        comp.to_csv(os.path.join(out_dir, "compound_summary.csv"), index=False)

    xlsx = os.path.join(out_dir, "NTS_results.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        pd.DataFrame({"parameter": list((params or {}).keys()), "value": [str(v) for v in (params or {}).values()]}).to_excel(xw, sheet_name="Run_info", index=False)
        lsum.to_excel(xw, sheet_name="Level_summary", index=False)
        if not comp.empty: _clean(comp).to_excel(xw, sheet_name="Compounds", index=False)
        if not ident.empty: _clean(ident).to_excel(xw, sheet_name="Identified_L1_2b", index=False)
        if not top.empty: _clean(top).head(200000).to_excel(xw, sheet_name="Top_per_feature", index=False)
        if not hits.empty:
            h5 = hits[hits["cand_rank"] < 5] if "cand_rank" in hits.columns else hits
            _clean(h5.drop(columns=["_ppm_abs"], errors="ignore")).head(500000).to_excel(xw, sheet_name="Candidates_top5", index=False)
    log.info("wrote %s", xlsx)

    # mirror plots for confirmed levels
    if ctx is not None and not ident.empty:
        pdir = os.path.join(out_dir, "plots"); os.makedirs(pdir, exist_ok=True)
        from . import insilico as IS
        n = 0
        for r in ident.sort_values("confidence_score", ascending=False).itertuples():
            if n >= plots_max: break
            q_mz, q_int = M2.string_to_peaks(getattr(r, "ms2_peaks", ""))
            if len(q_mz) == 0: continue
            lm = li = None
            if ctx.lib is not None and r.lib_score and r.lib_score > 0:
                try:
                    res = ctx.lib.search(float(r.mz_apex), r.polarity, q_mz, q_int, cfg["library"]["precursor_ppm"], cfg["library"]["frag_tol_da"],
                                         cfg["library"]["frag_tol_ppm"], 0.3, cfg["library"]["remove_precursor"])
                    m = res[res["ik1"] == r.ik1]
                    if not m.empty: lm, li = ctx.lib.spectrum(m.iloc[0]["lib_id"])
                except Exception:
                    pass
            pm = None
            if r.level in ("2b",) and isinstance(r.smiles, str):
                pm = IS.predict_product_ions(r.smiles, r.adduct if r.adduct != "direct" else "[M+H]+", top=60)["mz"].values
            fn = f"{r.sample}_{r.feature_id.split('_')[-1]}_{r.level}_{str(r.name)[:30]}".replace("/", "-").replace(" ", "_").replace(":", "")
            try:
                mirror_plot(os.path.join(pdir, fn + ".png"), f"{r.name}  L{r.level}  {r.adduct}  RT {r.rt_apex:.2f}  conf {r.confidence_score}",
                            q_mz, q_int, lm, li, pm, float(r.mz_apex))
                n += 1
            except Exception as e:
                log.debug("plot failed: %s", e)
    return dict(features=feats, hits=hits, top=top, compounds=comp, level_summary=lsum, xlsx=xlsx)
