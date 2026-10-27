"""Orchestration + outputs (Excel, CSV, chart) for `pec-screen`."""
from __future__ import annotations
import logging, os
import numpy as np
import pandas as pd

from . import pecscreen as PS
from .common import resolve
from .properties import read_property_tables, resolve_properties
from .fate import FateRunner

log = logging.getLogger("hrms_nts")


def _read_cf(path):
    if not path or not os.path.exists(path): return {}
    from .common import read_table
    t = read_table(path); t.columns = [c.lower().strip() for c in t.columns]
    key = "ik1" if "ik1" in t.columns else "name"
    return dict(zip(t[key].astype(str), t["cf"].astype(float)))


def _predict(cands, F, cfg, props_path, no_fate, backend):
    """Properties -> fate (only for chemicals detected in water) -> PEC vs MEC.  Shared by pec-screen (mzML) and cd-screen (Compound Discoverer tables)."""
    pc = cfg["pec_screen"]
    ppath = props_path or ([resolve(cfg, x.partition("::")[0]) + ("::" + x.partition("::")[2] if "::" in x else "") for x in ([pc["properties"]] if isinstance(pc.get("properties"), str) else (pc.get("properties") or []))] or None)
    table = read_property_tables(ppath, pc.get("properties_sheet")) if ppath else None
    props = resolve_properties(cands, table, pc)
    endo = None
    if pc.get("exclude_endogenous", True):
        from .endogenous import flag_table
        epath = resolve(cfg, pc.get("endogenous_extra"))
        endo = flag_table(cands, epath)[["cand_id", "endogenous", "endogenous_category", "endogenous_confidence"]]
        log.info("endogenous / biogenic metabolites flagged: %d of %d candidates (excluded from the PEC-vs-MEC ranking)", int(endo.endogenous.sum()), len(endo))
    poly = None
    if pc.get("exclude_polymers", True):
        from .polymers import flag_polymers
        poly = flag_polymers(cands)
        log.info("PEG/PPG-type oligomers flagged: %d of %d candidates (excluded)", int(poly.polymer.sum()), len(poly))
    fate = FateRunner(cfg, backend=backend)
    kr = {}
    need_fate = set(F[(F.matrix == "water") & (F.n_det > 0)].cand_id)           # PEC only where the chemical was detected in the irrigation water
    if endo is not None: need_fate -= set(endo[endo.endogenous].cand_id)
    if poly is not None: need_fate -= set(poly[poly.polymer].cand_id)
    l4 = set(cands[cands.level.astype(str) == "4"].cand_id)
    if l4:                                                # Level 4 = formula match, structure = top-ranked database isomer: only worth the (slow) fate model where the solid matrix also shows the signal
        solid = set(F[(F.matrix == "grain") & (F.n_det > 0)].cand_id)
        skipped = (need_fate & l4) - solid; need_fate -= skipped
        log.info("Level 4: fate model only for chemicals detected in water AND the solid matrix (%d of %d Level 4 candidates detected in water skipped)", len(skipped), len(need_fate & l4) + len(skipped))
    todo = [p for _, p in props.iterrows() if p.cand_id in need_fate]
    log.info("fate model for %d of %d candidates (detected in water)", len(todo), len(props))
    for _, p in props.iterrows():
        if p.cand_id not in need_fate:
            kr[p.cand_id] = dict(k_low=np.nan, k_high=np.nan, k_central=np.nan, n_grid=0, k_mode="not needed (no water detection)")
    for n, p in enumerate(todo, 1):
        if no_fate:
            kr[p.cand_id] = dict(k_low=np.nan, k_high=np.nan, k_central=np.nan, n_grid=0, k_mode="skipped (--no-fate)"); continue
        try: kr[p.cand_id] = fate.k_range(p)
        except Exception as e:
            log.error("fate failed for %s: %s", p.cand_id, e); kr[p.cand_id] = dict(k_low=np.nan, k_high=np.nan, k_central=np.nan, n_grid=0, k_mode=f"failed: {e}")
        if n % 10 == 0:
            log.info("fate %d/%d (%d Excel evaluations)", n, len(todo), fate.n_calls); fate.save()
    fate.close()
    cf = _read_cf(resolve(cfg, pc.get("corrections")))
    P, S = PS.predict_and_compare(cands, F, props, kr, cfg, cf)
    if endo is not None:
        S = S.merge(endo, on="cand_id", how="left")
        m = S.endogenous.fillna(False).astype(bool)
        S.loc[m, "class"] = ("endogenous / biogenic (" + S.endogenous_category.astype(str) + ", " + S.endogenous_confidence.astype(str) + ") - excluded")[m]
    if poly is not None:
        S = S.merge(poly, on="cand_id", how="left")
        m = S.polymer.fillna(False).astype(bool)
        S.loc[m, "class"] = ("polymer / contaminant series (" + S.polymer_category.astype(str) + ") - excluded")[m]
    return P, S, props, ppath


def run_pec_screen(cfg: dict, map_path: str, out: str, min_level=None, reuse=True, props_path=None, no_fate=False, backend=None, sink=None) -> dict:
    PS.normalise_config(cfg); pc = cfg["pec_screen"]; min_level = min_level or pc["min_level"]
    os.makedirs(out, exist_ok=True)
    sink = sink or pc.get("sink", "grain"); pc["sink"] = sink; pc["_sink_label"] = "soil" if sink == "soil" else "grain"
    smap = PS.read_sample_map(map_path, cfg)
    smap = smap[smap.matrix.isin(["water", sink])].copy()
    if sink == "soil":                                  # internally the solid matrix keeps the name 'grain'; soil changes factor, PEC quantity and labels
        pc["grain_factor"] = float(pc["soil_factor"]); pc["fruit_to_grain_dw"] = float(pc["soil_pec_to_ng_g"]); pc["_quantity"] = "soil"
        smap["matrix"] = smap.matrix.replace({"soil": "grain"})
    else:
        pc["_quantity"] = "fruit"
    ref = PS.reference_response(cfg)
    log.info("reference %s: area %.4g counts.min at %.3g ug/mL -> %.4g per ug/mL (RT %.2f)", ref["name"], ref["area"], ref["conc_ug_mL"], ref["area_per_ug_mL"], ref["rt"])
    pc["_grain_floor_ng_g"] = PS.grain_floor_ng_g(cfg, ref)
    ident = PS.identify_water(smap, cfg, out, min_level, reuse)
    cands = PS.build_candidates(ident, cfg, min_level)
    if cands.empty: raise SystemExit("no candidates at the requested level in the water samples")
    log.info("candidates (%s or better): %d  %s", min_level, len(cands), cands["level"].value_counts().to_dict())
    removed = pd.DataFrame()
    if pc["id_filter"].get("enabled", True):
        from .idfilter import apply_id_filter
        cands, removed, _ = apply_id_filter(cands, cfg, out)
        if cands.empty: raise SystemExit("no candidates left after the identification filter")
    D, shifts = PS.quantify_injections(cands, smap, cfg, ref)
    D = PS.apply_blank_check(D, float(pc["blank_ratio"]))
    F = PS.aggregate_fields(D)
    P, S, props, ppath = _predict(cands, F, cfg, props_path, no_fate, backend)
    extra = [c for c in ("rt_err_min", "rt_status", "mb_eu_listed", "mb_eu_matched", "mb_eu_ref", "mb_eu_support", "id_filter") if c in cands.columns]
    S = cands[["cand_id", "name", "formula", "level", "polarity", "adduct", "mz", "rt", "confidence", "source", "n_ident_samples", "cas", "ik1", "smiles"] + extra].merge(S, on="cand_id")
    S = S.sort_values(["class", "PEC_over_MEC_geomean"]).reset_index(drop=True)
    res = dict(max_width=pc["max_interval_width"], summary=S, by_field=P, candidates=cands, detections=D, fields=F, properties=props, shifts=shifts, ref=ref, floor=pc["_grain_floor_ng_g"], sink_label=pc["_sink_label"])
    write_outputs(res, cfg, out, map_path, smap, ppath, min_level)
    return res


def run_cd_screen(cfg: dict, workbook: str, out: str, props_path=None, no_fate=False, backend=None, combine_with: str | None = None, sink=None) -> dict:
    """PEC vs MEC from the Compound Discoverer semi-quant tables of the project workbook (water: one max value per chemical; grain: per field)."""
    from .cdimport import read_cd_sheets, build_cd_tables, cd_candidates
    PS.normalise_config(cfg); pc = cfg["pec_screen"]; os.makedirs(out, exist_ok=True)
    sink = sink or pc.get("sink", "grain"); pc["sink"] = sink; pc["_sink_label"] = "soil" if sink == "soil" else "grain"
    if sink == "soil": pc["fruit_to_grain_dw"] = float(pc["soil_pec_to_ng_g"]); pc["_quantity"] = "soil"
    else: pc["_quantity"] = "fruit"
    raw = read_cd_sheets(workbook, os.path.join(resolve(cfg, "work"), "cd_sheets_raw.pkl"))
    water, grain, st = build_cd_tables(raw, cfg, sink)
    cands, F, D = cd_candidates(water, grain, st, cfg)
    log.info("CD tables: %d water rows, %d grain rows -> %d chemicals (%d with structure; %d in water, %d in grain, %d in both)", len(water), len(grain), len(cands),
             int(cands.structure_found.sum()), int(cands.in_water.sum()), int(cands.in_grain.sum()), int((cands.in_water & cands.in_grain).sum()))
    pc["_grain_floor_ng_g"] = float(grain.ng_g.min())
    c2 = cands[cands.smiles != ""].reset_index(drop=True)                      # PEC needs a structure
    P, S, props, ppath = _predict(c2, F[F.cand_id.isin(c2.cand_id)], cfg, props_path, no_fate, backend)
    nos = cands[cands.smiles == ""]
    S = c2[["cand_id", "name", "formula", "level", "polarity", "adduct", "mz", "rt", "confidence", "source", "n_ident_samples", "cas", "ik1", "smiles"]].merge(S, on="cand_id")
    S = S.merge(c2[["cand_id", "rt_water", "rt_grain"]], on="cand_id", how="left")
    # CD pairs water and grain by NAME only: a retention-time mismatch means a different compound (e.g. 17b-estradiol 8.29 min in water, 6.18 min in grain)
    conf = S.n_pairs.gt(0) & (S.rt_water - S.rt_grain).abs().gt(float(pc.get("cd_rt_tol_min", 0.5)))
    S.loc[conf, "class"] = ("identity conflict: RT differs between water and grain (" + S.rt_water.round(2).astype(str) + " vs " + S.rt_grain.round(2).astype(str) + " min)")[conf]
    S = S.sort_values(["class", "PEC_over_MEC_geomean"]).reset_index(drop=True)
    ref = dict(name="Compound Discoverer (project formula)", conc_ug_mL=np.nan, area=np.nan, area_per_ug_mL=np.nan, rt=np.nan, apex=np.nan, source=workbook)
    smap = pd.DataFrame([dict(file=workbook, matrix="water+grain", field="1-9", role="sample", path="", note=f"CD semi-quant sheets; {len(nos)} chemicals without a structure skipped")])
    res = dict(max_width=pc["max_interval_width"], summary=S, by_field=P, candidates=c2, detections=D, fields=F, properties=props, shifts=pd.DataFrame(), ref=ref,
               floor=pc["_grain_floor_ng_g"], no_structure=nos, sink_label=pc["_sink_label"])
    write_outputs(res, cfg, out, workbook, smap, ppath, "CD")
    nos[["name", "in_water", "in_grain"]].to_csv(os.path.join(out, "chemicals_without_structure.csv"), index=False)
    if combine_with:
        combine_and_plot(res, combine_with, out)
    return res


def combine_and_plot(cd_res: dict, tk_dir: str, out: str):
    """One table + chart with the toolkit (mzML) results and the Compound Discoverer results; flags chemicals present in both."""
    tk = pd.read_csv(os.path.join(tk_dir, "pec_screen_summary.csv")); tkp = pd.read_csv(os.path.join(tk_dir, "pec_vs_mec_by_field.csv"))
    cd, cdp = cd_res["summary"].copy(), cd_res["by_field"].copy()
    tk["set"] = "toolkit (mzML)"; cd["set"] = "Compound Discoverer"
    tk["cand_id"] = "TK_" + tk["cand_id"]; tkp["cand_id"] = "TK_" + tkp["cand_id"]; cd["cand_id"] = "CD_" + cd["cand_id"]; cdp["cand_id"] = "CD_" + cdp["cand_id"]
    both = (set(tk.ik1.dropna()) & set(cd.ik1.dropna())) - {""}
    for d in (tk, cd): d["in_both_sets"] = d.ik1.isin(both)
    S = pd.concat([tk, cd], ignore_index=True); P = pd.concat([tkp, cdp], ignore_index=True)
    S.to_csv(os.path.join(out, "combined_summary.csv"), index=False); P.to_csv(os.path.join(out, "combined_by_field.csv"), index=False)
    png = plot_chart(dict(summary=S, by_field=P, max_width=float(cd_res.get("max_width", 100)), sink_label=cd_res.get("sink_label", "grain")), os.path.join(out, "pec_vs_mec_combined_chart.png"), max_rows=40)
    log.info("combined: %d toolkit + %d CD chemicals, %d in both -> %s", len(tk), len(cd), len(both), png)
    return S, P, png


def plot_chart(res: dict, path: str, max_rows=30, title="PEC vs MEC for all identifications (paired fields), ranked by PEC/MEC"):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Patch
    S, P = res["summary"], res["by_field"]
    S = S[~S["class"].str.startswith(("identity conflict", "endogenous", "polymer"))]
    S = S[(S.n_pairs > 0) & ((S.pec_interval_width <= float(res.get("max_width", 100))) | S["class"].str.startswith(("MEC >>", "PEC >>")))].copy()
    if S.empty: return None
    S["_o"] = np.log10(S.PEC_over_MEC_geomean); S = S.sort_values("_o", ascending=False).head(max_rows)
    TEAL = "#1b9aaa"; ORANGE = "#d9822b"
    rows = [(r, P[(P.cand_id == r["cand_id"]) & P.MEC_ng_g.notna() & P.PEC_central.notna()]) for r in S.to_dict("records")]
    fig, ax = plt.subplots(figsize=(12, max(3.5, 0.55 * len(rows) + 2.2)), dpi=150)
    for i, (r, g) in enumerate(rows):
        y = float(i)
        lo, hi = g.PEC_low.min(), g.PEC_high.max()
        col = ORANGE if str(r.get("set", "")).startswith("Compound") else TEAL
        ax.add_patch(Rectangle((lo, y - 0.34), max(hi - lo, lo * 0.02), 0.26, facecolor="none", edgecolor=col, hatch="////", lw=1.1, zorder=2))
        ax.scatter(g.PEC_central, [y - 0.21] * len(g), s=14, color=col, zorder=3)
        ml, mh = g.MEC_ng_g.min(), g.MEC_ng_g.max()
        ax.add_patch(Rectangle((ml, y + 0.04), max(mh - ml, ml * 0.02), 0.26, facecolor="white", edgecolor="#222222", lw=1.1, zorder=2))
        tent = g.grain_id.isin(["MS1 only", "MS2 weak"])
        ax.scatter(g.MEC_ng_g[~tent], [y + 0.17] * int((~tent).sum()), s=16, color="#222222", zorder=4)
        ax.scatter(g.MEC_ng_g[tent], [y + 0.17] * int(tent.sum()), s=22, facecolor="white", edgecolor="#222222", zorder=4)
        f = r["PEC_over_MEC_geomean"]
        lab = (f"{f:.2g}×" if f >= 1 else f"÷{1 / f:.2g}") + f"  (n={r['n_pairs']})"
        ax.text(1.03, y, lab, transform=ax.get_yaxis_transform(), va="center", fontsize=8)
    ax.set_xscale("log"); ax.set_ylim(len(rows) - 0.4, -0.9)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{r['name'][:36]}  [{'CD' if str(r.get('set', '')).startswith('Compound') else r['level']}]" for r, _ in rows], fontsize=8.5)
    lab = res.get("sink_label", "grain")
    ax.set_xlabel("Concentration in %s (ng/g dw, semi-quantified)" % ("maize grain" if lab == "grain" else "soil"), fontsize=10); ax.grid(axis="x", color="#d9d9d9", lw=0.6); ax.set_axisbelow(True)
    for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
    ax.text(1.03, -0.68, "Factor difference\n(geo-mean PEC/MEC)", transform=ax.get_yaxis_transform(), fontsize=9, fontweight="bold", va="center")
    ax.legend(handles=[Patch(facecolor="none", edgecolor=TEAL, hatch="////", label="PEC range over fields - toolkit mzML re-integration (bar = property bounds, dot = central)"),
                       *([Patch(facecolor="none", edgecolor=ORANGE, hatch="////", label="PEC range over fields - Compound Discoverer semi-quant (water = one max value)")] if "set" in S.columns else []),
                       Patch(facecolor="white", edgecolor="#222222", label="MEC range over fields (semi-quant %s; open dot = MS1-only / weak MS2 match)" % lab)],
              loc="upper center", bbox_to_anchor=(0.5, -0.12 if len(rows) > 6 else -0.25), fontsize=8, frameon=False)
    ax.set_title(title, fontsize=11, loc="left")
    fig.subplots_adjust(left=0.27, right=0.84, bottom=0.18, top=0.92)
    fig.savefig(path); plt.close(fig); return path


def write_outputs(res, cfg, out, map_path, smap, props_path, min_level):
    S, P, C, D, F, PR = res["summary"], res["by_field"], res["candidates"], res["detections"], res["fields"], res["properties"]
    names = C.set_index("cand_id")["name"].to_dict()
    for df in (P, D, F, PR):
        if "cand_id" in df.columns and "name" not in df.columns: df.insert(1, "name", df["cand_id"].map(names))
    S.to_csv(os.path.join(out, "pec_screen_summary.csv"), index=False)
    P.to_csv(os.path.join(out, "pec_vs_mec_by_field.csv"), index=False)
    D.to_csv(os.path.join(out, "detections_all_injections.csv"), index=False)
    C.drop(columns=["ms2_peaks"]).to_csv(os.path.join(out, "candidates.csv"), index=False)
    PR.to_csv(os.path.join(out, "properties_used.csv"), index=False)
    need = PR[PR.prop_status.isin(["table-partial", "estimated"])]["cand_id"]
    C[C.cand_id.isin(need)][["name", "cas", "smiles", "inchikey", "formula"]].to_csv(os.path.join(out, "candidates_for_comptox_batch.csv"), index=False)
    if "polymer" in S.columns:
        S[S.polymer.fillna(False).astype(bool)][["cand_id", "name", "level", "polymer_category", "water_ng_L_max", "MEC_ng_g_max"]].to_csv(os.path.join(out, "polymers_flagged.csv"), index=False)
    if "endogenous" in S.columns:
        S[S.endogenous.fillna(False).astype(bool)][["cand_id", "name", "level", "endogenous_category", "endogenous_confidence", "water_ng_L_max", "MEC_ng_g_max"]].to_csv(
            os.path.join(out, "endogenous_flagged.csv"), index=False)
    S["level_group"] = S.level.astype(str).map({"1": "Level 1", "2a": "Level 2", "2b": "Level 2", "3": "Level 3", "4": "Level 4"}).fillna("other")      # 2a + 2b reported together
    res["summary"] = S
    groups = {}
    if S.level_group.nunique() > 1:                                                    # confidence groups are reported separately (Level 2 = 2a + 2b, Level 3 = tentative, Level 4 = formula only)
        for gname, gdf in S.groupby("level_group"):
            slug = gname.lower().replace(" ", "").replace("-", "")
            gdf.to_csv(os.path.join(out, f"pec_screen_summary_{slug}.csv"), index=False)
            groups[gname] = (gdf, plot_chart(dict(res, summary=gdf), os.path.join(out, f"pec_vs_mec_chart_{slug}.png"), title=f"PEC vs MEC - {gname} identifications (paired fields)"))
    gi = S.grain_identity.fillna("").astype(str)                              # best MS2 evidence over the fields where the solid matrix detection was made
    S["grain_evidence"] = np.where(gi == "MS1 only", "MS1 only", np.where(gi.isin(["MS2 confirmed", "MS2 weak"]), "MS2 (confirmed / weak)", ""))
    res["summary"] = S
    lab = cfg["pec_screen"].get("_sink_label", "grain")
    for tag, sel, ttl in (("ms1only", "MS1 only", f"MS1-only {lab} detections (no MS2 in the sample - matrix suppression / low intensity)"), ("ms2", "MS2 (confirmed / weak)", f"{lab} detections with MS2 evidence")):
        gdf = S[S.grain_evidence == sel]
        if len(gdf):
            gdf.to_csv(os.path.join(out, f"pec_screen_summary_{lab}_{tag}.csv"), index=False)
            groups[f"{lab} {tag}"] = (gdf, plot_chart(dict(res, summary=gdf), os.path.join(out, f"pec_vs_mec_chart_{lab}_{tag}.png"), title=f"PEC vs MEC - {ttl}"))
    png = plot_chart(res, os.path.join(out, "pec_vs_mec_chart.png"))
    settings = pd.DataFrame([("min_level", min_level), ("reference", f"{res['ref']['name']} {res['ref']['conc_ug_mL']} ug/mL, area {res['ref']['area']:.4g} counts.min"),
                             ("water_factor (ng/L per ug/mL)", cfg["pec_screen"]["water_factor"]), ("grain_factor (ng/g dw per ug/mL)", cfg["pec_screen"]["grain_factor"]),
                             ("grain detection floor (ng/g dw)", f"{res['floor']:.3g}"), ("blank_ratio", cfg["pec_screen"]["blank_ratio"]),
                             ("fate", f"depth {cfg['pec_screen']['fate']['depth_m']} m, irrigation {cfg['pec_screen']['fate']['irrigation_L_m2_d']} L/m2/d"),
                             ("properties table", "; ".join(props_path) if isinstance(props_path, list) else props_path or "none (structure-based estimates + property grid)"), ("sample map", map_path)], columns=["setting", "value"])
    x = os.path.join(out, "pec_screen_results.xlsx")
    with pd.ExcelWriter(x, engine="openpyxl") as xw:
        S.to_excel(xw, sheet_name="Summary", index=False)
        for gname, (gdf, _) in groups.items(): gdf.to_excel(xw, sheet_name=("Summary_" + gname.replace(" ", ""))[:31], index=False)
        P.to_excel(xw, sheet_name="PEC_vs_MEC_by_field", index=False)
        C.drop(columns=["ms2_peaks"]).to_excel(xw, sheet_name="Candidates", index=False)
        F.to_excel(xw, sheet_name="Fields", index=False)
        D.to_excel(xw, sheet_name="Detections", index=False)
        PR.to_excel(xw, sheet_name="Properties_used", index=False)
        res["shifts"].to_excel(xw, sheet_name="Grain_RT_shifts", index=False)
        smap.drop(columns=["path"]).to_excel(xw, sheet_name="Sample_map", index=False)
        settings.to_excel(xw, sheet_name="Settings", index=False)
    if png:
        from openpyxl import load_workbook
        from openpyxl.drawing.image import Image
        wb = load_workbook(x); ws = wb.create_sheet("Chart", 0); im = Image(png); im.width, im.height = im.width * 0.55, im.height * 0.55; ws.add_image(im, "A1"); wb.save(x)
    print("\n=== PEC-SCREEN SUMMARY (classes) ===")
    print(S.groupby("class").size().to_string())
    if groups:
        for gname, (gdf, _) in groups.items():
            print(f"\n--- {gname}: {len(gdf)} chemicals ---"); print(gdf["class"].str.replace(r" \(.*", "", regex=True).str.replace(r" \[.*", "", regex=True).str[:40].value_counts().to_string())
    print(f"\nresults: {out}\n  Excel: {x}" + (f"\n  chart: {png}" if png else ""))
