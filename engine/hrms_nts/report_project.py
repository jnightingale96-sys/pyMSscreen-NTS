"""NTS_report.xlsx for a workflow.Project: identified compounds with adducts, RT, product-ion evidence (library / in-silico
predicted / MetFrag), peak rating, presence confidence, replicate presence, QC, blanks, calibration, quantification, semi-quant."""
from __future__ import annotations
import os, logging
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")
LEVEL_ORDER = ["1", "2a", "2b", "3", "4", "5"]
GREEN, AMBER, RED, BLUE = "C6EFCE", "FFEB9C", "FFC7CE", "DDEBF7"


def _adducts(hits, T, rt_tol=0.05):
    """Adducts of the same structure (InChIKey block 1) co-eluting in the same sample -> 'adduct (ppm)' list."""
    if hits is None or hits.empty:
        return pd.Series("", index=T.index)
    H = hits[["sample", "ik1", "adduct", "rt_apex", "ppm_error"]].dropna(subset=["ik1"])
    out = []
    for r in T.itertuples():
        h = H[(H["sample"] == r.sample) & (H.ik1 == r.ik1) & ((H.rt_apex - r.rt_apex).abs() <= rt_tol)]
        a = h.sort_values("ppm_error", key=abs).drop_duplicates("adduct")
        out.append(", ".join(f"{x.adduct} ({x.ppm_error:+.1f} ppm)" for x in a.itertuples()))
    return pd.Series(out, index=T.index)


def _ms2_evidence(r):
    ev = []
    if (r.get("lib_matched") or 0) >= 2 and (r.get("lib_score") or 0) > 0:
        ev.append(f"library {r['lib_score']:.2f} / {int(r['lib_matched'])} ions ({r.get('lib_name', '')}; {r.get('lib_source', '')})")
    if (r.get("bi_expl_adj_n") or 0) > 0:
        ev.append(f"in-silico predicted: {r['bi_expl_adj_n']:.1f} ions explained (score {r.get('bi_score', 0):.2f})")
    if pd.notna(r.get("mf_score")) and r.get("mf_score"):
        ev.append(f"MetFrag {r['mf_score']:.2f} ({int(r.get('mf_expl') or 0)} ions)")
    return "; ".join(ev) if ev else ("MS2 present, not explained" if (r.get("ms2_n_peaks") or 0) else "no MS2")


def write_project_report(P, path):
    S = P._load("screen"); Q = P._load("qa"); C = P._load("calibration"); TG = P._load("targeted"); SQ = P._load("semiquant")
    sheets = {}
    notes = [("Data", P.data_dir), ("Results", P.out_dir)]
    if P.sheet is not None:
        sheets["Samples"] = P.sheet
        notes.append(("Injections", ", ".join(f"{k}: {v}" for k, v in P.sheet.type.value_counts().items())))
    T = None
    if Q is not None:
        T = Q["identifications"].copy()
        notes.append(("Replicate rule", f"present in a sample group when detected in >= {Q['rule']} of its replicate injections (blank-flagged features excluded)"))
        notes.append(("Presence confidence", f"pass >= {Q['presence_threshold']:g} % (weighted evidence: peak rating, S/N, Q1 accuracy, isotopes, MS2, replicates, blank ratio)"))
    elif S is not None:
        T = S["top"][S["top"].level.isin(LEVEL_ORDER[:4])].copy()
    if T is not None and len(T):
        T["adducts_seen"] = _adducts(S["hits_all"] if S else None, T)
        T["ms2_evidence"] = [_ms2_evidence(r) for _, r in T.iterrows()]
        T["product_ions_predicted"] = [" + ".join(p for p, ok in (("built-in", (r.get("bi_expl_adj_n") or 0) > 0),
                                                                  ("MetFrag", pd.notna(r.get("mf_score")) and (r.get("mf_score") or 0) > 0)) if ok)
                                       for _, r in T.iterrows()]
        T["lv"] = T.level.map({k: i for i, k in enumerate(LEVEL_ORDER)})
        pres = Q["presence"] if Q is not None else pd.DataFrame()
        pres_txt = {}
        if len(pres):
            for cid, g in pres.groupby("consensus_id"):
                pres_txt[cid] = "; ".join(f"{x.group} {int(x.n_detected)}/{int(x.n_injections)}{' ✓' if x.present else ''}" for x in g.itertuples())
        qcp = Q["qc"].set_index("consensus_id") if Q is not None and len(Q["qc"]) else None
        cols = dict(name="Compound", formula="Formula", smiles="SMILES", inchikey="InChIKey", level="Level", level_basis="Level basis",
                    adducts_seen="Adducts (co-eluting)", rt_apex="RT (min)", ppm_error="Q1 error (ppm)", rt_pred_qsar="RT predicted QSAR",
                    rt_pred_rti="RT predicted RTI", rt_pred_global="RT predicted (SMRT transfer)", rt_agree="RT agrees", std_rt_diff="RT diff to standard", ms2_evidence="Product-ion evidence",
                    product_ions_predicted="Product ions predicted", bi_expl_peaks="Explained product ions", iso_score="Isotope score",
                    peak_rating="Peak rating (0-10)", presence_confidence_pct="Presence confidence (%)", int_area="Peak area", int_sn="S/N",
                    blank_flag="In blank", blank_ratio="Sample/blank ratio", origin_hint="Origin (reference databases)",
                    norman_lists="NORMAN SLE lists", known_tp="Known transformation product", tp_of="TP of (parent, literature)",
                    np_sources="Natural-product / lipid databases", chebi_roles="ChEBI roles", lipid_class="Lipid class",
                    sample="Sample", feature_id="Feature", consensus_id="Consensus feature")
        det = T[[c for c in cols if c in T.columns]].rename(columns=cols)
        det["Present (replicates)"] = T.get("consensus_id", pd.Series(index=T.index, dtype=object)).map(pres_txt).fillna("")
        if qcp is not None:
            det["QC pass"] = T.consensus_id.map(qcp["qc_pass"])
        sheets["Detections_per_sample"] = det.sort_values(["Compound", "Sample"])
        # one row per compound: best level, best evidence row
        best = T.sort_values(["lv", "presence_confidence_pct" if "presence_confidence_pct" in T.columns else "max_int"], ascending=[True, False])
        g = best.groupby(["name", "ik1"], dropna=False)
        comp = g.head(1).copy()
        comp["n_samples"] = comp.set_index(["name", "ik1"]).index.map(g["sample"].nunique()).values
        comp["samples"] = comp.set_index(["name", "ik1"]).index.map(g["sample"].apply(lambda s: ", ".join(sorted(set(s))))).values
        if "presence_confidence_pct" in T.columns:
            comp["max_conf"] = comp.set_index(["name", "ik1"]).index.map(g["presence_confidence_pct"].max()).values
        cc = comp[[c for c in cols if c in comp.columns and c not in ("sample", "feature_id")]].rename(columns=cols)
        cc.insert(cc.columns.get_loc("Level") + 1, "Detected in (n samples)", comp["n_samples"].values)
        cc["Samples"] = comp["samples"].values
        cc["Present (replicates)"] = comp.get("consensus_id", pd.Series(index=comp.index, dtype=object)).map(pres_txt).fillna("").values
        if qcp is not None:
            cc["QC pass"] = comp.consensus_id.map(qcp["qc_pass"]).values
        if SQ is not None:
            st = SQ["table"]
            agg = dict(conc=("sq_conc", "median"), low=("sq_low", "median"), high=("sq_high", "median"))
            if SQ.get("method") == "is":
                agg.update(nn=("sq_nearest_IS", "first"), dist=("sq_nn_distance", "median"), flags=("sq_flags", "first"))
            else:
                agg.update(ed=("elucidation_distance", "median"), cal=("sq_calibrant", "first"), ok=("sq_ed_ok", "first"))
            sq = st.groupby("name").agg(**{k: v for k, v in agg.items() if v[0] in st.columns})
            unit = st["sq_unit"].iloc[0] if len(st) else ""
            cc[f"Semi-quant median ({unit})"] = comp.name.map(sq["conc"]).values
            cc["Semi-quant range (uncertainty)"] = [f"{a:.3g} - {b:.3g}" if pd.notna(a) else "" for a, b in zip(comp.name.map(sq["low"]), comp.name.map(sq["high"]))]
            if SQ.get("method") == "is":
                cc["Nearest internal standard"] = comp.name.map(sq.get("nn", pd.Series(dtype=object))).values
                cc["Distance to nearest IS"] = comp.name.map(sq.get("dist", pd.Series(dtype=float))).values
                cc["Semi-quant flags"] = comp.name.map(sq.get("flags", pd.Series(dtype=object))).values
            else:
                cc["Elucidation distance"] = comp.name.map(sq["ed"]).values
                cc["ED < 0.5 (close calibrant)"] = comp.name.map(sq["ok"]).values
                cc["Calibrant used"] = comp.name.map(sq["cal"]).values
        sheets["Identified_compounds"] = cc
        if len(pres):
            sheets["Presence_by_group"] = pres
        if Q is not None and len(Q["qc"]):
            sheets["QC_features"] = Q["qc"]
        lv = T.groupby(["sample", "level"]).size().unstack(fill_value=0).reindex(columns=[l for l in LEVEL_ORDER if l in set(T.level)])
        sheets["Level_summary"] = lv.reset_index()
    if TG is not None:
        sheets["IS_QC"] = TG["is_qc"]; sheets["Injection_QC"] = TG["injection_qc"]
        notes.append(("Targeted / IS", TG["xlsx"]))
    if C is not None:
        sheets["Calibration"] = C["summary"]; sheets["Calibration_points"] = C["points"]
        notes.append(("Calibration", f"linear and quadratic weighted fits; quadratic only when curvature is significant (F-test, AICc, monotonic); plots: {C['pdf']}"))
        if TG is not None and len(C["fits"]):
            from .calibration import quantify
            a = TG["assigned"]; cal_files = set(P.sheet.loc[P.sheet.type == "cal", "file"]) if P.sheet is not None else set()
            lst = pd.read_csv(os.path.join(P.out_dir, "calibration_list.csv")) if os.path.exists(os.path.join(P.out_dir, "calibration_list.csv")) else None
            isn = dict(zip(lst["name"], lst["is_name"].fillna(""))) if lst is not None and "is_name" in lst.columns else {}
            ar = a.set_index(["name", "file"])["int_area"]
            rows = []
            for x in a[a.name.isin(C["fits"]) & ~a.file.isin(cal_files)].itertuples():
                i = isn.get(x.name, "")
                resp = x.int_area / ar[(i, x.file)] if i and (i, x.file) in ar.index and ar[(i, x.file)] > 0 else (np.nan if i else x.int_area)
                rows.append(dict(name=x.name, file=x.file, response=resp))
            if rows:
                sheets["Quantification"] = quantify(C["fits"], pd.DataFrame(rows)).assign(unit=C["units"])
    ME = P._load("matrix_effect")
    if ME is not None:
        sheets["Matrix_effect"] = ME["table"]
        notes.append(("Matrix effect", "ME % = 100 x (area in spiked control extract / area in solvent standard - 1); < -20 % suppression, > +20 % enhancement"))
    if SQ is not None:
        keep = [c for c in ("sample", "name", "level", "rt", "area", "sq_conc", "sq_low", "sq_high", "sq_unit", "sq_method", "sq_nearest_IS",
                            "sq_nn_distance", "sq_recovery", "sq_ME_pct_at_rt", "sq_matrix_correction", "sq_flags", "rt_pred_from_properties",
                            "rt_property_plausible", "elucidation_distance", "sq_ed_ok", "sq_calibrant", "sq_logie", "sq_logp", "sq_logp_source",
                            "smiles") if c in SQ["table"].columns]
        sheets["Semi_quant"] = SQ["table"][keep]
        sheets["Properties"] = SQ["properties"]
        if SQ.get("method") == "is":
            sheets["IS_response_factors"] = SQ["is_rf"]
            notes.append(("Semi-quant", "IS surrogate: molar response of the 5 most similar internal standards (inverse-distance weighted; logKow, TPSA, "
                                        "molar volume, MW, H-bond donors/acceptors, charge) from the IS-mix injection; corrected for matrix effect; "
                                        "uncertainty x/÷3 near the IS, x/÷5 / x/÷10 further away, x/÷10 for void-eluting compounds"))
        else:
            notes.append(("Semi-quant", "elucidation distance in z-scored (logIE, RT, MW, logP) space to the nearest calibrant; conc = area / RF_cal x 10^(logIE_cal - logIE); "
                                        "logP from CompTox (experimental > predicted) else RDKit; ED < 0.5 = close analogue; report the x/÷ uncertainty factor"))
    notes += [("Reference databases", "Origin labels from the free structure databases (matched on the InChIKey first block): NORMAN Suspect "
                                      "List Exchange collections (incl. transformation-product lists S38/S66/S68/S74/S78/S81/S121 with "
                                      "literature parent-TP pairs), COCONUT, LOTUS, NP Atlas, LIPID MAPS and ChEBI roles. Context only - "
                                      "they do not change the identification level."),
              ("Levels","1 standard RT + library MS2; 2a library MS2 (cosine >= 0.70); 2b >= 3 decoy-corrected in-silico explained product ions; "
                         "3 1-2 explained ions / weak library; 4 exact mass (<= 5 ppm); 5 exact mass, no candidate (Schymanski et al. 2014)"),
              ("Peak rating", "0-10 from points across the peak, FWHM/base shape, jaggedness, modality and S/N (Compound-Discoverer-style; components in the confidence JSON)")]
    import openpyxl
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils import get_column_letter
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        pd.DataFrame(notes, columns=["Topic", "Note"]).to_excel(xw, sheet_name="Read_me", index=False)
        for k in ("Identified_compounds", "Detections_per_sample", "Presence_by_group", "Level_summary", "QC_features", "IS_QC", "Injection_QC",
                  "Calibration", "Calibration_points", "Quantification", "Matrix_effect", "Semi_quant", "IS_response_factors", "Properties", "Samples"):
            if k in sheets and len(sheets[k]):
                sheets[k].to_excel(xw, sheet_name=k, index=False)
    wb = openpyxl.load_workbook(path)
    hdr = PatternFill("solid", fgColor="1F6F78"); white = Font(color="FFFFFF", bold=True); fill = lambda c: PatternFill("solid", fgColor=c)
    for ws in wb:
        for cell in ws[1]:
            cell.fill = hdr; cell.font = white; cell.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = "B2"
        for i, col in enumerate(ws.columns, 1):
            L = max((len(str(c.value)) for c in list(col)[:300] if c.value is not None), default=8)
            ws.column_dimensions[get_column_letter(i)].width = min(max(10, L * 0.9), 60)
        col = {c.value: c.column for c in ws[1]}
        for r in ws.iter_rows(min_row=2):
            if "Level" in col:
                v = str(r[col["Level"] - 1].value)
                r[col["Level"] - 1].fill = fill(GREEN if v in ("1", "2a") else BLUE if v == "2b" else AMBER if v == "3" else RED)
            if "Presence confidence (%)" in col and r[col["Presence confidence (%)"] - 1].value is not None:
                v = r[col["Presence confidence (%)"] - 1].value
                r[col["Presence confidence (%)"] - 1].fill = fill(GREEN if v >= (Q["presence_threshold"] if Q else 70) else AMBER if v >= 50 else RED)
            if "Peak rating (0-10)" in col and r[col["Peak rating (0-10)"] - 1].value is not None:
                v = r[col["Peak rating (0-10)"] - 1].value
                r[col["Peak rating (0-10)"] - 1].fill = fill(GREEN if v >= 7 else AMBER if v >= 5 else RED)
        if ws.max_row > 1 and ws.title != "Read_me":
            ws.auto_filter.ref = ws.dimensions
    ws = wb["Read_me"]; ws.title = "Methods_notes"; ws.column_dimensions["A"].width = 22; ws.column_dimensions["B"].width = 150
    for r in ws.iter_rows(min_row=2):
        for c in r: c.alignment = Alignment(wrap_text=True, vertical="top")
    # cover sheet
    import datetime
    from . import __version__
    cv = wb.create_sheet("Cover", 0)
    cv.sheet_view.showGridLines = False
    cv["B2"] = "HRMS non-target screening report"; cv["B2"].font = Font(size=18, bold=True, color="1F6F78")
    cv["B3"] = f"pyMSscreen-NTS engine {__version__}  ·  generated {datetime.datetime.now():%d %b %Y %H:%M}"
    cv["B3"].font = Font(size=10, color="5A6663")
    info = [("Data folder", P.data_dir), ("Results folder", P.out_dir)]
    if P.sheet is not None:
        info.append(("Injections", ", ".join(f"{v} {k}" for k, v in P.sheet.type.value_counts().items())))
    if "Identified_compounds" in sheets:
        lv = sheets["Identified_compounds"]["Level"].astype(str).value_counts()
        info.append(("Compounds identified", ", ".join(f"Level {k}: {lv[k]}" for k in LEVEL_ORDER if k in lv.index)))
    if TG is not None:
        info.append(("Internal standards found", f"{TG['assigned'].name.nunique()} (see IS_QC)"))
    if ME is not None and len(ME["table"]):
        m = ME["table"].ME_pct
        info.append(("Matrix effect", f"median {m.median():+.0f} % ({(m < -20).sum()} suppressed, {(m > 20).sum()} enhanced of {len(m)})"))
    if SQ is not None:
        info.append(("Semi-quantification", "IS surrogate (kNN)" if SQ.get("method") == "is" else "elucidation distance"))
    info.append(("Sheets", "Identified_compounds = one row per compound; Detections_per_sample = every detection; "
                           "Methods_notes = definitions and settings"))
    for i, (k, v) in enumerate(info, start=5):
        cv[f"B{i}"] = k; cv[f"B{i}"].font = Font(bold=True, color="1B2422")
        cv[f"C{i}"] = str(v); cv[f"C{i}"].alignment = Alignment(wrap_text=True, vertical="top")
    cv.column_dimensions["A"].width = 3; cv.column_dimensions["B"].width = 26; cv.column_dimensions["C"].width = 110
    wb.save(path)
    P._report_notes = notes
    log.info("project report -> %s", path)
    return path
