"""Excel workbook + EIC integration PDF for the targeted screen (hrms_nts.targeted)."""
from __future__ import annotations
import os, logging
import numpy as np
import pandas as pd

log = logging.getLogger("hrms_nts")
GREEN, RED, AMBER, BLUE = "C6EFCE", "FFC7CE", "FFEB9C", "DDEBF7"
SALTS = ("HCl", "hydrochloride", "oxalate", "perchlorate", "sodium salt", "potassium salt", "tartrate", "citrate", "maleate", "fumarate")


def _detected_table(ch, comp, blank, tc):
    cm = comp.set_index("name")
    from .targeted import likely_negative
    bl = blank.groupby("name").blank_area.max() if len(blank) else pd.Series(dtype=float)
    bl_sig = set(blank.loc[blank.verdict.str.startswith("Significant"), "name"]) if len(blank) else set()
    blank_files = set(blank.blank_file) if len(blank) else set()
    med_rt = ch[~ch.file.isin(blank_files)].groupby("name").int_apex_rt.median()
    rows = []
    for x in ch.sort_values(["name", "file"]).itertuples():
        c = cm.loc[x.name]
        flags = []
        if x.file in blank_files:
            flags.append("peak in a BLANK file")
        if x.adduct.endswith("+") and c.smiles and likely_negative(c.smiles) and not str(x.level).startswith(("1", "2a")):
            flags.append("acidic compound (ESI- preferred): weak ESI+ evidence only - tentative")
        if x.name in med_rt.index and abs(x.int_apex_rt - med_rt[x.name]) > 0.3:
            flags.append(f"RT differs from the other files by {x.int_apex_rt - med_rt[x.name]:+.2f} min - probably a different compound")
        if x.name in bl_sig and x.file not in blank_files:
            flags.append("same signal in the blank (> 5 %) - background, not the compound")
        q1ok = abs(x.ppm) <= tc["q1_ppm"]
        q4 = "yes" if getattr(x, "q4_top3_all", False) is True and x.n_expl else ("partly" if x.n_q4 else "NO / none")
        lvl = x.level
        basis = {"1": "library MS2 + RT of an authentic standard on this method", "2a": f"library MS2 {x.lib_score:.2f} ({int(x.lib_matched)} ions)",
                 "2b": f"{int(x.n_q4)} product ions explained by the formula within {tc['q4_ppm']:g} ppm", "3": "1-2 product ions / weak library match",
                 "4": "Q1 exact mass only"}.get(lvl.split(" ")[0], "precursor outside Q1 tolerance")
        r = {"File": x.file, "Compound": x.name, "Role": c.role, "Formula (neutral)": c.formula, "Monoisotopic mass (neutral)": round(c.mono, 4),
             "Adduct": x.adduct, "Q1 theoretical": round(x.mz_theo, 4), "Q1 measured": round(x.mz_meas, 4), "Q1 error (ppm)": round(x.ppm, 2),
             f"Q1 <= {tc['q1_ppm']:g} ppm": "yes" if q1ok else "NO",
             "RT apex (min)": round(x.int_apex_rt, 2), "Integration start (min)": round(x.int_start_rt, 2), "Integration end (min)": round(x.int_end_rt, 2),
             f"Peak area (±{tc['integrate_ppm']:g} ppm XIC)": int(x.int_area), "Peak height": int(x.int_height), "S/N": int(x.int_sn),
             "Points across peak": int(x.int_n_points), "FWHM (s)": round(x.int_fwhm_s, 1) if np.isfinite(x.int_fwhm_s) else None,
             "Peak rating": x.int_rating}
        if len(bl):
            ba = bl.get(x.name, np.nan)
            r["Blank area (same RT window)"] = int(ba) if np.isfinite(ba) else None
            r["Blank as % of this peak"] = round(ba / x.int_area * 100, 2) if np.isfinite(ba) and x.int_area else None
        r.update({"Confidence level": lvl, "Confidence basis": basis + ("; co-elutes with an isobaric compound" if "isobaric" in lvl else ""),
                  "Q4 product ions (theo, formula, measured, error)": x.q4_ions if isinstance(x.q4_ions, str) else "",
                  f"Q4 ions <= {tc['q4_ppm']:g} ppm": f"{int(x.n_q4)}/{int(x.n_expl)}",
                  f"Q4 <= {tc['q4_ppm']:g} ppm": q4,
                  "Library match": (f"{x.lib_score:.2f} / {int(x.lib_matched)} ions ({x.lib_kind}: {x.lib_name}; {x.lib_source})" if x.lib_matched else ""),
                  "MS2 source": x.ms2_source if isinstance(x.ms2_source, str) else "",
                  "Expected RT (list / standards.csv)": round(x.rt_expected, 2) if np.isfinite(x.rt_expected) else None,
                  "RT predicted (QSAR / RTI)": ", ".join(f"{v:.2f}" for v in (x.rt_pred_qsar, x.rt_pred_rti) if isinstance(v, float) and np.isfinite(v)),
                  "Isobaric group": c.isobaric_group, "Flags": "; ".join(flags)})
        rows.append(r)
    return pd.DataFrame(rows)


def _not_detected(comp, ch, P, info, tc):
    from .targeted import likely_negative, ADDUCTS, MAIN
    found = set(ch.name) if len(ch) else set()
    pols = set(info.polarity)
    rows = []
    for c in comp[~comp.name.isin(found)].itertuples():
        why = []
        if pols == {"positive"} and c.smiles and likely_negative(c.smiles):
            why.append("acidic, no basic nitrogen - ionises in ESI- (files are ESI+ only)")
        near = P[(P.name == c.name)] if len(P) else P
        if len(near):
            b = near.sort_values("height", ascending=False).iloc[0]
            why.append(f"nearest candidate {b.adduct} at {b.rt:.2f} min, {b.ppm:+.1f} ppm, height {int(b.height)} - outside acceptance")
        if not why:
            why.append(f"no peak >= {tc['min_apex']:g} counts within {tc['wide_ppm']:g} ppm in any file (not in the mix, below detection, or eluting in the void)")
        rows.append({"Compound": c.name, "Formula (neutral)": c.formula, **{f"{MAIN[p]}": round(c.mono + ADDUCTS[p][MAIN[p]], 4) for p in ("positive", "negative")},
                     "Reason": "; ".join(why)})
    return pd.DataFrame(rows)


def _mass_check(comp):
    from .targeted import ADDUCTS
    rows = []
    for c in comp.itertuples():
        theo = c.mono + ADDUCTS["positive"]["[M+H]+"]
        r = {"Compound": c.name, "Formula (neutral, free base/acid)": c.formula, "Q1 [M+H]+ calculated": round(theo, 4),
             "Neutral monoisotopic mass (PubChem 'exact mass')": round(c.mono, 4), "Average MW, free base (vendor 'MW')": round(c.avg_mw, 2),
             "[M-H]- calculated": round(c.mono + ADDUCTS["negative"]["[M-H]-"], 4), "Unlabelled parent [M+H]+": round(c.parent_mono + ADDUCTS["positive"]["[M+H]+"], 4),
             "Salt in name (vendor MW includes it)": next((s for s in SALTS if s.lower() in c.name.lower()), ""),
             "SMILES source": c.smiles_source}
        if np.isfinite(c.q1_list):
            e = (c.q1_list - theo) / theo * 1e6
            r.update({"List Q1": c.q1_list, "List vs calculated (ppm)": round(e, 1), "List Q1 OK?": "yes" if abs(e) < 1 else "WRONG"})
        rows.append(r)
    return pd.DataFrame(rows)


def _notes(info, ch, blank, comp, tc, cfg):
    n = []
    n.append(("Files", "; ".join(f"{r.file}: {r.instrument}, {r.acquisition}, {r.polarity}, {r.role}"
                                  + (f", conc {r.conc:g}" if np.isfinite(r.conc) else "") + (f", started {r.start_time}" if r.start_time else "")
                                  + (" (profile data centroided on reading)" if r.profile_centroided else "") + (f" - NO SCANS {r.gaps}" if r.gaps else "")
                                  + (" - TRUNCATED" if r.truncated else "") for r in info.itertuples())))
    n.append(("Q1 / Q4", f"Q1 = precursor ion (MS1), acceptance <= {tc['q1_ppm']:g} ppm. Q4 = product ions (DIA: fragments whose XIC correlates with the "
                         f"precursor; DDA: precursor-matched MS2), annotated as sub-formulae of the (labelled) formula, acceptance <= {tc['q4_ppm']:g} ppm. "
                         "Error = (measured - theoretical) / theoretical x 1e6. Search window is +-%g ppm so a calibration offset cannot hide a compound." % tc["wide_ppm"]))
    n.append(("Confidence", f"1 = library MS2 + RT of an authentic standard on this method (list rt / standards.csv, +-{cfg['levels']['level1_rt_tol']:g} min). "
                            f"2a = label-aware library cosine >= {tc['lib_cosine']:.2f} with >= {tc['lib_min_matched']} ions (installed spectral library; a library "
                            "fragment may match with k deuteriums / 13C retained). 2b = >= 3 product ions explained within Q4. 3 = 1-2 product ions or weak library match. "
                            "4 = Q1 only. 'isobaric, unresolved' = compounds with the same formula sharing one peak that MS2 could not separate."))
    n.append(("Peak integration", f"+-{tc['integrate_ppm']:g} ppm XIC (most intense centroid per scan); bounds at baseline + 5 % of apex or a valley; "
                                  "area = trapezoid above baseline (counts x s). S/N = height / max(noise, noise floor); profile exports omit zero points so the "
                                  "baseline is often empty and the noise floor (centroid threshold) is used. Good: S/N >= 100, >= 8 points, not valley-limited; "
                                  "Fair: S/N >= 10, >= 5 points; Poor: below. Every integration is drawn in EIC_integration_QC.pdf."))
    if len(ch):
        a = ch.int_area
        n.append(("Area scale", f"Assigned peaks: area {a.min():,.0f} - {a.max():,.0f} (median {a.median():,.0f}) counts x s; Q1 error median {ch.ppm.median():+.2f} ppm, "
                                f"{(ch.ppm.abs() <= tc['q1_ppm']).mean() * 100:.0f} % within {tc['q1_ppm']:g} ppm."))
    if len(blank):
        n.append(("Blank", "Blank XICs integrated over each compound's sample RT window without peak picking (any signal counts). Verdicts: "
                           + ", ".join(f"{k}: {v}" for k, v in blank.verdict.str.split(" (", regex=False).str[0].value_counts().items())))
    iso = comp.loc[comp.isobaric_group != "", "isobaric_group"].unique()
    if len(iso):
        n.append(("Isobaric compounds", "Same formula (cannot be told apart by mass): " + " | ".join(iso) +
                  ". A peak is given to one member only when its library MS2 match is clearly better; otherwise flagged unresolved."))
    n.append(("RT", "RT is reported, not used as an identification criterion unless an expected RT is given (list column rt or user_data/standards.csv), "
                    "which also enables Level 1. RT models (QSAR / NORMAN RTI, trained on standards.csv) give a plausibility prediction. "
                    "rt_anchors_for_standards.csv = RTs of the detected labelled standards for their native compounds; add them with --update-standards."))
    return n


def write_report(out, comp, info, P, Fr, ch, blank, cal, quant, natives, anchors, tc, cfg):
    import openpyxl
    from openpyxl.styles import Font, PatternFill, Alignment
    from openpyxl.utils import get_column_letter
    path = os.path.join(out, "targeted_results.xlsx")
    det = _detected_table(ch, comp, blank, tc) if len(ch) else pd.DataFrame()
    nd = _not_detected(comp, ch, P, info, tc)
    mc = _mass_check(comp)
    bk = blank.rename(columns={"name": "Compound"}) if len(blank) else pd.DataFrame()
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        pd.DataFrame(_notes(info, ch, blank, comp, tc, cfg), columns=["Topic", "Note"]).to_excel(xw, sheet_name="Notes", index=False)
        det.to_excel(xw, sheet_name="Detected", index=False)
        if len(bk): bk.round(3).to_excel(xw, sheet_name="Blank_check", index=False)
        nd.to_excel(xw, sheet_name="Not_detected", index=False)
        if len(cal): cal.round(6).to_excel(xw, sheet_name="Calibration", index=False)
        if len(quant): quant.round(5).to_excel(xw, sheet_name="Quantification", index=False)
        if len(natives): natives.round(4).to_excel(xw, sheet_name="Native_at_IS_RT", index=False)
        if len(anchors): anchors.to_excel(xw, sheet_name="RT_anchors", index=False)
        mc.to_excel(xw, sheet_name="Mass_check", index=False)
        info.to_excel(xw, sheet_name="Files", index=False)
        if len(P): P.round(4).to_excel(xw, sheet_name="All_candidate_peaks", index=False)
        if len(Fr): Fr.round(4).to_excel(xw, sheet_name="Product_ions", index=False)
    wb = openpyxl.load_workbook(path)
    hdr = PatternFill("solid", fgColor="1F6F78"); white = Font(color="FFFFFF", bold=True)
    fill = lambda c: PatternFill("solid", fgColor=c)
    for ws in wb:
        for cell in ws[1]:
            cell.fill = hdr; cell.font = white; cell.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = "C2" if ws.title == "Detected" else "B2"
        for i, col in enumerate(ws.columns, 1):
            L = max((len(str(c.value)) for c in list(col)[1:200] if c.value is not None), default=8)
            ws.column_dimensions[get_column_letter(i)].width = min(max(12, L * 0.95), 70)
        ws.row_dimensions[1].height = 45
    if "Detected" in wb.sheetnames and len(det):
        ws = wb["Detected"]; col = {c.value: c.column for c in ws[1]}
        for r in ws.iter_rows(min_row=2):
            for k in (f"Q1 <= {tc['q1_ppm']:g} ppm", f"Q4 <= {tc['q4_ppm']:g} ppm"):
                v = r[col[k] - 1].value; r[col[k] - 1].fill = fill(GREEN if v == "yes" else AMBER if v == "partly" else RED)
            v = str(r[col["Confidence level"] - 1].value)
            r[col["Confidence level"] - 1].fill = fill(GREEN if v[:1] == "1" or v.startswith("2a") else BLUE if v.startswith("2b") else AMBER if v.startswith("3") else RED)
            v = r[col["Peak rating"] - 1].value
            r[col["Peak rating"] - 1].fill = fill(GREEN if v == "Good" else AMBER if str(v).startswith("Fair") else RED)
            if r[col["Flags"] - 1].value:
                r[col["Flags"] - 1].fill = fill(AMBER)
            if "Blank as % of this peak" in col:
                v = r[col["Blank as % of this peak"] - 1].value
                if v is not None: r[col["Blank as % of this peak"] - 1].fill = fill(GREEN if v < tc["blank_trace_pct"] else AMBER if v < tc["blank_low_pct"] else RED)
    if "Blank_check" in wb.sheetnames:
        ws = wb["Blank_check"]; col = {c.value: c.column for c in ws[1]}
        for r in ws.iter_rows(min_row=2):
            v = str(r[col["verdict"] - 1].value)
            r[col["verdict"] - 1].fill = fill(GREEN if v.startswith(("Clean", "Minimal")) else AMBER if v.startswith("Low") else RED)
    ws = wb["Mass_check"]; col = {c.value: c.column for c in ws[1]}
    if "List Q1 OK?" in col:
        for r in ws.iter_rows(min_row=2):
            v = r[col["List Q1 OK?"] - 1].value
            if v: r[col["List Q1 OK?"] - 1].fill = fill(GREEN if v == "yes" else RED)
    ws = wb["Notes"]; ws.column_dimensions["A"].width = 18; ws.column_dimensions["B"].width = 150
    for r in ws.iter_rows(min_row=2):
        for c in r: c.alignment = Alignment(wrap_text=True, vertical="top")
    wb.save(path)
    log.info("targeted workbook -> %s", path)
    return path


def eic_pdf(path, ch, res, comp, blank_files, tc):
    """One panel per compound: +-integrate_ppm XIC of every file overlaid (blanks thick), integrated region shaded."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pal = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#7d4bc2", "#c43c3c"]
    runs = {r["info"]["file"]: r["run"] for r in res}
    files = list(runs)
    col = {f: pal[i % len(pal)] for i, f in enumerate(files)}
    names = sorted(ch.name.unique(), key=str.lower)
    plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#888"})
    with PdfPages(path) as pdf:
        for k in range(0, len(names), 12):
            fig, axs = plt.subplots(4, 3, figsize=(11.7, 8.3))
            for ax, n in zip(axs.flat, names[k:k + 12]):
                q = ch[ch.name == n]; mz0 = q.mz_theo.iloc[0]; c0 = q.int_apex_rt.median()
                for f in files:
                    rt, I, _ = runs[f].xic_max(mz0, tc["integrate_ppm"])
                    m = np.abs(rt - c0) <= 1.0
                    b = f in blank_files
                    ax.plot(rt[m], I[m], color=col[f], lw=1.6 if b else 1.1, label=f + (" (blank)" if b else ""), zorder=5 if b else 2)
                    x = q[q.file == f]
                    if len(x):
                        x = x.iloc[0]; s = (rt >= x.int_start_rt) & (rt <= x.int_end_rt)
                        ax.fill_between(rt[s], x.int_baseline, I[s], color=col[f], alpha=0.25, lw=0)
                ax.set_title(f"{n[:34]}  {q.adduct.iloc[0]} {mz0:.4f} ±{tc['integrate_ppm']:g} ppm", fontsize=7, loc="left")
                ax.set_xlabel("RT (min)"); ax.set_ylabel("counts"); ax.ticklabel_format(axis="y", style="plain")
            for ax in list(axs.flat)[len(names[k:k + 12]):]:
                ax.axis("off")
            h, l = axs.flat[0].get_legend_handles_labels()
            fig.legend(h, l, loc="upper center", ncol=min(4, len(files)), frameon=False, fontsize=8)
            fig.suptitle("Shaded = integrated region above baseline; ±1 min around the apex.", y=0.955, fontsize=8, color="#555")
            fig.tight_layout(rect=(0, 0, 1, 0.94)); pdf.savefig(fig); plt.close(fig)
    log.info("EIC QC -> %s", path)
