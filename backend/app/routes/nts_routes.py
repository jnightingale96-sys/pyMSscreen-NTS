"""
Routes for the HRMS NTS engine (``engine/hrms_nts``): non-target screening (Level 5-1),
semi-quantification and product-ion prediction.

Long runs (``run``, ``calibrate``) are started as ``python -m hrms_nts ...`` subprocesses so the
engine behaves exactly as it does from the command line and the web server stays responsive.
Short calls (semi-quant re-calculation, product-ion prediction) run in-process.
"""
import datetime
import glob
import os
import subprocess
import sys
import threading
import uuid

import numpy as np
import pandas as pd
from flask import Blueprint, jsonify, request, send_file

from app.config import BASE_WORKING_DIR
from app.utils import resolve_path

nts_bp = Blueprint("nts", __name__, url_prefix="/nts")

LEVEL_ORDER = ["1", "2a", "2b", "3", "4", "5"]
RESULT_TABLES = {
    "compounds": "compound_summary.csv",
    "identified": "identified_level1_2b.csv",
    "top": "top_candidates_per_feature.csv",
    "levels": "level_summary.csv",
    "features": "features.csv",
}
TOP_COLUMNS = ["sample", "feature_id", "level", "confidence_score", "name", "formula", "adduct", "mz_apex", "ppm_error",
               "rt_apex", "std_rt_diff", "rt_agree", "rt_conflict", "lib_score", "lib_matched", "lib_name", "bi_score",
               "bi_expl_adj_n", "mf_score", "iso_score", "max_int", "area", "background_flag", "blank_flag",
               "stereo_ambiguous", "level_basis", "smiles", "inchikey"]
DOWNLOADABLE = (".csv", ".xlsx", ".log", ".png", ".yaml", ".json", ".pdf")

_JOBS = {}
_JOBS_LOCK = threading.Lock()


def _engine_cfg(overrides=None):
    from hrms_nts.common import load_config
    return load_config(None, overrides)


def _clean_records(df, limit=None):
    """DataFrame -> JSON-safe list of dicts (NaN/inf -> None)."""
    if limit:
        df = df.head(limit)
    df = df.replace([np.inf, -np.inf], np.nan).astype(object).where(df.notna(), None)
    return df.to_dict(orient="records")


def _run_dir_arg(name="run_dir"):
    src = request.args if request.method == "GET" else (request.get_json(silent=True) or {})
    run_dir = resolve_path((src.get(name) or "").strip())
    if not run_dir or not os.path.isdir(run_dir):
        return None
    return os.path.abspath(run_dir)


def _start_job(kind, args, log_path, out_path, n_files=0):
    job_id = uuid.uuid4().hex[:12]
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logf = open(log_path, "w", encoding="utf-8", errors="replace")
    env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUNBUFFERED="1")
    proc = subprocess.Popen([sys.executable, "-m", "hrms_nts", *args], stdout=logf, stderr=subprocess.STDOUT,
                            cwd=os.path.dirname(log_path), env=env,
                            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    with _JOBS_LOCK:
        _JOBS[job_id] = dict(id=job_id, kind=kind, proc=proc, log=log_path, logf=logf, out=out_path, n_files=n_files,
                             started_at=datetime.datetime.now().isoformat(timespec="seconds"), cancelled=False)
    return job_id


def _tail(path, n=40):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            return f.read().splitlines()[-n:]
    except OSError:
        return []


# ------------------------------------------------------------------------------------------------ info
@nts_bp.route("/info", methods=["GET"])
def info():
    """Engine version, option lists and which data sources were found (the `doctor` check)."""
    import shutil
    from hrms_nts import __version__
    from hrms_nts.common import resolve
    try:
        cfg = _engine_cfg()
    except Exception as e:
        return jsonify({"error": f"Could not load engine config: {e}"}), 500
    checks = []
    for k in ("susdat", "pubchemlite", "user_compounds", "library_db", "library_sources_dir", "standards", "calibrants", "metfrag_jar", "inhouse_msp"):
        p = resolve(cfg, cfg["paths"].get(k))
        checks.append({"item": f"paths.{k}", "path": p, "found": bool(p and os.path.exists(p))})
    from hrms_nts.insilico import metfrag_available
    mf_ok, _, java = metfrag_available(cfg)          # resolves a bundled/relative JRE (e.g. data/jre/bin/java.exe)
    checks.append({"item": "java (MetFrag)", "path": shutil.which(java), "found": bool(mf_ok)})
    sq = cfg["semiquant"]
    library = None
    lib_path = resolve(cfg, cfg["paths"].get("library_db"))
    if lib_path and os.path.exists(lib_path):
        import json
        import sqlite3
        con = sqlite3.connect(lib_path)
        library = {"path": lib_path, "n_spectra": con.execute("select count(*) from spectra").fetchone()[0]}
        con.close()
        man = os.path.join(os.path.dirname(lib_path), "spectral_library_manifest.json")
        if os.path.exists(man):
            with open(man, encoding="utf-8") as f:
                m = json.load(f)
            library.update(built=m.get("built"), unique_structures=m.get("main", {}).get("unique_structures"),
                           sources=[s["source"] for s in m.get("main", {}).get("by_source", [])])
    return jsonify({
        "version": __version__,
        "library": library,
        "default_working_dir": BASE_WORKING_DIR,
        "config_layers": cfg.get("_config_layers", []),
        "checks": checks,
        "matrices": {k: v.get("unit", "") for k, v in sq["matrices"].items()},
        "semiquant_modes": ["single", "logie", "curve", "none"],
        "semiquant_default": sq.get("mode", "single"),
        "single_reference": sq.get("single_reference", {}),
        "uncertainty_factor": sq.get("uncertainty_factor", 10),
        "insilico_engines": ["builtin", "metfrag", "both", "none"],
        "insilico_default": cfg["insilico"].get("engine", "builtin"),
        "ppm_default": cfg["level4"]["ppm"],
        "min_intensity_default": cfg["ms1"]["min_intensity"],
        "adducts": {"positive": ["[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+"],
                    "negative": ["[M-H]-", "[M+Cl]-", "[M+HCOO]-", "[M+CH3COO]-"]},
    })


@nts_bp.route("/list_mzml", methods=["GET"])
def list_mzml():
    folder = resolve_path(request.args.get("folder", "").strip())
    if not folder or not os.path.isdir(folder):
        return jsonify({"error": f"Folder not found: {folder}"}), 400
    files = sorted(set(glob.glob(os.path.join(folder, "*.mzML")) + glob.glob(os.path.join(folder, "*.mzml"))))
    return jsonify({"folder": os.path.abspath(folder),
                    "files": [{"name": os.path.basename(f), "size_mb": round(os.path.getsize(f) / 1e6, 1)} for f in files]})


# ------------------------------------------------------------------------------------------------ screening run
@nts_bp.route("/start_run", methods=["POST"])
def start_run():
    body = request.get_json(force=True) or {}
    folder = resolve_path((body.get("folder") or "").strip())
    if not os.path.isdir(folder):
        return jsonify({"error": f"mzML folder not found: {folder}"}), 400
    names = body.get("files") or []
    files = [os.path.join(folder, n) for n in names if os.path.exists(os.path.join(folder, n))]
    if not files:
        return jsonify({"error": "Select at least one mzML file."}), 400
    out = resolve_path((body.get("out_dir") or "").strip()) or os.path.join(
        folder, "NTS_out_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    os.makedirs(out, exist_ok=True)

    args = []
    for kv in body.get("overrides") or []:          # free-form `--set section.key=value`
        if isinstance(kv, str) and "=" in kv:
            args += ["--set", kv.strip()]
    args += ["run", *files, "--out", out, "--matrix", body.get("matrix", "water"),
             "--polarity", body.get("polarity", "auto"), "--semiquant", body.get("semiquant", "single"),
             "--insilico", body.get("insilico", "builtin")]
    if body.get("ppm"):
        args += ["--ppm", str(float(body["ppm"]))]
    if body.get("min_intensity"):
        args += ["--min-intensity", str(float(body["min_intensity"]))]
    if body.get("no_library"):
        args.append("--no-library")
    job_id = _start_job("run", args, os.path.join(out, "gui_run.log"), out, n_files=len(files))
    return jsonify({"job_id": job_id, "out_dir": out})


@nts_bp.route("/start_calibrate", methods=["POST"])
def start_calibrate():
    body = request.get_json(force=True) or {}
    mzml = resolve_path((body.get("mzml") or "").strip())
    if not os.path.isfile(mzml):
        return jsonify({"error": f"mzML not found: {mzml}"}), 400
    out_csv = os.path.join(os.path.dirname(mzml), f"calibrant_areas_{os.path.splitext(os.path.basename(mzml))[0]}.csv")
    args = ["calibrate", mzml, "--out", out_csv, "--polarity", body.get("polarity", "auto")]
    job_id = _start_job("calibrate", args, out_csv[:-4] + ".log", out_csv)
    return jsonify({"job_id": job_id, "out": out_csv})


# ------------------------------------------------------------------------------------------------ targeted screen (IS / targets / calibrants)
@nts_bp.route("/start_targeted", methods=["POST"])
def start_targeted():
    """body: folder, files[], list (CSV/XLSX), blanks[] (file names), sample_sheet, q1_ppm, q4_ppm, min_apex, polarity,
    update_standards, out_dir."""
    body = request.get_json(force=True) or {}
    folder = resolve_path((body.get("folder") or "").strip())
    if not os.path.isdir(folder):
        return jsonify({"error": f"mzML folder not found: {folder}"}), 400
    files = [os.path.join(folder, n) for n in (body.get("files") or []) if os.path.exists(os.path.join(folder, n))]
    if not files:
        return jsonify({"error": "Select at least one mzML file."}), 400
    lst = resolve_path((body.get("list") or "").strip())
    if not lst or not os.path.isfile(lst):
        return jsonify({"error": f"Compound list not found: {lst}"}), 400
    out = resolve_path((body.get("out_dir") or "").strip()) or os.path.join(
        folder, "targeted_out_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    os.makedirs(out, exist_ok=True)
    args = ["targeted", *files, "--list", lst, "--out", out, "--polarity", body.get("polarity", "auto")]
    blanks = [os.path.splitext(b)[0] for b in (body.get("blanks") or [])]
    if blanks:
        args += ["--blank", *blanks]
    sheet = resolve_path((body.get("sample_sheet") or "").strip())
    roles = body.get("roles") or {}
    if roles and not (sheet and os.path.isfile(sheet)):          # roles chosen on the page -> sample sheet in the results folder
        sheet = os.path.join(out, "sample_sheet.csv")
        pd.DataFrame([{"file": os.path.splitext(n)[0], "type": r.get("type", "sample"), "conc": r.get("conc") or None}
                      for n, r in roles.items()]).to_csv(sheet, index=False)
    if sheet and os.path.isfile(sheet):
        args += ["--sample-sheet", sheet]
    for k in ("q1_ppm", "q4_ppm", "min_apex"):
        if body.get(k):
            args += ["--" + k.replace("_", "-"), str(float(body[k]))]
    if body.get("update_standards"):
        args.append("--update-standards")
    job_id = _start_job("targeted", args, os.path.join(out, "gui_targeted.log"), out, n_files=len(files))
    return jsonify({"job_id": job_id, "out_dir": out})


@nts_bp.route("/targeted_results", methods=["GET"])
def targeted_results():
    run_dir = _run_dir_arg()
    if not run_dir:
        return jsonify({"error": "results folder not found"}), 400
    x = os.path.join(run_dir, "targeted_results.xlsx")
    if not os.path.exists(x):
        return jsonify({"error": "targeted_results.xlsx not found"}), 404
    sheets = pd.read_excel(x, sheet_name=None)
    keep = ("Detected", "Blank_check", "Not_detected", "Calibration", "Quantification", "Native_at_IS_RT", "RT_anchors", "Notes", "Files")
    out = {k: _clean_records(v) for k, v in sheets.items() if k in keep}
    out["files"] = [f for f in ("targeted_results.xlsx", "EIC_integration_QC.pdf", "rt_anchors_for_standards.csv",
                                "calibrants_for_semiquant.csv", "run.log") if os.path.exists(os.path.join(run_dir, f))]
    return jsonify(out)


@nts_bp.route("/job", methods=["GET"])
def job():
    job_id = request.args.get("job_id", "")
    with _JOBS_LOCK:
        j = _JOBS.get(job_id)
    if not j:
        return jsonify({"id": job_id, "state": "unknown", "message": "Job not found (server restarted?)"}), 200
    rc = j["proc"].poll()
    lines = _tail(j["log"], 400)
    done_samples = sum(1 for ln in lines if " features, " in ln and " candidate rows" in ln)
    failed = sum(1 for ln in lines if "sample failed" in ln)
    if rc is None:
        state = "running"
    else:
        if not j["logf"].closed:
            j["logf"].close()
        state = "cancelled" if j["cancelled"] else ("done" if rc == 0 else "error")
        if state == "done" and j["kind"] == "run" and j["n_files"] and failed >= j["n_files"]:
            state = "error"          # the CLI exits 0 even when every sample failed
    progress = None
    if j["kind"] == "targeted" and j["n_files"]:
        screened = sum(1 for ln in lines if "targeted: " in ln and " | " in ln)
        progress = 100 if state == "done" else round(90 * min(screened, j["n_files"]) / j["n_files"])
    if j["kind"] == "run" and j["n_files"]:
        progress = 100 if state == "done" else round(95 * min(done_samples + failed, j["n_files"]) / j["n_files"])
    res = dict(id=job_id, kind=j["kind"], state=state, returncode=rc, out=j["out"], started_at=j["started_at"],
               progress=progress, samples_done=done_samples, samples_failed=failed, n_files=j["n_files"], log=lines[-40:],
               warnings=[ln for ln in lines if " WARNING " in ln or "sample failed" in ln][-20:])
    if j["kind"] == "calibrate" and state == "done" and os.path.exists(j["out"]):
        res["table"] = _clean_records(pd.read_csv(j["out"]))
    return jsonify(res)


@nts_bp.route("/cancel", methods=["POST"])
def cancel():
    job_id = (request.get_json(force=True) or {}).get("job_id", "")
    with _JOBS_LOCK:
        j = _JOBS.get(job_id)
    if not j:
        return jsonify({"error": "Job not found"}), 404
    if j["proc"].poll() is None:
        j["cancelled"] = True
        j["proc"].terminate()
    return jsonify({"status": "ok"})


# ------------------------------------------------------------------------------------------------ results
@nts_bp.route("/runs", methods=["GET"])
def runs():
    """Result folders (anything containing level_summary.csv) directly inside `folder` (or `folder` itself)."""
    folder = resolve_path(request.args.get("folder", "").strip())
    if not folder or not os.path.isdir(folder):
        return jsonify({"runs": []})
    cands = [folder] + [os.path.join(folder, d) for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]
    out = []
    for d in cands:
        ls = os.path.join(d, "level_summary.csv")
        if os.path.exists(ls):
            out.append({"run_dir": os.path.abspath(d), "name": os.path.basename(d),
                        "modified": datetime.datetime.fromtimestamp(os.path.getmtime(ls)).isoformat(timespec="minutes")})
    out.sort(key=lambda r: r["modified"], reverse=True)
    return jsonify({"runs": out})


@nts_bp.route("/results", methods=["GET"])
def results():
    run_dir = _run_dir_arg()
    if not run_dir:
        return jsonify({"error": "Results folder not found"}), 400
    table = request.args.get("table", "compounds")
    fn = RESULT_TABLES.get(table)
    if not fn:
        return jsonify({"error": f"unknown table {table}"}), 400
    path = os.path.join(run_dir, fn)
    if not os.path.exists(path):
        return jsonify({"rows": [], "columns": [], "total": 0})
    df = pd.read_csv(path, low_memory=False)
    if "level" in df.columns:
        df["level"] = df["level"].astype(str)
    if "best_level" in df.columns:
        df["best_level"] = df["best_level"].astype(str)
    levels = request.args.get("levels")
    lvl_col = "best_level" if "best_level" in df.columns else ("level" if "level" in df.columns else None)
    if levels and lvl_col:
        df = df[df[lvl_col].isin(levels.split(","))]
    q = (request.args.get("q") or "").strip().lower()
    if q and "name" in df.columns:
        df = df[df["name"].astype(str).str.lower().str.contains(q, regex=False)]
    if table in ("top", "identified"):
        keep = [c for c in TOP_COLUMNS if c in df.columns] + [c for c in df.columns if c.startswith("conc_") or c.startswith("sq_")]
        df = df[keep]
    total = len(df)
    limit = int(request.args.get("limit", 2000))
    return jsonify({"columns": list(df.columns), "rows": _clean_records(df, limit), "total": total})


@nts_bp.route("/plots", methods=["GET"])
def plots():
    run_dir = _run_dir_arg()
    if not run_dir:
        return jsonify({"plots": []})
    files = sorted(glob.glob(os.path.join(run_dir, "plots", "*.png")))
    return jsonify({"plots": [os.path.basename(f) for f in files]})


@nts_bp.route("/file", methods=["GET"])
def file():
    """Serve a result file that lives inside a results folder (no path escapes)."""
    run_dir = _run_dir_arg()
    name = request.args.get("name", "")
    if not run_dir or not name:
        return jsonify({"error": "missing run_dir or name"}), 400
    path = os.path.abspath(os.path.join(run_dir, name))
    if os.path.commonpath([path, run_dir]) != run_dir or not path.lower().endswith(DOWNLOADABLE) or not os.path.isfile(path):
        return jsonify({"error": "file not available"}), 404
    return send_file(path, as_attachment=request.args.get("download") == "1")


@nts_bp.route("/export_to_prescreen", methods=["POST"])
def export_to_prescreen():
    """Write identified compounds as a pyMSscreen compounds CSV (ID,Name,SMILES,mz) so they can be
    checked in the Extraction / Plotting pages."""
    body = request.get_json(force=True) or {}
    run_dir = _run_dir_arg()
    if not run_dir:
        return jsonify({"error": "Results folder not found"}), 400
    levels = body.get("levels") or ["1", "2a", "2b"]
    path = os.path.join(run_dir, "compound_summary.csv")
    if not os.path.exists(path):
        return jsonify({"error": "compound_summary.csv not found in results folder"}), 404
    df = pd.read_csv(path)
    df["best_level"] = df["best_level"].astype(str)
    df = df[df["best_level"].isin(levels) & df["smiles"].notna()].drop_duplicates("smiles")
    if df.empty:
        return jsonify({"error": f"No compounds at levels {', '.join(levels)}"}), 400
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    rows = []
    for i, r in enumerate(df.itertuples(), 1):
        m = Chem.MolFromSmiles(r.smiles)
        if m is None:
            continue
        rows.append({"ID": f"NTS{i:04d}", "Name": r.name, "SMILES": r.smiles, "mz": round(Descriptors.ExactMolWt(m), 6),
                     "nts_level": r.best_level})
    target_dir = resolve_path((body.get("target_dir") or "").strip()) or run_dir
    if not os.path.isdir(target_dir):
        return jsonify({"error": f"Target folder not found: {target_dir}"}), 400
    out = os.path.join(target_dir, f"compounds_from_NTS_L{'-'.join(levels)}.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    return jsonify({"status": "success", "path": out, "n": len(rows)})


# ------------------------------------------------------------------------------------------------ semi-quant
@nts_bp.route("/calibrants", methods=["GET", "POST"])
def calibrants():
    from hrms_nts.common import resolve
    cfg = _engine_cfg()
    path = resolve(cfg, cfg["paths"]["calibrants"])
    cols = ["name", "smiles", "rt", "conc_ug_ml", "area", "polarity"]
    if request.method == "GET":
        if not path or not os.path.exists(path):
            return jsonify({"path": path, "rows": []})
        from hrms_nts.semiquant import calibrants_table
        df = calibrants_table(path)
        return jsonify({"path": path, "rows": _clean_records(df.reindex(columns=cols))})
    rows = (request.get_json(force=True) or {}).get("rows") or []
    df = pd.DataFrame(rows).reindex(columns=cols)
    df = df[df["name"].astype(str).str.strip().ne("") & df["name"].notna()]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(",".join(cols) + "\n")
        f.write("# Semi-quant calibrants on THIS method (edited in pyMSscreen-NTS). >= 3 levels per compound enables curve mode.\n")
        df.to_csv(f, index=False, header=False)
    return jsonify({"status": "success", "path": path, "n": len(df)})


@nts_bp.route("/semiquant", methods=["POST"])
def semiquant():
    """Re-calculate concentrations for an existing results folder with another mode / matrix / reference,
    without re-running the screening."""
    from hrms_nts.common import resolve
    from hrms_nts.semiquant import calibrants_table, semiquantify
    body = request.get_json(force=True) or {}
    run_dir = _run_dir_arg()
    if not run_dir:
        return jsonify({"error": "Results folder not found"}), 400
    src = os.path.join(run_dir, "top_candidates_per_feature.csv")
    if not os.path.exists(src):
        return jsonify({"error": "top_candidates_per_feature.csv not found"}), 404
    mode = body.get("mode", "single")
    matrix = body.get("matrix", "water")
    ov = {"semiquant": {"mode": mode}}
    if body.get("single_reference"):
        sr = body["single_reference"]
        ov["semiquant"]["single_reference"] = {"name": sr.get("name"), "conc_ug_mL": float(sr["conc_ug_mL"]), "area": float(sr["area"])}
    if body.get("uncertainty_factor"):
        ov["semiquant"]["uncertainty_factor"] = float(body["uncertainty_factor"])
    cfg = _engine_cfg(ov)
    if matrix not in cfg["semiquant"]["matrices"]:
        return jsonify({"error": f"unknown matrix {matrix}"}), 400
    top = pd.read_csv(src, low_memory=False)
    top["level"] = top["level"].astype(str)
    levels = body.get("levels") or ["1", "2a", "2b", "3", "4"]
    top = top[top["level"].isin(levels)]
    top = top.drop(columns=[c for c in top.columns if c.startswith("conc_") or c.startswith("sq_") or c == "below_loq"])
    cal = None
    cp = resolve(cfg, cfg["paths"].get("calibrants"))
    if cp and os.path.exists(cp):
        cal = calibrants_table(cp)
    parts = []
    try:
        for pol, g in top.groupby("polarity"):
            parts.append(semiquantify(g, cfg, matrix, cal, polarity=pol, area_col=cfg["semiquant"].get("area_column", "area")))
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    res = pd.concat(parts, ignore_index=True) if parts else top
    out = os.path.join(run_dir, f"semiquant_{mode}_{matrix}.csv")
    res.to_csv(out, index=False)
    c = f"conc_{matrix}"
    show = [x for x in ["sample", "level", "name", "rt_apex", "area", "max_int", "sq_method", "sq_calibrant",
                        "conc_ug_mL_extract", c, f"{c}_low", f"{c}_high", "sq_unit", "smiles"] if x in res.columns]
    res = res.sort_values(["level", c], ascending=[True, False]) if c in res.columns else res
    return jsonify({"status": "success", "path": out, "unit": cfg["semiquant"]["matrices"][matrix].get("unit", ""),
                    "columns": show, "rows": _clean_records(res[show], 3000), "total": len(res)})


# ------------------------------------------------------------------------------------------------ product ions
@nts_bp.route("/structure.svg", methods=["GET"])
def structure_svg():
    """2D structure drawn locally with RDKit (works offline, unlike an online depiction service)."""
    from flask import Response
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    m = Chem.MolFromSmiles(request.args.get("smiles", ""))
    size = max(80, min(int(request.args.get("size", 260)), 1200))
    if m is None:
        return Response('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"/>', mimetype="image/svg+xml")
    d = rdMolDraw2D.MolDraw2DSVG(size, size)
    d.drawOptions().clearBackground = False
    d.DrawMolecule(m)
    d.FinishDrawing()
    return Response(d.GetDrawingText(), mimetype="image/svg+xml")


@nts_bp.route("/predict_ions", methods=["POST"])
def predict_ions():
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from hrms_nts.insilico import predict_product_ions, score_candidate, suggest_transitions
    body = request.get_json(force=True) or {}
    smiles = (body.get("smiles") or "").strip()
    m = Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:
        return jsonify({"error": "Invalid or empty SMILES"}), 400
    adducts = body.get("adducts") or ["[M+H]+"]
    top = int(body.get("top", 20))
    peaks = _parse_peaks(body.get("experimental_peaks") or "")
    out = {"smiles": smiles, "formula": rdMolDescriptors.CalcMolFormula(m), "exact_mass": Descriptors.ExactMolWt(m), "results": []}
    for ad in adducts:
        p = predict_product_ions(smiles, ad, top=top)
        if p.empty:
            out["results"].append({"adduct": ad, "error": "adduct not supported or no ions predicted"})
            continue
        r = {"adduct": ad, "precursor_mz": float(p["precursor_mz"].iloc[0]),
             "ions": _clean_records(p[["mz", "label", "kind", "rank_score"]]),
             "transitions": _clean_records(suggest_transitions(smiles, ad))}
        if peaks is not None:
            sc = score_candidate(smiles, ad, peaks[:, 0], peaks[:, 1], r["precursor_mz"])
            r["experimental"] = {"peaks": peaks.tolist(), **sc}
        out["results"].append(r)
    return jsonify(out)


def _parse_peaks(text):
    """'mz intensity' pairs, one per line (space, tab, comma or ':' separated)."""
    pts = []
    for ln in str(text).replace(";", "\n").splitlines():
        bits = ln.replace(",", " ").replace(":", " ").replace("\t", " ").split()
        if len(bits) >= 2:
            try:
                pts.append((float(bits[0]), float(bits[1])))
            except ValueError:
                continue
    return np.array(pts) if pts else None
