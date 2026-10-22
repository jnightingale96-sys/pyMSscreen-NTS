"""Command line:  python -m hrms_nts <command> [options]
  doctor | build-library | run | validate | prescreen | predict-ions | calibrate | convert
Global options: --config FILE   --set section.key=value (repeatable, YAML value)   --version"""
from __future__ import annotations
import argparse, glob, os, sys, logging, datetime, json, subprocess, shutil
import numpy as np
import pandas as pd

from . import __version__
from .common import load_config, resolve, setup_logging

log = logging.getLogger("hrms_nts")


def _files(inputs):
    out = []
    for p in inputs:
        if os.path.isdir(p):
            out += sorted(glob.glob(os.path.join(p, "*.mzML")) + glob.glob(os.path.join(p, "*.mzml")))
        elif any(ch in p for ch in "*?["):
            out += sorted(glob.glob(p))
        else:
            out.append(p)
    return list(dict.fromkeys(out))


def cmd_doctor(a, cfg):
    import importlib
    print(f"hrms_nts {__version__}  python {sys.version.split()[0]}")
    print("  config layers:", " <- ".join(cfg.get("_config_layers", [])))
    for m in ["numpy", "pandas", "scipy", "pyteomics", "lxml", "rdkit", "sklearn", "matplotlib", "openpyxl", "yaml"]:
        try:
            mod = importlib.import_module(m); print(f"  [ok]  {m} {getattr(mod, '__version__', '')}")
        except Exception as e:
            print(f"  [MISSING] {m}  ->  run install_nts.bat")
    for k in ("susdat", "pubchemlite", "user_compounds", "library_db", "standards", "calibrants", "metfrag_jar", "inhouse_msp"):
        p = resolve(cfg, cfg["paths"].get(k))
        print(f"  {'[ok]' if p and os.path.exists(p) else '[--]'}  paths.{k}: {p}")
    for p in cfg["paths"].get("library_msp") or []:
        print(f"  {'[ok]' if os.path.exists(resolve(cfg, p)) else '[--]'}  library_msp: {p}")
    print("  java:", shutil.which(cfg["paths"].get("java", "java")) or "NOT FOUND (MetFrag engine disabled; built-in predictor still works)")


def cmd_build_library(a, cfg):
    from .library import build_library
    msps = [resolve(cfg, p) for p in (cfg["paths"].get("library_msp") or [])]
    src = [os.path.splitext(os.path.basename(p))[0] for p in msps]
    if a.only:
        keep = [i for i, p in enumerate(msps) if any(o.lower() in p.lower() for o in a.only)]
        msps = [msps[i] for i in keep]; src = [src[i] for i in keep]
    ih = resolve(cfg, cfg["paths"].get("inhouse_msp"))
    if ih and os.path.exists(ih):
        msps.append(ih); src.append("inhouse")
    msps = [p for p in msps if os.path.exists(p)]
    db = resolve(cfg, cfg["paths"]["library_db"])
    n = build_library(msps, db, src, limit=a.limit, append=False)
    print(f"library built: {n} spectra -> {db}")


def cmd_run(a, cfg):
    from .pipeline import Context, process_sample
    from .report import write_all
    from .aggregate import blank_filter, align_features, presence_table
    files = _files(a.input)
    if not files:
        sys.exit("no .mzML files found (convert RAW files first: `convert` command or MSConvert/ThermoRawFileParser)")
    if a.semiquant: cfg["semiquant"]["mode"] = a.semiquant
    if a.insilico: cfg["insilico"]["engine"] = a.insilico
    if a.ppm: cfg["level4"]["ppm"] = a.ppm
    if a.min_intensity: cfg["ms1"]["min_intensity"] = a.min_intensity
    out = a.out or os.path.join(os.getcwd(), "NTS_out_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    os.makedirs(out, exist_ok=True)
    setup_logging(logfile=os.path.join(out, "run.log"))
    log.info("hrms_nts %s | %d file(s) | out=%s", __version__, len(files), out)
    ctx = Context(cfg, need_library=not a.no_library)
    results = []
    for f in files:
        try:
            r = process_sample(f, ctx, a.matrix, a.polarity, out, None, a.semiquant != "none")
            results.append(r)
        except Exception:
            log.exception("sample failed: %s", f)
    # blanks / alignment
    if results:
        allf = pd.concat([r["features"] for r in results if not r["features"].empty], ignore_index=True)
        blanks = cfg.get("blank", {}).get("samples") or []
        allf = blank_filter(allf, blanks, cfg["blank"]["ratio"])
        allf = align_features(allf)
        for r in results:
            sub = allf[allf["sample"] == r["sample"]][["feature_id", "consensus_id", "blank_flag", "blank_ratio"]]
            for k in ("features", "top"):
                if not r[k].empty:
                    r[k] = r[k].drop(columns=[c for c in ("consensus_id", "blank_flag", "blank_ratio") if c in r[k].columns]).merge(sub, on="feature_id", how="left")
        pres = presence_table(allf, "area")
        pres.reset_index().to_csv(os.path.join(out, "feature_presence_by_sample.csv"), index=False)
    from .provenance import write_provenance
    write_provenance(out, cfg, files, sys.argv)
    res = write_all(results, out, cfg, ctx, params=dict(version=__version__, files=";".join(os.path.basename(x) for x in files),
                    matrix=a.matrix, polarity=a.polarity, semiquant=cfg["semiquant"]["mode"], insilico=cfg["insilico"]["engine"],
                    ppm=cfg["level4"]["ppm"], rt_tol=cfg["rt"]["tolerance_min"], config=cfg["_config_path"]))
    print("\n=== LEVEL SUMMARY ===")
    print(res["level_summary"].to_string(index=False))
    print(f"\nresults: {out}\n  Excel: {res['xlsx']}")


def cmd_validate(a, cfg):
    from .validate import validate
    rep, summ = validate(a.results, a.expected, a.absent, a.rt_tol, a.out)
    pd.set_option("display.width", 200)
    print(rep[["name", "sample", "detected", "best_level", "ppm_error", "rt_diff", "top_candidate", "meets_min_level"]].round(3).to_string(index=False))
    print("\n=== VALIDATION SUMMARY ===")
    for k, v in summ.items():
        print(f"  {k:32s} {v}")
    print("\nwritten: validation_report.csv / validation_summary.csv in", a.out or a.results)


def cmd_prescreen(a, cfg):
    from .database import load_user_compounds
    from .prescreen import prescreen
    from . import mzml_io
    from .common import guess_polarity
    comp = load_user_compounds(a.compounds)
    files = _files(a.input); out = a.out or os.path.join(os.getcwd(), "prescreen_out"); os.makedirs(out, exist_ok=True)
    setup_logging(logfile=os.path.join(out, "run.log"))
    allr = []
    for f in files:
        run = mzml_io.read_mzml(f, None, cfg["ms1"]["store_min_intensity"], cache_dir=resolve(cfg, cfg["paths"].get("cache_dir")))
        run.polarity = a.polarity if a.polarity != "auto" else (guess_polarity(f, cfg.get("polarity_patterns")) or run.polarity or "positive")
        allr.append(prescreen(run, comp, a.adducts.split(",") if a.adducts else None, a.ppm, a.ms1_thresh, a.ms2_thresh, a.s2n, a.rt_tol))
    res = pd.concat(allr, ignore_index=True)
    res.to_csv(os.path.join(out, "summary_table.csv"), index=False)
    print(res.groupby("tag")["qa_pass"].agg(["sum", "size"]).rename(columns={"sum": "pass", "size": "tested"}))
    print("->", os.path.join(out, "summary_table.csv"))


def cmd_predict(a, cfg):
    from .insilico import predict_product_ions, suggest_transitions
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    smiles = a.smiles
    if a.csv:
        df = pd.read_csv(a.csv); smiles = None
    else:
        df = pd.DataFrame({"name": [a.name or smiles], "smiles": [smiles]})
    rows = []
    for r in df.itertuples():
        for ad in a.adduct.split(","):
            p = predict_product_ions(r.smiles, ad, top=a.top)
            if p.empty: continue
            p.insert(0, "name", r.name); p.insert(1, "adduct", ad); rows.append(p)
            print(f"\n{r.name}  {ad}  precursor m/z {p['precursor_mz'].iloc[0]:.4f}")
            print(p.head(a.top)[["mz", "label", "kind", "rank_score"]].to_string(index=False))
            print("suggested transitions:\n", suggest_transitions(r.smiles, ad).to_string(index=False))
    if rows and a.out:
        pd.concat(rows).to_csv(a.out, index=False); print("->", a.out)


def cmd_calibrate(a, cfg):
    from .semiquant import measure_calibrants, calibrants_table
    from . import mzml_io
    from .common import guess_polarity
    cal = calibrants_table(resolve(cfg, cfg["paths"]["calibrants"]) if not a.calibrants else a.calibrants)
    run = mzml_io.read_mzml(a.mzml, None, cfg["ms1"]["store_min_intensity"], cache_dir=resolve(cfg, cfg["paths"].get("cache_dir")))
    run.polarity = a.polarity if a.polarity != "auto" else (guess_polarity(a.mzml) or run.polarity or "positive")
    res = measure_calibrants(run, cal)
    print(res.to_string(index=False))
    if a.out:
        res.to_csv(a.out, index=False); print("->", a.out, "(copy the 'area' column into calibrants.csv / config single_reference.area)")


def cmd_convert(a, cfg):
    trfp = a.trfp or resolve(cfg, cfg["paths"].get("thermorawfileparser"))
    if not trfp or not os.path.exists(trfp):
        sys.exit(f"ThermoRawFileParser not found ({trfp}); set paths.thermorawfileparser, use --trfp, or convert with ProteoWizard msconvert (peak picking on)")
    os.makedirs(a.out, exist_ok=True)
    subprocess.run([trfp, "-d", a.input, "-o", a.out, "-f", "2", "-l", "1"], check=True)   # -f 2 = indexed mzML, centroid by default


def main(argv=None):
    ap = argparse.ArgumentParser(prog="hrms_nts", description="HRMS non-target screening: Level 5/4/3/2b/2a/1 + semi-quant")
    ap.add_argument("--config", help="config.yaml (default: next to the package; config.local.yaml is merged on top automatically)")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override a config value, e.g. --set level4.ppm=3 (repeatable)")
    ap.add_argument("--version", action="version", version=f"hrms_nts {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("doctor")
    p = sub.add_parser("build-library"); p.add_argument("--only", nargs="*"); p.add_argument("--limit", type=int)
    p = sub.add_parser("run")
    p.add_argument("input", nargs="+", help="mzML files / folders / globs")
    p.add_argument("--out"); p.add_argument("--matrix", default="water", help="water|grain|soil|extract (see config semiquant.matrices)")
    p.add_argument("--polarity", default="auto", choices=["auto", "positive", "negative"])
    p.add_argument("--semiquant", choices=["single", "logie", "curve", "none"]); p.add_argument("--insilico", choices=["builtin", "metfrag", "both", "none"])
    p.add_argument("--ppm", type=float); p.add_argument("--min-intensity", type=float); p.add_argument("--no-library", action="store_true")
    p = sub.add_parser("validate", help="score a results folder against expected.csv (and optionally absent.csv)")
    p.add_argument("results"); p.add_argument("expected"); p.add_argument("--absent"); p.add_argument("--rt-tol", type=float, default=0.3); p.add_argument("--out")
    p = sub.add_parser("prescreen"); p.add_argument("compounds"); p.add_argument("input", nargs="+"); p.add_argument("--out")
    p.add_argument("--polarity", default="auto"); p.add_argument("--adducts"); p.add_argument("--ppm", type=float, default=5)
    p.add_argument("--ms1-thresh", type=float, default=1e5); p.add_argument("--ms2-thresh", type=float, default=2500)
    p.add_argument("--s2n", type=float, default=3); p.add_argument("--rt-tol", type=float, default=0.5)
    p = sub.add_parser("predict-ions"); p.add_argument("--smiles"); p.add_argument("--name"); p.add_argument("--csv")
    p.add_argument("--adduct", default="[M+H]+"); p.add_argument("--top", type=int, default=15); p.add_argument("--out")
    p = sub.add_parser("calibrate"); p.add_argument("mzml"); p.add_argument("--calibrants"); p.add_argument("--polarity", default="auto"); p.add_argument("--out")
    p = sub.add_parser("convert"); p.add_argument("input"); p.add_argument("out"); p.add_argument("--trfp")
    a = ap.parse_args(argv)
    import yaml
    ov = {}
    for kv in a.set:
        k, _, v = kv.partition("=")
        d = ov
        keys = k.split(".")
        for kk in keys[:-1]:
            d = d.setdefault(kk, {})
        d[keys[-1]] = yaml.safe_load(v)
    cfg = load_config(a.config, ov)
    setup_logging()
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")      # database SMILES that RDKit rejects are skipped; per-structure parse errors only flood run.log
    {"doctor": cmd_doctor, "build-library": cmd_build_library, "run": cmd_run, "validate": cmd_validate, "prescreen": cmd_prescreen,
     "predict-ions": cmd_predict, "calibrate": cmd_calibrate, "convert": cmd_convert}[a.cmd](a, cfg)


if __name__ == "__main__":
    main()
