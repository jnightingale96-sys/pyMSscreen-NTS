"""Builds HRMS_NTS_Colab.ipynb (run: python colab/make_notebook.py). Edit the cells here, not in the .ipynb."""
import json, os

cells = []
md = lambda s: cells.append({"cell_type": "markdown", "metadata": {}, "source": s.strip("\n").splitlines(keepends=True)})
code = lambda s: cells.append({"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None, "outputs": [],
                               "source": s.strip("\n").splitlines(keepends=True)})

md("""
# HRMS non-target screening (NTS) + semi-quantification — pyMSscreen-NTS engine 0.5
Runs in Google Colab or locally in Jupyter / VS Code (box 0 picks the source).
Step by step, one box per step. Run the boxes in order; each step saves its results to the project folder, so you can stop,
close Colab and continue later from any step.

| Step | What it does |
|---|---|
| 0 | Install the engine + databases (NORMAN SusDat, PubChemLite, 970k-spectrum MS/MS library, NORMAN RTI) |
| 1 | Point to your data folder (.raw / .mzML) and a results folder on Google Drive (1b: or copy the data from a personal OneDrive) |
| 2 | Convert Thermo .raw to centroided mzML (SCIEX .wiff/.wiff2: convert with msconvert on Windows first) |
| 3 | Sample sheet: sample / blank / QC / calibration, replicate groups, calibration levels |
| 4 | QA/QC settings: replicate rule (e.g. 2/3), blank ratio, QC RSD, presence-confidence threshold, Q1/Q3 tolerances |
| 5 | Calibration: parent compounds, their internal standards, concentrations → linear or quadratic curves |
| 6 | Internal standards / targets in every injection: Q1 & Q3 errors, levels, integration, blanks |
| 7 | Non-target screening Level 5 → 1 (adducts, RT models, library + in-silico product ions) |
| 8 | QA: replicate presence, QC stability, blank rule, peak rating (0–10), presence confidence (%) |
| 9 | Semi-quantification by elucidation distance (logIE, RT, MW, logP from CompTox) |
| 10 | Chromatogram viewer: any compound × injection, Q1 and Q3 traces, integration, rating |
| 11 | Report (NTS_report.xlsx) and downloads |
""")

md("""## 0 · Setup (once per session, ~3 min in Colab)
Runs in **Google Colab** or **locally** (Jupyter / VS Code on this PC) - the same boxes either way.
* **GitHub** (Colab default): the engine is cloned from your repository; the databases + MetFrag jar (`hrms_nts_data.zip`) and,
  optionally, MetFrag's offline spectra (`metfrag_spectra.zip`) are downloaded from a Release. For a **private** repository add a
  GitHub token (fine-grained, read-only *Contents*) as a Colab secret named `GH_TOKEN` (key icon on the left).
* **Google Drive zip**: the complete package zip (`HRMS_NTS_Ultimate_*.zip`) on your Drive.
* **Local folder (this PC)**: the unzipped package (`...\\HRMS_NTS_Ultimate`) - uses its bundled Java + MetFrag; no session time limit,
  best for MetFrag runs (~30 min per sample).

**MetFrag** needs Java: in Colab tick `INSTALL_METFRAG` (installs OpenJDK, ~30 s). Free Colab sessions stop after ~12 h or when idle.""")
code(r'''
#@title 0 · Setup { display-mode: "form" }
SOURCE = "GitHub"  #@param ["GitHub", "Google Drive zip", "Local folder (this PC)"]
GITHUB_REPO = "jnightingale96-sys/pyMSscreen-NTS"  #@param {type:"string"}
GITHUB_BRANCH = "main"  #@param {type:"string"}
DATA_RELEASE_TAG = "data-v1"  #@param {type:"string"}
PACKAGE_ZIP = "/content/drive/MyDrive/HRMS_NTS/HRMS_NTS_Ultimate_20261102_v0.5.zip"  #@param {type:"string"}
LOCAL_PACKAGE = r"C:\Users\jnigh\Desktop\HRMS_NTS_Ultimate_build\HRMS_NTS_Ultimate"  #@param {type:"string"}
INSTALL_METFRAG = True  #@param {type:"boolean"}
import os, sys, zipfile, subprocess, glob, json, shutil, urllib.request
try:
    import google.colab
    IN_COLAB = True
except ImportError:
    IN_COLAB = False
if IN_COLAB:
    from google.colab import drive, output
    drive.mount("/content/drive")
    output.enable_custom_widget_manager()


def fetch_release_asset(repo, tag, name, dest, token=""):
    """Download one asset of a GitHub Release (works for private repositories with a token). Returns False if absent."""
    hdr = {"Authorization": f"Bearer {token}"} if token else {}
    rel = json.load(urllib.request.urlopen(urllib.request.Request(f"https://api.github.com/repos/{repo}/releases/tags/{tag}", headers=hdr)))
    asset = next((a for a in rel["assets"] if a["name"] == name), None)
    if asset is None:
        return False
    print(f"downloading {name} ({asset['size'] / 1e9:.2f} GB) ...")
    req = urllib.request.Request(asset["url"], headers={**hdr, "Accept": "application/octet-stream"})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 24)
    return True


OVERRIDES = {}
if SOURCE.startswith("Local"):
    ROOT = LOCAL_PACKAGE                       # its engine/config.local.yaml points to the bundled data, Java and MetFrag
elif SOURCE == "GitHub":
    try:
        from google.colab import userdata
        TOKEN = userdata.get("GH_TOKEN")
    except Exception:
        TOKEN = os.environ.get("GH_TOKEN", "")    # empty for a public repository
    ROOT = "/content/pyMSscreen-NTS" if IN_COLAB else os.path.abspath("pyMSscreen-NTS")
    if not os.path.exists(ROOT):
        auth = f"{TOKEN}@" if TOKEN else ""
        subprocess.run(["git", "clone", "-q", "--depth", "1", "-b", GITHUB_BRANCH, f"https://{auth}github.com/{GITHUB_REPO}.git", ROOT], check=True)
    data_dir = os.path.join(ROOT, "engine", "data"); os.makedirs(data_dir, exist_ok=True)
    tmp = os.path.join(ROOT, "_download.zip")
    if not os.path.exists(os.path.join(data_dir, "spectral_library.sqlite")):
        fetch_release_asset(GITHUB_REPO, DATA_RELEASE_TAG, "hrms_nts_data.zip", tmp, TOKEN)
        zipfile.ZipFile(tmp).extractall(data_dir); os.remove(tmp)
    if INSTALL_METFRAG and not os.path.isdir(os.path.join(data_dir, "mb_converted")):
        if fetch_release_asset(GITHUB_REPO, DATA_RELEASE_TAG, "metfrag_spectra.zip", tmp, TOKEN):
            zipfile.ZipFile(tmp).extractall(data_dir); os.remove(tmp)
    with open(os.path.join(ROOT, "engine", "config.local.yaml"), "w") as f:      # paths of the downloaded databases
        f.write('paths:\n  susdat: "data/susdat.csv"\n  pubchemlite: "data/PubChemLite_exposomics.csv"\n  library_msp: []\n'
                '  library_db: "data/spectral_library.sqlite"\n  metfrag_jar: "data/MetFragCommandLine-2.6.11.jar"\n'
                '  metfrag_mb_folder: "data/mb_converted"\npec_screen:\n  id_filter:\n    massbank_eu: "data/MassBankEU_compounds.csv"\n'
                '    norman_workbook: "data/Updated_NORMAN_with_TPSA.xlsx"\n')
else:
    ROOT = "/content/HRMS_NTS_Ultimate"
    if not os.path.exists(ROOT):
        print("unzipping the package (1.2 GB) ...")
        zipfile.ZipFile(PACKAGE_ZIP).extractall("/content")
if IN_COLAB:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "rdkit", "pyteomics", "psims", "plotly", "ipywidgets", "requests"], check=True)
    if INSTALL_METFRAG and shutil.which("java") is None:
        print("installing Java for MetFrag ...")
        subprocess.run("apt-get -qq update > /dev/null && apt-get -qq install -y openjdk-17-jre-headless > /dev/null", shell=True, check=True)
    # Linux: system Java instead of the bundled Windows runtime; RAW converter installed in box 2; caches on the local disk
    OVERRIDES = {"paths": {"java": "java", "thermorawfileparser": None, "cache_dir": "/content/nts_cache"}}
ENGINE = os.path.join(ROOT, "engine")
if ENGINE not in sys.path:
    sys.path.insert(0, ENGINE)
from hrms_nts import __version__
from hrms_nts.workflow import Project, install_thermorawfileparser
from hrms_nts import colab_ui as ui
from hrms_nts.insilico import metfrag_available
from hrms_nts.common import load_config
mf_ok = metfrag_available(load_config(None, OVERRIDES))[0]
print(f"hrms_nts {__version__} ready ({'Colab' if IN_COLAB else 'local'}); MetFrag {'available' if mf_ok else 'NOT available - use IN_SILICO = builtin'}")
''')

md("## 1 · Data and results folders\n`DATA_DIR` holds the .raw and/or .mzML files (sub-folders are searched). Results go to `OUT_DIR` — keep it on Drive so nothing is lost when Colab disconnects.")
code(r'''
#@title 1 · Folders { display-mode: "form" }
DATA_DIR = "/content/drive/MyDrive/HRMS_NTS/my_batch"  #@param {type:"string"}
OUT_DIR = "/content/drive/MyDrive/HRMS_NTS/my_batch_results"  #@param {type:"string"}
P = Project(DATA_DIR, OUT_DIR, overrides=OVERRIDES)
ins = P.inputs()
print(f"{len(ins['raw'])} Thermo .raw, {len(ins['sciex'])} SCIEX .wiff/.wiff2, {len(ins['mzml'])} mzML")
for f in ins["mzml"][:30]: print("  ", os.path.basename(f))
''')

md("""## 1b · Data from a personal OneDrive (optional - skip if your data are on Google Drive)
Copies a OneDrive folder into this Colab machine with **rclone** and makes it the data folder (results still go to `OUT_DIR`).
One-off on your PC: download rclone (rclone.org/downloads, Windows zip), open a command prompt in its folder and run
`rclone authorize "onedrive"` - a browser opens, sign in with your **personal** Microsoft account, then copy everything between
`--->` and `<---` (one line starting `{"access_token"`). In Colab add it as a secret named `RCLONE_ONEDRIVE_TOKEN` (key icon).
University / work OneDrive accounts usually need IT approval for rclone - use a personal account or Google Drive.""")
code(r'''
#@title 1b · Copy data from OneDrive { display-mode: "form" }
USE_ONEDRIVE = False  #@param {type:"boolean"}
ONEDRIVE_FOLDER = "NTS_CoLab_Test"  #@param {type:"string"}
if USE_ONEDRIVE:
    import json, urllib.request, shutil, subprocess
    if shutil.which("rclone") is None:
        subprocess.run("curl -s https://rclone.org/install.sh | sudo bash > /dev/null", shell=True, check=True)
    try:
        from google.colab import userdata
        tok = userdata.get("RCLONE_ONEDRIVE_TOKEN")
    except Exception:
        tok = os.environ.get("RCLONE_ONEDRIVE_TOKEN", "")
    t = json.loads(tok)
    me = json.load(urllib.request.urlopen(urllib.request.Request("https://graph.microsoft.com/v1.0/me/drive",
                                                                 headers={"Authorization": f"Bearer {t['access_token']}"})))
    subprocess.run(["rclone", "config", "create", "onedrive", "onedrive", "token", tok, "drive_id", me["id"],
                    "drive_type", me["driveType"], "--non-interactive"], check=True, stdout=subprocess.DEVNULL)
    LOCAL_DATA = os.path.join("/content/onedrive_data", os.path.basename(ONEDRIVE_FOLDER.rstrip("/")))
    print(f"copying onedrive:{ONEDRIVE_FOLDER} -> {LOCAL_DATA} ...")
    subprocess.run(["rclone", "copy", f"onedrive:{ONEDRIVE_FOLDER}", LOCAL_DATA, "--include", "*.{mzML,mzml,raw,RAW}",
                    "--transfers", "4", "--stats-one-line", "--stats", "30s", "-v"], check=True)
    DATA_DIR = LOCAL_DATA
    P = Project(DATA_DIR, OUT_DIR, overrides=OVERRIDES)
    ins = P.inputs()
    print(f"{len(ins['raw'])} Thermo .raw, {len(ins['mzml'])} mzML now in {DATA_DIR}")
else:
    print("OneDrive not used - data come from DATA_DIR in box 1")
''')

md("## 2 · Convert .raw → mzML (skip if you already have mzML)\nThermo .raw files are converted with ThermoRawFileParser (centroided, indexed mzML). SCIEX .wiff/.wiff2 cannot be read on Linux: convert them on Windows with ProteoWizard msconvert (*Peak picking: Vendor, MS levels 1-*) and put the mzML in `DATA_DIR`. Profile mzML is fine — the engine centroids it.")
code(r'''
#@title 2 · Convert Thermo .raw { display-mode: "form" }
if P.inputs()["raw"]:
    if IN_COLAB:
        subprocess.run("apt-get -qq install -y mono-runtime > /dev/null", shell=True)    # ThermoRawFileParser .NET Framework builds
        new = P.convert_raw(install_thermorawfileparser("/content/trfp"))
    else:
        new = P.convert_raw()                                                                # bundled ThermoRawFileParser.exe
    print(f"converted {len(new)} file(s)")
else:
    print("no .raw files - nothing to convert")
print(len(P.mzml_files()), "mzML files available")
''')

md("## 3 · Sample sheet\nFor every injection choose the **type** (`sample`, `blank`, `qc` = pooled QC, `cal` = calibration standard, `exclude`), the **group** (replicate injections of the same sample share a group — the replicate rule is applied per group) and, for `cal`, the **concentration**. The table is pre-filled from the file names; check it and press *Save*. Alternatively upload a CSV with columns `file,type,group,replicate,conc` and set `SHEET_CSV`.")
code(r'''
#@title 3 · Sample sheet { display-mode: "form" }
SHEET_CSV = ""  #@param {type:"string"}
if SHEET_CSV:
    print(P.set_sample_sheet(SHEET_CSV))
else:
    ui.sample_sheet_editor(P)
''')

md("## 4 · QA/QC settings")
code(r'''
#@title 4 · QA/QC settings { display-mode: "form" }
REPLICATE_RULE = "2/3"  #@param {type:"string"}
BLANK_RATIO = 3.0  #@param {type:"number"}
QC_MAX_RSD_PCT = 30.0  #@param {type:"number"}
PRESENCE_CONFIDENCE_PASS_PCT = 70.0  #@param {type:"number"}
Q1_PPM = 5.0  #@param {type:"number"}
Q3_PPM = 10.0  #@param {type:"number"}
MIN_LEVEL_REPORTED = "3"  #@param ["1", "2a", "2b", "3", "4"]
P.cfg["blank"]["ratio"] = BLANK_RATIO
P.cfg["level4"]["ppm"] = Q1_PPM
print(f"replicate rule {REPLICATE_RULE}: a compound is present in a sample group when found in at least that share of its injections")
''')

md("## 5 · Calibration (needs `cal` injections in the sample sheet)\nEnter the **parent compounds**, their **internal standards** (role `IS`, with the label e.g. `D10`) and for each parent the IS it is normalised to. Concentrations come from the sample sheet. Each curve is fitted **linear and quadratic** (1/x weighted); the quadratic is used only when the curvature is significant (F-test, AICc) and monotonic. The response factors also feed the semi-quantification (step 9).")
code(r'''
#@title 5a · Calibration compounds { display-mode: "form" }
CAL_CSV = ""  #@param {type:"string"}
CAL_UNITS = "ng/mL"  #@param {type:"string"}
import pandas as pd
get_cal_table = (lambda: pd.read_csv(CAL_CSV)) if CAL_CSV else ui.calibration_editor(P)
''')
code(r'''
#@title 5b · Fit calibration curves { display-mode: "form" }
WEIGHTING = "1/x"  #@param ["1/x", "1/x2", "none"]
MIN_LEVELS_QUADRATIC = 5  #@param {type:"integer"}
cal = P.calibrate(get_cal_table(), WEIGHTING, MIN_LEVELS_QUADRATIC, CAL_UNITS)
display(cal["summary"][[c for c in ["name", "model", "reason", "n_levels", "r2", "r2_linear", "r2_quadratic", "lod", "loq", "accuracy_pct", "warning"] if c in cal["summary"]]])
ui.show_calibration(cal)
print("PDF:", cal["pdf"])
''')

md("## 6 · Internal standards / targets in every injection (optional)\nA list with `name` + `formula` (labels as `[2H]`, `[13C]`) or `smiles` + `label`; vendor IS sheets (Name, Chemical formula, Q1) work as they are. Gives Q1 and Q3 errors, level, integration, blank check and IS QC (area RSD, response drift per injection).")
code(r'''
#@title 6 · Targeted IS screen { display-mode: "form" }
IS_LIST = ""  #@param {type:"string"}
if IS_LIST:
    tg = P.targeted(IS_LIST, Q1_PPM, Q3_PPM)
    display(tg["is_qc"]); display(tg["injection_qc"])
    print("workbook:", tg["xlsx"])
else:
    print("no IS list given - skipped")
''')

md("## 7 · Non-target screening (Level 5 → 1)\nEvery sample, QC and blank: features, adducts ([M+H]+, [M+Na]+, [M+NH4]+, [M+K]+ / [M-H]-, [M+Cl]-, [M+HCOO]-), NORMAN SusDat + PubChemLite candidates, RT plausibility (RT models from `standards.csv` + NORMAN RTI), library MS2 (970k spectra; DIA/SWATH via correlated pseudo-MS2), built-in in-silico product ions. ~2–5 min per injection.")
code(r'''
#@title 7 · Screening { display-mode: "form" }
MATRIX = "water"  #@param ["water", "grain", "soil", "extract"]
IN_SILICO = "builtin"  #@param ["builtin", "both", "metfrag", "none"]
# both = built-in + MetFrag (~30 min per injection). MetFrag-only Level 2b counts as tentative unless the built-in scorer agrees (decoy-corrected).
sc = P.screen(MATRIX, IN_SILICO)
display(sc["level_summary"])
''')

md("## 8 · QA: replicates, QC, blanks, peak rating, presence confidence\n**Peak rating** 0–10 (points across the peak, FWHM/base shape, jaggedness, modality, S/N — Compound-Discoverer style). **Presence confidence** % combines the peak rating, S/N, Q1 accuracy, isotope pattern, MS2 evidence, replicate frequency and blank ratio; the pass threshold is set in step 4.")
code(r'''
#@title 8 · QA { display-mode: "form" }
q = P.qa(REPLICATE_RULE, QC_MAX_RSD_PCT, MIN_LEVEL_REPORTED, Q1_PPM, PRESENCE_CONFIDENCE_PASS_PCT)
I = q["identifications"]
print(f"{len(I)} detections at Level <= {MIN_LEVEL_REPORTED}; {int(I.presence_pass.sum())} pass {PRESENCE_CONFIDENCE_PASS_PCT:g} % presence confidence")
cols = ["sample", "name", "level", "adduct", "rt_apex", "ppm_error", "lib_score", "bi_expl_adj_n", "peak_rating", "presence_confidence_pct", "blank_flag"]
display(I[[c for c in cols if c in I.columns]].sort_values(["level", "presence_confidence_pct"], ascending=[True, False]).head(50))
if len(q["presence"]): display(q["presence"].head(30))
''')

md("## 9 · Semi-quantification (elucidation distance)\nEach identified compound is matched to the nearest calibrant in z-scored (logIE, RT, MW, logP) space; conc = area / RF × 10^(logIE_cal − logIE). logP comes from the **CompTox Dashboard** (experimental median, else predicted) — get a free CTX API key from the EPA CTX API team and store it as a Colab secret named `CTX_API_KEY` (key icon on the left), or upload a Dashboard batch-search export. ED < 0.5 = close analogue. Concentrations carry a ×/÷ factor (default 10).")
code(r'''
#@title 9 · Semi-quantification { display-mode: "form" }
MIN_LEVEL_SEMIQUANT = "2b"  #@param ["1", "2a", "2b", "3"]
COMPTOX_BATCH_EXPORT = ""  #@param {type:"string"}
MATRIX_FACTOR = 1.0  #@param {type:"number"}
UNIT = "ng/mL extract"  #@param {type:"string"}
try:
    from google.colab import userdata
    KEY = userdata.get("CTX_API_KEY")
except Exception:
    KEY = os.environ.get("CTX_API_KEY", "")
sq = P.semiquant(KEY, COMPTOX_BATCH_EXPORT or None, 0.5, MIN_LEVEL_SEMIQUANT, MATRIX_FACTOR, UNIT)
t = sq["table"]
display(t[[c for c in ["sample", "name", "level", "sq_conc", "sq_low", "sq_high", "elucidation_distance", "sq_ed_ok", "sq_calibrant", "sq_logp_source"] if c in t.columns]].head(50))
''')

md("## 10 · Chromatogram viewer\nPick a compound and an injection: **Q1** = precursor ion trace (integrated region shaded), **Q3** = product-ion traces (observed in the MS2 of this peak; dotted = in-silico predicted when no MS2 is available). Shows area, S/N, peak rating, Q1 error, presence confidence and ion ratios.")
code(r'''
#@title 10 · Chromatograms { display-mode: "form" }
ui.chromatogram_explorer(P, MIN_LEVEL_REPORTED)
''')

md("## 11 · Report\n`NTS_report.xlsx`: identified compounds (level, adducts, RT measured/predicted, product-ion evidence: library / in-silico predicted / MetFrag, peak rating, presence confidence, replicate presence, QC, blanks, semi-quant + ED), detections per sample, presence by group, QC, IS QC, calibration, quantification, properties.")
code(r'''
#@title 11 · Report { display-mode: "form" }
path = P.report()
print("report:", path)
DOWNLOAD = True  #@param {type:"boolean"}
if DOWNLOAD and IN_COLAB:
    from google.colab import files
    files.download(path)
''')

md("## 11b · Copy the results to OneDrive (optional)\nNeeds box 1b to have been run in this session (rclone set up).")
code(r'''
#@title 11b · Results to OneDrive { display-mode: "form" }
ONEDRIVE_RESULTS_FOLDER = "NTS_CoLab_Test_results"  #@param {type:"string"}
if shutil.which("rclone") and "onedrive" in subprocess.run(["rclone", "listremotes"], capture_output=True, text=True).stdout:
    subprocess.run(["rclone", "copy", OUT_DIR, f"onedrive:{ONEDRIVE_RESULTS_FOLDER}", "--exclude", "_*.pkl",
                    "--stats-one-line", "--stats", "30s", "-v"], check=True)
    print(f"results copied to onedrive:{ONEDRIVE_RESULTS_FOLDER}")
else:
    print("rclone / OneDrive not set up - run box 1b first")
''')


nb = {"nbformat": 4, "nbformat_minor": 0, "metadata": {"colab": {"provenance": [], "toc_visible": True},
                                                        "kernelspec": {"name": "python3", "display_name": "Python 3"},
                                                        "language_info": {"name": "python"}}, "cells": cells}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "HRMS_NTS_Colab.ipynb")
json.dump(nb, open(out, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print("wrote", out, len(cells), "cells")
