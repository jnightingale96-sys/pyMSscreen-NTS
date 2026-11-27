"""Builds HRMS_NTS_Colab.ipynb (run: python colab/make_notebook.py). Edit the cells here, not in the .ipynb."""
import json, os

cells = []
md = lambda s: cells.append({"cell_type": "markdown", "metadata": {}, "source": s.strip("\n").splitlines(keepends=True)})
code = lambda s: cells.append({"cell_type": "code", "metadata": {"cellView": "form"}, "execution_count": None, "outputs": [],
                               "source": s.strip("\n").splitlines(keepends=True)})

md("""
# HRMS non-target screening (NTS) + semi-quantification — pyMSscreen-NTS engine 0.6
Runs in Google Colab or locally in Jupyter / VS Code (box 0 picks the source).
Step by step, one box per step. Run the boxes in order; each step saves its results to the project folder, so you can stop,
close Colab and continue later from any step.

| Step | What it does |
|---|---|
| 0 | Install the engine + databases (NORMAN SusDat, PubChemLite, 970k-spectrum MS/MS library, NORMAN RTI) |
| 1 | Point to your data folder (.raw / .mzML) and a results folder on Google Drive (1c: upload files from this computer) |
| 2 | Convert Thermo .raw to centroided mzML (SCIEX .wiff/.wiff2: convert with msconvert on Windows first) |
| 3 | Templates (sample sheet, calibration list, IS list) to download and fill in; then the sample sheet |
| 4 | QA/QC settings: replicate rule (e.g. 2/3), blank ratio, QC RSD, presence-confidence threshold, Q1/Q3 tolerances |
| 5 | Calibration: parent compounds, their internal standards, concentrations → linear or quadratic curves |
| 6 | Internal standards / targets in every injection: Q1 & Q3 errors, levels, integration, blanks; 6b matrix effect (solvent vs spiked extract) |
| 7 | Non-target screening Level 5 → 1 (adducts, RT models, library + in-silico product ions) |
| 8 | QA: replicate presence, QC stability, blank rule, peak rating (0–10), presence confidence (%) |
| 9 | Physico-chemical properties from CompTox (automatic batch search), then semi-quantification from the internal standards |
| 10 | Chromatogram viewer: any compound × injection, Q1 and Q3 traces, integration, rating |
| 11 | Report: NTS_report.html (summary with figures) + NTS_report.xlsx (all tables), downloads |
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
DATA_RELEASE_TAG = "data-v2"  #@param {type:"string"}
ZENODO_RECORD = ""  #@param {type:"string"}
PACKAGE_ZIP = "/content/drive/MyDrive/HRMS_NTS/HRMS_NTS_Ultimate_20261102_v0.5.zip"  #@param {type:"string"}
LOCAL_PACKAGE = r"C:\path\to\HRMS_NTS_Ultimate"  #@param {type:"string"}
INSTALL_METFRAG = True  #@param {type:"boolean"}
import os, sys, zipfile, subprocess, glob, json, shutil, urllib.request, urllib.error
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
    """Download one asset of a GitHub Release (works for private repositories with a token). Returns False if absent.
    Falls back to the older release 'data-v1' when the requested tag or asset does not exist."""
    hdr = {"Authorization": f"Bearer {token}"} if token else {}
    try:
        rel = json.load(urllib.request.urlopen(urllib.request.Request(f"https://api.github.com/repos/{repo}/releases/tags/{tag}", headers=hdr)))
    except urllib.error.HTTPError as e:
        if e.code == 404 and tag != "data-v1":
            print(f"release {tag} not found - using data-v1")
            return fetch_release_asset(repo, "data-v1", name, dest, token)
        raise
    if tag != "data-v1" and not any(a["name"] == name for a in rel["assets"]):
        return fetch_release_asset(repo, "data-v1", name, dest, token)
    asset = next((a for a in rel["assets"] if a["name"] == name), None)
    if asset is None:
        return False
    print(f"downloading {name} ({asset['size'] / 1e9:.2f} GB) ...")
    req = urllib.request.Request(asset["url"], headers={**hdr, "Accept": "application/octet-stream"})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 24)
    return True


def fetch_zenodo_file(record, name, dest):
    """Download one file of a public Zenodo record and check its MD5. Returns False if the record has no such file."""
    import hashlib
    rec = json.load(urllib.request.urlopen(f"https://zenodo.org/api/records/{record}"))
    f = next((x for x in rec.get("files", []) if x["key"] == name), None)
    if f is None:
        return False
    print(f"downloading {name} ({f['size'] / 1e9:.2f} GB) from Zenodo record {record} ...")
    md5 = hashlib.md5()
    with urllib.request.urlopen(f["links"]["self"]) as r, open(dest, "wb") as out:
        for chunk in iter(lambda: r.read(1 << 24), b""):
            md5.update(chunk); out.write(chunk)
    if f.get("checksum", "md5:").split(":", 1)[1] not in ("", md5.hexdigest()):
        os.remove(dest)
        raise SystemExit(f"{name}: checksum mismatch - the download is incomplete; run box 0 again")
    return True


def fetch_data(name, dest, token=""):
    """Databases come from the Zenodo record when ZENODO_RECORD is set, otherwise from the GitHub Release."""
    record = ZENODO_RECORD.strip().rstrip("/").split("/")[-1]
    if record:
        return fetch_zenodo_file(record, name, dest)
    return fetch_release_asset(GITHUB_REPO, DATA_RELEASE_TAG, name, dest, token)


OVERRIDES = {}
if SOURCE.startswith("Local"):
    ROOT = LOCAL_PACKAGE                       # its engine/config.local.yaml points to the bundled data, Java and MetFrag
elif SOURCE == "GitHub":
    try:
        from google.colab import userdata
        TOKEN = userdata.get("GH_TOKEN")
    except Exception:
        TOKEN = os.environ.get("GH_TOKEN", "")    # empty for a public repository
    TOKEN = "".join((TOKEN or "").split())           # pasted secrets often carry a trailing line break / spaces
    ROOT = "/content/pyMSscreen-NTS" if IN_COLAB else os.path.abspath("pyMSscreen-NTS")
    if not os.path.exists(os.path.join(ROOT, "engine")):
        shutil.rmtree(ROOT, ignore_errors=True)
        auth = f"x-access-token:{TOKEN}@" if TOKEN else ""
        r = subprocess.run(["git", "clone", "-q", "--depth", "1", "-b", GITHUB_BRANCH, f"https://{auth}github.com/{GITHUB_REPO}.git", ROOT],
                           capture_output=True, text=True)
        if r.returncode != 0:                         # never echo the token
            raise SystemExit("Cloning the GitHub repository failed: " + r.stderr.replace(TOKEN, "***").strip() +
                             "\nCheck GITHUB_REPO, that the GH_TOKEN secret has Notebook access, and that the token can read this repository "
                             "(Contents: Read-only) and has not expired.")
        subprocess.run(["git", "-C", ROOT, "remote", "set-url", "origin", f"https://github.com/{GITHUB_REPO}.git"])   # token not kept on disk
    else:                                             # already cloned in this runtime: update to the latest code (databases are kept)
        auth = f"x-access-token:{TOKEN}@" if TOKEN else ""
        r = subprocess.run(["git", "-C", ROOT, "fetch", "-q", "--depth", "1", f"https://{auth}github.com/{GITHUB_REPO}.git", GITHUB_BRANCH],
                           capture_output=True, text=True)
        if r.returncode == 0:
            subprocess.run(["git", "-C", ROOT, "reset", "-q", "--hard", "FETCH_HEAD"], check=True)
            print("engine code updated from GitHub")
        else:
            print("could not update the engine code (" + r.stderr.replace(TOKEN, "***").strip()[:200] + ") - using the existing copy")
        for m in [m for m in list(sys.modules) if m == "hrms_nts" or m.startswith("hrms_nts.")]:
            del sys.modules[m]                        # re-import the updated engine below
    data_dir = os.path.join(ROOT, "engine", "data"); os.makedirs(data_dir, exist_ok=True)
    tmp = os.path.join(ROOT, "_download.zip")
    # data-v2 adds the free-database reference index, TP pairs, NORMAN SLE candidates and the SMRT RT model
    if not os.path.exists(os.path.join(data_dir, "spectral_library.sqlite")) or             (DATA_RELEASE_TAG != "data-v1" and not os.path.exists(os.path.join(data_dir, "reference_index.csv.gz"))):
        fetch_data("hrms_nts_data.zip", tmp, TOKEN)
        zipfile.ZipFile(tmp).extractall(data_dir); os.remove(tmp)
    if INSTALL_METFRAG and not os.path.isdir(os.path.join(data_dir, "mb_converted")):
        if fetch_data("metfrag_spectra.zip", tmp, TOKEN):
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
if "Project" not in globals():
    raise SystemExit("Box 0 has not finished in this session - run box 0 first (it must end with 'hrms_nts ... ready').")
P = Project(DATA_DIR, OUT_DIR, overrides=OVERRIDES)
ins = P.inputs()
print(f"{len(ins['raw'])} Thermo .raw, {len(ins['sciex'])} SCIEX .wiff/.wiff2, {len(ins['mzml'])} mzML")
for f in ins["mzml"][:30]: print("  ", os.path.basename(f))
''')

md("""## 1b · Data from a personal OneDrive (optional, advanced - most users: Google Drive or box 1c)
**Note:** Colab has shown a policy warning ("code that is disallowed") when rclone runs. Prefer Google Drive or box 1c.
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
        subprocess.run("apt-get -qq install -y rclone > /dev/null", shell=True, check=True)
    import re
    try:
        from google.colab import userdata
        tok = userdata.get("RCLONE_ONEDRIVE_TOKEN")
    except Exception as e:
        raise SystemExit(f"Cannot read the Colab secret RCLONE_ONEDRIVE_TOKEN ({type(e).__name__}): add it under the key icon with exactly "
                         "that name and switch 'Notebook access' on, then run this box again.")
    m = re.search(r"\{.*\}", tok or "", re.S)               # tolerate the ---> / <--- arrows or extra spaces around the token
    if not m:
        raise SystemExit("RCLONE_ONEDRIVE_TOKEN is empty or is not the {access_token: ...} line printed by rclone authorize onedrive - "
                         "copy that whole line (from { to }) into the secret.")
    tok = "".join(m.group(0).split())     # copying from a Windows console inserts line breaks where it wrapped; the token has none
    try:
        t = json.loads(tok)
    except json.JSONDecodeError:
        t = None
    if t is None:
        raise SystemExit("RCLONE_ONEDRIVE_TOKEN is incomplete or damaged. Copying from the black console window often breaks it: on your PC run\n"
                         "    rclone authorize \"onedrive\" > %USERPROFILE%\\Desktop\\onedrive_token.txt\n"
                         "open that file in Notepad, Ctrl+A, Ctrl+C, paste it as the secret value, then delete the file.")
    ONEDRIVE_FOLDER = ONEDRIVE_FOLDER.split("=", 1)[-1].strip().strip("/")      # tolerate 'ONEDRIVE_FOLDER = ...' typed into the box
    me = json.load(urllib.request.urlopen(urllib.request.Request("https://graph.microsoft.com/v1.0/me/drive",
                                                                 headers={"Authorization": f"Bearer {t['access_token']}"})))
    subprocess.run(["rclone", "config", "create", "onedrive", "onedrive", "token", tok, "drive_id", me["id"],
                    "drive_type", me["driveType"], "--non-interactive"], check=True, stdout=subprocess.DEVNULL)
    # OneDrive web links show the drive itself as '/Documents/', so 'Documents/X' may really be 'X' at the top level
    exists = lambda f: subprocess.run(["rclone", "lsf", "--max-depth", "1", f"onedrive:{f}"], capture_output=True).returncode == 0
    candidates = [ONEDRIVE_FOLDER] + ([ONEDRIVE_FOLDER.split("/", 1)[1]] if ONEDRIVE_FOLDER.lower().startswith("documents/") else [])
    found = next((f for f in candidates if exists(f)), None)
    if found is None:
        top = subprocess.run(["rclone", "lsd", "onedrive:"], capture_output=True, text=True).stdout
        names = [ln.split(None, 4)[-1] for ln in top.splitlines() if ln.strip()]
        raise SystemExit(f"OneDrive folder '{ONEDRIVE_FOLDER}' not found. Folders at the top of your OneDrive:\n  " + "\n  ".join(names) +
                         "\nPut the folder path (as in this list, sub-folders with /) into ONEDRIVE_FOLDER and run this box again.")
    ONEDRIVE_FOLDER = found
    LOCAL_DATA = os.path.join("/content/onedrive_data", os.path.basename(ONEDRIVE_FOLDER.rstrip("/")))
    print(f"copying onedrive:{ONEDRIVE_FOLDER} -> {LOCAL_DATA} ...")
    r = subprocess.run(["rclone", "copy", f"onedrive:{ONEDRIVE_FOLDER}", LOCAL_DATA, "--include", "*.{mzML,mzml,raw,RAW}",
                        "--transfers", "4", "--stats-one-line", "--stats", "30s", "-v"])
    if r.returncode != 0:
        raise SystemExit(f"rclone copy failed (exit {r.returncode}) - see the rclone messages above")
    DATA_DIR = LOCAL_DATA
    P = Project(DATA_DIR, OUT_DIR, overrides=OVERRIDES)
    ins = P.inputs()
    print(f"{len(ins['raw'])} Thermo .raw, {len(ins['mzml'])} mzML now in {DATA_DIR}")
else:
    print("OneDrive not used - data come from DATA_DIR in box 1")
''')

md("""## 1c · Upload files from this computer (optional)
Opens your browser's file picker; the chosen .mzML / .raw files are saved into `DATA_DIR` on Google Drive. Fine for a few files
(< ~2 GB in total); for large batches copy them into the Drive folder with the Google Drive app or drive.google.com.""")
code(r'''
#@title 1c · Upload mzML / raw files { display-mode: "form" }
UPLOAD_NOW = False  #@param {type:"boolean"}
if UPLOAD_NOW:
    ui.upload_files(DATA_DIR, accept=(".mzml", ".raw"))
    P = Project(DATA_DIR, OUT_DIR, overrides=OVERRIDES)
ins = P.inputs()
print(f"{len(ins['raw'])} Thermo .raw, {len(ins['mzml'])} mzML in {DATA_DIR}")
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

md("""## 3 · Templates and sample sheet
**3a** writes three fill-in templates (with example rows and a README) to `OUT_DIR/templates` and downloads them:
`sample_sheet_template.csv`, `calibration_list_template.csv` (one column per calibration level, L1 ... L8, so every compound can
have its own concentrations) and `is_list_template.csv` (name, formula, monoisotopic mass to 5 decimals, label, concentration).
**Sample types:** `sample`, `blank`, `qc` (pooled QC), `cal` (calibration level - give `cal_level` L1 ... L8), `is_mix` (IS mix in
solvent: response factors for semi-quantification), `solvent_std` and `matrix_spike` (the same mix in 10 % MeOH / water and spiked
into a control extract: matrix effect, box 6b), `exclude`. Replicate injections of one sample share a `group`.""")
code(r'''
#@title 3a · Download templates { display-mode: "form" }
paths = P.templates()
print(open(paths["README_templates.txt"], encoding="utf-8").read())
DOWNLOAD_TEMPLATES = True  #@param {type:"boolean"}
if DOWNLOAD_TEMPLATES:
    ui.download(*[v for k, v in paths.items() if k.endswith(".csv")])
''')
code(r'''
#@title 3b · Sample sheet { display-mode: "form" }
#@markdown Upload your filled-in sample sheet (tick UPLOAD_SHEET), give its path, or leave both empty to edit the table below.
UPLOAD_SHEET = False  #@param {type:"boolean"}
SHEET_CSV = ""  #@param {type:"string"}
if UPLOAD_SHEET:
    SHEET_CSV = (ui.upload_files(os.path.join(OUT_DIR, "inputs"), accept=(".csv", ".xlsx")) or [""])[0]
if SHEET_CSV:
    s = P.set_sample_sheet(SHEET_CSV)
    print(s.type.value_counts().to_string()); display(s)
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

md("## 5 · Calibration (needs `cal` injections in the sample sheet)\nUse `calibration_list_template.csv` (box 3a): one row per **target** and per **internal standard** (role `IS`, label e.g. `D10`), the IS each target is normalised to (`is_name`), the formula, the monoisotopic mass (5 decimals; checked against the formula) and the **concentration of that compound at every level** (columns L1 ... L8; the sample sheet says which injection is which level). Upload it below or use the editor. Each curve is fitted **linear and quadratic** (1/x weighted); the quadratic is used only when the curvature is significant (F-test, AICc) and monotonic. The response factors also feed the semi-quantification (step 9).")
code(r'''
#@title 5a · Calibration compounds { display-mode: "form" }
UPLOAD_CAL = False  #@param {type:"boolean"}
CAL_CSV = ""  #@param {type:"string"}
CAL_UNITS = "ng/mL"  #@param {type:"string"}
import pandas as pd
if UPLOAD_CAL:
    CAL_CSV = (ui.upload_files(os.path.join(OUT_DIR, "inputs"), accept=(".csv", ".xlsx")) or [""])[0]
if CAL_CSV:
    chk = P.check_list(CAL_CSV, "cal")
    display(chk[[c for c in ["name", "formula", "monoisotopic_mass", "mass_from_formula", "role", "is_name", "check"] if c in chk]])
    if (chk.check != "OK").any(): print("Fix the rows whose check is not OK, then upload again.")
    get_cal_table = lambda: pd.read_csv(CAL_CSV) if CAL_CSV.endswith(".csv") else pd.read_excel(CAL_CSV)
else:
    get_cal_table = ui.calibration_editor(P)
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

md("## 6 · Internal standards / targets in every injection\nUse `is_list_template.csv` (box 3a): name, formula (labels as `[2H]`, `[13C]`), monoisotopic mass (5 decimals, checked), SMILES, label and the **concentration in the IS mix** (needed for semi-quantification). Vendor IS sheets (Name, Chemical formula, Q1) also work. Gives Q1 and Q3 errors, level, integration, blank check and IS QC (area RSD, response drift per injection).")
code(r'''
#@title 6 · Targeted IS screen { display-mode: "form" }
UPLOAD_IS_LIST = False  #@param {type:"boolean"}
IS_LIST = ""  #@param {type:"string"}
if UPLOAD_IS_LIST:
    IS_LIST = (ui.upload_files(os.path.join(OUT_DIR, "inputs"), accept=(".csv", ".xlsx")) or [""])[0]
if IS_LIST:
    chk = P.check_list(IS_LIST, "is")
    bad = chk[chk.check != "OK"]
    print(f"list check: {len(chk) - len(bad)} of {len(chk)} rows OK")
    if len(bad): display(bad[[c for c in ["name", "formula", "monoisotopic_mass", "mass_from_formula", "check"] if c in bad]])
    tg = P.targeted(IS_LIST, Q1_PPM, Q3_PPM)
    display(tg["is_qc"]); display(tg["injection_qc"])
    print("workbook:", tg["xlsx"])
else:
    print("no IS list given - skipped")
''')

md("""## 6b · Matrix effect (after the run)
Needs injections of type `solvent_std` (the standard / IS mix in 10 % MeOH / water) and `matrix_spike` (the same mix spiked into
a control / blank matrix extract) in the sample sheet, and box 6. ME % = 100 × (area in spiked extract / area in solvent − 1):
below −20 % = ion suppression, above +20 % = enhancement. The map across the chromatogram can correct semi-quantification (box 9b).""")
code(r'''
#@title 6b · Matrix effect { display-mode: "form" }
if P.sheet is not None and (P.sheet.type == "solvent_std").any() and (P.sheet.type == "matrix_spike").any():
    me = P.matrix_effect()
    from IPython.display import Image
    display(Image(me["png"]))
    display(me["table"][["name", "rt", "ME_pct", "effect"]].sort_values("ME_pct"))
else:
    print("no solvent_std + matrix_spike injections in the sample sheet - matrix effect not calculated")
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

md("""## 9a · Physico-chemical properties (CompTox Chemicals Dashboard)
Looks up every identified compound, internal standard and calibrant in the US EPA **CompTox Dashboard batch search** -
automatically, no key needed. Each match is checked against the structure (name searches can return a different compound);
logKow and pKa are the experimental median where measured, otherwise the OPERA prediction. If the automatic search fails
(website down) the box prints the identifier list and the manual steps: run the batch search at
https://comptox.epa.gov/dashboard/batch-search and upload the Excel export here (tick UPLOAD_EXPORT).""")
code(r'''
#@title 9a · CompTox properties { display-mode: "form" }
AUTOMATIC = True  #@param {type:"boolean"}
UPLOAD_EXPORT = False  #@param {type:"boolean"}
exports = ui.upload_files(os.path.join(OUT_DIR, "comptox"), accept=(".xlsx",)) if UPLOAD_EXPORT else None
try:
    ct = P.comptox(auto=AUTOMATIC and not exports, exports=exports)
    t = ct["table"]
    print(t.comptox_status.value_counts().to_string())
    display(t[[c for c in ["name", "DTXSID", "preferred_name", "logKow", "logKow_source", "pKa_a", "pKa_b", "comptox_status"] if c in t]])
except RuntimeError as e:
    print(e)
''')

md("""## 9b · Semi-quantification
**IS surrogate (recommended)** - each identified compound takes the molar response of its 5 most similar internal standards
(logKow, TPSA, molar volume, MW, H-bond donors/acceptors, charge; Good, resolved, non-void IS only) measured in the IS-mix
injection (`is_mix`, else `solvent_std`, else `qc`). Matrix correction: `is_recovery` (IS response in each sample near the
compound's RT) or `me_map` (the matrix-effect map from box 6b). Uncertainty ×/÷3 near the IS, ×/÷5 / ×/÷10 further away,
×/÷10 for void-eluting compounds. **Elucidation distance** (`ed`) - the earlier method (nearest calibrant), kept for comparison.

**Units:** concentrations are first in the injected extract (µg/mL). To report per sample, give the final extract volume and the
amount of sample extracted (g or L); they are then converted to ng/g or ng/L.""")
code(r'''
#@title 9b · Semi-quantification { display-mode: "form" }
METHOD = "is"  #@param ["is", "ed"]
MIN_LEVEL_SEMIQUANT = "2b"  #@param ["1", "2a", "2b", "3"]
MATRIX_CORRECTION = "is_recovery"  #@param ["is_recovery", "me_map", "none"]
IS_MIX_FILES = ""  #@param {type:"string"}
IS_CONC_UG_ML = 0.1  #@param {type:"number"}
REPORT_PER = "extract (µg/mL)"  #@param ["extract (µg/mL)", "sample mass (ng/g)", "sample volume (ng/L)"]
EXTRACT_VOLUME_ML = 1.0  #@param {type:"number"}
SAMPLE_AMOUNT = 1.0  #@param {type:"number"}
#@markdown SAMPLE_AMOUNT = grams extracted (for ng/g) or litres extracted (for ng/L).
if REPORT_PER.startswith("extract"):
    factor, unit = 1.0, "µg/mL extract"
else:
    factor = 1000.0 * EXTRACT_VOLUME_ML / SAMPLE_AMOUNT          # µg/mL x mL / (g or L) x 1000 ng/µg
    unit = "ng/g" if "mass" in REPORT_PER else "ng/L"
sq = P.semiquant(min_level=MIN_LEVEL_SEMIQUANT, matrix_factor=factor, unit=unit, method=METHOD,
                 is_mix=[f.strip() for f in IS_MIX_FILES.split(",") if f.strip()] or None,
                 is_conc_ug_ml=IS_CONC_UG_ML or None,
                 matrix_correction={"none": False}.get(MATRIX_CORRECTION, MATRIX_CORRECTION))
t = sq["table"]
if sq.get("method") == "is":
    rf = sq["is_rf"]; print(f"{int(rf.use.sum())} of {len(rf)} internal standards used for response factors")
    if (~rf.use).any(): display(rf[~rf.use][["name", "rt", "rating", "reason"]])
cols = ["sample", "name", "level", "sq_conc", "sq_low", "sq_high", "sq_unit", "sq_nearest_IS", "sq_nn_distance", "sq_recovery",
        "sq_ME_pct_at_rt", "sq_flags", "rt", "rt_pred_from_properties", "rt_property_plausible", "elucidation_distance", "sq_calibrant"]
display(t[[c for c in cols if c in t.columns]].head(100))
''')

md("## 10 · Chromatogram viewer\nPick a compound and an injection: **Q1** = precursor ion trace (integrated region shaded), **Q3** = product-ion traces (observed in the MS2 of this peak; dotted = in-silico predicted when no MS2 is available). Shows area, S/N, peak rating, Q1 error, presence confidence and ion ratios.")
code(r'''
#@title 10 · Chromatograms { display-mode: "form" }
ui.chromatogram_explorer(P, MIN_LEVEL_REPORTED)
''')

md("## 11 · Report\n**NTS_report.html** - one file to open in any browser or e-mail: summary, injections, QC, calibration curves, matrix-effect map, identified compounds (level, evidence, confidence) and semi-quantitative concentrations with uncertainty, methods and definitions. **NTS_report.xlsx** - every table (cover sheet, identified compounds, detections per sample, presence by group, QC, IS QC, calibration, matrix effect, semi-quant, properties).")
code(r'''
#@title 11 · Report { display-mode: "form" }
rep = P.report()
print("HTML report: ", rep["html"]); print("Excel report:", rep["xlsx"])
SHOW_HERE = True  #@param {type:"boolean"}
DOWNLOAD = True  #@param {type:"boolean"}
if SHOW_HERE:
    from IPython.display import HTML
    import html as _h
    display(HTML(f'<iframe srcdoc="{_h.escape(open(rep["html"], encoding="utf-8").read())}" style="width:100%;height:700px;border:1px solid #ccc"></iframe>'))
if DOWNLOAD:
    ui.download(rep["html"], rep["xlsx"])
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
