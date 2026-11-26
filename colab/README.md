# HRMS NTS in Google Colab

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jnightingale96-sys/pyMSscreen-NTS/blob/main/colab/HRMS_NTS_Colab.ipynb)

`HRMS_NTS_Colab.ipynb` runs the whole non-target screening workflow in Google Colab, one box per step. Nothing is installed on
your computer; the results are written to your Google Drive.

| Box | What it does |
|---|---|
| 0 | Installs the engine (from this repository) and downloads the databases (Release `data-v1`) |
| 1 · 1c | Data folder on Google Drive (`.mzML` / Thermo `.raw`); 1c uploads files from your computer |
| 2 | Converts Thermo `.raw` to mzML (SCIEX `.wiff` must be converted with msconvert on Windows first) |
| 3a · 3b | Downloads **fill-in templates** (sample sheet, calibration list, IS list) and reads your sample sheet |
| 4 | QA/QC settings: replicate rule (e.g. 2/3), blank ratio, QC RSD, presence-confidence threshold, Q1/Q3 tolerances |
| 5 | Calibration: one curve per compound, each level with its own concentration (columns L1 … L8) |
| 6 · 6b | Internal standards / targets in every injection; matrix effect (spiked control extract vs solvent standard) |
| 7 | Non-target screening Level 5 → 1 (adducts, RT models, library + in-silico product ions, MetFrag optional) |
| 8 | QA: replicate presence, QC stability, blank rule, peak rating (0–10), presence confidence (%) |
| 9a · 9b | Physico-chemical properties from the CompTox Dashboard (automatic), then semi-quantification from the internal standards |
| 10 | Chromatogram viewer: any compound × injection, Q1 and Q3 traces, integration, rating |
| 11 | Report: `NTS_report.html` (one file with figures) + `NTS_report.xlsx` (all tables) |

## Quick start

1. Click the badge above (or in Colab: *File → Open notebook → GitHub*, tick *Include private repos* while the repository is private).
2. Run box 0. It asks for permission to use your Google Drive (results are saved there) and takes ~3 min.
   *Private repository:* add a Colab secret named `GH_TOKEN` (key icon in the left bar) holding a fine-grained GitHub token with
   read-only *Contents* access to the repository, and switch *Notebook access* on.
3. Put your mzML files in `My Drive/HRMS_NTS/my_batch` (or any folder; set `DATA_DIR` in box 1), or upload them in box 1c.
4. Run box 3a, download the templates, fill them in and upload them in boxes 3b, 5a and 6. The boxes check every row
   (formula vs monoisotopic mass, roles, calibration levels) before anything runs.
5. Run the remaining boxes in order. Every box saves its result to the results folder, so you can close Colab and continue
   later from any step.

**Sample types** (sample sheet): `sample`, `blank`, `qc` (pooled QC), `cal` (calibration level, with `cal_level` L1 … L8),
`is_mix` (internal-standard mix in solvent, used for semi-quantification response factors), `solvent_std` / `matrix_spike`
(the same mix in 10 % MeOH/water and spiked into a control extract: matrix effect), `exclude`.

## Databases (all free)

Box 0 downloads them from the Release (`data-v2`, falls back to `data-v1`). Details, versions and licences: `engine/DATA_SOURCES.md`.

| Purpose | Databases |
|---|---|
| Candidates (Level 4/3) | NORMAN SusDat · PubChemLite for Exposomics (2026-09, with predicted CCS) · all 137 NORMAN Suspect List Exchange collections |
| MS/MS spectra (Level 2a) | GNPS · MassBank (2026.03) · MoNA · MS-DIAL public libraries incl. lipids (970k spectra, 75k structures) |
| In-silico product ions (Level 2b) | built-in fragmenter (decoy-corrected) · MetFrag with offline MassBank/GNPS/MoNA spectra |
| Retention time | NORMAN RTI · METLIN SMRT global model (80k compounds) transferred to your method with the internal standards / standards as anchors |
| Origin labels | NORMAN SLE list membership · literature parent → TP pairs (7,898) · COCONUT · LOTUS · NP Atlas · LIPID MAPS · ChEBI roles |
| Physico-chemical properties | US EPA CompTox Dashboard (live batch search) |

Commercial libraries (NIST, Wiley, mzCloud) cannot be redistributed; a lab that owns one can add its MSP export as `inhouse_msp`.

## Semi-quantification

Box 9b estimates concentrations for identified compounds without their own standards. Each compound takes the molar response
of its five most similar internal standards (logKow, TPSA, molar volume, MW, H-bond donors/acceptors, charge), measured in the
IS-mix injection, corrected for matrix effect (IS recovery near the compound's RT, or the matrix-effect map from box 6b).
The reported range is ×/÷3 close to the IS, ×/÷5 or ×/÷10 further away. Concentrations are converted to ng/g or ng/L from the
extract volume and sample amount you enter. The earlier elucidation-distance method is kept as an option.

## MetFrag

Box 0 option `INSTALL_METFRAG` (default on) installs OpenJDK and downloads `metfrag_spectra.zip` from the Release; choose
`IN_SILICO = both` in box 7. MetFrag runs one Java call per feature (~30 min per injection): fine for a few injections, but
free Colab sessions stop after ~12 h or when idle, so run large MetFrag batches locally.

## Locally (Jupyter / VS Code) – same notebook

Set box 0 `SOURCE = Local folder (this PC)` and `LOCAL_PACKAGE` to the unzipped package; the bundled Java, MetFrag and
ThermoRawFileParser are used and there is no session limit.

## The same steps without the notebook

```python
import sys; sys.path.insert(0, r"...\engine")
from hrms_nts.workflow import Project
P = Project(r"D:\batch", r"D:\batch_results")
P.templates()                                    # fill-in templates + README
P.set_sample_sheet("sample_sheet.csv")           # file,type,group,replicate,cal_level
P.calibrate("calibration_list.csv")             # name,formula,monoisotopic_mass,smiles,label,role,is_name,L1..L8
P.targeted("is_list.csv"); P.matrix_effect()
P.screen(); P.qa("2/3"); P.comptox(); P.semiquant(method="is"); P.report()
P.chromatogram("Carbamazepine", "WW1_rep1")
```

## Maintaining the notebook

Edit `colab/make_notebook.py` and run `python colab/make_notebook.py` to regenerate `HRMS_NTS_Colab.ipynb`. Databases for the
Release are built with `python colab/make_data_zip.py <engine/data> hrms_nts_data.zip`. Box 0 updates an existing Colab
clone from GitHub, so a pushed change is picked up by re-running box 0.

Licences: the databases have their own licences (`engine/DATA_SOURCES.md`; some library records are non-commercial). Keep the
Release private unless you have checked that redistribution is allowed.
