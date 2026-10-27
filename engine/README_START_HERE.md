# HRMS non-target screening toolkit - start here

Screens centroided LC-HRMS (DDA) data in **mzML** format for suspects and unknowns: accurate-mass matching (Level 4/5), MS/MS library matching (Level 2a), in-silico product-ion
matching (Level 2b/3), retention-time plausibility checks, optional semi-quantification, Excel/CSV reports and MS2 mirror plots.
Everything runs **offline on your own machine** (internet is needed once, to install the Python packages).

## What is in the zip
| Folder / file | What it is |
|---|---|
| `hrms_nts/` | the program (Python package) |
| `run_nts.bat`, `install_nts.bat` | Windows launchers (menu), one-time installer |
| `data/susdat.csv` | NORMAN SusDat suspect list (~100k chemicals, predicted RT index) |
| `data/PubChemLite_exposomics.csv` | PubChemLite for Exposomics (~400k chemicals) |
| `data/spectral_library.sqlite` | ready-built MS/MS library, 970k experimental spectra (GNPS, MassBank, MoNA, MS-DIAL libraries) |
| `data/Updated_NORMAN_with_TPSA.xlsx` | NORMAN table with retention-time index (RTI), applicability domain, TPSA - used for the RT prediction |
| `data/MassBankEU_compounds.csv` | measured MassBank-EU fragments, used as extra support for in-silico identifications |
| `data/MetFragCommandLine-2.6.11.jar`, `data/jre/`, `data/mb_converted/` | MetFrag product-ion prediction, its Java runtime and offline spectra |
| `data/ThermoRawFileParser/` | converts Thermo `.raw` to mzML |
| `config.yaml`, `config.local.yaml` | settings (the local file already points at `data/`) |
| `user_data/` | your own standards / calibrants (blank templates) |
| `DATA_SOURCES.md` | where each database comes from and its licence |

## 1. Install (once)
1. Unzip to a **short path without spaces**, e.g. `C:\NTS\HRMS_NTS_Toolkit` (about 5 GB unpacked; keep 10 GB free for outputs and caches).
2. Install Python 3.10 or newer from python.org (tick "Add python.exe to PATH") if you do not have it.
3. Double-click **`install_nts.bat`** (needs internet; installs numpy, pandas, RDKit, pyteomics ... into a local `.venv`). It ends with a `doctor` check.
4. Double-click **`run_nts.bat`**, choose **1 (doctor)**: every data path should say found / OK. (Java: found because the bundled runtime is used.)

## 2. Prepare your data
* Input must be **centroided mzML with MS1 and DDA MS2**. For Thermo `.raw`: menu **7** in `run_nts.bat` (or `hrms-nts convert <raw folder> <mzml folder>`).
  If you convert with ProteoWizard MSConvert, switch on *peak picking (vendor, MS levels 1-)*.
* Polarity is read from the mzML; if that fails put `POS` or `NEG` in the file name, or pass `--polarity positive|negative`.
* Blanks: include your solvent / procedural blank mzML files in the same run and list their file names (without `.mzML`) under `blank: samples: [...]` in `config.local.yaml`
  (or `--set blank.samples=[blank1,blank2]`) so features present in the blanks are flagged (`blank_flag`).

## 3. Screen
Menu **3**, or on the command line (from the toolkit folder):
```
run_nts.bat run "D:\data\mzml\*.mzML" --matrix water --semiquant none --insilico builtin
run_nts.bat run "D:\data\mzml\*.mzML" --matrix water --semiquant none --insilico both     # adds MetFrag (slower, usually more 2b candidates)
```
`--matrix` is `water`, `grain`, `soil` or `extract` (only used for semi-quant factors). Expect roughly 1-6 minutes per sample depending on feature count and in-silico engine.
Outputs go to `work/<run name>/`: a multi-sheet Excel workbook, CSV tables, mirror plots for identified compounds, `provenance.json` (versions, settings, data fingerprints) and `config_used.yaml`.

Other useful menu items: **5** predict product ions for any SMILES, **6** measure calibrant areas (semi-quant), **9** validate against a list of expected compounds.

## 4. How to read the confidence levels
* **Level 1** authentic standard on your method (standard RT + MS2 match; needs `user_data/standards.csv`).
* **Level 2a** experimental library spectrum matches (cosine >= 0.70, >= 2 fragments), RT not contradicting.
* **Level 2b** no library spectrum but >= 3 explained product ions from the in-silico predictors and a clear top candidate, RT not contradicting.
* **Level 3** tentative: MS2 with >= 2 explained ions, RT plausible.
* **Level 4** accurate mass / formula only (the structure is just the best-ranked database isomer). **Level 5** exact mass, no candidate.
Evidence columns in the tables show exactly why each level was given. The 0-100 `confidence_score` is a ranking aid, not a probability.

## 5. Retention-time prediction (what it does)
Each candidate can get a predicted RT from the NORMAN retention-time index (RTI, `susdat.csv` / `Updated_NORMAN_with_TPSA.xlsx`) mapped onto *your* gradient, plus a QSAR model on the structure.
**Both are calibrated with authentic standards you ran on your own method**: fill at least **6** compounds with a retention time per polarity in `user_data/standards.csv`
(columns: name, smiles, inchikey, rt_pos, rt_neg - minutes). Without them the RT check is skipped (the log says "RT model ... skipped") and Level 2b/3 rely on MS/MS evidence alone.
A good set is a mix spanning the gradient (polar to non-polar). Candidates whose RT differs from the prediction by more than the tolerance (`rt:` section of `config.yaml`) are flagged `rt_conflict` and cannot reach Level 2.
On our own wastewater/maize data the RTI-to-RT error was about 1 min, so treat RT as a plausibility filter, not proof.

## 6. Semi-quantification
`--semiquant none` is the safe default. `single` / `logie` / `curve` need calibrant areas measured on your own instrument (`user_data/calibrants.csv`, menu 6); the shipped file is only a template.
Our semi-quantitative values carried an uncertainty of about an order of magnitude or more, so use them for ranking, not for absolute concentrations.

## 7. Please be aware
* The tool is a research prototype: tests pass and it was used on real wastewater, maize and soil data, but it is **not validated** for regulatory or forensic use. Check important identifications by hand (mirror plots) and with a standard.
* Not included: the project-specific fate-model workbook and the PEC-vs-MEC part (`pec-screen`) needs it - ignore that menu entry; no sample data of ours.
* Level 4/3 lists are long and contain many false positives; use the blank flag, RT check and `endogenous` filters, and sort by level and confidence.
* If you use results in a publication, cite the databases listed in `DATA_SOURCES.md` and the toolkit (`CITATION.cff`).

## Troubleshooting
* `doctor` shows NOT FOUND: you probably moved or renamed the `data` folder - edit `config.local.yaml` or put it back.
* "No suitable Python found": run `install_nts.bat` again (needs internet) or install Python and tick "Add to PATH".
* Out of memory on big files: raise `ms1.store_min_intensity` (e.g. `--set ms1.store_min_intensity=30000`) and process fewer files at a time.
* Questions about the method: see `README.md` (full documentation), `VALIDATION.md` and `CHANGELOG.md`.
