# CLAUDE.md — HRMS non-target screening package (pyMSscreen-NTS)

You are working inside a self-contained LC-HRMS analysis package. Use it to screen mzML data for suspects and unknowns
(Level 5 → 1), semi-quantify them, predict product ions and report results. Everything runs offline once installed.
Read this file fully before running anything; `engine/README.md` and `engine/VALIDATION.md` have the method details.

## Layout
| Path | What |
|---|---|
| `engine/` | the analysis engine, Python package `hrms_nts` (CLI: `python -m hrms_nts ...`). **Run engine commands from this folder** (relative paths in its config resolve here). |
| `engine/data/` | bundled databases: NORMAN SusDat, PubChemLite, `spectral_library.sqlite` (970,218 MS/MS spectra, 75,561 structures, 28 libraries), NORMAN RTI workbook, MassBank-EU fragments, MetFrag jar + offline spectra, Java runtime (`jre/`), ThermoRawFileParser. Sources/licences: `engine/DATA_SOURCES.md`. |
| `engine/config.yaml` | all thresholds (defaults, never hard-code numbers in code). `engine/config.local.yaml` = this machine's paths (already relative to `engine/data`). Override per run with `--set section.key=value`. |
| `engine/user_data/` | `standards.csv` (authentic standards + RTs on the user's LC method), `calibrants.csv` (semi-quant calibrants) — **blank templates; must be filled by the user**. |
| `backend/`, `frontend/` | pyMSscreen web app (Flask + React) with extra NTS / Semi-quant / Product-ion pages. For people; you normally use the CLI. |
| `START_pyMSscreen-NTS.bat` | installs `.venv` (first run) and starts the web app on http://127.0.0.1:8765. `START_pyMSscreen-NTS.bat install` = install only. |

Python to use: `.venv\Scripts\python.exe` in the package root (created by the installer). From `engine/` that is `..\.venv\Scripts\python`.

## First-time setup (do this once, check each step)
1. Python ≥ 3.10 must exist (`py -3 --version`). Then run `START_pyMSscreen-NTS.bat install` from the package root (needs internet for pip, 5–10 min).
   If you cannot run .bat files: `py -3 -m venv .venv` then `.venv\Scripts\python -m pip install -r requirements-windows.txt`.
2. `cd engine` and run `..\.venv\Scripts\python -m hrms_nts doctor`. Every data path must be `[ok]` and the library must report 970,218 spectra.
   `user_compounds`, `inhouse_msp`, `library_sources_dir` may be `[--]` (optional).
3. `..\.venv\Scripts\python -m pytest tests -q` — all tests must pass before you analyse real data. Report failures; do not "fix" thresholds to make tests pass.

## Standard analysis workflow
1. **Inventory the data.** List the mzML files (size, polarity from name/metadata). Ask the user which files are samples, blanks, standards, QCs and which matrix each is.
2. **Check the files are usable** before screening:
   - MS1 + MS2, **DDA or DIA/SWATH**, centroided **or profile** (since 0.4.0: profile spectra are centroided on reading, DIA windows are
     detected and a correlation pseudo-MS2 is built per feature; the instrument family is read from the mzML header and SCIEX files get TOF
     thresholds - `config.yaml` `instrument.profiles`). The log line `read <file>: <instrument>, <DDA|DIA>, ...` tells you what was detected.
   - A warning `has NO MS1 scans between X and Y min` means an incomplete export - ask for a re-export from the raw file.
   - Check the acquisition order (`startTimeStamp` in the mzML header / `Files` sheet of `targeted`) before interpreting blanks (carry-over).
   - Must be **complete**: a valid mzML ends with `</indexedmzML>` or `</mzML>`. The engine reads truncated files only up to the break and logs `is TRUNCATED` — tell the user and ask for the `.raw` file instead of reporting partial results.
   - Thermo `.raw` → mzML: `python -m hrms_nts convert <raw folder> <mzml folder>` (bundled ThermoRawFileParser, centroided).
3. **RT calibration (strongly recommended):** fill `engine/user_data/standards.csv` (`name,smiles,inchikey,rt_pos,rt_neg`, minutes) with ≥ 8 standards per polarity
   spanning the gradient, measured on the user's method. Without them the RT models are skipped (log: "RT model ... skipped") and the same compound
   can appear at several RTs — say so in the report.
4. **Blanks:** pass blank file stems with `--set "blank.samples=[blank1,blank2]"` so features in blanks get `blank_flag`.
5. **Screen:**
   ```
   ..\.venv\Scripts\python -m hrms_nts run "D:\data\mzml" --out "D:\data\NTS_out_<name>" --matrix water --semiquant none --insilico builtin
   ```
   - Default to `--insilico builtin` (~2 min per sample). `--insilico both` adds MetFrag (~30 min per sample: one Java call per feature) —
     **see the MetFrag caveat under Reporting rules before using its Level 2b hits.** `--polarity positive|negative` if it cannot be detected.
   - For many files run it in the background and watch `<out>/run.log`.
6. **Read the results** (`<out>/`):
   - `level_summary.csv` features per level per sample · `compound_summary.csv` best level per compound across samples (the paper table)
   - `identified_level1_2b.csv` Level 1/2a/2b with evidence · `top_candidates_per_feature.csv` every feature · `all_candidates.csv` every candidate
   - `NTS_results.xlsx` all of the above · `plots/` MS2 mirror plots for identifications · `run.log`, `provenance.json`, `config_used.yaml` (always keep these with results)
   - Evidence columns: `ppm_error`, `rt_agree`/`rt_conflict`, `std_rt_diff`, `lib_score`/`lib_matched`/`lib_name`/`lib_source`, `bi_score`/`bi_expl_adj_n` (built-in), `mf_score` (MetFrag),
     `iso_score`, `background_flag`, `blank_flag`, `stereo_ambiguous`, `level_basis`, `confidence_score` (ranking aid, **not** a probability).
7. **Semi-quantification** only if calibrant areas measured on the same instrument exist: `python -m hrms_nts calibrate <standard.mzML> --out calibrant_areas.csv`,
   copy areas into `user_data/calibrants.csv`, then `run ... --semiquant logie` (elucidation-distance / ionisation-efficiency) or `single` / `curve`.
   Matrices and unit factors are in `config.yaml` `semiquant.matrices`. Always report the ×/÷ uncertainty factor (default 10).
8. **Product ions / transitions** for any structure: `python -m hrms_nts predict-ions --smiles "<SMILES>" --adduct "[M+H]+"`.
9. **Targeted screen of internal standards / targets / calibrants** (Q1 precursor and Q4 product-ion errors, levels, integration QC, blanks,
   calibration curves, RT anchors) - use this first whenever the user has an IS list or wants to check mass accuracy / calibration:
   ```
   ..\.venv\Scripts\python -m hrms_nts targeted <mzML...> --list <IS list .csv/.xlsx> --blank <blank stem or pattern> [--sample-sheet file,type,conc] [--q1-ppm 5 --q4-ppm 10]
   ```
   - The list may be a vendor sheet (Name, Chemical formula with `[2H]`/`[13C]`, Q1) or name + SMILES + label; `examples/drug_IS_pos_75.csv` is a validated example.
     Do not trust list masses: the `Mass_check` sheet recomputes every Q1 from the formula (explain neutral vs [M+H]+, average MW, salt forms if they disagree).
   - Output `targeted_results.xlsx` (Detected with Q1/Q4 pass/fail, peak area/S/N/rating, blank %, level, flags; Blank_check; Not_detected with reasons;
     Calibration/Quantification when `cal` files are given; Native_at_IS_RT; RT_anchors; Mass_check) + `EIC_integration_QC.pdf` (look at it).
   - Report the IS RTs as RT anchors and offer `--update-standards` (adds them to `user_data/standards.csv` -> RT models + Level 1 in `run`).
   - Negative-mode compounds are flagged when only ESI+ files exist; isobaric IS (same formula) are split by library MS2 / batch RT or flagged unresolved.
10. **Targeted prescreen** of a compound list (pyMSscreen-style QA): `python -m hrms_nts prescreen compounds.csv <mzML...>` (`ID,Name,SMILES,mz`).
11. **Validation** against known spiked/absent compounds: `python -m hrms_nts validate <out> expected.csv --absent absent.csv` (see `engine/VALIDATION.md`).

## Confidence levels (Schymanski et al. 2014 scheme, as implemented)
- **1** authentic standard on this method: standard RT within ±0.2 min and MS2 library match.
- **2a** library spectrum: cosine ≥ 0.70, ≥ 2 matched fragments, same adduct, RT not contradicting.
- **2b** no library spectrum: ≥ 3 decoy-adjusted explained product ions (≥ 60 m/z), score ≥ 0.20, clear top candidate, RT not in conflict; never from neutral mass < 120 Da or a `background_flag` ion.
- **3** MS2 with ≥ 2 explained ions or weak library match, RT plausible. **4** exact mass/formula only (structure = best-ranked isomer). **5** exact mass, no candidate.

## Reporting rules (important)
- Report identifications with their level and the evidence; never present Level 3/4 as identified compounds. Level 4/3 lists contain many false positives.
- Same name at several RTs, in-source fragments, adduct partners and stereoisomers (`stereo_ambiguous`) are common — group and flag them, do not double-count.
- **MetFrag caveat (`--insilico both`/`metfrag`):** the Level 2b rule accepts MetFrag's raw explained-peak count (`mf_expl`, not decoy-corrected) and
  its relative score (`mf_score`, the top candidate of each feature is always 1.0). In a soil test this raised Level 2b from 3 (built-in only) to 62,
  many of them implausible ([M+Cl]⁻/[M+HCOO]⁻/[M+CH₃COO]⁻ adducts, RT < 1 min, `bi_expl_adj_n` = 0). Report a 2b as such only if `bi_expl_adj_n` ≥ 3
  (built-in, decoy-corrected) or it is otherwise supported (MassBank-EU fragments, RT agreement); call MetFrag-only 2b "tentative (Level 3)" and say why.
- Ubiquitous contaminants (phthalates, PEG/ethoxylates, fatty-acid amides, siloxanes, surfactants) are often in blanks — check `blank_flag` before claiming them.
- Semi-quantitative concentrations are order-of-magnitude estimates (validation pending): give the ×/÷ factor and never call them measured concentrations.
- Say which files were truncated, profile-mode or missing MS2, and whether RT models were active.
- Keep `provenance.json` + `config_used.yaml` with every result you hand over. Cite the databases (`engine/DATA_SOURCES.md`) and pyMSscreen.

## Figure style (the user's preference)
- PEC vs MEC / concentration comparisons: box-and-whisker per chemical, chemicals on the Y axis with dashed separators, log concentration on the X axis,
  turquoise boxes = PEC, white boxes = MEC, median line, min–max whiskers, every field/sample as a dot; no gridlines; plain-number ticks (0.001 not 10⁻³);
  axis fitted tightly to the data. Use `hrms_nts/boxchart.py` (`plot_box_chart`). State in the caption what box/line/whiskers/dots mean.

## Spectral library maintenance
- `engine/data/spectral_library.sqlite` is ready to use. To rebuild from downloaded libraries: put them in one folder (`massbank/ msdial/ mona/ gnps/ lipidblast/`),
  set `paths.library_sources_dir` in `config.local.yaml`, run `python -m hrms_nts build-library` (hours; keeps the previous library as `*_prev.sqlite`, writes a manifest
  with a self-check). `python tools/package_library.py` zips it.

## Not included / not for general use
- `pec-screen` / `cd-screen` (PEC-vs-MEC for one irrigation project) need a fate-model workbook that is not shipped — do not use them unless the user provides it.
- `engine/tools/` contains project analysis scripts (RT anchor fits, MetFrag comparisons, figures); read before running, several expect project-specific inputs.

## Working rules for you
- Never write results into the user's raw-data folders unless asked; use `--out` to a results folder.
- Do not change thresholds in `config.yaml` silently: use `--set` for a run, and record any change (it is saved in `config_used.yaml`). Threshold/level-rule changes need validation data (`engine/CONTRIBUTING.md`).
- Run long jobs in the background and report progress from `run.log`; do not poll in tight loops.
- If a database is missing (`doctor` shows `[--]`), tell the user — do not download replacements without asking.
