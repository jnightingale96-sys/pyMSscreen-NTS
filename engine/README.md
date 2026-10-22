# HRMS NTS toolkit (Level 5 → 4 → 3 → 2b → 2a → 1, semi-quantification, product-ion prediction)

Python + `.bat` launchers. Input: **centroided mzML** (DDA: MS1 + MS2). Output: CSV tables, a multi-sheet Excel workbook,
MS2 mirror plots for identified compounds.

**Runs entirely offline on your own machine - no Claude, no internet and no API credits are needed at run time.**
It is a normal Python package meant to live in a git repo and be developed/validated by a team (see `CONTRIBUTING.md`, `VALIDATION.md`).

## Install and first run
```
install_nts.bat                        # creates .venv and runs `pip install -e .[dev]`   (or use your own Python >= 3.10)
copy config.local.example.yaml config.local.yaml     # point it at YOUR copies of SusDat / PubChemLite / MSP libraries (never committed)
hrms-nts doctor                        # shows config layers and what is found / missing
hrms-nts build-library                 # one-off: MassBank/GNPS/in-house MSP -> work/spectral_library.sqlite
run_tests.bat                          # pytest, 22 tests (also `pytest tests`)
python examples/make_demo_mzml.py      # tiny synthetic dataset to try the pipeline without instrument data
```
Data sources are **not** shipped (size/licence): NORMAN SusDat (norman-network.com/nds/susdat), PubChemLite for Exposomics, MassBank NIST-format MSP
(github.com/MassBank/MassBank-data), GNPS libraries (gnps-external.ucsd.edu/gnpslibrary), optional MetFrag jar and ThermoRawFileParser. Check the current download pages.

## Commands (`run_nts.bat` = interactive menu or direct mode; `hrms-nts` / `python -m hrms_nts` are equivalent)
```
hrms-nts run  <mzML folder|files> --matrix water|grain|soil|extract --semiquant single|logie|curve|none --insilico builtin|metfrag|both|none
hrms-nts validate <results folder> expected.csv [--absent absent.csv]    # recovery, ppm/RT accuracy, top-1, false-ID rate, semi-quant fold error
hrms-nts prescreen compounds.csv <mzML...>                               # targeted pyMSscreen-style QA table
hrms-nts predict-ions --smiles "NC(=O)N1c2ccccc2C=Cc2ccccc12" --adduct "[M+H]+"
hrms-nts calibrate <standard.mzML>                                       # measure calibrant areas with THIS toolkit (needed for semi-quant)
hrms-nts convert <raw folder> <mzML folder>                              # ThermoRawFileParser, centroided
hrms-nts --set level4.ppm=3 --set ms1.bg_occupancy=0.5 run ...           # any config value can be overridden per run
python tools/benchmark_insilico.py 200 7 positive                        # target/decoy benchmark of the built-in product-ion scorer
```
Configuration is layered: `config.yaml` (portable defaults, committed) <- `config.local.yaml` (this machine, git-ignored) <- `$HRMS_NTS_CONFIG` <- `--set`.
Every output folder contains `provenance.json` (toolkit version, git commit, package versions, input/data fingerprints, command) and `config_used.yaml`.

## Workflow and what each level means

| Step | Module | What it does |
|---|---|---|
| 1 | `mzml_io` | pyteomics reader, cached; MS1 peak arrays + DDA MS2 index |
| 2 | `features` | **all** centroids ≥ threshold are clustered in m/z (ppm) and split on RT gaps and chromatographic valleys (so continuous background does not become one run-length feature); apex, width, area, continuity, **m/z occupancy → `background_flag`**; isotopologue and adduct-partner flags |
| 3 | `level4` | accurate-mass match of every feature against every adduct of its polarity (pos: H, Na, K, NH4; neg: −H, Cl, HCOO, CH3COO) in NORMAN SusDat + PubChemLite + your list; isotope-pattern score (M+1/M+2 from the formula); co-eluting adducts of the same compound counted |
| 4 | `rt_prediction` | RT plausibility (±0.5 min): (a) SusDat predicted RTI mapped to your run, (b) RDKit-descriptor QSAR trained on your `standards.csv`; `rt_agree` / `rt_conflict`; authentic-standard RT difference |
| 5 | `ms2` | consensus MS2 per feature (≤3 scans closest to apex), pyMSscreen-style QA flags |
| 6 | `library` | cosine match against experimental spectra (MassBank/GNPS/MoNA/in-house); adduct-consistent; also *library-only* identifications |
| 7 | `insilico` | built-in product-ion predictor (neutral-loss rules + 1–2 bond cleavages, decoy-corrected scoring) and optional **MetFrag** (Java + jar) |
| 8 | `levels` | confidence assignment with auditable evidence columns |
| 9 | `semiquant` | single-reference, logIE/elucidation-distance corrected, or own calibration curve; matrix factors; ×/÷ uncertainty factor |
| 10 | `aggregate`, `report` | blank flag, cross-sample alignment/presence, Excel + CSV + mirror plots |

* **Level 1** – authentic standard on this method: standard RT within ±0.2 min **and** MS2 matches a library spectrum (`level1_require_inhouse_ms2: true` to demand your own standard's spectrum).
* **Level 2a** – library spectrum: cosine ≥ 0.70 and ≥ 2 matched fragments, same adduct; not contradicted by RT.
* **Level 2b** – no library spectrum, but ≥ 3 *decoy-adjusted* explained product ions (≥ 60 m/z), score ≥ 0.20, a **clear top** planar structure (margin ≥ 0.15, or ≥ 0.5 if single candidate), RT not in conflict. Never from neutral mass < 120 Da or a `background_flag` ion.
* **Level 3** – MS2 present with ≥ 2 explained ions (or weak library match) and RT plausible.
* **Level 4** – mass (≤ 5 ppm) / formula only. **Level 5** – exact mass only, no candidate.
* Stereoisomers share an InChIKey first block: the toolkit flags `stereo_ambiguous` (e.g. artemotil vs alpha-arteether cannot be separated by MS/MS).
* `confidence_score` (0–100) is a ranking aid, **not** a probability.

## Validation status (be honest in the paper)
* `pytest tests` (22 tests): feature picking, valley-splitting, ppm matching, isotope maths, cosine, level rules, config layering, semi-quant arithmetic, RT models, blank filter/alignment, MetFrag wrapper, and an end-to-end CLI run + `validate` on a generated mzML - **pass**.
* Real data, `WW_Additional_maize_Negative11.mzML` (2,823 features; ~5.5 min/sample with the built-in scorer, ~1 min without in-silico):
  4-nitrophenol recovered at RT 6.105 min, apex 1.11e9 counts (identical to the earlier Level-4 analysis) -> **Level 1** with the standard RT in `standards.csv`.
  17 Level-2a hits (mostly WW surfactants: lauryl/tetradecyl sulfate, LAS, cumene-sulfonate, sulisobenzone, valsartan acid), 3 Level-2b, 63 Level-3.
  57 % of features carry `background_flag` (ion present in >30 % of all scans); these cannot be promoted to 2b/3. Tune `ms1.bg_occupancy` / `bg_peak_to_median`.
* Built-in product-ion scorer, target/decoy benchmark on real MassBank spectra (7 mass-matched decoys each, decoys are NOT isomers):
  positive [M+H]+ (n=200): AUC 0.84, top-1 62 % (chance 12 %); 2b-style criterion passes 29 % of true structures vs 1 % of decoys.
  negative [M-H]- (n=150): AUC 0.80, top-1 57 %; criterion passes 15 % of true vs 0 % of decoys.
  => conservative and informative, but screening-grade: use `--insilico both` (MetFrag) and review mirror plots before reporting Level 2b.
* **Not validated:** absolute semi-quantification. `single` mode is only meaningful when the reference area was measured with *this* toolkit (`calibrate`);
  the area in `config.yaml` is a Compound Discoverer number (placeholder - the toolkit's trapezoid areas are on a different scale, e.g. 4-nitrophenol
  gives 5.4e7 here). Keep the x/÷10 uncertainty.
* RT models need >= 8 standards with SMILES+RT (QSAR) or >= 6 (RTI map); without them RT is "unknown" and never promotes a level.
* Library search uses MS2 precursors from DDA; DIA is not supported. Library-only identifications require <= 5 ppm to the library precursor.
* Known limitation: the built-in fragmenter predicts even-electron product ions + neutral losses, not rearrangement chemistry; MetFrag shares that limitation.

## pyMSscreen
`prescreen` re-implements the pyMSscreen idea (https://gitlab.com/uniluxembourg/lcsb/eci/pymsscreen, Artistic License 2.0):
compound list (`ID, Name, SMILES, mz, precursor_mz`) → EIC + MS2 extraction → QA columns `qa_ms1_exists, qa_ms1_good_int, qa_ms1_above_noise,
qa_ms2_exists, qa_ms2_good_int, qa_ms2_near, qa_pass` with the same defaults (MS1 1e5, MS2 2500, S/N 3, ±0.5 min, 5 ppm).
Differences: noise from EIC median/MAD (pyMSscreen uses the mean), ppm error reported, optional expected RT. No code was copied; cite the original paper if you use the concept.

## Files you maintain (git-ignored copies live in `user_data/`; templates in `examples/`)
* `user_data/standards.csv` - `name,smiles,inchikey,rt_pos,rt_neg` (authentic standards on this LC method: RT models, Level 1)
* `user_data/calibrants.csv` - `name,smiles,rt,conc_ug_ml,area,polarity` (semi-quant; >= 3 levels per compound enables `curve`)
* optional `compounds.csv` (`paths.user_compounds`) - extra suspects / targets, pyMSscreen format
* optional in-house MS2 spectra (MSP) via `paths.inhouse_msp` -> tagged `inhouse` -> Level 1
* `expected.csv` / `absent.csv` - ground truth for `validate` (see `examples/`)

## Licence
Not yet chosen - the team should agree one before the first release (`CITATION.cff` has a placeholder).
