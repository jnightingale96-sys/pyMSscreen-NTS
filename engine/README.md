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
hrms-nts build-library                 # one-off: all MS/MS libraries -> work/spectral_library.sqlite (or use a ready-made one, below)
run_tests.bat                          # pytest (also `pytest tests`)
python examples/make_demo_mzml.py      # tiny synthetic dataset to try the pipeline without instrument data
```

### MS/MS spectral library
* **Ready-made:** the 2026-10-23 build (970,218 spectra, 75,561 structures, 28 libraries: MassBank, MS-DIAL public libraries, MoNA LC-MS/MS, GNPS)
  ships as `spectral_library.sqlite` + manifest; point `paths.library_db` at it.
* **Rebuild:** put the downloaded libraries in one folder (`massbank/ msdial/ mona/ gnps/ lipidblast/`), set `paths.library_sources_dir`, run
  `hrms-nts build-library`. MSP and MoNA/MetFrag `.mb` files are read; GC-EI, MALDI and DART spectra are dropped; cross-source duplicates are removed
  (curated libraries win over GNPS); missing InChIKeys are derived from SMILES; the previous library is kept as `*_prev.sqlite`; a manifest with
  counts, MD5s and a self-check search is written. `python tools/package_library.py` zips it.
Data sources are **not** shipped (size/licence): NORMAN SusDat (norman-network.com/nds/susdat), PubChemLite for Exposomics, MassBank NIST-format MSP
(github.com/MassBank/MassBank-data), GNPS libraries (gnps-external.ucsd.edu/gnpslibrary), optional MetFrag jar and ThermoRawFileParser. Check the current download pages.

## Commands (`run_nts.bat` = interactive menu or direct mode; `hrms-nts` / `python -m hrms_nts` are equivalent)
```
hrms-nts run  <mzML folder|files> --matrix water|grain|soil|extract --semiquant single|logie|curve|none --insilico builtin|metfrag|both|none
hrms-nts validate <results folder> expected.csv [--absent absent.csv]    # recovery, ppm/RT accuracy, top-1, false-ID rate, semi-quant fold error
hrms-nts prescreen compounds.csv <mzML...>                               # targeted pyMSscreen-style QA table
hrms-nts predict-ions --smiles "NC(=O)N1c2ccccc2C=Cc2ccccc12" --adduct "[M+H]+"
hrms-nts targeted <mzML...> --list IS.csv --blank eq5 [--sample-sheet sheet.csv]   # internal standards / targets / calibrants: Q1+Q4 errors, levels, blanks, curves
hrms-nts calibrate <standard.mzML>                                       # measure calibrant areas with THIS toolkit (needed for semi-quant)
hrms-nts convert <raw folder> <mzML folder>                              # ThermoRawFileParser, centroided
hrms-nts --set level4.ppm=3 --set ms1.bg_occupancy=0.5 run ...           # any config value can be overridden per run
python tools/benchmark_insilico.py 200 7 positive                        # target/decoy benchmark of the built-in product-ion scorer
```
Configuration is layered: `config.yaml` (portable defaults, committed) <- `config.local.yaml` (this machine, git-ignored) <- `$HRMS_NTS_CONFIG` <- `--set`.
Every output folder contains `provenance.json` (toolkit version, git commit, package versions, input/data fingerprints, command) and `config_used.yaml`.

## SCIEX, profile and DIA (SWATH) data (0.4.0)
* **Instrument profiles** (`config.yaml` `instrument.profiles`, `instrument.type: auto`): the mzML header is read and SCIEX files get TOF-appropriate
  thresholds (MS1 store 150 / feature 3 000 counts, product-ion tolerance 15 ppm, targeted peak 600 counts) instead of the Orbitrap defaults
  (10 000 / 100 000). `--set` always wins over a profile. Add your own profile (e.g. `bruker_qtof`) the same way.
* **Profile spectra** are centroided on reading (`hrms_nts/centroid.py`: 5-point smoothing for maxima, intensity-weighted centroid of the points
  >= 50 % of the apex). On X500R data this removes the 3-4 split centroids per peak that raw local maxima give at low counts. Centroided files are untouched.
* **DIA** (isolation windows wider than `dia.min_window_da` = 5 Da, e.g. SWATH 25-200 Da) is detected automatically and stored per window. For every
  feature a **pseudo-MS2** is built from the fragments of its window whose XIC correlates with the precursor MS1 XIC (Pearson r >= 0.7 on smoothed
  traces, apex within 3 cycles) - `hrms_nts/dia.py`. It fills the same `ms2_*` columns as DDA, so library search, in-silico scoring and the Level
  rules are unchanged. Co-eluting precursors in the same window with identical elution profiles cannot be separated (true for any DIA software).
* **Scan gaps** (> 1 min without MS1, e.g. an interrupted export) are reported as warnings and in the targeted `Files` sheet.

## Step-by-step project workflow and Google Colab (0.5.0)
`colab/HRMS_NTS_Colab.ipynb` (one box per step) or from Python:
```python
from hrms_nts.workflow import Project
P = Project(r"D:\batch", r"D:\batch_results")
P.convert_raw()                         # Thermo .raw -> mzML (ThermoRawFileParser); SCIEX .wiff: msconvert on Windows
P.set_sample_sheet("sheet.csv")         # file,type(sample|blank|qc|cal),group,replicate,conc   (P.guess_sample_sheet() to start)
P.calibrate("cal_list.csv")             # name,smiles|formula,label,role(target|IS),is_name -> linear / quadratic curves
P.targeted("IS_list.csv")               # IS / targets: Q1/Q3, levels, integration, blanks, IS QC
P.screen("water")                       # NTS Level 5 -> 1
P.qa("2/3")                             # replicate rule, QC, blanks, peak rating, presence confidence
P.semiquant(method="is", is_mix="ISMIX", is_conc_ug_ml=0.1)   # 0.6.0: kNN of labelled-IS response factors + matrix correction
P.semiquant(method="ed", api_key=...)   # 0.5.0: elucidation distance to calibration standards + CompTox logP
P.report(); P.chromatogram("Carbamazepine", "S1_rep1")
```
Each step saves its results in the results folder. Details and test status: `colab/README.md`, VALIDATION.md items 8-11.

## Targeted screen: internal standards, targets, calibrants (`targeted`)
```
hrms-nts targeted <mzML...> --list examples/drug_IS_pos_75.csv --sample-sheet examples/sample_sheet_example.csv
hrms-nts targeted <folder> --list DrugIS_pos.xlsx --blank eq5 "*blank*" --q1-ppm 5 --q4-ppm 10 --update-standards
```
* **Compound list** (CSV/XLSX): `name` + `formula` (labels as `[2H]`, `[13C]`, `[15N]` or `D`) **or** `smiles` (+ `label` e.g. `D4`, `13C1,D3`, or
  `nD`/`n13C`; the label is also read from the name, e.g. `Fentanyl-D5`). Optional: `role` (IS/target/calibrant), `Q1` (checked against the
  formula -> `Mass_check`), `rt` (expected RT -> Level 1), `is_name` (IS used for area ratios), `polarity`. Vendor IS sheets (Name, Chemical formula,
  Adduct, Precursor (Q1) Mass) load as they are; missing SMILES are looked up by parent name in SusDat/PubChemLite (needed for library matching).
* **Sample sheet** (optional): `file,type,conc` with type `sample | blank | qc | cal`; `cal` files with a concentration give weighted (1/x)
  calibration curves (R2, LOD/LOQ, back-calculated accuracy) and quantification of the other files; `calibrants_for_semiquant.csv`
  can be pasted into `user_data/calibrants.csv` for the NTS `curve` / `logie` semi-quantification.
* Per compound and file: **Q1** precursor error (search +-40 ppm so a calibration offset cannot hide a peak; acceptance `targeted.q1_ppm` = 5),
  candidate peaks, MS2 (DIA pseudo-MS2 or DDA spectrum), **Q4** product ions annotated as sub-formulae of the labelled formula (acceptance
  `targeted.q4_ppm` = 10), **label-aware library cosine** against the installed library (a library fragment of the unlabelled parent may match
  with k deuteriums / 13C retained) -> Level 1 / 2a / 2b / 3 / 4. Isobaric compounds (same formula) are split by their library match, or by the
  RT at which they were library-identified in another file of the batch; otherwise flagged `isobaric, unresolved`.
* **Integration**: +-10 ppm XIC, explicit start/end, baseline, S/N against the centroid noise floor, rating Good/Fair/Poor; every integration is
  drawn in `EIC_integration_QC.pdf` (all files overlaid, blanks thick).
* **Blank check**: each blank XIC is integrated over the compound's sample RT window without peak picking; verdict clean / minimal (< 1 %) /
  low (1-5 %) / significant (> 5 %).
* **Flags** (Detected sheet): acidic compound with weak ESI+ evidence only, RT inconsistent with the other files, peak in a blank file,
  same signal in the blank.
* **RT**: reported, not an identification criterion unless an expected RT is given. RT models (QSAR / NORMAN RTI on `standards.csv`) add a
  plausibility prediction; `rt_anchors_for_standards.csv` holds the RTs of the detected labelled standards for their native compounds
  (`--update-standards` appends them to `user_data/standards.csv`, which then trains the RT models and enables Level 1 in `run`).
* **Natives**: the unlabelled analogue of every labelled IS is integrated at the IS RT (exact mass + co-elution + library MS2) -> `Native_at_IS_RT`.
* Outputs: `targeted_results.xlsx` (Notes, Detected, Blank_check, Not_detected, Calibration, Quantification, Native_at_IS_RT, RT_anchors,
  Mass_check, Files, All_candidate_peaks, Product_ions), `EIC_integration_QC.pdf`, CSVs, `provenance.json`.
* Web app: page **IS / targets**. Menu: `run_nts.bat` option **T**.

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
* DDA: library search uses the precursor-matched MS2. DIA/SWATH (since 0.4.0): correlation pseudo-MS2 per feature (see "SCIEX, profile and DIA data"). Library-only identifications require <= 5 ppm to the library precursor.
* Known limitation: the built-in fragmenter predicts even-electron product ions + neutral losses, not rearrangement chemistry; MetFrag shares that limitation.

## PEC-vs-MEC screening of all identified chemicals (`pec-screen`)
Compares what the fate model predicts for **maize grain** from the irrigation water (PEC) with what is measured in the grain (MEC), for every chemical the NTS run identified at Level 2b or better - **both sides semi-quantified with the same single-reference method** (no literature concentrations):

    conc (ug/mL extract) = ref_conc x area / area_ref        area_ref: carbamazepine standard, integrated by THIS integrator
    water ng/L = conc x 1000 x 120        grain ng/g dw = conc x 0.5/3 x 1000        PEC (ng/g dw) = k_fruit x C_water(ng/L) x 0.005

    python -m hrms_nts pec-screen my_map.csv --out PEC_out --properties user_data/properties_project.csv --properties openera_export.xlsx::Sheet1

1. **Sample map** (`examples/pec_map_example.csv`): one row per injection - matrix (water|grain), field label, role (sample|blank|exclude). Field labels are stored inside Thermo `.raw` headers (`tools/raw_sample_names.py`).
2. **Candidates** = Level >= `--min-level` identifications in the water samples (identification is cached in `<out>/_ident`; rerun is fast).
3. **Water**: EIC re-integration in every injection, blank check (apex must exceed 3x the highest blank), mean over detected injections per field. **Grain**: same compound found by m/z + RT (auto RT-shift per injection from MS2-confirmed anchors) and **MS2 cosine vs the water spectrum** (confirmed / weak / MS1-only / conflict = different compound, not counted).
4. **Properties** (`--properties`, first match wins; InChIKey > CAS > name): toolkit CSV and/or a CompTox batch export (OPERA). Missing logKow/Kaw/DT50 -> structure-based estimates (flagged) and a **property grid**, so the PEC is a *range*, never an invented point. Zwitterions and permanent cations are reported as "fate model not applicable". `candidates_for_comptox_batch.csv` lists what to look up.
5. **Fate model**: the project's Excel Trapp/SimpleTreat workbook via COM (`paths.fate_model`; Windows + Excel + pywin32). Every k_fruit is cached (`work/kfruit_cache.csv`), so colleagues without Excel can rerun; or put `k_fruit` in the properties table.
6. **Classes** (per chemical, geometric mean PEC/MEC over paired fields, interval overlap with the x10 window): *water-explained* | *MEC >> PEC (other source / contamination?)* | *PEC >> MEC* | *predicted but not found in grain* | *PEC below grain detection floor* | *grain only* | *water only* | *fate model not applicable*.

Limits to keep in mind: semi-quant has a ~x10 uncertainty per compound (ionisation efficiency); only fields with both water and grain data are paired; grain MS2 is DDA, so low-level grain hits are often MS1-only (shown with open markers); OPERA biodegradation half-lives are used as soil DT50 when no BIOWIN4 / DT50 is in the table (flagged in `properties_used.csv`).

### Same comparison from Compound Discoverer tables (`cd-screen`)
    python -m hrms_nts cd-screen Israel_Modelling_Data.xlsx --out CD_out --properties ... --combine-with PEC_out

Reads the workbook sheets 'Semi quant_Modelling_WW' (one max value per chemical, ng/L) and 'Semi quant modelling_Plants' (ng/g per field 1-9) and runs the same properties -> fate -> PEC-vs-MEC code. Because the water sheet has no per-field values, the one water value is used for every grain field. Structures are taken only from rows where the CompTox block matches the name (it is misaligned on ~45 % of the Plants rows), then from the suspect database by exact name; several names with one InChIKey are merged. CD pairs water and grain by name only, so a retention-time difference > 0.5 min between the two sheets is reported as **identity conflict** (e.g. 17b-estradiol 8.29 vs 6.18 min). `--combine-with` writes one table and chart (`pec_vs_mec_combined_chart.png`) with the mzML-based results (teal) and the CD results (orange).

### Endogenous metabolites, RT plausibility and MetFrag second opinion
* **Endogenous filter** (`pec_screen.exclude_endogenous`, on by default; `--keep-endogenous` to disable): chemicals a plant makes itself (amino acids/dipeptides, TCA and other small acids, sugars, nucleosides, cytokinins, oxylipins and plant fatty acids, phenolic acids, flavonoids, glycosides/glycolipids, chlorophyll catabolites) are taken out of the ranking and listed in `endogenous_flagged.csv` with the reason and confidence (*high* / *probable*). Curated names + SMARTS rules + the CD 'Class' text; add your own in `user_data/endogenous_extra.csv`.
* **RT plausibility without standards** (`tools/rt_anchor_fit.py`, `tools/rt_check_levels.py`): robust Level 1/2a identifications of the run are used as RT anchors for the NORMAN `Pred_RTI` -> RT map. On the Israel water run the anchors give a leave-one-out RMSE of ~1 min (73 % within +-1 min), so the +-0.5 min criterion cannot be applied; a +-2 x RMSE window still exposes gross conflicts: **44 % of the Level 2b candidates fall outside it** (7 % of the anchors), i.e. a large share of 2b assignments is doubtful. Real standards on the method would tighten this - add them to `user_data/standards.csv`.
* **Identification filter** (`pec_screen.id_filter`, on by default; `--no-id-filter` to disable): Level 1/2a are kept; a Level 2b is kept only if it is retention-time plausible (RTI->RT map fitted on the run's own robust 2a anchors, tolerance 2 x leave-one-out RMSE, floor 1 min) **or** its MS2 contains the fragments MassBank-EU recorded for that structure. 2b chemicals with no NORMAN Pred_RTI and no MassBank-EU support are dropped. Written: `candidates_removed_by_id_filter.csv`, `rt_map.json`, and `rt_err_min / rt_status / mb_eu_*` columns in the summary.
* **MassBank-EU table** (`pec_screen.id_filter.massbank_eu`, `hrms_nts/massbank_eu.py`): NORMAN MassBankEU_Cmpds_*_wMS_DTXSIDs_*.csv. Synonyms in NAME are joined by `|` and some rows start with a record-code token (`LID_235.1805_10.1`) - both are handled; rows are merged per InChIKey. On the Israel run 20 of the 23 Level 2a identifications listed there have their recorded fragments in the feature MS2 (87 %).
* **Applicability domain and the papers' RT models** (`id_filter.use_applicability_domain`, `tools/rt_models_compare.py`): following the RTI paper (Aalizadeh 2021) RT only counts inside the RTI model's domain (`Uncertainty_RTI` = 'Covered by Model', from SusDat / `Updated_NORMAN_with_TPSA.xlsx`), and the anchors inside it set the tolerance (leave-one-out RMSE 0.95 min, tolerance +-1.9 min). We also tested the improvements suggested by the MARS paper (Celma 2022: model sqrt(RT)) and a TPSA / logP descriptor correction: with ~30 anchors none beats the plain RTI->RT map (LOO RMSE 1.04-1.22 vs 1.05 min), so the simple map is kept. MultiConditionRT (Souihi 2022) needs training data for other columns/eluents and PaDEL descriptors; its gradient-mapping step is what the anchor map already does.
* **MetFrag product ions and candidate databases** (`tools/metfrag_compare.py`): scores each candidate's MS2 against the NORMAN workbook, NORMAN SusDat and PubChemLite structures at its mass and writes `work/metfrag_compare.csv` plus `work/metfrag_product_ions.csv` (every observed MS2 peak with the fragment formula MetFrag assigns to the assigned structure, next to the built-in predictor's label). MetFrag explains ~47 % of the MS2 peaks of the assigned structures vs ~18 % for the built-in predictor (`insilico.run_metfrag` now returns `mf_expl_formulas`).
* **MetFrag second opinion** (`tools/metfrag_verify.py`): scores every SusDat + PubChemLite structure at the candidate's mass (FragmenterScore + offline MassBank) and a drug-rich NORMAN collection subset (source code S1). Fragmenter scores alone do not separate isomers (confirmed carbamazepine/lidocaine/venlafaxine rank 2/2/6 among isomers), so treat it as one more flag, not as proof. `tools/id_flags.py` merges the RT and MetFrag flags into `pec_screen_summary_flagged.csv`.

### Soil instead of grain (`pec-screen --sink soil`)
Same screen with soil as the solid matrix: MEC = soil semi-quant (project soil formula, ng/g dw = ug/mL x 4500), PEC = k_soil (fate model steady-state soil concentration per mg/L irrigation water, Excel cell Input!G39) x C_water x 0.001. A water row with field `*` in the sample map is the mean of all usable irrigation-water injections and is paired with every soil field - use it when the water injections cannot be matched to soil fields one by one. Always check the scan polarity first: the Israel soil files in the folder named 'Soil Positive' are NEGATIVE-mode acquisitions, so the soil screen has to be run on negative water identifications (see `work/pec_map_israel_soil_neg.csv`). Soil injections run ~0.3 min earlier than the water run (automatic RT shift handles it).

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
