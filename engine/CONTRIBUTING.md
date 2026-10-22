# Working on the toolkit as a team

## Set up
```
git clone <repo>  &&  cd HRMS_NTS_Toolkit
install_nts.bat                       # or: python -m venv .venv && .venv\Scripts\activate && pip install -e .[dev]
copy config.local.example.yaml config.local.yaml    # edit: paths to SusDat, PubChemLite, MSP libraries, (Java + MetFrag jar optional)
hrms-nts doctor                       # shows what is found / missing
hrms-nts build-library                # one-off MS/MS library (SQLite)
pytest tests                          # must be green
```
Data (SusDat, PubChemLite, libraries, mzML) is never committed (`.gitignore`); each person points `config.local.yaml` at their copy.

## Code map
`hrms_nts/` - `mzml_io` (reader) -> `features` (MS1 picking, background flag) -> `level4` (mass/isotope match) -> `rt_prediction` -> `ms2` ->
`library` (spectral library) -> `insilico` (product ions, MetFrag) -> `levels` (confidence rules) -> `semiquant` -> `aggregate`/`report`.
`pipeline.py` orchestrates one sample; `cli.py` is the command line; `validate.py` scores a result folder against ground truth;
`provenance.py` stamps every result folder. All thresholds live in `config.yaml` - **no magic numbers in code**.

## Rules
* One topic per branch/merge request; `pytest tests` green; add or update a test with every behaviour change.
* Changing a **threshold default or level rule** needs: (a) a validation result on real data (`VALIDATION.md`), (b) an entry in `CHANGELOG.md`.
* Never hard-code personal paths; use `config.local.yaml`. Never commit data or credentials.
* Use `--set section.key=value` for parameter sweeps; the merged config is saved as `config_used.yaml` in every run.
* Bump `hrms_nts/__init__.py::__version__` for releases; results carry the version + git commit in `provenance.json`.

## Good first tasks (open)
* Replace the logIE proxy in `semiquant.logie_estimate` with a trained ionisation-efficiency model.
* Add DIA (SWATH / MS-all) MS2 deconvolution.
* Isomer-aware in-silico benchmark (`tools/benchmark_insilico.py`).
* Parallelise samples (`multiprocessing`) and the in-silico stage (currently ~5 min / sample / 2,800 features).
* GUI front-end over `pipeline.process_sample`.
