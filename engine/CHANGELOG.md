# Changelog

## 0.2.0 (2026-10-20) - team-ready
* Layered configuration (`config.yaml` portable defaults <- `config.local.yaml` <- `$HRMS_NTS_CONFIG` <- `--set key=value`); no personal paths in the repo.
* `validate` command (recovery, ppm/RT accuracy, top-1 rate, false-identification rate with `absent.csv`, semi-quant fold error vs `true_conc_ug_ml`).
* Provenance (`provenance.json`, `config_used.yaml`) written with every run.
* pytest suite (22 tests incl. end-to-end CLI run on a generated mzML), `pyproject.toml` with `hrms-nts` entry point, CI definitions, `.gitignore`.
* Docs: README, VALIDATION.md (protocol + open items), CONTRIBUTING.md.

## 0.1.0 (2026-10-20) - first working version
* Levels 5-4-3-2b-2a-1, features/Level 4/RT/MS2/library/in-silico (built-in + MetFrag)/semi-quant/prescreen/reporting.
* Background-occupancy flag; adduct-consistent library matching; low-mass and background caps on in-silico levels; decoy-adjusted explained-ion counts.
