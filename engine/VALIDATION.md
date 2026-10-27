# Validation protocol (for the team)

The tool reports *evidence* for each level; it does not prove identifications. Before results are used in a paper, the team should run
the checks below on **your own data** and record the outcome in `validation/results/` (commit the CSV summaries + `provenance.json`).

## 0. Automated checks (every change)
`pytest tests` (or `run_tests.bat`) - 22 tests: feature picking, ppm matching, isotopes, cosine, level rules, config layering,
semi-quant arithmetic, MetFrag wrapper (skipped if Java/jar absent) and an end-to-end CLI run on a generated mzML.
A change that alters a **level rule or threshold** must update/extend `tests/test_core.py::test_level_logic` in the same commit.

## 1. Ground-truth dataset (needed once per LC-MS method / matrix)
| Sample set | Purpose |
|---|---|
| Matrix blank (unspiked, same extraction) | false-identification rate; background behaviour |
| Matrix spiked with 50-100 authentic standards across classes/polarity, >= 3 levels (e.g. 0.5/5/50 ng/mL) | recovery, ppm/RT accuracy, level reached, semi-quant fold error |
| Solvent blank / procedural blank | blank filter |
| Real samples with confirmed compounds (e.g. carbamazepine, 4-nitrophenol) | sanity check |
Acquire in DDA with the same method as the real samples. Put the standards' RTs in `user_data/standards.csv` **from a separate injection or a
separate half of the standards** (do not train the RT model on the compounds you then test it on).

## 2. Run and score
```
hrms-nts run <mzML folder> --matrix <water|grain|soil|extract> --insilico both --out validation/results/<date>_<set>
hrms-nts validate validation/results/<date>_<set> expected.csv --absent absent.csv
```
`expected.csv`: `name,smiles|inchikey,sample,adduct,rt,min_level,true_conc_ug_ml` (see `examples/expected_example.csv`).
`absent.csv`: compounds known **not** to be in the sample (not spiked / in blank) - gives the false-identification rate per level.

## 3. Metrics reported by `validate` and suggested acceptance criteria (PROPOSALS - the team must agree and fix them before testing)
| Metric | Meaning | Suggested target |
|---|---|---|
| `recovery_pct` | expected compounds found at >= Level 4 | >= 90 % at the spike level of interest |
| `reached_level<=2a/1` | share reaching library / standard confirmation | report, no target |
| `median_abs_ppm` | mass accuracy | <= 2 ppm |
| `median_abs_rt_diff_min` | RT agreement (needs `rt` in expected.csv) | report |
| `top1_rate_pct` | correct structure ranked first within its feature | >= 80 % |
| `absent_flagged_at_<=L2b` / `<=L2a` / `<=L1` | false identifications | 0 at L1 and L2a; <= 1-2 % at L2b |
| `semiquant_median_fold_error`, `within_3x_pct`, `within_10x_pct` | semi-quant accuracy vs spiked concentration | report; literature value for single-reference semi-quant is typically within ~1 order of magnitude |

## 4. What still needs validating (known open items)
1. **Absolute semi-quant**: `config.yaml` `single_reference.area` came from Compound Discoverer; re-measure with `hrms-nts calibrate` on this toolkit, then check fold error on spiked compounds (class-specific bias is expected: ionisation efficiency varies by orders of magnitude).
2. **Background flag** (`ms1.bg_occupancy`, `ms1.bg_peak_to_median`): 57 % of features were flagged in the first wastewater test; sweep the thresholds against the blank and spiked data and record the recovery/FP trade-off (`--set ms1.bg_occupancy=0.5`).
3. **Library thresholds** (`levels.lib_cosine`, `lib_min_matched`, `library.frag_tol_*`): instrument/collision-energy dependent; check false 2a rate on `absent.csv`. GNPS contains lower-quality spectra than MassBank - compare libraries (`build-library --only MassBank`).
4. **Level 2b**: built-in scorer benchmark (`tools/benchmark_insilico.py`): AUC 0.84 (pos) / 0.80 (neg), decoys were not isomers. Repeat with isomeric decoys and with `--insilico both`.
5. **RT models** (`rt.tolerance_min` = 0.5 min): need >= 8 standards; report the cross-validated RMSE printed in `run.log` and whether +-0.5 min is appropriate (`rt.adaptive_k`).
6. **Isotope scoring**, **adduct handling** ([M+HCOO]-, [M+Cl]-, [M+CH3COO]- produce many chance candidates for small molecules - see the 120 Da rule) and **Level 1 definition** (`levels.level1_require_inhouse_ms2`).
7. DIA (SWATH) support (0.4.0) is validated for the targeted screen (below); the NTS `run` on DIA data still needs a spiked ground-truth
   set (pseudo-MS2 recovery and false 2a/2b rate vs DDA on the same samples).

## 4b. Targeted screen on SCIEX X500R SWATH data (2026-10-27)
Data: 75 deuterated / 13C drug internal standards (`examples/drug_IS_pos_75.csv`), SCIEX X500R, ESI+, profile mzML, SWATH 11 windows,
files BLK1 (IS-spiked QC), TWB_2010pool_1, TWB_2024pool_1 (wastewater pools), eq5 (solvent blank, injected last).
Reference: an independent hand-run analysis of the same files (same algorithms, different code; DrugIS_simple_Q1_Q4.xlsx).
* `targeted` reproduced it exactly: 166/166 compound-file detections, identical RT (0 differences > 0.05 min), identical areas (all within 1 %),
  identical level for 165/166 (Nitrazepam-D5 in one pool 3 vs 4).
* Q1: median -0.80 ppm, 98 % within 5 ppm (max 5.6 ppm, the most intense peaks - TOF detector saturation); Q4: SD ~6 ppm at SWATH ion counts.
* Levels (all files): 2a 63, 2b 49, 3 29, 4 12, isobaric-unresolved 13 (morphine/hydromorphone/norcodeine-D3, ephedrine/pseudoephedrine-D3:
  co-elute at the void), Q1 > 5 ppm 3. Amitriptyline-D3 vs EDDP-D3 and oxymorphone-D3 vs noroxycodone-D3 split correctly by library match /
  batch RT. Not found: the 7 ESI- compounds (flagged by the acid rule), nicotine/anabasine/anatabine-D4, atenolol-D7, cotinine-D3, OH-cotinine-D3.
* Blank (eq5): 58/62 no signal, 4 single-scan spikes; the only "significant" is the weak 17-min carboxy-ibuprofen-D3 ESI+ peak, which is
  background (same signal in the blank) - flagged in Detected.
* Scan gap 7.20-16.25 min in TWB_2024pool_1 reported automatically.
* Profile centroiding: synthetic test (`tests/test_targeted.py`) < 2 ppm, one centroid per peak; DIA pseudo-MS2 keeps co-eluting fragments and
  rejects a same-window interference eluting 0.4 min later.

## 5. Record every validation run
Copy `validation_summary.csv`, `provenance.json`, `config_used.yaml` into `validation/results/<date>_<set>/` and add one line to `validation/LOG.md`
(date, who, dataset, toolkit commit, parameters changed, headline metrics, decision).
