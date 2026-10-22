"""End-to-end: synthetic mzML on disk -> CLI `run` -> outputs -> CLI `validate`.  Needs no external databases or instrument data."""
import os, json, shutil, yaml, pandas as pd, pytest
from mzml_writer import write_mzml, synthetic_scans
from hrms_nts.cli import main

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture()
def workspace(tmp_path):
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config.yaml"), encoding="utf-8"))
    cfg["paths"].update(susdat=None, pubchemlite=None, user_compounds="compounds.csv", library_msp=[], library_db="work/none.sqlite",
                        standards=None, calibrants=None, cache_dir="work/cache", metfrag_jar=None)
    (tmp_path / "config.yaml").write_text(yaml.safe_dump(cfg))
    (tmp_path / "compounds.csv").write_text("ID,Name,SMILES\nC1,Carbamazepine,NC(=O)N1c2ccccc2C=Cc2ccccc12\nC2,Decoy,CC\n")
    write_mzml(str(tmp_path / "demo_POS_01.mzML"), synthetic_scans())
    return tmp_path


def test_run_and_validate(workspace):
    out = workspace / "out"
    main(["--config", str(workspace / "config.yaml"), "run", str(workspace / "demo_POS_01.mzML"), "--out", str(out), "--no-library",
          "--semiquant", "single", "--matrix", "grain", "--insilico", "builtin"])
    for f in ("features.csv", "all_candidates.csv", "top_candidates_per_feature.csv", "compound_summary.csv", "NTS_results.xlsx", "provenance.json",
              "config_used.yaml", "level_summary.csv", "run.log"):
        assert (out / f).exists(), f
    prov = json.load(open(out / "provenance.json"))
    assert prov["toolkit_version"] and prov["inputs"][0]["sha256_partial"]
    top = pd.read_csv(out / "top_candidates_per_feature.csv", low_memory=False)
    cb = top[top.name.astype(str) == "Carbamazepine"]
    assert len(cb) == 1 and abs(cb.rt_apex.iloc[0] - 8.0) < 0.1 and abs(cb.ppm_error.iloc[0]) < 2
    assert str(cb.level.iloc[0]) in ("4", "3", "2b"), "no library/standard supplied -> cannot exceed 2b"
    assert cb.ms2_n_peaks.iloc[0] >= 3 and cb.bi_expl_n.iloc[0] >= 2
    assert (top.name.astype(str) == "Decoy").sum() == 0
    assert "conc_grain" in top.columns and cb.conc_grain.iloc[0] > 0
    # validate
    (workspace / "expected.csv").write_text("name,smiles,sample,adduct,rt,min_level\nCarbamazepine,NC(=O)N1c2ccccc2C=Cc2ccccc12,*,[M+H]+,8.0,4\n")
    main(["--config", str(workspace / "config.yaml"), "validate", str(out), str(workspace / "expected.csv")])
    summ = pd.read_csv(out / "validation_summary.csv", index_col=0)["value"]
    assert float(summ["recovery_pct"]) == 100.0
