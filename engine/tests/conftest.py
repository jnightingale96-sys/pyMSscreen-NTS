import os, sys
import numpy as np, pandas as pd, pytest
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hrms_nts.mzml_io import Run
from hrms_nts.database import CompoundDB
from mzml_writer import synthetic_scans, CBZ


@pytest.fixture(scope="session")
def fake_run():
    scans = synthetic_scans()
    ms1 = [s for s in scans if s["level"] == 1]
    rt = np.array([s["rt_min"] for s in ms1]); mzs = [np.asarray(s["mz"])[np.argsort(s["mz"])] for s in ms1]
    its = [np.asarray(s["int"], np.float32)[np.argsort(s["mz"])] for s in ms1]
    off = np.concatenate([[0], np.cumsum([len(m) for m in mzs])])
    ms2 = [dict(scan_id=f"s{i}", rt=s["rt_min"], precursor_mz=s["precursor_mz"], precursor_int=5e7, iso_width=1.0, ce=30.0,
                mz=np.asarray(s["mz"]), inten=np.asarray(s["int"], np.float32), tic=float(np.sum(s["int"]))) for i, s in enumerate(scans) if s["level"] == 2]
    return Run("fake.mzML", "fake_POS", "positive", rt, off, np.concatenate(mzs), np.concatenate(its).astype(np.float32), ms2)


@pytest.fixture()
def tiny_db():
    df = pd.DataFrame({"db_id": ["X1", "X2"], "name": ["Carbamazepine", "Decoy"], "formula": ["C15H12N2O", "C15H12N2O"],
                       "exact_mass": [236.094963, 236.094963 + 0.01], "smiles": ["NC(=O)N1c2ccccc2C=Cc2ccccc12", "CC"],
                       "inchikey": ["FFGPTBGBLSHEPO-UHFFFAOYSA-N", "AAAAAAAAAAAAAA-UHFFFAOYSA-N"]})
    df = df.assign(cas="", source="test", pred_rti_pos=np.nan, pred_rti_neg=np.nan, direct_precursor_mz=np.nan, db_level="",
                   ik1=lambda d: d.inchikey.str.split("-").str[0])
    return CompoundDB(df)


@pytest.fixture(scope="session")
def cfg():
    from hrms_nts.common import load_config
    return load_config(os.path.join(ROOT, "config.yaml"))     # portable defaults only (config.local.yaml is merged, harmless)
