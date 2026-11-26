import pickle
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("rdkit")
from hrms_nts.rt_prediction import GlobalRT, RTModel, global_descriptors, GLOBAL_DESC_NAMES


class LogPIndex:
    """Stand-in for the SMRT model: global index = 3 x logP + 10 (logP is the first descriptor)."""
    def predict(self, X):
        return 3.0 * np.asarray(X)[:, 0] + 10.0


SMI = ["NC(=O)N1c2ccccc2C=Cc2ccccc12", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", "COc1ccc(cc1)C(CN(C)C)C1(O)CCCCC1",
       "CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c1ccccc1", "CNC1(CCCCC1=O)c1ccccc1Cl", "CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1",
       "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"]


def test_global_descriptors_shape():
    d = global_descriptors(SMI[0])
    assert d is not None and len(d) == len(GLOBAL_DESC_NAMES)
    assert global_descriptors("not a smiles") is None


def test_transfer_recovers_a_linear_method(tmp_path):
    p = tmp_path / "g.pkl"
    pickle.dump(dict(model=LogPIndex(), cv="test"), open(p, "wb"))
    g = GlobalRT(str(p))
    idx = g.index(SMI)
    st = pd.DataFrame(dict(name=[f"c{i}" for i in range(len(SMI))], smiles=SMI, rt=0.5 * idx - 2.0))
    m = RTModel("positive", tol=0.5).fit_global(st, g)
    assert m.n_anchor == len(SMI) and m.glob_rmse < 0.3
    h = pd.DataFrame(dict(smiles=SMI[:3], rt_apex=(0.5 * idx[:3] - 2.0) + np.array([0.0, 0.1, 5.0])))
    a = m.annotate(h)
    assert np.allclose(a.rt_pred_global, 0.5 * idx[:3] - 2.0, atol=0.3)
    assert list(a.rt_agree) == [True, True, False] and bool(a.rt_conflict.iloc[2])


def test_too_few_anchors_leaves_model_unset(tmp_path):
    p = tmp_path / "g.pkl"
    pickle.dump(dict(model=LogPIndex()), open(p, "wb"))
    st = pd.DataFrame(dict(name=["a", "b", "c"], smiles=SMI[:3], rt=[1.0, 2.0, 3.0]))
    m = RTModel("positive").fit_global(st, GlobalRT(str(p)))
    assert m.glob_map is None and np.isnan(m.predict_global(SMI[:2])).all()
