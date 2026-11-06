import numpy as np
import pandas as pd
import pytest

pytest.importorskip("rdkit")
from hrms_nts import surrogate as SG  # noqa: E402

SMI = {
    "Carbamazepine-D10": "NC(=O)N1c2ccccc2C=Cc2ccccc12",
    "Diazepam-D5": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "Venlafaxine-D6": "COc1ccc(cc1)C(CN(C)C)C1(O)CCCCC1",
    "Methadone-D9": "CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c1ccccc1",
    "Cocaine-D3": "COC(=O)C1C(CC2CCC1N2C)OC(=O)c1ccccc1",
    "Ketamine-D4": "CNC1(CCCCC1=O)c1ccccc1Cl",
    "Fentanyl-D5": "CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1",
    "Indomethacin-D4": "COc1ccc2n(C(=O)c3ccc(Cl)cc3)c(C)c(CC(=O)O)c2c1",
}
RT = {"Carbamazepine-D10": 8.2, "Diazepam-D5": 10.4, "Venlafaxine-D6": 4.95, "Methadone-D9": 7.45, "Cocaine-D3": 3.6,
      "Ketamine-D4": 2.83, "Fentanyl-D5": 6.1, "Indomethacin-D4": 12.4}


def _comp():
    return pd.DataFrame([dict(name=n, smiles=s, avg_mw=300.0) for n, s in SMI.items()])


def _assigned(scale_sample=0.5):
    rows = []
    for i, n in enumerate(SMI):
        a = 1e5 * (1 + 0.1 * i)
        rows.append(dict(name=n, file="ISMIX", int_area=a, int_rating="Good", int_apex_rt=RT[n], level="2a"))
        rows.append(dict(name=n, file="S1", int_area=a * scale_sample, int_rating="Good", int_apex_rt=RT[n], level="2a"))
    return pd.DataFrame(rows)


def test_physchem_classes_and_curated_override():
    p = SG.physchem_table(_comp()).set_index("name")
    assert p.loc["Methadone-D9", "ionisable_class"] == "aliphatic_amine"
    assert p.loc["Indomethacin-D4", "ionisable_class"] == "acid"
    assert p.loc["Carbamazepine-D10", "ionisable_class"] == "neutral"
    assert p.loc["Methadone-D9", "f_cation"] > 0.99 and p.loc["Carbamazepine-D10", "f_cation"] == 0
    cur = pd.DataFrame(dict(smiles=[SMI["Methadone-D9"]], logKow=[3.93], pKa_b_final=[9.07], pKa_a_final=[np.nan]))
    q = SG.physchem_table(_comp(), cur).set_index("name")
    assert q.loc["Methadone-D9", "logKow"] == pytest.approx(3.93)
    assert q.loc["Methadone-D9", "logKow_source"] == "curated"


def test_is_response_factors_flags_void_and_conc():
    a = _assigned()
    a.loc[a.name == "Ketamine-D4", "int_apex_rt"] = 1.2          # void
    rf = SG.is_response_factors(a, _comp(), ["ISMIX"], conc_ug_ml=0.1, t0=1.0, void_margin=0.5).set_index("name")
    assert not rf.loc["Ketamine-D4", "use"] and "void" in rf.loc["Ketamine-D4", "reason"]
    assert rf.loc["Carbamazepine-D10", "use"]
    assert rf.loc["Carbamazepine-D10", "logRF"] == pytest.approx(np.log10(1e5 / (0.1 / 300.0)))


def test_semiquant_is_recovers_known_concentration_and_matrix_correction():
    a = _assigned(scale_sample=0.5)                              # 50 % suppression in sample S1
    comp = _comp()
    rf = SG.is_response_factors(a, comp, ["ISMIX"], conc_ug_ml=0.1)
    props = SG.physchem_table(comp)
    rec = SG.is_recovery(a, rf, ["ISMIX"])
    assert np.allclose(rec.recovery, 0.5)
    # target = carbamazepine (unlabelled) at 0.1 ug/mL in S1 with the same molar RF as the IS median, suppressed 50 %
    med = 10 ** np.median(rf.logRF)
    mw = props.set_index("name").loc["Carbamazepine-D10", "MW"]
    tgt = pd.DataFrame([dict(name="Carbamazepine-D10", smiles=SMI["Carbamazepine-D10"], rt=8.2, sample="S1",
                             area=0.1 / mw * med * 0.5)])
    out = SG.semiquantify_is(tgt, rf, props, method="median", recovery=rec)
    assert out.sq_conc_ug_ml.iloc[0] == pytest.approx(0.1, rel=1e-6)
    assert out.sq_recovery.iloc[0] == pytest.approx(0.5)
    assert round(out.sq_high.iloc[0] / out.sq_conc_ug_ml.iloc[0], 6) in (3.0, 5.0, 10.0)


def test_void_and_unrepresented_class_widen_uncertainty():
    comp = _comp()
    rf = SG.is_response_factors(_assigned(), comp, ["ISMIX"], conc_ug_ml=0.1)
    rf.loc[rf.name == "Indomethacin-D4", "use"] = False          # no acid among the usable IS
    props = SG.physchem_table(comp)
    tgt = pd.DataFrame([dict(name="Indomethacin-D4", smiles=SMI["Indomethacin-D4"], rt=1.1, sample="S1", area=1e5)])
    out = SG.semiquantify_is(tgt, rf, props, k=3)
    assert "void-eluting" in out.sq_flags.iloc[0] and "not represented" in out.sq_flags.iloc[0]
    assert out.sq_high.iloc[0] / out.sq_conc_ug_ml.iloc[0] == pytest.approx(10.0)


def test_rt_property_model_interval_contains_training_points():
    props = SG.physchem_table(_comp())
    rt = props.name.map(RT).to_numpy()
    m = SG.RTPropertyModel(t0=1.0).fit(rt, props)
    ok = m.plausible(rt, props, level=0.99)
    assert ok.mean() >= 0.75
    iv = m.predict_interval(props)
    assert (iv.rt_lo < iv.rt_hi).all()


def test_project_semiquant_is_end_to_end(tmp_path):
    from hrms_nts.workflow import Project
    lst = tmp_path / "is_list.csv"
    pd.DataFrame([dict(name=n, smiles=s, label="D5", role="IS", conc=0.1) for n, s in SMI.items()]).to_csv(lst, index=False)
    P = Project(str(tmp_path), str(tmp_path / "out"))
    a = _assigned(0.5)
    a["level"] = "2a"
    P._save("targeted", dict(assigned=a, list_path=str(lst)))
    ident = pd.DataFrame([dict(sample="S1", name="Carbamazepine", smiles=SMI["Carbamazepine-D10"], level="2a",
                               rt_apex=8.3, int_area=5e4)])
    P._save("qa", dict(identifications=ident))
    res = P.semiquant(method="is", is_mix="ISMIX")
    t = res["table"]
    assert res["method"] == "is"
    assert res["is_rf"].use.sum() == len(SMI)
    assert np.isfinite(t.sq_conc.iloc[0]) and t.sq_recovery.iloc[0] == pytest.approx(0.5)
    assert "rt_property_plausible" in t.columns
