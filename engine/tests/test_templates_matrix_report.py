import os
import numpy as np
import pandas as pd
import pytest

from hrms_nts import templates as T
from hrms_nts import matrix as MX


def test_templates_written_and_examples_pass_the_checker(tmp_path):
    p = T.write_templates(str(tmp_path))
    assert set(p) == {"sample_sheet_template.csv", "calibration_list_template.csv", "is_list_template.csv", "README_templates.txt"}
    assert (T.check_compound_list(p["is_list_template.csv"], "is").check == "OK").all()
    assert (T.check_compound_list(p["calibration_list_template.csv"], "cal").check == "OK").all()
    s = pd.read_csv(p["sample_sheet_template.csv"])
    assert set(s.type) <= set(T.SAMPLE_TYPES)


def test_checker_flags_wrong_mass_missing_conc_and_bad_is_name():
    df = pd.DataFrame([dict(name="X-D5", formula="C16H8[2H]5ClN2O", monoisotopic_mass=289.1, role="IS", conc=np.nan)])
    c = T.check_compound_list(df, "is").check.iloc[0]
    assert "ppm" in c and "no concentration" in c
    cal = pd.DataFrame([dict(name="A", formula="C15H12N2O", monoisotopic_mass=236.09496, role="target", is_name="nope", L1=1, L2=2)])
    c = T.check_compound_list(cal, "cal").check.iloc[0]
    assert "is_name" in c and "levels" in c


def test_level_concentrations_long_table():
    lv = T.level_concentrations(T.CAL_LIST)
    assert lv[(lv.name == "Carbamazepine") & (lv.cal_level == "L3")].conc.iloc[0] == 5


def test_matrix_effect_and_map(tmp_path):
    sheet = pd.DataFrame(dict(file=["S1", "M1", "M2"], type=["solvent_std", "matrix_spike", "matrix_spike"]))
    rows = []
    for i, (n, rt, me) in enumerate([("a", 2.0, -0.5), ("b", 5.0, 0.0), ("c", 9.0, 0.3)]):
        rows.append(dict(name=n, file="S1", int_area=1000.0, int_rating="Good", int_apex_rt=rt))
        rows += [dict(name=n, file=f, int_area=1000.0 * (1 + me), int_rating="Good", int_apex_rt=rt) for f in ("M1", "M2")]
    me = MX.matrix_effect(pd.DataFrame(rows), sheet)
    assert me.set_index("name").ME_pct.round(1).to_dict() == {"a": -50.0, "b": 0.0, "c": 30.0}
    assert me.set_index("name").effect["a"] == "suppression"
    mp = MX.me_map(me, window_min=1.0)
    assert MX.me_at_rt(mp, 2.0) == pytest.approx(-50.0)
    assert np.isnan(MX.me_at_rt(mp, 30.0))
    assert os.path.exists(MX.plot_me(me, mp, str(tmp_path / "me.png")))


def test_project_matrix_effect_and_reports(tmp_path):
    from hrms_nts.workflow import Project
    for f in ("S1", "M1"):
        open(tmp_path / f"{f}.mzML", "w").write("<mzML/>")
    P = Project(str(tmp_path), str(tmp_path / "out"))
    P.set_sample_sheet(pd.DataFrame(dict(file=["S1", "M1"], type=["solvent_std", "matrix_spike"])))
    a = pd.DataFrame([dict(name="a", file="S1", int_area=1000.0, int_rating="Good", int_apex_rt=3.0, level="2a"),
                      dict(name="a", file="M1", int_area=400.0, int_rating="Good", int_apex_rt=3.0, level="2a")])
    P._save("targeted", dict(assigned=a, list_path="", is_qc=pd.DataFrame(dict(name=["a"], area_rsd_pct=[5.0])),
                             injection_qc=pd.DataFrame(dict(file=["S1", "M1"], median_IS_response_rel=[1.0, 0.4])), xlsx=""))
    me = P.matrix_effect()
    assert me["table"].ME_pct.iloc[0] == pytest.approx(-60.0)
    rep = P.report()
    assert os.path.exists(rep["xlsx"]) and os.path.exists(rep["html"])
    h = open(rep["html"], encoding="utf-8").read()
    assert "Matrix effect" in h and "data:image/png;base64" in h
    import openpyxl
    wb = openpyxl.load_workbook(rep["xlsx"])
    assert wb.sheetnames[0] == "Cover" and "Matrix_effect" in wb.sheetnames
