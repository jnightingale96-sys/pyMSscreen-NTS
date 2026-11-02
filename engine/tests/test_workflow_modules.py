"""Calibration model choice, peak rating / presence confidence, replicate rule, sample-sheet guessing, CompTox parsing."""
import numpy as np
import pandas as pd

from hrms_nts import calibration as CAL, quality as QL, qa as QA


def test_quadratic_chosen_only_when_curved():
    x = np.array([1, 2, 5, 10, 20, 50, 100.0])
    curved = 1000 * x - 4 * x ** 2 + 50
    ch, lin, quad, why = CAL.choose(x, curved)
    assert ch["model"] == "quadratic" and "significant" in why
    assert abs(CAL.invert(ch, [CAL.predict(ch, 30.0)])[0] - 30.0) < 1e-6
    straight = 1000 * x + 50 + np.random.default_rng(1).normal(0, 1, len(x))
    ch2, *_ , why2 = CAL.choose(x, straight)
    assert ch2["model"] == "linear"
    few, *_ , why3 = CAL.choose(x[:4], curved[:4])
    assert few["model"] == "linear" and "needs" in why3


def test_calibrate_table_and_quantify():
    x = np.array([1, 2, 5, 10, 20.0]); y = 2000 * x + 100
    cal = pd.DataFrame(dict(name="A", conc=x, response=y))
    summ, fits = CAL.calibrate_table(cal)
    assert summ.loc[0, "model"] == "linear" and summ.loc[0, "r2"] > 0.9999 and not summ.loc[0, "warning"]
    q = CAL.quantify(fits, pd.DataFrame(dict(name=["A", "B"], file=["s1", "s1"], response=[2000 * 7 + 100, 5.0])))
    assert abs(q.loc[0, "conc"] - 7) < 1e-6 and q.loc[1, "flag"] == "no calibration curve"


def test_peak_rating_gaussian_beats_noise():
    rt = np.arange(0, 2, 0.013)
    g = 2e4 * np.exp(-0.5 * ((rt - 1) / 0.03) ** 2)
    good = QL.peak_rating(rt, g, 0.9, 1.1, sn=130)
    noisy = QL.peak_rating(rt, np.random.default_rng(0).uniform(0, 400, len(rt)), 0.9, 1.1, sn=2)
    assert good["peak_rating"] >= 8 and noisy["peak_rating"] < 4


def test_presence_confidence_threshold():
    hi, c = QL.presence_confidence(9, 200, 1.0, 5, 0.9, 0.85, 5, 4, 1.0, 50)
    lo, _ = QL.presence_confidence(2, 4, 7.5, 5, None, 0.2, 1, 0, 0.33, 1.2)
    assert hi >= 70 > lo and set(c) == {"rating", "sn", "ppm", "iso", "ms2", "replicates", "blank"}


def test_replicate_rule():
    assert QA.required(3, "2/3") == 2 and QA.required(2, "2/3") == 2 and QA.required(4, "0.5") == 2 and QA.required(3, "1") == 1
    sheet = pd.DataFrame(dict(file=["a1", "a2", "a3", "b1"], type=["sample"] * 3 + ["blank"], group=["A", "A", "A", "B"]))
    al = pd.DataFrame(dict(sample=["a1", "a2", "a1", "b1"], consensus_id=["C1", "C1", "C2", "C1"], area=[10, 12, 5, 1],
                           blank_flag=[False, False, False, False]))
    p = QA.replicate_presence(al, QA.read_sample_sheet(sheet), "2/3").set_index("consensus_id")
    assert bool(p.loc["C1", "present"]) and not bool(p.loc["C2", "present"]) and p.loc["C1", "n_injections"] == 3


def test_guess_sample_sheet():
    s = QA.guess_sample_sheet(["WW1_rep1.mzML", "WW1_rep2.mzML", "Blank_1.mzML", "QC_pool_3.mzML", "cal_5ngmL.mzML", "eq5.mzML"]).set_index("file")
    assert s.loc["WW1_rep1", "group"] == s.loc["WW1_rep2", "group"] == "WW1"
    assert s.loc["Blank_1", "type"] == "blank" and s.loc["eq5", "type"] == "blank" and s.loc["QC_pool_3", "type"] == "qc"
    assert s.loc["cal_5ngmL", "type"] == "cal" and s.loc["cal_5ngmL", "conc"] == 5


def test_comptox_prefers_experimental(monkeypatch, tmp_path):
    from hrms_nts import comptox
    ct = comptox.CompTox(api_key="x", cache_path=str(tmp_path / "c.json"), delay=0)
    pages = {"/chemical/search/equal/": [{"dtxsid": "DTXSID1", "preferredName": "Carbamazepine"}],
             "/chemical/property/summary/": [{"propName": "LogKow: Octanol-Water", "experimentalCount": 3, "experimentalMedian": 2.45, "predictedMedian": 2.2},
                                             {"propName": "pKa Basic Apparent", "experimentalCount": 0, "predictedMedian": 0.4, "predictedCount": 2}]}
    monkeypatch.setattr(ct, "_get", lambda path: next(v for k, v in pages.items() if path.startswith(k)))
    p = ct.properties("FFGPTBGBLSHEPO-UHFFFAOYSA-N")
    assert p["logkow"] == 2.45 and "experimental" in p["logkow_source"] and p["pka_base"] == 0.4 and "predicted" in p["pka_base_source"]
    ct.save(); assert (tmp_path / "c.json").exists()
