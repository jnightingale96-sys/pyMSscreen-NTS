"""Profile centroiding, DIA pseudo-MS2, isotope-labelled masses and the targeted-screen helpers."""
import numpy as np
import pandas as pd

from hrms_nts.centroid import centroid_profile
from hrms_nts.mzml_io import Run
from hrms_nts import targeted as T
from hrms_nts.dia import pseudo_ms2


def _profile(centres, heights, spacing_ppm=7.0, fwhm_ppm=12.0, floor=20.0, rng=None):
    """Sparse TOF-like profile (zero points omitted) with Poisson-ish noise on the counts."""
    rng = rng or np.random.default_rng(0)
    mz, it = [], []
    for c, h in zip(centres, heights):
        sd = c * fwhm_ppm * 1e-6 / 2.355
        x = c + np.arange(-12, 13) * c * spacing_ppm * 1e-6
        y = h * np.exp(-0.5 * ((x - c) / sd) ** 2)
        y = rng.poisson(np.clip(y, 0, None)).astype(float)
        k = y >= floor
        mz.append(x[k]); it.append(y[k])
    return np.concatenate(mz), np.concatenate(it)


def test_centroid_profile_accuracy_and_no_split():
    centres = [118.08626, 342.25876, 794.89592]
    mz, it = _profile(centres, [3000, 800, 400])
    m, i = centroid_profile(mz, it, min_height=150)
    assert len(m) == 3                                   # one centroid per peak even at low counts
    for c in centres:
        assert np.min(np.abs(m - c) / c * 1e6) < 2.0     # within 2 ppm


def test_labelled_masses_match_vendor_q1():
    # values from a vendor IS sheet (all correct except the hydroxycotinine typo, see VALIDATION.md)
    for f, q1 in (("C8H5[2H]4NO2", 156.0957), ("C5[13C]3H10N4O2", 198.0977), ("C6[13C]H6[2H]2ClN3O4S2", 300.9877),
                  ("C10H6D9N", 159.1842), ("C18H21[2H]3I3N3O8", 794.8959)):
        d = T.parse_formula(f)
        assert abs(T.mono(d) + T.ADDUCTS["positive"]["[M+H]+"] - q1) < 1e-4
    assert T.parse_formula("C10H9[2H]3N2O2") == {"C": 10, "H": 9, "D": 3, "N": 2, "O": 2}
    assert T.unlabel(T.parse_formula("C5[13C]3H10N4O2")) == {"C": 8, "H": 10, "N": 4, "O": 2}


def test_names_and_labels():
    assert T.parent_name("Fentanyl-D5 oxalate") == "Fentanyl"
    assert T.parent_name("6-Acetylmorphine-D3 (MAM-D3)") == "6-Acetylmorphine"
    assert T.parent_name("Anatabine-2,4,5,6-D4") == "Anatabine"
    assert T.label_from_text("cis-Tramadol-¹³C,D3") == (3, 1, 0)
    assert T.label_from_text("Carbamazepine-D10") == (10, 0, 0)
    assert "ritalinic acid" in T._name_variants(T.parent_name("(±)-threo-Ritalinic acid-D10 HCl"))


def test_load_vendor_sheet_and_mass_check(tmp_path):
    p = tmp_path / "is.csv"
    pd.DataFrame({"Name": ["Acetaminophen-D4", "trans-3'-Hydroxycotinine-D3", "Morphine-D3", "Hydromorphone-D3"],
                  "Chemical formula": ["C8H5[2H]4NO2", "C10H9[2H]3N2O2", "C17H16[2H]3NO3", "C17H16[2H]3NO3"],
                  "Precursor (Q1) Mass (Da)": [156.0957, 330.111977, 289.1626, 289.1626]}).to_csv(p, index=False)
    c = T.load_compounds(str(p))
    assert list(c.role.unique()) == ["IS"]
    err = (c.q1_list - (c.mono + T.ADDUCTS["positive"]["[M+H]+"])).abs()
    assert err.iloc[0] < 1e-4 and err.iloc[1] > 100          # the hydroxycotinine Q1 typo is caught
    assert (c.isobaric_group != "").sum() == 2                # morphine-D3 / hydromorphone-D3 share a formula


def test_label_aware_cosine_matches_deuterium_shifted_fragments():
    lib = dict(mz=np.array([105.0699, 188.1434, 281.2011]), it=np.array([30.0, 100.0, 10.0]), pmz=337.2274)
    obs = [dict(mz=105.0699, inten=30.0), dict(mz=193.1748, inten=100.0), dict(mz=281.2011, inten=10.0)]   # 188 kept 5 D
    s_shift = T.label_cosine(obs, lib, 5, 0, 342.2588, shift=True)
    s_plain = T.label_cosine(obs, lib, 5, 0, 342.2588, shift=False)
    assert s_shift[0] > 0.99 and s_shift[1] == 3 and "+5D" in s_shift[2]
    assert s_plain[0] < 0.5


def test_fit_calibration_weighted():
    x = np.array([1, 2, 5, 10, 20, 50.0]); y = 1000 * x + 50
    f = T.fit_calibration(x, y, "1/x")
    assert abs(f["slope"] - 1000) < 1e-6 and f["r2"] > 0.999999


def _dia_run():
    """Synthetic DIA run: precursor 300.1 (MS1) co-eluting with fragments 150.05 and 200.07 in window 275-325; an interfering
    fragment 180.0 in the same window elutes 0.4 min later and must not be assigned."""
    rt = np.arange(0, 3, 0.01)
    gauss = lambda c, h: h * np.exp(-0.5 * ((rt - c) / 0.03) ** 2)
    p = gauss(1.5, 5e4); fa = gauss(1.5, 2e4); fb = gauss(1.5, 8e3); fi = gauss(1.9, 3e4)
    mz1, in1, off1 = [], [], [0]
    for k in range(len(rt)):
        mz1.append(np.array([300.1])); in1.append(np.array([p[k] + 1], np.float32)); off1.append(off1[-1] + 1)
    dmz, dint, doff, drt, dwin = [], [], [0], [], []
    for k in range(len(rt)):
        m = np.array([150.05, 180.0, 200.07]); i = np.array([fa[k] + 1, fi[k] + 1, fb[k] + 1], np.float32)
        dmz.append(m); dint.append(i); doff.append(doff[-1] + 3); drt.append(rt[k] + 0.002); dwin.append(300.0)
    run = Run("x.mzML", "x", "positive", rt, np.array(off1), np.concatenate(mz1), np.concatenate(in1))
    run.acquisition = "DIA"
    run.dia = dict(rt=np.array(drt), win=np.array(dwin), off=np.array(doff), mz=np.concatenate(dmz), inten=np.concatenate(dint),
                   windows=[(300.0, 275.0, 325.0, 30.0)])
    return run


def test_dia_pseudo_ms2_keeps_coeluting_rejects_interference():
    run = _dia_run()
    w, fr = pseudo_ms2(run, 300.1, 1.4, 1.6)
    got = sorted(round(f["mz"], 2) for f in fr)
    assert w[0] == 300.0 and got == [150.05, 200.07]


def test_integrate_and_rating():
    rt = np.arange(0, 3, 0.013)
    I = 2e4 * np.exp(-0.5 * ((rt - 1.5) / 0.03) ** 2)
    g = T.integrate(rt, I, 1.5, 150)
    assert abs(g["apex_rt"] - 1.5) < 0.02 and g["n_points"] >= 8 and T.integration_rating(g) == "Good"
    area_true = 2e4 * 0.03 * np.sqrt(2 * np.pi) * 60
    assert abs(g["area"] - area_true) / area_true < 0.05
