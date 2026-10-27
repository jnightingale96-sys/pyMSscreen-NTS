"""pec-screen: integration, semi-quant arithmetic, blank logic, grain RT shift / MS2 identity, properties, fate grid, classification."""
import copy, os
import numpy as np, pandas as pd, pytest
from conftest import ROOT
from mzml_writer import write_mzml, synthetic_scans, CBZ
from hrms_nts import pecscreen as PS
from hrms_nts.ms2 import peaks_to_string
from hrms_nts.properties import resolve_properties, read_property_table, species_type, structure_flags
from hrms_nts.fate import FateRunner

FRAG_MZ = [167.0730, 180.0808, 192.0808, 194.0964, 237.1022]; FRAG_INT = [1e4, 2e4, 3e4, 1e5, 5e3]


def _scaled(scale=1.0, dt=0.0, seed=0):
    sc = synthetic_scans(seed=seed)
    for s in sc:
        s["rt_min"] = s["rt_min"] + dt
        if s["level"] == 1:
            mz = np.asarray(s["mz"]); it = np.asarray(s["int"], float).copy()
            it[(np.abs(mz - CBZ) < 0.01) | (np.abs(mz - CBZ - 1.00335) < 0.01)] *= scale
            s["int"] = it
    return sc


@pytest.fixture()
def pcfg(cfg, tmp_path):
    c = copy.deepcopy(cfg)
    c["paths"]["cache_dir"] = str(tmp_path / "cache")
    c["pec_screen"]["reference"].update(mzml=str(tmp_path / "std.mzML"), rt=None, mz=CBZ, conc_ug_mL=0.01, name="Carbamazepine")
    c["pec_screen"]["rt_shift_min_anchors"] = 1; c["pec_screen"]["rt_shift_min_apex"] = 3e5
    write_mzml(str(tmp_path / "std.mzML"), _scaled(1.0))
    return c


def _cand(ms2=True):
    peaks = peaks_to_string(FRAG_MZ, FRAG_INT) if ms2 else peaks_to_string([100.0, 120.0, 130.0, 140.0], [1e4, 1e4, 1e4, 1e4])
    return pd.DataFrame([dict(cand_id="FFGPTBGBLSHEPO_pos", name="carbamazepine", mz=CBZ, rt=8.0, polarity="positive", ms2_peaks=peaks, ik1="FFGPTBGBLSHEPO")])


def _map(tmp_path, rows):
    out = []
    for stem, matrix, field, role, scans in rows:
        p = str(tmp_path / f"{stem}.mzML"); write_mzml(p, scans)
        out.append(dict(path=p, stem=stem, matrix=matrix, field=field, role=role, polarity="positive", note="", scale=1.0))
    return pd.DataFrame(out)


def test_integrate_eic_gaussian_area_and_floor():
    rt = np.linspace(0, 10, 2001); y = 1e6 * np.exp(-0.5 * ((rt - 5) / 0.05) ** 2)
    a, r, apex = PS.integrate_eic(rt, y, 5.0, 0.25, 0.1, floor=0.0)
    assert abs(a - 1e6 * 0.05 * np.sqrt(2 * np.pi)) / a < 0.08 and abs(r - 5.0) < 0.01     # ~5 % low by design (2 % apex cut-off); cancels vs the standard
    assert PS.integrate_eic(rt, 0.01 * y, 5.0, 0.25, 0.1, floor=3e4)[0] == 0.0


def test_semiquant_arithmetic_water_and_grain(pcfg, tmp_path):
    ref = PS.reference_response(pcfg)
    assert ref["conc_ug_mL"] == 0.01 and ref["area"] > 1e5
    smap = _map(tmp_path, [("w1", "water", "A", "sample", _scaled(0.5, seed=1)), ("g1", "grain", "A", "sample", _scaled(0.25, seed=2))])
    D, _ = PS.quantify_injections(_cand(), smap, pcfg, ref)
    w = D[D["sample"] == "w1"].iloc[0]; g = D[D["sample"] == "g1"].iloc[0]
    assert w.conc == pytest.approx(0.005 * 120000, rel=0.05)           # 0.5 x std = 0.005 ug/mL -> x1000 x120 = 600 ng/L
    assert g.conc == pytest.approx(0.0025 * 166.6667, rel=0.05)         # 0.25 x std = 0.0025 ug/mL -> x0.5/3 x1000 = 0.417 ng/g
    assert w.id_status == "MS2 confirmed" and w.ms2_matched >= 4


def test_blank_check_and_identity_rules(pcfg, tmp_path):
    ref = PS.reference_response(pcfg)
    smap = _map(tmp_path, [("blank", "water", "", "blank", _scaled(0.2, seed=3)), ("hi", "water", "A", "sample", _scaled(1.0, seed=4)),
                           ("lo", "water", "B", "sample", _scaled(0.3, seed=5))])
    D, _ = PS.quantify_injections(_cand(), smap, pcfg, ref); D = PS.apply_blank_check(D, 3.0)
    det = D[D.role == "sample"].set_index("sample").detected
    assert bool(det["hi"]) and not bool(det["lo"])                       # 1.0/0.2 = 5x blank passes; 0.3/0.2 = 1.5x does not
    # a spectrum that contradicts the water identification is not counted
    smap2 = _map(tmp_path, [("hi2", "water", "A", "sample", _scaled(1.0, seed=6))])
    D2, _ = PS.quantify_injections(_cand(ms2=False), smap2, pcfg, ref); D2 = PS.apply_blank_check(D2, 3.0)
    assert D2.iloc[0].id_status == "MS2 conflict" and not bool(D2.iloc[0].detected)


def test_grain_rt_shift_is_estimated_from_ms2_anchors(pcfg, tmp_path):
    ref = PS.reference_response(pcfg)
    smap = _map(tmp_path, [("gshift", "grain", "A", "sample", _scaled(0.5, dt=0.2, seed=7))])
    D, sh = PS.quantify_injections(_cand(), smap, pcfg, ref)
    assert sh.iloc[0].rt_shift_min == pytest.approx(0.2, abs=0.03) and D.iloc[0].area > 0 and D.iloc[0].id_status == "MS2 confirmed"
    pcfg["pec_screen"]["rt_shift"] = 0.0                                # no correction: the peak is 0.2 min away, outside +-0.15
    D0, _ = PS.quantify_injections(_cand(), smap, pcfg, ref)
    assert D0.iloc[0].area == 0.0


def test_injection_scale_multiplies_concentration(pcfg, tmp_path):
    ref = PS.reference_response(pcfg)
    smap = _map(tmp_path, [("g1", "grain", "A", "sample", _scaled(0.5, seed=1))]); smap["scale"] = 0.1
    D, _ = PS.quantify_injections(_cand(), smap, pcfg, ref)
    assert D.iloc[0].conc == pytest.approx(0.1 * 0.005 * 166.6667, rel=0.05)


def test_aggregate_fields_uses_detected_injections_only():
    D = pd.DataFrame(dict(cand_id=["a"] * 3, matrix=["water"] * 3, field=["1"] * 3, role=["sample"] * 3, detected=[True, True, False], blank_level=[False, False, True],
                          conc=[100.0, 300.0, 5.0], id_status=["MS2 confirmed", "MS2 confirmed", "MS1 only"]))
    F = PS.aggregate_fields(D).iloc[0]
    assert F.n_inj == 3 and F.n_det == 2 and F.conc == pytest.approx(200.0) and F.blank_level_n == 1


# ---------------- properties -------------------------------------------------------------------------------------------------------
def test_species_and_structure_flags():
    assert species_type(4.5, np.nan) == "Acid" and species_type(np.nan, 9.5) == "Base" and species_type(np.nan, np.nan) == "Neutral"
    assert species_type(4.6, 10.8) == "Zwitterion"
    assert structure_flags("NCC1(CC(O)=O)CCCCC1")["has_acid"] and structure_flags("NCC1(CC(O)=O)CCCCC1")["has_base"]       # gabapentin
    assert structure_flags("CCCCCCCCCCCC[N+](C)(C)C")["quaternary"]


def test_comptox_export_units_and_matching(tmp_path, cfg):
    t = pd.DataFrame({"INPUT": ["Diphenylguanidine"], "INCHIKEY": ["OWRCNXZUPFZXOS-UHFFFAOYSA-N"], "AVERAGE_MASS": [211.268],
                      "VAPOR_PRESSURE_MMHG_OPERA_PRED": [2.53e-07], "WATER_SOLUBILITY_MOL/L_OPERA_PRED": [0.004759],
                      "HENRYS_LAW_ATM-M3/MOLE_OPERA_PRED": [3.71e-09], "OCTANOL_WATER_PARTITION_COEFF_LOGKOW_OPERA_PRED": [2.959],
                      "OPERA_PKAB_OPERA_PRED": [9.84], "BIODEGRADATION_HALF_LIFE_DAYS_OPERA_PRED": [19.8]})
    p = tmp_path / "cx.csv"; t.to_csv(p, index=False)
    tab = read_property_table(str(p))
    assert tab.vp_pa.iloc[0] == pytest.approx(2.53e-07 * 133.322) and tab.solubility_mg_L.iloc[0] == pytest.approx(0.004759 * 211.268 * 1000)
    c = pd.DataFrame([dict(cand_id="x_pos", name="dpg", smiles="N=C(Nc1ccccc1)Nc1ccccc1", ik1="OWRCNXZUPFZXOS", cas="")])
    r = resolve_properties(c, tab, cfg["pec_screen"]).iloc[0]
    assert r.prop_status == "table" and r.species == "Base" and r.kaw == pytest.approx(3.71e-09 / (8.2057e-5 * 298.15), rel=1e-3) and "OPERA" in r.dt50_src
    r2 = resolve_properties(c.assign(ik1="UNKNOWN", name="zzz"), tab, cfg["pec_screen"]).iloc[0]       # not in table -> estimated, flagged, gridded
    assert r2.prop_status == "estimated" and r2.logkow_src.startswith("estimated") and ";" in r2.grid_dt50


def test_zwitterion_and_quaternary_are_not_applicable(cfg):
    c = pd.DataFrame([dict(cand_id="g", name="gabapentin", smiles="NCC1(CC(O)=O)CCCCC1", ik1="UGJMXCAKCUNAIE", cas=""),
                      dict(cand_id="q", name="DTMA", smiles="CCCCCCCCCCCC[N+](C)(C)C", ik1="VICYBMUVWHJEFT", cas="")])
    r = resolve_properties(c, None, cfg["pec_screen"])
    assert (r.prop_status == "not-applicable").all()


# ---------------- fate grid / cache -------------------------------------------------------------------------------------------------
def _props(**kw):
    d = dict(cand_id="c", mw=250.0, logkow=3.0, logkow_src="table", pka_acid=np.nan, pka_base=np.nan, species="Neutral", kaw=1e-8, biowin4=np.nan, dt50_soil_d=30.0,
             k_fruit_user=np.nan, not_applicable="", grid_logkow="2.5;3.0;3.5", grid_kaw="1e-9;1e-6", grid_dt50="10;100", prop_status="table")
    d.update(kw); return pd.Series(d)


def test_fate_point_grid_and_cache(cfg, tmp_path):
    calls = []
    def backend(mw, logkow, pka, base, kaw, dt50_soil_d=None, biowin4=None):
        calls.append(1); return 0.1 / (1 + logkow) * (dt50_soil_d or 30) / 100 / (1 + kaw * 1e4)
    cfg2 = copy.deepcopy(cfg); fr = FateRunner(cfg2, backend=backend, cache_path=str(tmp_path / "k.csv"))
    pt = fr.k_range(_props()); assert pt["k_mode"] == "point" and pt["k_low"] == pt["k_high"] and len(calls) == 1
    fr.k_range(_props()); assert len(calls) == 1                                    # cached
    gr = fr.k_range(_props(kaw=np.nan, dt50_soil_d=np.nan, logkow_src="estimated (Crippen)"))
    assert gr["k_mode"] == "grid" and gr["n_grid"] == 3 * 2 * 2 and gr["k_low"] < gr["k_central"] < gr["k_high"]
    assert fr.k_range(_props(k_fruit_user=0.07))["k_central"] == 0.07
    assert np.isnan(fr.k_range(_props(not_applicable="zwitterion"))["k_central"])
    fr.close(); assert os.path.exists(tmp_path / "k.csv")


# ---------------- classification --------------------------------------------------------------------------------------------------
def _compare(cfg, water, grain, k, applicable=True):
    cfg2 = copy.deepcopy(cfg); cfg2["pec_screen"]["_grain_floor_ng_g"] = 0.05
    cands = pd.DataFrame([dict(cand_id="c", ik1="X")])
    props = pd.DataFrame([dict(cand_id="c", prop_status="table" if applicable else "not-applicable", not_applicable="" if applicable else "zwitterion")])
    F = pd.DataFrame([dict(cand_id="c", matrix="water", field=f, n_inj=1, n_det=1, conc=v, id_status="MS2 confirmed") for f, v in water.items()] +
                     [dict(cand_id="c", matrix="grain", field=f, n_inj=1, n_det=1 if np.isfinite(v) else 0, conc=v, id_status="MS2 confirmed") for f, v in grain.items()])
    return PS.predict_and_compare(cands, F, props, {"c": k}, cfg2)[1].iloc[0]


K = dict(k_low=0.01, k_high=0.01, k_central=0.01, n_grid=1, k_mode="point")


def test_classification_classes(cfg):
    # PEC = 0.01 x 1000 ng/L x 0.005 = 0.05 ng/g
    assert _compare(cfg, {"1": 1000.0}, {"1": 0.1}, K)["class"].startswith("water-explained")          # ratio 0.5
    assert _compare(cfg, {"1": 1000.0}, {"1": 20.0}, K)["class"].startswith("MEC >> PEC")                # ratio 0.0025
    assert _compare(cfg, {"1": 1000.0}, {"1": 0.001}, K)["class"].startswith("PEC >> MEC")               # ratio 50
    assert _compare(cfg, {"1": 1000.0}, {"1": np.nan}, dict(K, k_low=1.0, k_high=1.0, k_central=1.0))["class"] == "predicted but not found in grain"
    assert _compare(cfg, {"1": 1000.0}, {"1": np.nan}, dict(K, k_low=0.005, k_high=0.005, k_central=0.005))["class"].startswith("PEC below grain detection floor")
    assert _compare(cfg, {"1": 1000.0}, {"1": 0.1}, K, applicable=False)["class"] == "fate model not applicable"
    s = _compare(cfg, {"1": 1000.0}, {"1": 20.0}, dict(K, k_low=0.001, k_high=100.0, k_central=0.3, k_mode="grid"))
    assert s["class"].startswith("unconstrained")                                                      # x1e5-wide property range overlaps the window but proves nothing
    s = _compare(cfg, {"1": 1000.0}, {"1": 0.1}, dict(K, k_low=0.005, k_high=0.05, k_central=0.016, k_mode="grid"))
    assert s["class"].startswith("water-explained") and "bounded" in s["class"]                        # narrow bounds inside the window
    assert _compare(cfg, {"1": 1000.0}, {"1": 500.0}, dict(K, k_low=0.001, k_high=1.0, k_central=0.03, k_mode="grid"))["class"].startswith("MEC >> PEC")   # conclusive even with wide bounds
    assert np.isclose(_compare(cfg, {"1": 1000.0, "2": 1000.0}, {"1": 0.05, "2": 0.5}, K).PEC_over_MEC_geomean, np.sqrt(1.0 * 0.1))


def test_soil_sink_quantity_and_broadcast_water(cfg, tmp_path):
    """soil: the backend is asked for quantity='soil' (own cache key); a water row with field '*' pairs with every soil field."""
    seen = []
    def backend(mw, logkow, pka, base, kaw, dt50_soil_d=None, biowin4=None, quantity="fruit"):
        seen.append(quantity); return 0.5 if quantity == "soil" else 0.01
    c = copy.deepcopy(cfg); c["pec_screen"]["_quantity"] = "soil"; c["pec_screen"]["_sink_label"] = "soil"; c["pec_screen"]["fruit_to_grain_dw"] = 0.001
    c["pec_screen"]["_grain_floor_ng_g"] = 0.05
    fr = FateRunner(c, backend=backend, cache_path=str(tmp_path / "k.csv"))
    k = fr.k_range(_props()); assert k["k_central"] == 0.5 and set(seen) == {"soil"}
    cands = pd.DataFrame([dict(cand_id="c", ik1="X")])
    props = pd.DataFrame([dict(cand_id="c", prop_status="table", not_applicable="")])
    F = pd.DataFrame([dict(cand_id="c", matrix="water", field="*", n_inj=4, n_det=4, conc=2000.0, id_status="MS2 confirmed")] +
                     [dict(cand_id="c", matrix="grain", field=f, n_inj=1, n_det=1, conc=v, id_status="MS2 confirmed") for f, v in (("7", 1.0), ("8", 2.0), ("9", 0.5))])
    P, S = PS.predict_and_compare(cands, F, props, {"c": k}, c)
    assert set(P.field) == {"7", "8", "9"} and (P.PEC_central == 0.5 * 2000.0 * 0.001).all() and S.iloc[0].n_pairs == 3
    assert S.iloc[0]["class"].startswith("water-explained")                                   # PEC 1.0 ng/g vs MEC 0.5-2 ng/g
    F2 = F[F.matrix == "water"]; P2, S2 = PS.predict_and_compare(cands, F2, props, {"c": k}, c)
    assert S2.iloc[0]["class"] == "water only (no soil data for these fields)"
